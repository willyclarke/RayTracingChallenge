//! Watches a scene file for `--window --watch`: [`FileEvents`] sleeps on
//! GCD vnode sources until the file (or, while it is missing, its folder)
//! changes, and [`Watch::changed`] then says whether its modification time
//! is new, so the window reloads the scene once per change.

use block2::RcBlock;
use dispatch2::{
    _dispatch_source_type_vnode, DispatchObject, DispatchQueue, DispatchRetained, DispatchSource,
    dispatch_source_type_t, dispatch_source_vnode_flags_t as Vnode,
};
use std::ffi::CString;
use std::fs;
use std::os::unix::ffi::OsStrExt;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, Weak};
use std::time::SystemTime;

/// A file and the modification time last seen.
pub struct Watch {
    path: PathBuf,
    modified: Option<SystemTime>,
}

impl Watch {
    /// Watch `path`, taking its current state as seen.
    pub fn new(path: impl Into<PathBuf>) -> Watch {
        let path = path.into();
        let modified = modified(&path);
        Watch { path, modified }
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    /// True once for each new modification time. A missing or unreadable
    /// file is unchanged, so a file being replaced is picked up once it is
    /// back.
    pub fn changed(&mut self) -> bool {
        match modified(&self.path) {
            Some(t) if self.modified != Some(t) => {
                self.modified = Some(t);
                true
            }
            _ => false,
        }
    }
}

fn modified(path: &Path) -> Option<SystemTime> {
    fs::metadata(path).and_then(|m| m.modified()).ok()
}

/// Calls `notify` (on a queue of its own) whenever the file at a path may
/// have changed: written in place, replaced by a rename (an atomic save),
/// or made again after a delete. Dropping it stops the watching.
pub struct FileEvents {
    inner: Arc<Mutex<Inner>>,
}

struct Inner {
    path: PathBuf,
    notify: Arc<dyn Fn() + Send + Sync>,
    queue: DispatchRetained<DispatchQueue>,
    /// The source on the file, or on its folder while the file is missing.
    source: Option<DispatchRetained<DispatchSource>>,
}

/// What a source watches.
#[derive(PartialEq)]
enum Watching {
    File,
    Folder,
}

impl FileEvents {
    /// Watch `path`; `Err` when neither it nor its folder can be opened.
    pub fn start(
        path: &Path,
        notify: impl Fn() + Send + Sync + 'static,
    ) -> Result<FileEvents, String> {
        let inner = Arc::new(Mutex::new(Inner {
            path: path.to_path_buf(),
            notify: Arc::new(notify),
            queue: DispatchQueue::new("rtc-metal.watch", None),
            source: None,
        }));
        arm(&inner, &mut inner.lock().unwrap())?;
        Ok(FileEvents { inner })
    }
}

impl Drop for FileEvents {
    fn drop(&mut self) {
        if let Some(source) = self.inner.lock().unwrap().source.take() {
            source.cancel();
        }
    }
}

/// Watch the file, or its folder while the file can't be opened.
fn arm(this: &Arc<Mutex<Inner>>, inner: &mut Inner) -> Result<Watching, String> {
    if let Some(source) = inner.source.take() {
        source.cancel();
    }
    let weak = Arc::downgrade(this);
    if let Some(fd) = open(&inner.path) {
        let mask = Vnode::DISPATCH_VNODE_WRITE.0
            | Vnode::DISPATCH_VNODE_EXTEND.0
            | Vnode::DISPATCH_VNODE_ATTRIB.0
            | Vnode::DISPATCH_VNODE_DELETE.0
            | Vnode::DISPATCH_VNODE_RENAME.0;
        inner.source = Some(vnode(fd, mask as usize, &inner.queue, move || {
            on_file_event(&weak)
        }));
        return Ok(Watching::File);
    }
    let folder = folder_of(&inner.path);
    let fd = open(folder).ok_or_else(|| format!("{}: can't watch it", folder.display()))?;
    inner.source = Some(vnode(
        fd,
        Vnode::DISPATCH_VNODE_WRITE.0 as usize,
        &inner.queue,
        move || on_folder_event(&weak),
    ));
    // The file may have come back before the folder was watched.
    if inner.path.exists() {
        return arm(this, inner);
    }
    Ok(Watching::Folder)
}

/// The file changed: notify, and after a delete or rename watch whatever
/// now has its path.
fn on_file_event(weak: &Weak<Mutex<Inner>>) {
    let Some(this) = weak.upgrade() else {
        return;
    };
    let mut inner = this.lock().unwrap();
    let data = inner.source.as_ref().map_or(0, |s| s.data());
    let gone = (Vnode::DISPATCH_VNODE_DELETE.0 | Vnode::DISPATCH_VNODE_RENAME.0) as usize;
    let changed = if data & gone != 0 {
        rearm(&this, &mut inner)
    } else {
        true
    };
    let notify = inner.notify.clone();
    drop(inner);
    if changed {
        notify();
    }
}

/// The folder changed: once the file is back, watch it and notify.
fn on_folder_event(weak: &Weak<Mutex<Inner>>) {
    let Some(this) = weak.upgrade() else {
        return;
    };
    let mut inner = this.lock().unwrap();
    if !inner.path.exists() {
        return;
    }
    let changed = rearm(&this, &mut inner);
    let notify = inner.notify.clone();
    drop(inner);
    if changed {
        notify();
    }
}

/// Watch the path anew; true when the file is there to reload.
fn rearm(this: &Arc<Mutex<Inner>>, inner: &mut Inner) -> bool {
    match arm(this, inner) {
        Ok(watching) => watching == Watching::File,
        Err(e) => {
            eprintln!("rtc-metal: --watch stopped: {e}");
            false
        }
    }
}

/// The folder `path` is in: `.` for a bare file name, whose parent is "".
fn folder_of(path: &Path) -> &Path {
    match path.parent() {
        Some(p) if !p.as_os_str().is_empty() => p,
        _ => Path::new("."),
    }
}

/// `path` opened for events only (it doesn't keep its volume from
/// unmounting).
fn open(path: &Path) -> Option<i32> {
    let c = CString::new(path.as_os_str().as_bytes()).ok()?;
    // SAFETY: a valid C string; the descriptor is closed by its source.
    let fd = unsafe { libc::open(c.as_ptr(), libc::O_EVTONLY) };
    (fd >= 0).then_some(fd)
}

/// An active vnode source on `fd` for `mask`, calling `handler` on
/// `queue`; it closes `fd` when cancelled.
fn vnode(
    fd: i32,
    mask: usize,
    queue: &DispatchQueue,
    handler: impl Fn() + Send + Sync + 'static,
) -> DispatchRetained<DispatchSource> {
    let kind: dispatch_source_type_t = (&raw const _dispatch_source_type_vnode).cast_mut();
    // SAFETY: a vnode source on an open descriptor with vnode flags; the
    // handlers are `Send + Sync` and dispatch copies the blocks.
    unsafe {
        let source = DispatchSource::new(kind, fd as usize, mask, Some(queue));
        let event: RcBlock<dyn Fn()> = RcBlock::new(handler);
        source.set_event_handler_with_block(RcBlock::as_ptr(&event));
        let cancel: RcBlock<dyn Fn()> = RcBlock::new(move || {
            libc::close(fd);
        });
        source.set_cancel_handler_with_block(RcBlock::as_ptr(&cancel));
        source.activate();
        source
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::{self, File};
    use std::time::{Duration, SystemTime};

    fn touch(path: &Path, t: SystemTime) -> std::io::Result<()> {
        File::options().write(true).open(path)?.set_modified(t)
    }

    #[test]
    fn changed_reports_each_new_modification_once() -> std::io::Result<()> {
        let path = std::env::temp_dir().join("rtc_metal_watch_changed.json");
        fs::write(&path, "{}")?;
        let t0 = SystemTime::now();
        touch(&path, t0)?;
        let mut watch = Watch::new(&path);
        assert!(!watch.changed());
        touch(&path, t0 + Duration::from_secs(1))?;
        assert!(watch.changed());
        assert!(!watch.changed());
        touch(&path, t0 + Duration::from_secs(1) + Duration::from_nanos(1))?;
        assert!(watch.changed());
        assert_eq!(watch.path(), path);
        Ok(())
    }

    #[test]
    fn a_missing_file_is_unchanged_until_it_reappears() -> std::io::Result<()> {
        let path = std::env::temp_dir().join("rtc_metal_watch_missing.json");
        fs::write(&path, "{}")?;
        let t0 = SystemTime::now();
        touch(&path, t0)?;
        let mut watch = Watch::new(&path);
        fs::remove_file(&path)?;
        assert!(!watch.changed());
        assert!(!watch.changed());
        fs::write(&path, "{}")?;
        touch(&path, t0 + Duration::from_secs(1))?;
        assert!(watch.changed());
        Ok(())
    }

    #[test]
    fn watching_a_file_that_does_not_exist_yet() -> std::io::Result<()> {
        let path = std::env::temp_dir().join("rtc_metal_watch_later.json");
        let _ = fs::remove_file(&path);
        let mut watch = Watch::new(&path);
        assert!(!watch.changed());
        fs::write(&path, "{}")?;
        assert!(watch.changed());
        Ok(())
    }

    /// A fresh folder for `name`, and `FileEvents` on `scene.json` in it
    /// whose notifications arrive on the receiver.
    fn events(name: &str) -> (PathBuf, FileEvents, std::sync::mpsc::Receiver<()>) {
        let dir =
            std::env::temp_dir().join(format!("rtc_metal_events_{name}_{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join("scene.json");
        fs::write(&path, "{}").unwrap();
        let (tx, rx) = std::sync::mpsc::channel();
        let events = FileEvents::start(&path, move || {
            let _ = tx.send(());
        })
        .unwrap();
        (path, events, rx)
    }

    /// Whether a notification comes within two seconds.
    fn next(rx: &std::sync::mpsc::Receiver<()>) -> bool {
        rx.recv_timeout(Duration::from_secs(2)).is_ok()
    }

    /// Drop notifications already sent.
    fn drain(rx: &std::sync::mpsc::Receiver<()>) {
        std::thread::sleep(Duration::from_millis(200));
        while rx.try_recv().is_ok() {}
    }

    /// Write `path` the way an editor's atomic save does.
    fn save_atomically(path: &Path, text: &str) {
        let tmp = path.with_extension("tmp");
        fs::write(&tmp, text).unwrap();
        fs::rename(&tmp, path).unwrap();
    }

    #[test]
    fn a_write_in_place_is_reported() {
        let (path, _events, rx) = events("in_place");
        fs::write(&path, "{ }").unwrap();
        assert!(next(&rx));
    }

    #[test]
    fn a_new_file_renamed_over_the_scene_is_reported_and_watched() {
        let (path, _events, rx) = events("renamed");
        save_atomically(&path, "{ }");
        assert!(next(&rx));
        drain(&rx);
        save_atomically(&path, "{  }");
        assert!(next(&rx));
        drain(&rx);
        fs::write(&path, "{   }").unwrap();
        assert!(next(&rx));
    }

    #[test]
    fn a_deleted_file_made_again_is_reported() {
        let (path, _events, rx) = events("deleted");
        fs::remove_file(&path).unwrap();
        drain(&rx);
        assert!(rx.recv_timeout(Duration::from_millis(300)).is_err());
        fs::write(&path, "{}").unwrap();
        assert!(next(&rx));
        drain(&rx);
        fs::write(&path, "{ }").unwrap();
        assert!(next(&rx));
    }

    #[test]
    fn a_bare_file_name_is_in_the_working_folder() {
        assert_eq!(folder_of(Path::new("scene.json")), Path::new("."));
        assert_eq!(folder_of(Path::new("scenes/a.json")), Path::new("scenes"));
        assert_eq!(folder_of(Path::new("/a.json")), Path::new("/"));
    }

    #[test]
    fn start_fails_for_a_missing_folder() {
        assert!(FileEvents::start(Path::new("/nonexistent/rtc/scene.json"), || {}).is_err());
    }
}
