//! The viewer's camera and the editor's load commands, shared between
//! rtc-metal and rtc-edit through a small memory-mapped file.
//!
//! The editor creates the file and starts the viewer with its path; the
//! viewer writes its camera whenever it changes, and the editor sends it a
//! scene to load. Each block has one writer and is a seqlock: its `seq` is
//! odd while being written, and a reader copies the block and retries
//! while `seq` is odd or changed. After a write the writer wakes the other
//! side, which sleeps in `os_sync_wait_on_address` (macOS 14.4+) until the
//! `seq` it saw changes. The file starts with [`SIGNATURE`], a hash of
//! this source file, so a viewer built from another version of this crate
//! refuses the file instead of misreading it. Keep this crate to this one
//! file, so the hash covers the whole layout.

use std::ffi::c_void;
use std::fs::{File, OpenOptions};
use std::os::fd::AsRawFd;
use std::os::unix::ffi::OsStrExt;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU32, AtomicU64, Ordering, fence};

/// A seqlock's counter: read it, and move it on by one (each write moves
/// it twice, odd while the block is written).
trait SeqWord {
    fn load(&self, order: Ordering) -> u32;
    fn bump(&self, order: Ordering);
}

impl SeqWord for AtomicU32 {
    fn load(&self, order: Ordering) -> u32 {
        AtomicU32::load(self, order)
    }

    /// Only the block's one writer bumps it.
    fn bump(&self, order: Ordering) {
        let s = AtomicU32::load(self, Ordering::Relaxed);
        self.store(s.wrapping_add(1), order);
    }
}

/// The camera seqlock, the low half of the camera word, whose high half
/// `end_viewer` may change at the same time: every access is to the whole
/// word, as atomics of two sizes on one word are undefined in Rust.
struct CameraSeq<'a>(&'a AtomicU64);

impl SeqWord for CameraSeq<'_> {
    fn load(&self, order: Ordering) -> u32 {
        self.0.load(order) as u32
    }

    /// The low half wraps without carrying into the generation.
    fn bump(&self, order: Ordering) {
        let next = |w: u64| (w & !0xffff_ffff) | (w as u32).wrapping_add(1) as u64;
        let _ = self
            .0
            .fetch_update(order, Ordering::Relaxed, |w| Some(next(w)));
    }
}

/// FNV-1a, 64 bits: a quick hash for noticing a change, not for security.
pub const fn fnv1a(s: &str) -> u64 {
    let bytes = s.as_bytes();
    let mut hash = 0xcbf29ce484222325u64;
    let mut i = 0;
    while i < bytes.len() {
        hash ^= bytes[i] as u64;
        hash = hash.wrapping_mul(0x100000001b3);
        i += 1;
    }
    hash
}

/// The hash of this file; it changes whenever the file does.
pub const SIGNATURE: u64 = fnv1a(include_str!("lib.rs"));

/// The longest load path, in bytes.
const PATH_MAX: usize = 1024;

/// The file's layout.
#[repr(C)]
struct Layout {
    signature: u64,
    size: u64,
    /// The camera's seqlock (low half) and the viewer generation (high
    /// half), in one word so a camera waiter wakes on either.
    camera_word: u64,
    from: [f64; 3],
    to: [f64; 3],
    up: [f64; 3],
    focal_distance: f64,
    /// The load's seqlock.
    load_seq: u32,
    fresh: u32,
    len: u32,
    path: [u8; PATH_MAX],
    /// One more than the generation of the viewer that last said it takes
    /// loads (0: none has); written by the viewer.
    ready: u32,
}

const SIZE: usize = size_of::<Layout>();

/// The camera as last written, with the write's number.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Snapshot {
    pub seq: u32,
    pub from: [f64; 3],
    pub to: [f64; 3],
    pub up: [f64; 3],
    pub focal_distance: f64,
}

/// A scene the editor asks the viewer to load: from its own camera when
/// `fresh` (a new scene), else keeping the view (edits).
#[derive(Clone, Debug, PartialEq)]
pub struct Load {
    pub seq: u32,
    pub path: PathBuf,
    pub fresh: bool,
}

const PROT_READ: i32 = 1;
const PROT_WRITE: i32 = 2;
const MAP_SHARED: i32 = 1;
const OS_SYNC_WAIT_ON_ADDRESS_SHARED: u32 = 1;
const OS_SYNC_WAKE_BY_ADDRESS_SHARED: u32 = 1;
const EINTR: i32 = 4;

unsafe extern "C" {
    fn mmap(addr: *mut c_void, len: usize, prot: i32, flags: i32, fd: i32, off: i64)
    -> *mut c_void;
    fn munmap(addr: *mut c_void, len: usize) -> i32;
    fn os_sync_wait_on_address(addr: *mut c_void, value: u64, size: usize, flags: u32) -> i32;
    fn os_sync_wake_by_address_all(addr: *mut c_void, size: usize, flags: u32) -> i32;
}

/// Sleep while the `size` bytes at `addr` hold `value` (or until a wake).
fn wait(addr: *mut c_void, value: u64, size: usize) -> Result<(), String> {
    // SAFETY: `addr` is a word of our live shared mapping.
    if unsafe { os_sync_wait_on_address(addr, value, size, OS_SYNC_WAIT_ON_ADDRESS_SHARED) } >= 0 {
        return Ok(());
    }
    let e = std::io::Error::last_os_error();
    if e.raw_os_error() == Some(EINTR) {
        Ok(())
    } else {
        Err(format!("waiting on the shared file: {e}"))
    }
}

/// Wake every process sleeping on the `size` bytes at `addr`.
fn wake(addr: *mut c_void, size: usize) {
    // SAFETY: as `wait`. A -1 (no one waiting) is fine.
    unsafe { os_sync_wake_by_address_all(addr, size, OS_SYNC_WAKE_BY_ADDRESS_SHARED) };
}

/// The shared file, mapped.
pub struct Shared {
    _file: File,
    map: *mut Layout,
}

// SAFETY: the mapping is plain memory owned by this value; every access is
// atomic or goes through a seqlock.
unsafe impl Send for Shared {}
unsafe impl Sync for Shared {}

impl Shared {
    /// Make the file at `path`, replacing any old one (a viewer still
    /// quitting keeps that), with the header written and nothing sent yet.
    pub fn create(path: &Path) -> Result<Shared, String> {
        let err = |e: std::io::Error| format!("{}: {e}", path.display());
        match std::fs::remove_file(path) {
            Err(e) if e.kind() != std::io::ErrorKind::NotFound => return Err(err(e)),
            _ => {}
        }
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create_new(true)
            .open(path)
            .map_err(err)?;
        // A new file reads as zeros: every word 0, nothing written.
        file.set_len(SIZE as u64).map_err(err)?;
        let shared = Shared::map(file, path)?;
        // SAFETY: our fresh mapping; no one else has the file's path yet.
        unsafe {
            (*shared.map).signature = SIGNATURE;
            (*shared.map).size = SIZE as u64;
        }
        Ok(shared)
    }

    /// Map the file at `path`, made by `create` from the same version of
    /// this crate.
    pub fn open(path: &Path) -> Result<Shared, String> {
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .open(path)
            .map_err(|e| format!("{}: {e}", path.display()))?;
        let len = file
            .metadata()
            .map_err(|e| format!("{}: {e}", path.display()))?
            .len();
        let stale = || {
            format!(
                "{}: made by another version; rebuild rtc-metal",
                path.display()
            )
        };
        if len < SIZE as u64 {
            return Err(stale());
        }
        let shared = Shared::map(file, path)?;
        // SAFETY: the header is written once, by `create`, before the path
        // is handed on.
        let (signature, size) = unsafe { ((*shared.map).signature, (*shared.map).size) };
        if signature != SIGNATURE || size != SIZE as u64 {
            return Err(stale());
        }
        Ok(shared)
    }

    fn map(file: File, path: &Path) -> Result<Shared, String> {
        // SAFETY: a fresh shared mapping of `SIZE` bytes of a file that has
        // at least that many; unmapped in `drop`.
        let p = unsafe {
            mmap(
                std::ptr::null_mut(),
                SIZE,
                PROT_READ | PROT_WRITE,
                MAP_SHARED,
                file.as_raw_fd(),
                0,
            )
        };
        if p as usize == usize::MAX {
            return Err(format!(
                "{}: {}",
                path.display(),
                std::io::Error::last_os_error()
            ));
        }
        Ok(Shared {
            _file: file,
            map: p.cast(),
        })
    }

    /// The camera seqlock and generation word.
    fn camera_word(&self) -> &AtomicU64 {
        // SAFETY: an aligned u64 of our live mapping, only accessed atomically.
        unsafe { AtomicU64::from_ptr(&raw mut (*self.map).camera_word) }
    }

    /// The camera seqlock: the camera word's low half.
    fn camera_seq(&self) -> CameraSeq<'_> {
        CameraSeq(self.camera_word())
    }

    fn ready_word(&self) -> &AtomicU32 {
        // SAFETY: an aligned u32 of our live mapping, only accessed atomically.
        unsafe { AtomicU32::from_ptr(&raw mut (*self.map).ready) }
    }

    fn load_seq_word(&self) -> &AtomicU32 {
        // SAFETY: an aligned u32 of our live mapping, only accessed atomically.
        unsafe { AtomicU32::from_ptr(&raw mut (*self.map).load_seq) }
    }

    /// Write a block under the seqlock `seq`: odd while `write` runs.
    fn seq_write(seq: &impl SeqWord, write: impl FnOnce()) {
        seq.bump(Ordering::Relaxed);
        fence(Ordering::Release);
        write();
        seq.bump(Ordering::Release);
    }

    /// Copy a block under the seqlock `seq`, with the `seq` it was copied at.
    fn seq_read<T>(seq: &impl SeqWord, read: impl Fn() -> T) -> (u32, T) {
        loop {
            let s = seq.load(Ordering::Acquire);
            if !s.is_multiple_of(2) {
                std::hint::spin_loop();
                continue;
            }
            let t = read();
            fence(Ordering::Acquire);
            if seq.load(Ordering::Relaxed) == s {
                return (s, t);
            }
        }
    }

    /// Write the camera, as the next write, and wake the editor.
    pub fn write_camera(&self, from: [f64; 3], to: [f64; 3], up: [f64; 3], focal_distance: f64) {
        let m = self.map;
        Shared::seq_write(&self.camera_seq(), || {
            // SAFETY: the viewer is the camera's one writer; readers retry.
            unsafe {
                (&raw mut (*m).from).write_volatile(from);
                (&raw mut (*m).to).write_volatile(to);
                (&raw mut (*m).up).write_volatile(up);
                (&raw mut (*m).focal_distance).write_volatile(focal_distance);
            }
        });
        wake(self.camera_word().as_ptr().cast(), 8);
    }

    /// The camera last written; `None` before the first write.
    pub fn read_camera(&self) -> Option<Snapshot> {
        let m = self.map;
        let (seq, (from, to, up, focal_distance)) = Shared::seq_read(&self.camera_seq(), || {
            // SAFETY: a copy, checked against the seqlock.
            unsafe {
                (
                    (&raw const (*m).from).read_volatile(),
                    (&raw const (*m).to).read_volatile(),
                    (&raw const (*m).up).read_volatile(),
                    (&raw const (*m).focal_distance).read_volatile(),
                )
            }
        });
        (seq != 0).then_some(Snapshot {
            seq,
            from,
            to,
            up,
            focal_distance,
        })
    }

    /// Sleep until a camera other than write `last` is written, and return
    /// it; `None` once the viewer of `generation` has ended.
    pub fn wait_camera(&self, last: u32, generation: u32) -> Result<Option<Snapshot>, String> {
        let word = self.camera_word();
        loop {
            let w = word.load(Ordering::Acquire);
            if (w >> 32) as u32 != generation {
                return Ok(None);
            }
            let seq = w as u32;
            if seq.is_multiple_of(2) && seq != last {
                if let Some(s) = self.read_camera().filter(|s| s.seq != last) {
                    return Ok(Some(s));
                }
                continue;
            }
            wait(word.as_ptr().cast(), w, 8)?;
        }
    }

    /// The current viewer's number; `end_viewer` moves it on.
    pub fn generation(&self) -> u32 {
        (self.camera_word().load(Ordering::Acquire) >> 32) as u32
    }

    /// The viewer says it takes loads, from the `load_seq` it has seen on:
    /// call after reading that.
    pub fn viewer_ready(&self) {
        self.ready_word()
            .store(self.generation().wrapping_add(1), Ordering::Release);
    }

    /// Whether the current viewer has said it takes loads: one sent before
    /// then may be missed.
    pub fn is_ready(&self) -> bool {
        self.ready_word().load(Ordering::Acquire) == self.generation().wrapping_add(1)
    }

    /// The viewer has ended: wake its camera waiters, which then return
    /// `None`.
    pub fn end_viewer(&self) {
        self.camera_word().fetch_add(1 << 32, Ordering::Release);
        wake(self.camera_word().as_ptr().cast(), 8);
    }

    /// Ask the viewer to load the scene at `path`, `fresh` for a new
    /// scene; `Err` for a path over 1024 bytes.
    pub fn send_load(&self, path: &Path, fresh: bool) -> Result<(), String> {
        let bytes = path.as_os_str().as_bytes();
        if bytes.len() > PATH_MAX {
            return Err(format!("{}: longer than {PATH_MAX} bytes", path.display()));
        }
        let m = self.map;
        Shared::seq_write(self.load_seq_word(), || {
            let mut buf = [0u8; PATH_MAX];
            buf[..bytes.len()].copy_from_slice(bytes);
            // SAFETY: the editor is the load's one writer; readers retry.
            unsafe {
                (&raw mut (*m).fresh).write_volatile(fresh as u32);
                (&raw mut (*m).len).write_volatile(bytes.len() as u32);
                (&raw mut (*m).path).write_volatile(buf);
            }
        });
        wake(self.load_seq_word().as_ptr().cast(), 4);
        Ok(())
    }

    /// The number of the last load sent; 0 before the first.
    pub fn load_seq(&self) -> u32 {
        self.load_seq_word().load(Ordering::Acquire)
    }

    /// Sleep until a load other than `last` is sent, and return it.
    pub fn wait_load(&self, last: u32) -> Result<Load, String> {
        let seq = self.load_seq_word();
        let m = self.map;
        loop {
            let s = seq.load(Ordering::Acquire);
            if s.is_multiple_of(2) && s != last {
                let (s, (fresh, len, path)) = Shared::seq_read(seq, || {
                    // SAFETY: a copy, checked against the seqlock.
                    unsafe {
                        (
                            (&raw const (*m).fresh).read_volatile(),
                            (&raw const (*m).len).read_volatile(),
                            (&raw const (*m).path).read_volatile(),
                        )
                    }
                });
                if s == last {
                    continue;
                }
                let len = (len as usize).min(PATH_MAX);
                return Ok(Load {
                    seq: s,
                    path: PathBuf::from(std::ffi::OsStr::from_bytes(&path[..len])),
                    fresh: fresh != 0,
                });
            }
            wait(seq.as_ptr().cast(), s as u64, 4)?;
        }
    }
}

impl Drop for Shared {
    fn drop(&mut self) {
        // SAFETY: the mapping `map` made, unmapped once.
        unsafe { munmap(self.map.cast(), SIZE) };
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    /// A fresh path for `name` in this process's test folder.
    fn path(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("rtc-share-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join(name);
        let _ = std::fs::remove_file(&p);
        p
    }

    #[test]
    fn fnv1a_matches_known_answers() {
        assert_eq!(fnv1a(""), 0xcbf29ce484222325);
        assert_eq!(fnv1a("a"), 0xaf63dc4c8601ec8c);
    }

    #[test]
    fn signature_is_the_hash_of_this_file() {
        assert_eq!(SIGNATURE, fnv1a(include_str!("lib.rs")));
    }

    #[test]
    fn read_is_none_until_the_first_write() -> Result<(), String> {
        let p = path("none");
        assert_eq!(Shared::create(&p)?.read_camera(), None);
        std::fs::remove_file(&p).unwrap();
        Ok(())
    }

    #[test]
    fn a_write_through_one_mapping_is_read_through_another() -> Result<(), String> {
        let p = path("two");
        let a = Shared::create(&p)?;
        let b = Shared::open(&p)?;
        a.write_camera([1., 2., 3.], [0., 1., 0.], [0., 1., 0.], 5.);
        assert_eq!(
            b.read_camera(),
            Some(Snapshot {
                seq: 2,
                from: [1., 2., 3.],
                to: [0., 1., 0.],
                up: [0., 1., 0.],
                focal_distance: 5.,
            })
        );
        a.write_camera([4., 5., 6.], [0., 1., 0.], [0., 1., 0.], 5.);
        assert_eq!(
            b.read_camera().map(|s| (s.seq, s.from)),
            Some((4, [4., 5., 6.]))
        );
        std::fs::remove_file(&p).unwrap();
        Ok(())
    }

    #[test]
    fn open_refuses_another_signature() -> Result<(), String> {
        use std::io::{Seek, Write};
        let p = path("signature");
        let _keep = Shared::create(&p)?;
        let mut f = std::fs::OpenOptions::new().write(true).open(&p).unwrap();
        f.seek(std::io::SeekFrom::Start(0)).unwrap();
        f.write_all(&(!SIGNATURE).to_ne_bytes()).unwrap();
        let e = Shared::open(&p).err().expect("opened");
        assert!(e.contains("rebuild rtc-metal"), "{e}");
        std::fs::remove_file(&p).unwrap();
        Ok(())
    }

    #[test]
    fn open_refuses_a_short_or_missing_file() {
        let p = path("short");
        std::fs::write(&p, [0u8; 10]).unwrap();
        assert!(Shared::open(&p).is_err());
        std::fs::remove_file(&p).unwrap();
        assert!(Shared::open(&path("missing")).is_err());
    }

    #[test]
    fn create_replaces_an_old_file() -> Result<(), String> {
        let p = path("replace");
        let old = Shared::create(&p)?;
        old.write_camera([1., 2., 3.], [0., 0., 0.], [0., 1., 0.], 1.);
        assert_eq!(Shared::create(&p)?.read_camera(), None);
        std::fs::remove_file(&p).unwrap();
        Ok(())
    }

    #[test]
    fn a_sent_load_reads_back_with_its_path_and_fresh() -> Result<(), String> {
        let p = path("load");
        let a = Shared::create(&p)?;
        let b = Shared::open(&p)?;
        a.send_load(Path::new("/x/y.json"), true)?;
        assert_eq!(
            b.wait_load(0)?,
            Load {
                seq: 2,
                path: PathBuf::from("/x/y.json"),
                fresh: true,
            }
        );
        a.send_load(Path::new("/z.json"), false)?;
        let l = b.wait_load(2)?;
        assert_eq!((l.path, l.fresh), (PathBuf::from("/z.json"), false));
        assert_eq!(b.load_seq(), 4);
        std::fs::remove_file(&p).unwrap();
        Ok(())
    }

    #[test]
    fn wait_load_sleeps_until_a_send() -> Result<(), String> {
        let p = path("wait-load");
        let a = Shared::create(&p)?;
        let b = Shared::open(&p)?;
        let last = b.load_seq();
        let waiter = std::thread::spawn(move || b.wait_load(last));
        std::thread::sleep(std::time::Duration::from_millis(50));
        assert!(!waiter.is_finished());
        a.send_load(Path::new("/s.json"), false)?;
        let l = waiter.join().unwrap()?;
        assert_eq!(l.path, PathBuf::from("/s.json"));
        std::fs::remove_file(&p).unwrap();
        Ok(())
    }

    #[test]
    fn wait_camera_returns_a_new_write() -> Result<(), String> {
        let p = path("wait-camera");
        let a = Shared::create(&p)?;
        let b = Shared::open(&p)?;
        let g = b.generation();
        let waiter = std::thread::spawn(move || b.wait_camera(0, g));
        std::thread::sleep(std::time::Duration::from_millis(50));
        a.write_camera([1., 2., 3.], [0., 0., 0.], [0., 1., 0.], 4.);
        let s = waiter.join().unwrap()?.expect("a camera");
        assert_eq!((s.seq, s.from, s.focal_distance), (2, [1., 2., 3.], 4.));
        std::fs::remove_file(&p).unwrap();
        Ok(())
    }

    #[test]
    fn wait_camera_ends_when_the_viewer_ends() -> Result<(), String> {
        let p = path("end");
        let a = Shared::create(&p)?;
        let b = Shared::open(&p)?;
        let c = Shared::open(&p)?;
        let g = b.generation();
        let waiter = std::thread::spawn(move || b.wait_camera(0, g));
        std::thread::sleep(std::time::Duration::from_millis(50));
        assert!(!waiter.is_finished());
        a.end_viewer();
        assert_eq!(waiter.join().unwrap()?, None);
        assert_eq!(c.wait_camera(0, g)?, None);
        std::fs::remove_file(&p).unwrap();
        Ok(())
    }

    #[test]
    fn a_viewer_is_ready_from_saying_so_until_it_ends() -> Result<(), String> {
        let p = path("ready");
        let editor = Shared::create(&p)?;
        let viewer = Shared::open(&p)?;
        assert!(!editor.is_ready());
        viewer.viewer_ready();
        assert!(editor.is_ready());
        editor.end_viewer();
        assert!(!editor.is_ready());
        Shared::open(&p)?.viewer_ready();
        assert!(editor.is_ready());
        std::fs::remove_file(&p).unwrap();
        Ok(())
    }

    #[test]
    fn a_reader_never_sees_a_half_written_camera() -> Result<(), String> {
        let p = path("torn");
        let a = Shared::create(&p)?;
        let b = Shared::open(&p)?;
        let writer = std::thread::spawn(move || {
            for k in 1..=20_000 {
                let k = k as f64;
                a.write_camera([k; 3], [k; 3], [k; 3], k);
            }
        });
        loop {
            let Some(s) = b.read_camera() else { continue };
            let k = s.focal_distance;
            assert!(
                [s.from, s.to, s.up].iter().flatten().all(|&v| v == k),
                "{s:?}"
            );
            if k == 20_000. {
                break;
            }
        }
        writer.join().unwrap();
        std::fs::remove_file(&p).unwrap();
        Ok(())
    }

    #[test]
    fn a_path_over_1024_bytes_is_refused() -> Result<(), String> {
        let p = path("long");
        let a = Shared::create(&p)?;
        let long = format!("/{}", "a".repeat(1024));
        let e = a.send_load(Path::new(&long), true).unwrap_err();
        assert!(e.contains("1024"), "{e}");
        assert_eq!(a.load_seq(), 0);
        a.send_load(Path::new(&long[..1024]), true)?;
        std::fs::remove_file(&p).unwrap();
        Ok(())
    }
}
