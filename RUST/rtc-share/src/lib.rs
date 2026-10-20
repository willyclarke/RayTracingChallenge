//! The viewer's camera, shared between rtc-metal and rtc-edit through a
//! small memory-mapped file.
//!
//! The editor creates the file and starts the viewer with its path; the
//! viewer writes its camera whenever it changes, and the editor reads it
//! when it likes. Every read and write holds an `flock` on the file, which
//! the kernel releases if a process dies. The file starts with
//! [`SIGNATURE`], a hash of this source file, so a viewer built from
//! another version of this crate refuses the file instead of misreading it.
//! Keep this crate to this one file, so the hash covers the whole layout.

use std::ffi::c_void;
use std::fs::{File, OpenOptions};
use std::os::fd::AsRawFd;
use std::path::Path;

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

/// The file's layout.
#[repr(C)]
#[derive(Clone, Copy)]
struct View {
    signature: u64,
    size: u64,
    /// +1 on every write; 0 until the first.
    seq: u64,
    from: [f64; 3],
    to: [f64; 3],
    up: [f64; 3],
    focal_distance: f64,
}

const SIZE: usize = size_of::<View>();

/// The camera as last written, with the write's number.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Snapshot {
    pub seq: u64,
    pub from: [f64; 3],
    pub to: [f64; 3],
    pub up: [f64; 3],
    pub focal_distance: f64,
}

const PROT_READ: i32 = 1;
const PROT_WRITE: i32 = 2;
const MAP_SHARED: i32 = 1;
const LOCK_SH: i32 = 1;
const LOCK_EX: i32 = 2;
const LOCK_UN: i32 = 8;

unsafe extern "C" {
    fn mmap(addr: *mut c_void, len: usize, prot: i32, flags: i32, fd: i32, off: i64)
    -> *mut c_void;
    fn munmap(addr: *mut c_void, len: usize) -> i32;
    fn flock(fd: i32, op: i32) -> i32;
}

/// The shared file, mapped.
pub struct Shared {
    file: File,
    view: *mut View,
}

// SAFETY: the mapping is plain memory owned by this value, and every access
// goes through `flock`.
unsafe impl Send for Shared {}

impl Shared {
    /// Make the file at `path`, replacing any old one (a viewer still
    /// quitting keeps that), with the header written and no camera yet.
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
        file.set_len(SIZE as u64).map_err(err)?;
        let shared = Shared::map(file, path)?;
        shared.locked(LOCK_EX, |v| {
            *v = View {
                signature: SIGNATURE,
                size: SIZE as u64,
                seq: 0,
                from: [0.0; 3],
                to: [0.0; 3],
                up: [0.0; 3],
                focal_distance: 0.0,
            }
        });
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
        let (signature, size) = shared.locked(LOCK_SH, |v| (v.signature, v.size));
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
            file,
            view: p.cast(),
        })
    }

    /// Run `f` on a copy of the view under the lock `op`, writing the copy
    /// back under `LOCK_EX`.
    fn locked<T>(&self, op: i32, f: impl FnOnce(&mut View) -> T) -> T {
        let fd = self.file.as_raw_fd();
        // SAFETY: `flock` on our own open file; `view` points at our live
        // mapping, accessed only under the lock and with volatile copies,
        // as another process writes it.
        unsafe {
            flock(fd, op);
            let mut v = std::ptr::read_volatile(self.view);
            let out = f(&mut v);
            if op == LOCK_EX {
                std::ptr::write_volatile(self.view, v);
            }
            flock(fd, LOCK_UN);
            out
        }
    }

    /// Write the camera, as the next write.
    pub fn write(&self, from: [f64; 3], to: [f64; 3], up: [f64; 3], focal_distance: f64) {
        self.locked(LOCK_EX, |v| {
            v.seq += 1;
            (v.from, v.to, v.up, v.focal_distance) = (from, to, up, focal_distance);
        });
    }

    /// The camera last written; `None` before the first write.
    pub fn read(&self) -> Option<Snapshot> {
        let v = self.locked(LOCK_SH, |v| *v);
        (v.seq != 0).then_some(Snapshot {
            seq: v.seq,
            from: v.from,
            to: v.to,
            up: v.up,
            focal_distance: v.focal_distance,
        })
    }
}

impl Drop for Shared {
    fn drop(&mut self) {
        // SAFETY: the mapping `map` made, unmapped once.
        unsafe { munmap(self.view.cast(), SIZE) };
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
        assert_eq!(Shared::create(&p)?.read(), None);
        std::fs::remove_file(&p).unwrap();
        Ok(())
    }

    #[test]
    fn a_write_through_one_mapping_is_read_through_another() -> Result<(), String> {
        let p = path("two");
        let a = Shared::create(&p)?;
        let b = Shared::open(&p)?;
        a.write([1., 2., 3.], [0., 1., 0.], [0., 1., 0.], 5.);
        assert_eq!(
            b.read(),
            Some(Snapshot {
                seq: 1,
                from: [1., 2., 3.],
                to: [0., 1., 0.],
                up: [0., 1., 0.],
                focal_distance: 5.,
            })
        );
        a.write([4., 5., 6.], [0., 1., 0.], [0., 1., 0.], 5.);
        assert_eq!(b.read().map(|s| (s.seq, s.from)), Some((2, [4., 5., 6.])));
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
        old.write([1., 2., 3.], [0., 0., 0.], [0., 1., 0.], 1.);
        assert_eq!(Shared::create(&p)?.read(), None);
        std::fs::remove_file(&p).unwrap();
        Ok(())
    }
}
