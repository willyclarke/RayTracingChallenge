//! Polls a file's modification time, for `--window --watch`: the window
//! asks [`Watch::changed`] a few times a second and reloads the scene when
//! it says yes.

use std::fs;
use std::path::{Path, PathBuf};
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
}
