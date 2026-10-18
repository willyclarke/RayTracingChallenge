//! The Metal viewer, `rtc-metal <preview> --window --watch`, run as a child
//! process that ends with the editor.

use std::os::raw::c_int;
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};

/// The viewer program: `env` (`$RTC_METAL`) when set and non-empty,
/// otherwise rtc-metal's release build next to this crate.
pub fn find(env: Option<&str>, crate_dir: &Path) -> Result<PathBuf, String> {
    if let Some(exe) = env.filter(|e| !e.is_empty()) {
        return Ok(PathBuf::from(exe));
    }
    let exe = crate_dir.join("../rtc-metal/target/release/rtc-metal");
    if exe.is_file() {
        Ok(exe)
    } else {
        Err(format!(
            "rtc-metal not found at {}; build it or set RTC_METAL",
            exe.display()
        ))
    }
}

pub struct Viewer {
    child: Child,
}

impl Viewer {
    /// Start `exe <preview> --window --watch`.
    pub fn start(exe: &Path, preview: &Path) -> Result<Viewer, String> {
        Command::new(exe)
            .arg(preview)
            .args(["--window", "--watch"])
            .stdin(Stdio::null())
            .spawn()
            .map(|child| Viewer { child })
            .map_err(|e| format!("starting {}: {e}", exe.display()))
    }

    /// Whether the viewer is still open.
    pub fn running(&mut self) -> bool {
        matches!(self.child.try_wait(), Ok(None))
    }

    pub fn id(&self) -> u32 {
        self.child.id()
    }
}

const SIGTERM: c_int = 15;

unsafe extern "C" {
    fn kill(pid: c_int, sig: c_int) -> c_int;
}

impl Drop for Viewer {
    /// Ask the viewer to quit, so it can offer to save its image; it is not
    /// waited for, as the user may take a while to answer.
    fn drop(&mut self) {
        // SAFETY: `kill` only sends a signal; the pid is our own child's.
        unsafe {
            kill(self.child.id() as c_int, SIGTERM);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::os::unix::fs::PermissionsExt;
    use std::time::{Duration, Instant};

    fn fresh(name: &str) -> PathBuf {
        let root = std::env::temp_dir().join(format!("rtc_edit_viewer_{name}"));
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(&root).unwrap();
        root
    }

    /// An executable shell script at `dir/name`.
    fn script(dir: &Path, name: &str, body: &str) -> PathBuf {
        let path = dir.join(name);
        fs::write(&path, format!("#!/bin/sh\n{body}\n")).unwrap();
        fs::set_permissions(&path, fs::Permissions::from_mode(0o755)).unwrap();
        path
    }

    #[test]
    fn find_prefers_the_environment() {
        let got = find(Some("/x/rtc-metal"), Path::new("/nowhere")).unwrap();
        assert_eq!(got, Path::new("/x/rtc-metal"));
    }

    #[test]
    fn find_falls_back_to_the_sibling_crate() {
        let root = fresh("find");
        let exe = root.join("rtc-metal/target/release/rtc-metal");
        fs::create_dir_all(exe.parent().unwrap()).unwrap();
        fs::write(&exe, "").unwrap();
        let crate_dir = root.join("rtc-edit");
        fs::create_dir_all(&crate_dir).unwrap();
        for env in [None, Some("")] {
            let got = find(env, &crate_dir).unwrap();
            assert_eq!(got.canonicalize().unwrap(), exe.canonicalize().unwrap());
        }
    }

    #[test]
    fn find_reports_a_missing_viewer() {
        let root = fresh("missing");
        let err = find(None, &root.join("rtc-edit")).unwrap_err();
        assert!(err.contains("RTC_METAL"), "{err}");
    }

    #[test]
    fn running_is_true_while_the_child_runs() {
        let root = fresh("runs");
        let exe = script(&root, "viewer", "exec sleep 30");
        let mut v = Viewer::start(&exe, Path::new("p.json")).unwrap();
        assert!(v.running());
    }

    #[test]
    fn running_is_false_after_the_child_exits() {
        let root = fresh("exits");
        let exe = script(&root, "viewer", "exit 0");
        let mut v = Viewer::start(&exe, Path::new("p.json")).unwrap();
        let start = Instant::now();
        while v.running() && start.elapsed() < Duration::from_secs(2) {
            std::thread::sleep(Duration::from_millis(20));
        }
        assert!(!v.running());
    }

    #[test]
    fn dropping_the_viewer_asks_it_to_quit() {
        let root = fresh("quit");
        let got = root.join("got");
        let ready = root.join("ready");
        let body = format!(
            "trap 'echo term > {}; kill $!; exit 0' TERM; : > {}; sleep 30 & wait",
            got.display(),
            ready.display()
        );
        let exe = script(&root, "viewer", &body);
        let v = Viewer::start(&exe, Path::new("p.json")).unwrap();
        // Wait for the shell to install its trap before the signal arrives.
        let start = Instant::now();
        while !ready.exists() {
            assert!(
                start.elapsed() < Duration::from_secs(5),
                "viewer script never installed its trap"
            );
            std::thread::sleep(Duration::from_millis(10));
        }
        drop(v);
        let start = Instant::now();
        while !got.exists() && start.elapsed() < Duration::from_secs(2) {
            std::thread::sleep(Duration::from_millis(20));
        }
        assert!(got.exists(), "viewer never received SIGTERM");
    }

    #[test]
    fn start_reports_a_missing_program() {
        let err = Viewer::start(Path::new("/no/such/rtc-metal"), Path::new("p.json"))
            .err()
            .expect("start should fail");
        assert!(err.contains("/no/such/rtc-metal"), "{err}");
    }
}
