//! The viewer run as the scene editor runs it.

use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

#[test]
fn a_viewer_that_cannot_use_the_shared_file_exits() {
    let share =
        std::env::temp_dir().join(format!("rtc-metal-missing-{}.share", std::process::id()));
    let mut viewer = Command::new(env!("CARGO_BIN_EXE_rtc-metal"))
        .args(["scenes/chapter7.json", "--window", "--share"])
        .arg(&share)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
        .unwrap();
    let start = Instant::now();
    let status = loop {
        if let Some(status) = viewer.try_wait().unwrap() {
            break Some(status);
        }
        if start.elapsed() > Duration::from_secs(20) {
            let _ = viewer.kill();
            break None;
        }
        std::thread::sleep(Duration::from_millis(50));
    };
    let status = status.expect("the viewer kept running without its shared file");
    assert_eq!(status.code(), Some(1));
}
