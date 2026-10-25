//! Tells the compiler which commit rtc-edit is built from, for About
//! rtc-edit and Report an Issue…: `RTC_EDIT_COMMIT`, such as
//! `7f1756a 2026-10-25`, with ` dirty` when tracked files have changes.

use std::process::Command;

/// `git args…` in this crate's folder: its trimmed output, if it succeeded.
fn git(args: &[&str]) -> Option<String> {
    let out = Command::new("git")
        .args(args)
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .output()
        .ok()?;
    out.status
        .success()
        .then(|| String::from_utf8_lossy(&out.stdout).trim().to_string())
}

fn main() {
    let commit = match (
        git(&["rev-parse", "--short", "HEAD"]),
        git(&["log", "-1", "--format=%cs"]),
    ) {
        (Some(hash), Some(date)) => {
            // Untracked files, like the example scenes, don't count.
            let clean = git(&["diff", "--quiet", "HEAD"]).is_some();
            format!("{hash} {date}{}", if clean { "" } else { " dirty" })
        }
        _ => "unknown".to_string(),
    };
    println!("cargo:rustc-env=RTC_EDIT_COMMIT={commit}");
    // A commit moves HEAD's branch and rewrites the index; staging or
    // unstaging rewrites the index.
    for name in ["HEAD", "index"] {
        if let Some(path) = git(&["rev-parse", "--path-format=absolute", "--git-path", name]) {
            println!("cargo:rerun-if-changed={path}");
        }
    }
    // An edit not yet staged changes neither: watch what rtc-edit is built
    // from, so the build says dirty. Not `.`, which holds `target/`.
    for path in [
        "src",
        "assets",
        "Cargo.toml",
        "../rtc-rust/src",
        "../rtc-share/src",
        "../rtc-config/src",
    ] {
        println!("cargo:rerun-if-changed={path}");
    }
    println!("cargo:rerun-if-changed=build.rs");
}
