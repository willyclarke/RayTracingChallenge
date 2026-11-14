//! `cargo xtask bundle [--install]`, run from `RUST/`: builds
//! `build/RTC Edit.app` with `RTC Viewer.app` in its `Contents/Helpers`,
//! signed ad hoc for this Mac, and with `--install` copies it to
//! `~/Applications`.

mod bundle;

use bundle::{App, crate_version, icon_sizes, info_plist, layout};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode};

const USAGE: &str = "usage: cargo xtask bundle [--install]";

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let install = match args.iter().map(String::as_str).collect::<Vec<_>>()[..] {
        ["bundle"] => false,
        ["bundle", "--install"] => true,
        _ => {
            eprintln!("{USAGE}");
            return ExitCode::from(2);
        }
    };
    // `RUST/`, from this crate's own folder rather than the current one.
    let root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("xtask lives in RUST/")
        .to_path_buf();
    let app_dir = root.join("build/RTC Edit.app");
    if let Err(e) = build(&root, &app_dir) {
        // A half-built bundle must not be mistaken for a good one.
        let _ = fs::remove_dir_all(&app_dir);
        eprintln!("xtask: {e}");
        return ExitCode::from(1);
    }
    println!("built {}", app_dir.display());
    if install {
        // The built bundle stays: a failed install can be tried again.
        match self::install(&app_dir) {
            Ok(at) => println!("installed {}", at.display()),
            Err(e) => {
                eprintln!("xtask: {e}");
                return ExitCode::from(1);
            }
        }
    }
    ExitCode::SUCCESS
}

/// Run `cmd`; its stderr names the failure of `step`.
fn run(step: &str, cmd: &mut Command) -> Result<(), String> {
    let out = cmd.output().map_err(|e| format!("{step} failed: {e}"))?;
    if out.status.success() {
        Ok(())
    } else {
        let stderr = String::from_utf8_lossy(&out.stderr);
        Err(format!("{step} failed: {}", stderr.trim()))
    }
}

fn io(step: &str) -> impl Fn(std::io::Error) -> String + '_ {
    move |e| format!("{step} failed: {e}")
}

fn version(root: &Path, krate: &str) -> Result<String, String> {
    let toml = root.join(krate).join("Cargo.toml");
    let text = fs::read_to_string(&toml).map_err(io("reading Cargo.toml"))?;
    crate_version(&text).ok_or(format!("{} has no version", toml.display()))
}

fn build(root: &Path, app_dir: &Path) -> Result<(), String> {
    for krate in ["rtc-metal", "rtc-edit"] {
        run(
            &format!("cargo build --release in {krate}"),
            Command::new("cargo")
                .args(["build", "--release"])
                .current_dir(root.join(krate)),
        )?;
    }
    let p = layout(app_dir);
    if app_dir.exists() {
        fs::remove_dir_all(app_dir).map_err(io("removing the old bundle"))?;
    }
    let apps = [
        (
            App {
                id: "no.sb21.rtc-edit".into(),
                name: "rtc-edit".into(),
                display: "RTC Edit".into(),
                version: version(root, "rtc-edit")?,
                documents: true,
            },
            &p.plist,
            &p.exe,
            &p.icns,
        ),
        (
            App {
                id: "no.sb21.rtc-metal".into(),
                name: "rtc-metal".into(),
                display: "RTC Viewer".into(),
                version: version(root, "rtc-metal")?,
                documents: false,
            },
            &p.helper_plist,
            &p.helper_exe,
            &p.helper_icns,
        ),
    ];
    for (app, plist, exe, icns) in &apps {
        let krate = root.join(&app.name);
        for dir in [exe.parent(), icns.parent()].into_iter().flatten() {
            fs::create_dir_all(dir).map_err(io("making the bundle's folders"))?;
        }
        fs::write(plist, info_plist(app)).map_err(io("writing Info.plist"))?;
        fs::copy(krate.join("target/release").join(&app.name), exe)
            .map_err(io(&format!("copying {}", app.name)))?;
        icon(&krate.join("assets/icon.ppm"), icns)?;
    }
    // The helper first: signing the outer app seals what it holds.
    for contents in [&p.helper_contents, &p.contents] {
        let bundle = contents.parent().expect("Contents is in a bundle");
        run(
            "codesign",
            Command::new("codesign")
                .args(["--force", "--sign", "-"])
                .arg(bundle),
        )?;
    }
    Ok(())
}

/// `ppm` (1024 px) as an `.icns` at `icns`, through an `.iconset` of PNGs.
fn icon(ppm: &Path, icns: &Path) -> Result<(), String> {
    if !ppm.is_file() {
        return Err(format!("icon failed: no {}", ppm.display()));
    }
    let set = std::env::temp_dir()
        .join(format!("xtask-{}", std::process::id()))
        .join("icon.iconset");
    let _ = fs::remove_dir_all(&set);
    fs::create_dir_all(&set).map_err(io("making the iconset"))?;
    let made = icon_sizes()
        .into_iter()
        .try_for_each(|(name, px)| {
            run(
                "sips",
                Command::new("sips")
                    .args([
                        "-s",
                        "format",
                        "png",
                        "-z",
                        &px.to_string(),
                        &px.to_string(),
                    ])
                    .arg(ppm)
                    .arg("--out")
                    .arg(set.join(name)),
            )
        })
        .and_then(|()| {
            run(
                "iconutil",
                Command::new("iconutil")
                    .args(["-c", "icns", "-o"])
                    .arg(icns)
                    .arg(&set),
            )
        });
    let _ = fs::remove_dir_all(set.parent().expect("the iconset is in a folder"));
    made
}

/// Replace `~/Applications/RTC Edit.app` with `app_dir`; where it went.
fn install(app_dir: &Path) -> Result<PathBuf, String> {
    let home = std::env::var_os("HOME").ok_or("install failed: no $HOME")?;
    let to = PathBuf::from(home).join("Applications/RTC Edit.app");
    install_to(app_dir, &to)?;
    Ok(to)
}

/// Copy `app_dir` beside `to` first and swap it in only once the copy is
/// whole, so a failure leaves the installed app as it was.
fn install_to(app_dir: &Path, to: &Path) -> Result<(), String> {
    let staged = to.with_extension("app.new");
    let old = to.with_extension("app.old");
    for leftover in [&staged, &old] {
        let _ = fs::remove_dir_all(leftover);
    }
    if let Err(e) = run("ditto", Command::new("ditto").arg(app_dir).arg(&staged)) {
        let _ = fs::remove_dir_all(&staged);
        return Err(e);
    }
    if to.exists() {
        fs::rename(to, &old).map_err(io("moving the installed app aside"))?;
    }
    if let Err(e) = fs::rename(&staged, to) {
        let _ = fs::rename(&old, to);
        let _ = fs::remove_dir_all(&staged);
        return Err(format!("installing failed: {e}"));
    }
    let _ = fs::remove_dir_all(&old);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scratch(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("xtask-test-{}-{name}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn installing_replaces_the_old_app_wholly() {
        let dir = scratch("replace");
        let (new, to) = (dir.join("new.app"), dir.join("Applications/rtc-edit.app"));
        fs::create_dir_all(new.join("Contents")).unwrap();
        fs::write(new.join("Contents/new"), "").unwrap();
        fs::create_dir_all(to.join("Contents")).unwrap();
        fs::write(to.join("Contents/stale"), "").unwrap();
        install_to(&new, &to).unwrap();
        assert!(to.join("Contents/new").is_file());
        assert!(!to.join("Contents/stale").exists());
        let names: Vec<_> = fs::read_dir(to.parent().unwrap())
            .unwrap()
            .map(|e| e.unwrap().file_name())
            .collect();
        assert_eq!(names, ["rtc-edit.app"]);
    }

    #[test]
    fn a_failed_install_keeps_the_old_app() {
        let dir = scratch("keep");
        let to = dir.join("Applications/rtc-edit.app");
        fs::create_dir_all(to.join("Contents")).unwrap();
        fs::write(to.join("Contents/old"), "").unwrap();
        assert!(install_to(&dir.join("missing.app"), &to).is_err());
        assert!(to.join("Contents/old").is_file());
        let names: Vec<_> = fs::read_dir(to.parent().unwrap())
            .unwrap()
            .map(|e| e.unwrap().file_name())
            .collect();
        assert_eq!(names, ["rtc-edit.app"]);
    }
}
