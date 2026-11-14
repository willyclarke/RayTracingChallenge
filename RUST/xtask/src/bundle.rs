//! What goes in the app bundle, and where: pure, so it can be tested.

use std::path::{Path, PathBuf};

/// One app of the bundle: RTC Edit, or RTC Viewer inside it. Its
/// executable and icon are named after `name` (the command); Finder, the
/// Dock and the menu bar show `display`.
pub struct App {
    pub id: String,
    pub name: String,
    pub display: String,
    pub version: String,
    /// Whether Finder opens scenes with this app: the `.rtc` type and the
    /// document types.
    pub documents: bool,
}

/// The app's `Info.plist`.
pub fn info_plist(app: &App) -> String {
    let icon = format!("{}.icns", app.name);
    let strings = [
        ("CFBundleIdentifier", app.id.as_str()),
        ("CFBundleName", &app.display),
        ("CFBundleDisplayName", &app.display),
        ("CFBundleExecutable", &app.name),
        ("CFBundleIconFile", &icon),
        ("CFBundleShortVersionString", &app.version),
        ("CFBundleVersion", &app.version),
        ("CFBundlePackageType", "APPL"),
        ("CFBundleInfoDictionaryVersion", "6.0"),
        ("LSMinimumSystemVersion", "14.4"),
    ];
    let mut text = String::from(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
         <!DOCTYPE plist PUBLIC \"-//Apple//DTD PLIST 1.0//EN\" \"http://www.apple.com/DTDs/PropertyList-1.0.dtd\">\n\
         <plist version=\"1.0\">\n<dict>\n",
    );
    for (k, v) in strings {
        text += &format!("\t<key>{k}</key>\n\t<string>{}</string>\n", escape(v));
    }
    if app.documents {
        text += &documents();
    }
    text += "\t<key>NSHighResolutionCapable</key>\n\t<true/>\n</dict>\n</plist>\n";
    text
}

/// The `.rtc` scene type RTC Edit declares and owns, and the document
/// types it opens: its own scenes (Owner) and any JSON (Alternate, so
/// Open With lists it without taking over double-clicks).
fn documents() -> String {
    let s = |v: &str| format!("<string>{}</string>", escape(v));
    let doc = |name: &str, rank: &str, ty: &str| {
        format!(
            "\t\t<dict>\n\
             \t\t\t<key>CFBundleTypeName</key>{}\n\
             \t\t\t<key>CFBundleTypeRole</key>{}\n\
             \t\t\t<key>LSHandlerRank</key>{}\n\
             \t\t\t<key>LSItemContentTypes</key><array>{}</array>\n\
             \t\t</dict>\n",
            s(name),
            s("Editor"),
            s(rank),
            s(ty)
        )
    };
    format!(
        "\t<key>UTExportedTypeDeclarations</key>\n\
         \t<array>\n\
         \t\t<dict>\n\
         \t\t\t<key>UTTypeIdentifier</key>{}\n\
         \t\t\t<key>UTTypeDescription</key>{}\n\
         \t\t\t<key>UTTypeConformsTo</key><array>{}</array>\n\
         \t\t\t<key>UTTypeTagSpecification</key>\n\
         \t\t\t<dict>\n\
         \t\t\t\t<key>public.filename-extension</key><array>{}</array>\n\
         \t\t\t</dict>\n\
         \t\t</dict>\n\
         \t</array>\n\
         \t<key>CFBundleDocumentTypes</key>\n\
         \t<array>\n{}{}\t</array>\n",
        s(SCENE_TYPE),
        s("RTC Scene"),
        s("public.json"),
        s("rtc"),
        doc("RTC Scene", "Owner", SCENE_TYPE),
        doc("JSON", "Alternate", "public.json"),
    )
}

/// The uniform type identifier of RTC Edit's scene files.
const SCENE_TYPE: &str = "no.sb21.rtc.scene";

fn escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

/// Where each file of the bundle goes.
pub struct Paths {
    pub contents: PathBuf,
    pub exe: PathBuf,
    pub plist: PathBuf,
    pub icns: PathBuf,
    pub helper_contents: PathBuf,
    pub helper_exe: PathBuf,
    pub helper_plist: PathBuf,
    pub helper_icns: PathBuf,
}

/// The files of `app_dir` (`…/RTC Edit.app`), with the viewer's app in
/// its `Contents/Helpers`.
pub fn layout(app_dir: &Path) -> Paths {
    let contents = app_dir.join("Contents");
    let helper_contents = contents.join("Helpers/RTC Viewer.app/Contents");
    Paths {
        exe: contents.join("MacOS/rtc-edit"),
        plist: contents.join("Info.plist"),
        icns: contents.join("Resources/rtc-edit.icns"),
        helper_exe: helper_contents.join("MacOS/rtc-metal"),
        helper_plist: helper_contents.join("Info.plist"),
        helper_icns: helper_contents.join("Resources/rtc-metal.icns"),
        contents,
        helper_contents,
    }
}

/// The files of an `.iconset` that `iconutil` turns into an `.icns`, with
/// their size in pixels.
pub fn icon_sizes() -> Vec<(String, u32)> {
    [16, 32, 128, 256, 512]
        .into_iter()
        .flat_map(|px| {
            [
                (format!("icon_{px}x{px}.png"), px),
                (format!("icon_{px}x{px}@2x.png"), px * 2),
            ]
        })
        .collect()
}

/// The `version` of a `Cargo.toml`'s `[package]`.
pub fn crate_version(cargo_toml: &str) -> Option<String> {
    let mut in_package = false;
    for line in cargo_toml.lines().map(str::trim) {
        if line.starts_with('[') {
            in_package = line == "[package]";
        } else if in_package
            && let Some(rest) = line.strip_prefix("version")
            && let Some(v) = rest.trim_start().strip_prefix('=')
        {
            return Some(v.trim().trim_matches('"').to_string());
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    fn app(display: &str) -> App {
        App {
            id: "no.sb21.rtc-edit".into(),
            name: "rtc-edit".into(),
            display: display.into(),
            version: "0.1.0".into(),
            documents: true,
        }
    }

    #[test]
    fn the_editor_plist_declares_the_scene_type() {
        let text = info_plist(&app("RTC Edit"));
        for want in [
            "<key>UTExportedTypeDeclarations</key>",
            "<string>no.sb21.rtc.scene</string>",
            "<string>RTC Scene</string>",
            "<key>UTTypeConformsTo</key>",
            "<key>public.filename-extension</key>",
            "<string>rtc</string>",
            "<key>CFBundleDocumentTypes</key>",
            "<string>Owner</string>",
            "<string>Alternate</string>",
            "<string>public.json</string>",
            "<string>Editor</string>",
        ] {
            assert!(text.contains(want), "{want} missing:\n{text}");
        }
        assert!(text.ends_with("</dict>\n</plist>\n"));
    }

    #[test]
    fn the_viewer_plist_has_no_document_types() {
        let text = info_plist(&App {
            documents: false,
            ..app("RTC Viewer")
        });
        assert!(!text.contains("UTExportedTypeDeclarations"));
        assert!(!text.contains("CFBundleDocumentTypes"));
    }

    #[test]
    fn the_plist_has_every_key() {
        let text = info_plist(&app("RTC Edit"));
        assert!(text.starts_with(
            "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<!DOCTYPE plist PUBLIC \"-//Apple//DTD PLIST 1.0//EN\" \"http://www.apple.com/DTDs/PropertyList-1.0.dtd\">"
        ));
        for (k, v) in [
            ("CFBundleIdentifier", "no.sb21.rtc-edit"),
            ("CFBundleName", "RTC Edit"),
            ("CFBundleDisplayName", "RTC Edit"),
            ("CFBundleExecutable", "rtc-edit"),
            ("CFBundleIconFile", "rtc-edit.icns"),
            ("CFBundleShortVersionString", "0.1.0"),
            ("CFBundleVersion", "0.1.0"),
            ("CFBundlePackageType", "APPL"),
            ("CFBundleInfoDictionaryVersion", "6.0"),
            ("LSMinimumSystemVersion", "14.4"),
        ] {
            let want = format!("<key>{k}</key>\n\t<string>{v}</string>");
            assert!(text.contains(&want), "{want}\n{text}");
        }
        assert!(text.contains("<key>NSHighResolutionCapable</key>\n\t<true/>"));
        assert!(text.ends_with("</dict>\n</plist>\n"));
    }

    #[test]
    fn plist_values_are_escaped() {
        let text = info_plist(&app("a&b<c"));
        assert!(text.contains("<string>a&amp;b&lt;c</string>"), "{text}");
        assert!(!text.contains("a&b<c"));
    }

    #[test]
    fn the_helper_sits_inside_the_editor() {
        let p = layout(Path::new("/b/RTC Edit.app"));
        let c = "/b/RTC Edit.app/Contents";
        let h = "/b/RTC Edit.app/Contents/Helpers/RTC Viewer.app/Contents";
        assert_eq!(p.contents, Path::new(c));
        assert_eq!(p.exe, Path::new(&format!("{c}/MacOS/rtc-edit")));
        assert_eq!(p.plist, Path::new(&format!("{c}/Info.plist")));
        assert_eq!(p.icns, Path::new(&format!("{c}/Resources/rtc-edit.icns")));
        assert_eq!(p.helper_contents, Path::new(h));
        assert_eq!(p.helper_exe, Path::new(&format!("{h}/MacOS/rtc-metal")));
        assert_eq!(p.helper_plist, Path::new(&format!("{h}/Info.plist")));
        assert_eq!(
            p.helper_icns,
            Path::new(&format!("{h}/Resources/rtc-metal.icns"))
        );
    }

    #[test]
    fn icon_sizes_are_what_iconutil_wants() {
        let want: Vec<(String, u32)> = [
            ("icon_16x16.png", 16),
            ("icon_16x16@2x.png", 32),
            ("icon_32x32.png", 32),
            ("icon_32x32@2x.png", 64),
            ("icon_128x128.png", 128),
            ("icon_128x128@2x.png", 256),
            ("icon_256x256.png", 256),
            ("icon_256x256@2x.png", 512),
            ("icon_512x512.png", 512),
            ("icon_512x512@2x.png", 1024),
        ]
        .into_iter()
        .map(|(n, px)| (n.to_string(), px))
        .collect();
        assert_eq!(icon_sizes(), want);
    }

    #[test]
    fn crate_version_reads_the_package() {
        let toml = "[package]\nname = \"x\"\nversion = \"0.1.0\"\n\n[dependencies]\nfoo = { version = \"9\" }\n";
        assert_eq!(crate_version(toml).as_deref(), Some("0.1.0"));
        assert_eq!(crate_version("[dependencies]\nversion = \"9\"\n"), None);
        assert_eq!(crate_version(""), None);
    }
}
