# Milestone 31: an app bundle

Branch `metal/31-app-bundle`, merged into develop with `--no-ff` when done.
Milestone 30 added the settings file; after 31 come Open Recent, a scene
file type, and the viewer's settings.

## Goal

rtc-edit starts only from a terminal, and the Dock shows its icon only
because it sets one at runtime. A `rtc-edit.app` built by one command
should start from Spotlight, the Dock or `~/Applications`, carry its own
viewer inside, and show two distinct icons: rtc-edit's red ball and the
viewer's blue one.

For this Mac only: signed ad hoc, no Apple account. Anyone else builds and
signs their own the same way.

Success:
- `cargo xtask bundle --install` builds both crates and installs
  `~/Applications/rtc-edit.app`; `codesign --verify --strict` passes.
- Started from Spotlight, rtc-edit shows its icon in the Dock; Start
  Viewer adds a second, blue icon named rtc-metal, with its own menu bar
  name.
- The bundled editor finds the bundled viewer wherever the app is put.
- `cargo run` and `rtc-edit` from a terminal work as today.

## Out of scope

- Opening a scene from Finder (double-click, drop on the Dock icon,
  `open -a rtc-edit x.json`): an Apple Event iced/winit don't pass on; the
  app opens with a new scene.
- A scene file type and extension (`.rtc` or `.rtcscene`, conforming to
  `public.json`): planned with the Apple Event work, after 31.
- Developer ID signing, notarisation, a `.dmg`, a universal binary.
- An update mechanism.

## Decisions

### xtask

- `RUST/xtask`: a binary crate; `RUST/.cargo/config.toml` has
  `[alias] xtask = "run --quiet --manifest-path xtask/Cargo.toml --"`.
  Cargo reads an alias's `--manifest-path` against the current folder, so
  `cargo xtask` runs from `RUST/` (from a crate's folder it fails with
  cargo's own "manifest path does not exist"). The paths xtask uses are
  found from its own `CARGO_MANIFEST_DIR`, not the current folder.
- `cargo xtask bundle` builds `RUST/build/rtc-edit.app` (git-ignored);
  `cargo xtask bundle --install` also replaces `~/Applications/rtc-edit.app`.
  Anything else prints a usage line and exits 2.
- Steps, in order:
  1. `cargo build --release` in `rtc-metal`, then `rtc-edit`.
  2. Each icon: `sips -z N N -s format png` from the 1024 px PPM into a
     temporary `.iconset`, then `iconutil -c icns`.
  3. Remove any old `build/rtc-edit.app`; write the two `Info.plist` files;
     copy the binaries and the `.icns` files in.
  4. `codesign --force --sign -` on the helper, then on the outer app (no
     `--deep`).
  5. With `--install`: `ditto` (which keeps the signature) to
     `~/Applications/rtc-edit.app.new`, move the old app aside to
     `.app.old`, rename the new one into place, remove the old one; a
     failure leaves the installed app as it was, and keeps `build/` (review
     fix: the first draft removed the installed app first, and a failed
     install also removed the build).
- A failed step prints `xtask: <step> failed: <its stderr>`, removes the
  partial bundle and exits 1.

### The layout

```
RTC Edit.app/Contents/
  Info.plist
  MacOS/rtc-edit
  Resources/rtc-edit.icns
  Helpers/RTC Viewer.app/Contents/
    Info.plist
    MacOS/rtc-metal
    Resources/rtc-metal.icns
```

### Names (eye-check change, 2026-10-26)

Mac apps show a capitalised name, so the bundles are `RTC Edit.app` and
`RTC Viewer.app`, with `CFBundleName` and `CFBundleDisplayName` set to
those; the executables, icon files and ids stay `rtc-edit` / `rtc-metal`.
`about::NAME = "RTC Edit"` names the app in rtc-edit's app menu (About,
Hide, Quit), About's name, "RTC Edit Help", the window titles and Report
an Issue's title and version line; the help page's command line and the
stderr prefixes keep the commands. rtc-metal's app menu says Hide RTC
Viewer. Everything below that says `rtc-edit.app` or `rtc-metal.app` as a
folder means these.

### Info.plist

| key | rtc-edit | rtc-metal |
|---|---|---|
| `CFBundleIdentifier` | `no.sb21.rtc-edit` | `no.sb21.rtc-metal` |
| `CFBundleName`, `CFBundleDisplayName` | `RTC Edit` | `RTC Viewer` |
| `CFBundleExecutable` | `rtc-edit` | `rtc-metal` |
| `CFBundleIconFile` | `rtc-edit.icns` | `rtc-metal.icns` |
| `CFBundleShortVersionString`, `CFBundleVersion` | rtc-edit's crate version | rtc-metal's crate version |
| `CFBundlePackageType` | `APPL` | `APPL` |
| `CFBundleInfoDictionaryVersion` | `6.0` | `6.0` |
| `NSHighResolutionCapable` | true | true |
| `LSMinimumSystemVersion` | `14.4` | `14.4` |

14.4 is what `os_sync_wait_on_address` in rtc-share needs. The versions are
read from each crate's `Cargo.toml` `version` line.

### The icons

- `rtc-edit/assets/icon.json`: About's scene (`about.json`) with a
  1024×1024 camera, rendered once to `rtc-edit/assets/icon.ppm` and
  committed. About keeps its embedded 256 px `about.ppm`.
- `rtc-metal/assets/icon.json`: the same layout with a blue glass sphere,
  rendered to `rtc-metal/assets/icon.ppm`.
- Each file's `comment` gives its re-render command:
  `rtc-metal assets/icon.json --samples 1024 -o assets/icon.ppm`.
- About 3 MB each in git.

### Finding the viewer (rtc-edit)

- `viewer::find(env: Option<&str>, exe: Option<&Path>, crate_dir: &Path)`
  tries, in order:
  1. `RTC_METAL`, when set and non-empty (as today, no file check).
  2. `<exe's folder>/../Helpers/rtc-metal.app/Contents/MacOS/rtc-metal`,
     when it is a file; `exe` is `std::env::current_exe()`.
  3. `<crate_dir>/../rtc-metal/target/release/rtc-metal`, when it is a
     file.
- Neither file: `rtc-metal not found at <bundled> or <development>; build
  it or set RTC_METAL`.
- The viewer is still started directly by its path, so the child, its pid
  and the exit watcher stay as they are; macOS shows it as `rtc-metal.app`
  because the binary sits in that bundle.

### Started from Finder, the Dock or Spotlight

- The current folder is `/`; milestone 30's panels already start in the
  remembered or `~/rtc` folders. An untitled scene's relative asset paths
  still resolve against it; the editor's Add buttons make none.
- No arguments: a new scene.
- stderr goes nowhere, so a settings write that fails at quit is lost;
  running errors still go to the status line.
- macOS runs one copy per bundle: another Dock click raises it. `rtc-edit`
  from a terminal still starts its own.
- `set_dock_icon` stays: the same picture in the bundle, the only icon
  from `cargo run`.

### Docs

- rtc-edit's README, "Building the app": the command, where the app
  lands, ad-hoc signing and how others build their own, quit the app
  before `--install`, what doesn't work yet (opening scenes from Finder).
- The milestone list: 31, then Open Recent, the scene file type with
  Finder opening, the viewer's settings.
- The layout table: `../xtask`, `assets/icon.*`. rtc-metal's README: its
  icon and re-render command.
- `.gitignore`: `/RUST/build/`, `/RUST/xtask/target`.

## Code

- `xtask/src/main.rs`: argument parsing and the steps (shells out).
- `xtask/src/bundle.rs` (pure):
  - `pub struct App { pub id, pub name, pub version: String }` (the icon
    and executable are named after `name`).
  - `pub fn info_plist(app: &App) -> String`.
  - `pub struct Paths { contents, exe, plist, icns, helper_contents,
    helper_exe, helper_plist, helper_icns: PathBuf }`;
    `pub fn layout(app_dir: &Path) -> Paths`.
  - `pub fn icon_sizes() -> Vec<(String, u32)>`: file name and pixels.
  - `pub fn crate_version(cargo_toml: &str) -> Option<String>`.
- `rtc-edit/src/viewer.rs`: `find` as above; `main.rs` passes
  `std::env::current_exe().ok().as_deref()`.

## Tests

- xtask:
  - `info_plist`: starts with the XML declaration and plist doctype; has
    each key and value in the table; `&` and `<` in a value are escaped.
  - `layout`: every path under `<app>/Contents` as in the layout.
  - `icon_sizes`: the ten names `icon_16x16.png`, `icon_16x16@2x.png`
    (32 px), … `icon_512x512@2x.png` (1024 px), with their pixels.
  - `crate_version`: the `version = "0.1.0"` line of a `[package]`; none
    when missing; a dependency's `version` is not taken.
- rtc-edit `find`:
  - `RTC_METAL` wins over a bundled viewer.
  - A bundled viewer (a file in a temporary `x.app/Contents/Helpers/...`)
    wins over the development path.
  - The development path when there is no bundle.
  - Neither: the error names both paths and `RTC_METAL`.
- By eye: install and start from Spotlight; the two Dock icons and names;
  Open… starts in `~/rtc` or the remembered folder; quit and relaunch from
  the Dock restores the frame; `codesign --verify --strict
  ~/Applications/rtc-edit.app`; Finder shows both icons (the helper's via
  Show Package Contents).

## Errors

- A build, `sips`, `iconutil`, `codesign` or copy fails: the message,
  the partial bundle removed, exit 1.
- `--install` while the installed app runs: it is replaced; the running
  copy keeps its files until it quits (the README says to quit first).
- No viewer found: the status line names both places.
