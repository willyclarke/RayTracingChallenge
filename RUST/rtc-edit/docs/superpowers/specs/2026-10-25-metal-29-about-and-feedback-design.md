# Milestone 29: About rtc-edit, a ray-traced icon, and Report an Issue

Branch `metal/29-about-and-feedback`, merged into develop with `--no-ff` when
done. Milestone 28 added the Help menu; 30 is planned as a config file.

## Goal

About rtc-edit shows only a name and a version, and there is no way to
report a problem to the project. The About panel should show what the
project is (an icon the ray tracer rendered, the exact build, credits), the
Dock should show the same icon, and the Help menu should open a filled-in
GitHub issue.

Success:
- rtc-edit ▸ About rtc-edit shows the ray-traced icon, `Version 0.1.0
  (<commit> <date>)` and the credits.
- The Dock shows the same icon while rtc-edit runs.
- Help ▸ Report an Issue… opens a new GitHub issue in the browser, with a
  title prefix and a body outline holding the rtc-edit version and the macOS
  version.

## Background

The Help menu also holds "Send rtc-edit Feedback". macOS adds it to every
app for users in Apple's beta program; it opens Feedback Assistant and
reports to Apple, not to this project. There is no supported way to change
it, so it stays, and Report an Issue… sits beside it.

## Out of scope

- An email feedback item.
- An app bundle or `.icns` (the binary's Finder icon stays generic).
- An update check.
- Removing or redirecting Apple's Send rtc-edit Feedback.
- Rendering the icon in `build.rs` or at start-up.

## Decisions

### The icon

- `assets/about.json`: a reflective sphere on a checkered floor under one
  light, the camera 256×256. It is the source; the image is rendered from it
  once and committed:
  `rtc-metal assets/about.json --samples 1024 -o assets/about.ppm`
  (the README gives this line).
- `assets/about.ppm` (binary P6, 256×256) is embedded with
  `include_bytes!` and decoded with rtc-rust's `Canvas::from_ppm`.
- About panel: `AboutMetadata.icon` (muda `Icon::from_rgba`).
- Dock: on the editor window's first Opened (where the menus install), an
  `NSImage` built from the same RGBA (`NSBitmapImageRep`) is set as
  `NSApplication::applicationIconImage`. Needs `objc2`, `objc2-foundation`
  and `objc2-app-kit` at rtc-metal's versions (0.6.5 / 0.3.2 / 0.3.2).
- A PPM that doesn't decode leaves both without the icon; nothing else
  changes.

### The version

- `build.rs` (new) sets `RTC_EDIT_COMMIT` for the compiler: `git rev-parse
  --short HEAD`, a space, `git log -1 --format=%cs`, plus ` dirty` when
  `git diff --quiet HEAD` reports changes to tracked files (untracked files,
  like the example scenes, don't count). Any git failure gives `unknown`.
- It reruns when the repository's `HEAD` or `index` changes (paths from
  `git rev-parse --git-path`), so a commit or a staged change updates it,
  and when `src`, `assets`, `Cargo.toml` or the path crates' sources change,
  so an unstaged edit says dirty too (review fix).
- About panel: `version` = `CARGO_PKG_VERSION`, `short_version` = the
  commit string; AppKit shows `Version <version> (<short_version>)`, so
  `Version 0.1.0 (7f1756a 2026-10-25)` (review fix: the first draft had
  them swapped).

### Credits and copyright

- `copyright`: `© 2026 Willy Clarke`.
- `credits`:
  > Built while working through The Ray Tracer Challenge by Jamis Buck.
  > Made with Rust, iced, muda and Metal.
  > github.com/willyclarke/RayTracingChallenge

### Report an Issue…

- `Command::ReportIssue`, label `Report an Issue…`, in the Help menu after
  `rtc-edit Help`; no key; never greyed. The help page lists it, being built
  from `ALL`.
- Choosing it runs `open <url>`. A failure to start `open` puts `couldn't
  open the browser: <error>` in the status line.
- The URL:
  `https://github.com/willyclarke/RayTracingChallenge/issues/new?title=<t>&body=<b>`
  with `t` = `rtc-edit: ` and `b` =

  ```
  **What happened**

  **What I expected**

  **Steps to reproduce**

  ---
  rtc-edit 0.1.0 (7f1756a 2026-10-25), macOS 26.0
  ```

  both percent-encoded. The macOS version is `sw_vers -productVersion` at the
  time of the click, `unknown` if that fails.
- Percent-encoding keeps `A–Z a–z 0–9 - _ . ~` and writes every other UTF-8
  byte as `%XX`; no new crate.

## Code

- `build.rs` (new): the version, as above.
- `src/about.rs` (new, pure):
  - `pub fn version() -> &'static str`: `env!("RTC_EDIT_COMMIT")`.
  - `pub fn icon(ppm: &[u8]) -> Option<(Vec<u8>, u32, u32)>`: RGBA, width,
    height; alpha 255.
  - `pub const CREDITS: &str`, `pub const COPYRIGHT: &str`.
  - `pub fn percent_encode(s: &str) -> String`.
  - `pub fn issue_url(version: &str, macos: &str) -> String`.
- `src/menu.rs`: `Command::ReportIssue` (label, Help menu, no key, always
  enabled, id `report-issue`); `AboutMetadata` filled from about.rs.
- `src/main.rs`: `Message::ReportIssue` runs `open`; the Dock icon is set on
  the editor's first Opened.
- `Cargo.toml`: the objc2 crates.
- `README.md`: the About panel, Report an Issue…, the re-render command.

## Tests

- about.rs:
  - `icon(include_bytes!("../assets/about.ppm"))` is 256×256 with
    256·256·4 bytes, every fourth byte 255.
  - `icon(b"not a ppm")` is `None`.
  - `percent_encode`: `"a b*\n:?"` → `"a%20b%2A%0A%3A%3F"`; `"AZaz09-_.~"`
    unchanged; `"é"` → `"%C3%A9"`.
  - `issue_url("0.1.0 (abc 2026-10-25)", "26.0")` starts with
    `https://github.com/willyclarke/RayTracingChallenge/issues/new?title=rtc-edit%3A%20&body=`,
    and the body contains the encoded `rtc-edit 0.1.0 (abc 2026-10-25), macOS
    26.0`.
  - `version()` is not empty and has no newline.
- menu.rs: `label(ReportIssue)` is `Report an Issue…`, `menu_of` is Help,
  `key` is `None`, `enabled` is true; ids still map back (the existing
  test).
- help.rs: the existing "every menu command appears once" test covers the
  new item.
- By eye: the About panel's icon, version line and credits; the Dock icon;
  Report an Issue… opens a filled-in issue page; the icon looks right in
  light and dark.

## Errors

- Bad icon bytes: no icon in the panel or the Dock.
- No git at build time: version `unknown`.
- `open` fails to start: the status line says so.
