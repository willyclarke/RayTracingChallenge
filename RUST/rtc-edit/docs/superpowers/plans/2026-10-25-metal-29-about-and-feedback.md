# Milestone 29: About, Dock icon, Report an Issue Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** About rtc-edit shows a ray-traced icon, the exact build and credits; the Dock shows the icon; Help ▸ Report an Issue… opens a filled-in GitHub issue.

**Architecture:** A new pure `src/about.rs` holds the version (from a new `build.rs`), credits, the icon decoder and the issue URL. menu.rs gains `Command::ReportIssue` and fills `AboutMetadata`; main.rs opens the URL and sets the Dock icon through objc2-app-kit.

**Tech Stack:** Rust 2024, muda 0.21 (`AboutMetadata`, `Icon::from_rgba`), rtc-rust `Canvas::from_ppm_binary`, objc2 0.6.5 / objc2-foundation 0.3.2 / objc2-app-kit 0.3.2.

**Spec:** `RUST/rtc-edit/docs/superpowers/specs/2026-10-25-metal-29-about-and-feedback-design.md`

Paths are relative to `RUST/rtc-edit/`. Baseline: 149 lib tests pass (`cargo test --lib`).

## Global Constraints

- Issue URL prefix: `https://github.com/willyclarke/RayTracingChallenge/issues/new?title=`; title `rtc-edit: `.
- Body (before encoding), exactly: `**What happened**\n\n**What I expected**\n\n**Steps to reproduce**\n\n---\nrtc-edit {version}, macOS {macos}\n`, where `{version}` is `0.1.0 (<commit string>)`.
- Percent-encoding keeps `A–Z a–z 0–9 - _ . ~`, every other UTF-8 byte → `%XX` (upper-case hex).
- Commit string: `<short hash> <YYYY-MM-DD>` plus ` dirty` when `git diff --quiet HEAD` fails; `unknown` on any git failure.
- About: `short_version` = `CARGO_PKG_VERSION`, `version` = commit string, `copyright` = `© 2026 Willy Clarke`, credits = three lines: `Built while working through The Ray Tracer Challenge by Jamis Buck.` / `Made with Rust, iced, muda and Metal.` / `github.com/willyclarke/RayTracingChallenge`.
- `Report an Issue…`: Help menu after `rtc-edit Help`, id `report-issue`, no key, always enabled.
- Icon 256×256, alpha 255; a PPM that doesn't decode → no icon anywhere, nothing else changes.
- Commit messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`; do not push.

## Review Focus

1. A build in a checkout with staged-but-uncommitted changes: the version must say ` dirty` (rerun on `index`), and after the commit drop it.
2. `sw_vers` missing or failing: the URL still opens, with `macOS unknown`.
3. `open` itself missing / failing to spawn: status line `couldn't open the browser: <error>`, no panic.
4. The Dock icon set before NSApp finished launching may be reset by AppKit: set it on the editor's first `Opened`, not in `boot`.
5. Report an Issue… chosen while a pending `a`/`d` waits: like any menu command it cancels the wait (accepted behaviour, milestone 28).

---

### Task 1: build.rs version and about.rs text, encoding and URL

**Files:**
- Create: `build.rs`, `src/about.rs`
- Modify: `src/lib.rs` (`pub mod about;`, alphabetical first)

**Interfaces:**
- Produces:
  - `pub fn version() -> &'static str` — `env!("RTC_EDIT_COMMIT")`
  - `pub const CREDITS: &str`, `pub const COPYRIGHT: &str`
  - `pub fn percent_encode(s: &str) -> String`
  - `pub fn issue_url(version: &str, macos: &str) -> String` (`version` is the whole `0.1.0 (…)` string)

- [ ] **Step 1: Write the failing tests** in `src/about.rs` `mod tests`:

```rust
#[test]
fn percent_encode_keeps_unreserved_and_encodes_the_rest() {
    assert_eq!(percent_encode("a b*\n:?"), "a%20b%2A%0A%3A%3F");
    assert_eq!(percent_encode("AZaz09-_.~"), "AZaz09-_.~");
    assert_eq!(percent_encode("é"), "%C3%A9");
}

#[test]
fn the_issue_url_has_the_title_and_the_versions() {
    let url = issue_url("0.1.0 (abc 2026-10-25)", "26.0");
    let prefix = "https://github.com/willyclarke/RayTracingChallenge/issues/new?title=rtc-edit%3A%20&body=";
    assert!(url.starts_with(prefix), "{url}");
    let body = &url[prefix.len()..];
    assert!(body.starts_with(&percent_encode("**What happened**\n\n")), "{body}");
    assert!(body.ends_with(&percent_encode("---\nrtc-edit 0.1.0 (abc 2026-10-25), macOS 26.0\n")), "{body}");
}

#[test]
fn the_version_is_one_line() {
    assert!(!version().is_empty());
    assert!(!version().contains('\n'));
}
```

- [ ] **Step 2: Run** `cargo test --lib about::` — Expected: FAIL to compile (functions missing).

- [ ] **Step 3: Implement**
  - `build.rs`: run git with `std::process::Command` from `CARGO_MANIFEST_DIR`; `rev-parse --short HEAD` + `log -1 --format=%cs` trimmed and joined with a space, ` dirty` appended when `diff --quiet HEAD` exits non-zero; any spawn failure or non-zero exit of the first two → `unknown`. `println!("cargo:rustc-env=RTC_EDIT_COMMIT={s}")`. Rerun: for `HEAD` and `index`, `git rev-parse --git-path <name>` (relative to the manifest dir) → `cargo:rerun-if-changed=<that path>`; also `cargo:rerun-if-changed=build.rs`.
  - `about.rs` (module doc: what About rtc-edit and Report an Issue… show): the constants and functions above, values from Global Constraints.

- [ ] **Step 4: Run** `cargo test --lib` — Expected: 152 passed. `cargo clippy --all-targets` clean.

- [ ] **Step 5: Commit**

```bash
git add build.rs src/about.rs src/lib.rs
git commit -m "ENHANCEMENT: rtc-edit knows its commit, and builds a GitHub issue URL (milestone 29)"
```

### Task 2: the ray-traced icon

**Files:**
- Create: `assets/about.json`, `assets/about.ppm` (rendered)
- Modify: `src/about.rs`

**Interfaces:**
- Produces: `pub fn icon(ppm: &[u8]) -> Option<(Vec<u8>, u32, u32)>` (RGBA rows top to bottom, width, height); `pub const ICON_PPM: &[u8] = include_bytes!("../assets/about.ppm");`

- [ ] **Step 1: Write `assets/about.json`** (rtc-rust scene format, see `../rtc-rust/README.md`):

```json
{
  "comment": "About rtc-edit's icon: the book's first real scene, a reflective sphere on a checkered floor. Re-render: rtc-metal assets/about.json --samples 1024 -o assets/about.ppm",
  "camera": { "width": 256, "height": 256, "field_of_view": 0.8,
              "from": [0, 1.6, -4.2], "to": [0, 0.9, 0], "up": [0, 1, 0] },
  "light": { "type": "point", "position": [-8, 10, -8], "intensity": [1, 1, 1] },
  "shapes": [
    { "type": "plane", "comment": "floor",
      "material": { "specular": 0, "reflective": 0.1,
                    "pattern": { "type": "checkers", "a": [0.95, 0.95, 0.95], "b": [0.15, 0.15, 0.2] } } },
    { "type": "sphere", "comment": "ball",
      "transform": [ { "translate": [0, 1, 0] } ],
      "material": { "color": [0.85, 0.2, 0.1], "diffuse": 0.7, "specular": 0.9,
                    "shininess": 300, "reflective": 0.25 } }
  ]
}
```

- [ ] **Step 2: Render** from `../rtc-metal`: `cargo run --release -- ../rtc-edit/assets/about.json --samples 1024 -o ../rtc-edit/assets/about.ppm`. Expected: the file exists, `head -c 2` is `P6`, and it opens in Preview as a 256×256 red ball on checkers. Fix the scene by eye if the ball is cut off or the floor is empty; keep it square.

- [ ] **Step 3: Write the failing tests** in `src/about.rs` `mod tests`:

```rust
#[test]
fn the_icon_decodes_to_opaque_rgba() {
    let (rgba, w, h) = icon(ICON_PPM).expect("icon");
    assert_eq!((w, h), (256, 256));
    assert_eq!(rgba.len(), 256 * 256 * 4);
    assert!(rgba.chunks(4).all(|p| p[3] == 255));
}

#[test]
fn bad_bytes_give_no_icon() {
    assert_eq!(icon(b"not a ppm"), None);
}
```

- [ ] **Step 4: Run** `cargo test --lib about::` — Expected: FAIL to compile (`icon`, `ICON_PPM` missing).

- [ ] **Step 5: Implement** `icon`: `Canvas::from_ppm_binary(ppm).ok()?`, then per pixel row-major `pixel_at(x, y)` → `(c.clamp(0,1) * 255 + 0.5) as u8` for x/y/z, then 255.

- [ ] **Step 6: Run** `cargo test --lib` — Expected: 154 passed; clippy clean.

- [ ] **Step 7: Commit**

```bash
git add assets/about.json assets/about.ppm src/about.rs
git commit -m "ENHANCEMENT: rtc-edit's icon, rendered by the ray tracer (milestone 29)"
```

### Task 3: menu — Report an Issue… and the About panel

**Files:**
- Modify: `src/menu.rs`, `src/main.rs` (`Message`, `menu_message`, `handle`)

**Interfaces:**
- Consumes: `about::{version, CREDITS, COPYRIGHT, icon, ICON_PPM, issue_url}` (Tasks 1–2).
- Produces: `Command::ReportIssue` (in `ALL` after `Help`, 19 commands); `Message::ReportIssue`.

- [ ] **Step 1: Write the failing test** in `src/menu.rs` `mod tests`:

```rust
#[test]
fn report_an_issue_is_in_help_and_always_enabled() {
    assert_eq!(label(Command::ReportIssue), "Report an Issue…");
    assert_eq!(menu_of(Command::ReportIssue), "Help");
    assert_eq!(key(Command::ReportIssue), None);
    assert_eq!(id(Command::ReportIssue), "report-issue");
    assert!(enabled(Command::ReportIssue, &State::default()));
    assert_eq!(ALL.last(), Some(&Command::ReportIssue));
}
```

- [ ] **Step 2: Run** `cargo test --lib menu::` — Expected: FAIL to compile.

- [ ] **Step 3: Implement**
  - menu.rs: the variant in `ALL`, `id`, `label`, `key` (`None` arm), `menu_of` (Help), `enabled` (always arm); the Help submenu gets `item(Command::ReportIssue)` after Help.
  - `AboutMetadata { name, short_version: Some(CARGO_PKG_VERSION), version: Some(about::version()), copyright: Some(COPYRIGHT), credits: Some(CREDITS), icon: about::icon(ICON_PPM).and_then(|(p, w, h)| Icon::from_rgba(p, w, h).ok()), ..Default }` (the old `version: CARGO_PKG_VERSION` goes).
  - main.rs: `Message::ReportIssue` (doc `/// Help ▸ Report an Issue….`), mapped in `menu_message`. In `handle`: `macos` = `sw_vers -productVersion` stdout trimmed, or `"unknown"` on any failure; `let version = format!("{} ({})", env!("CARGO_PKG_VERSION"), about::version());` then `std::process::Command::new("open").arg(issue_url(&version, &macos)).spawn()`; on `Err(err)` set `e.status = format!("couldn't open the browser: {err}")`.

- [ ] **Step 4: Run** `cargo test --lib` — Expected: 155 passed (help.rs's "every menu command appears once" covers the new item); clippy clean.

- [ ] **Step 5: Commit**

```bash
git add src/menu.rs src/main.rs
git commit -m "ENHANCEMENT: rtc-edit's About panel shows the icon, build and credits; Help ▸ Report an Issue… (milestone 29)"
```

### Task 4: the Dock icon, README, eye check

**Files:**
- Modify: `Cargo.toml`, `src/main.rs` (the main window's `Opened` arm), `README.md`

**Interfaces:**
- Consumes: `about::{icon, ICON_PPM}`.

- [ ] **Step 1: Dependencies.** `Cargo.toml`, with a comment `# The Dock icon.`: `objc2 = "0.6.5"`, `objc2-foundation = "0.3.2"`, `objc2-app-kit = { version = "0.3.2", default-features = false, features = ["NSApplication", "NSResponder", "NSImage", "NSImageRep", "NSBitmapImageRep", "NSGraphics"] }`.

- [ ] **Step 2: Implement** `fn set_dock_icon()` in main.rs, called once beside `menus.install()` in the main window's first `Opened`:
  - `let Some((rgba, w, h)) = about::icon(about::ICON_PPM) else { return };`, `MainThreadMarker::new()` (return if `None`).
  - `NSBitmapImageRep::initWithBitmapDataPlanes_…_bytesPerRow_bitsPerPixel(alloc, null_mut(), w, h, 8, 4, true, false, NSDeviceRGBColorSpace, w * 4, 32)`; copy `rgba` into `bitmapData()` (`std::ptr::copy_nonoverlapping`, length `rgba.len()`); `NSImage::initWithSize(NSSize::new(w, h))`, `addRepresentation`, `NSApplication::sharedApplication(mtm).setApplicationIconImage(Some(&image))`.
  - Every `unsafe` block gets a `// SAFETY:` line (null planes: AppKit allocates; the copy length equals the rep's `w*h*4` bytes).

- [ ] **Step 3: Run** `cargo test --lib` and `cargo clippy --all-targets` — Expected: 155 passed, clean; `cargo build` succeeds.

- [ ] **Step 4: README.** In the menu table: rtc-edit row mentions About's icon, build and credits; a Help row `Report an Issue…` — "open a new GitHub issue with rtc-edit's and macOS's versions filled in". A short "The icon" paragraph: rendered by rtc-metal from `assets/about.json`, shown in About and the Dock, re-rendered by the Step 2 command of Task 2. The milestone list gets `- **29:** About with a ray-traced icon and the build; Report an Issue…`. Layout table: `build.rs`, `src/about.rs`, `assets/`.

- [ ] **Step 5: Check by eye** (`cargo run -- examples/someboxes.json`):
  - the Dock shows the ball while rtc-edit runs;
  - About rtc-edit: icon, `Version 0.1.0 (<hash> <date>)` (with ` dirty` if tracked files differ), copyright, the three credit lines;
  - Help ▸ Report an Issue… opens the browser on a new issue titled `rtc-edit: ` with the outline and the version line;
  - Apple's Send rtc-edit Feedback is still there and untouched.

- [ ] **Step 6: Commit**

```bash
git add Cargo.toml Cargo.lock src/main.rs README.md
git commit -m "ENHANCEMENT: rtc-edit's Dock icon is the ray-traced ball (milestone 29)"
```
