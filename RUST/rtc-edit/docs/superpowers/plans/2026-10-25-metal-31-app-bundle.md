# App Bundle Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `cargo xtask bundle [--install]` builds `RUST/build/rtc-edit.app` with `rtc-metal.app` inside, each with its own icon, signed ad hoc; the bundled editor finds the bundled viewer.

**Architecture:** A new binary crate `RUST/xtask`: a pure `bundle.rs` (plist text, layout, icon sizes, crate version) and a `main.rs` that shells out to cargo, `sips`, `iconutil`, `codesign` and `ditto`. rtc-edit's `viewer::find` gains a "next to my executable" step. Two 1024 px icon scenes are rendered once by rtc-metal and committed.

**Tech Stack:** Rust 2024, std only in xtask; macOS `sips`, `iconutil`, `codesign`, `ditto`.

**Spec:** `rtc-edit/docs/superpowers/specs/2026-10-25-metal-31-app-bundle-design.md`

## Global Constraints

- Bundle ids `no.sb21.rtc-edit`, `no.sb21.rtc-metal`; `LSMinimumSystemVersion` `14.4`; `CFBundleInfoDictionaryVersion` `6.0`; `CFBundlePackageType` `APPL`; `NSHighResolutionCapable` true.
- Layout: `rtc-edit.app/Contents/{Info.plist, MacOS/rtc-edit, Resources/rtc-edit.icns, Helpers/rtc-metal.app/Contents/{Info.plist, MacOS/rtc-metal, Resources/rtc-metal.icns}}`.
- Output `RUST/build/rtc-edit.app`; `--install` → `~/Applications/rtc-edit.app` via `ditto`.
- Sign: `codesign --force --sign -`, helper first, outer second, no `--deep`.
- Errors: `xtask: <step> failed: <stderr>`, partial bundle removed, exit 1; bad arguments: usage, exit 2.
- xtask finds paths from its own `CARGO_MANIFEST_DIR`, never the current folder.
- `cargo run` / terminal `rtc-edit` behave as today.
- Commit prefixes ENHANCEMENT/FIX/DOC; messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. Do not push.

## Review Focus

1. A bundle copied elsewhere (e.g. `/Applications`, or Desktop) still finds its own viewer: `find` must use `current_exe`, not a build-time path (Task 1 test with a temp bundle).
2. Running `cargo xtask bundle` twice: the second run replaces the old bundle cleanly, with no stale files and a valid signature (by eye in Task 4: run twice, `codesign --verify --strict`).
3. A failed step (e.g. `iconutil` on a missing PPM) leaves no half-built `build/rtc-edit.app` that the user might install (Task 4: by eye with the PPM renamed).
4. A value with `&` or `<` in the plist breaks the XML (Task 3 test).
5. `--install` over an existing `~/Applications/rtc-edit.app` leaves no mix of old and new files (Task 4: remove then `ditto`).

---

### Task 1: rtc-edit finds a viewer beside its executable

**Files:**
- Modify: `RUST/rtc-edit/src/viewer.rs` (`find` and its tests), `RUST/rtc-edit/src/main.rs` (the `find` call)

**Interfaces:**
- Produces: `pub fn find(env: Option<&str>, exe: Option<&Path>, crate_dir: &Path) -> Result<PathBuf, String>`; `pub fn bundled(exe: &Path) -> PathBuf` = `<exe's folder>/../Helpers/rtc-metal.app/Contents/MacOS/rtc-metal`.

- [ ] **Step 1: Update the three existing tests to the new signature and add two**

Existing calls gain `None` for `exe`. New:

```rust
#[test]
fn find_prefers_a_viewer_in_its_own_bundle() {
    let root = fresh("bundled");
    let app = root.join("rtc-edit.app/Contents");
    let exe = app.join("MacOS/rtc-edit");
    let viewer = app.join("Helpers/rtc-metal.app/Contents/MacOS/rtc-metal");
    for f in [&exe, &viewer] { fs::create_dir_all(f.parent().unwrap()).unwrap(); fs::write(f, "").unwrap(); }
    // A development build beside it too: the bundle's own wins.
    let dev = root.join("rtc-metal/target/release/rtc-metal");
    fs::create_dir_all(dev.parent().unwrap()).unwrap(); fs::write(&dev, "").unwrap();
    let got = find(None, Some(&exe), &root.join("rtc-edit")).unwrap();
    assert_eq!(got.canonicalize().unwrap(), viewer.canonicalize().unwrap());
    assert_eq!(find(Some("/x/rtc-metal"), Some(&exe), &root.join("rtc-edit")).unwrap(), Path::new("/x/rtc-metal"));
}
```

and `find_reports_a_missing_viewer` becomes: with `exe = root.join("a.app/Contents/MacOS/rtc-edit")` (no files), the error contains `"RTC_METAL"`, `"Helpers/rtc-metal.app/Contents/MacOS/rtc-metal"` and `"rtc-metal/target/release/rtc-metal"`.

- [ ] **Step 2: Run them to verify they fail**

Run: `cd RUST/rtc-edit && cargo test --lib find_`
Expected: compile errors (`find` takes 2 arguments).

- [ ] **Step 3: Implement**

`find`: env (non-empty) → as is; `exe.map(bundled)` when `is_file()`; dev path when `is_file()`; else `Err(format!("rtc-metal not found at {} or {}; build it or set RTC_METAL", bundled.display(), dev.display()))` — with no `exe`, only the dev path is named. Doc comment lists the order. `main.rs`: `find(env.as_deref(), std::env::current_exe().ok().as_deref(), Path::new(env!("CARGO_MANIFEST_DIR")))`.

- [ ] **Step 4: Run the suite**

Run: `cd RUST/rtc-edit && cargo test --lib && cargo clippy --all-targets -- -D warnings`
Expected: all pass (159 + 1); no warnings.

- [ ] **Step 5: Commit**

`ENHANCEMENT: rtc-edit finds the viewer inside its own app bundle first (milestone 31)`

### Task 2: the two icons

**Files:**
- Create: `RUST/rtc-edit/assets/icon.json`, `RUST/rtc-edit/assets/icon.ppm`, `RUST/rtc-metal/assets/icon.json`, `RUST/rtc-metal/assets/icon.ppm`

- [ ] **Step 1: Write the scenes**

- `rtc-edit/assets/icon.json`: `about.json` with `"width": 1024, "height": 1024` and the comment `rtc-edit's app icon: About's scene at 1024 px. Re-render: rtc-metal assets/icon.json --samples 1024 -o assets/icon.ppm`.
- `rtc-metal/assets/icon.json`: the same camera, light and floor; the ball `"color": [0.1, 0.3, 0.9], "diffuse": 0.3, "specular": 0.9, "shininess": 300, "reflective": 0.6, "transparency": 0.8, "refractive_index": 1.5`; comment names it rtc-metal's icon with the same re-render line.

- [ ] **Step 2: Render**

Run, from `RUST/rtc-metal`: `cargo run --release -- ../rtc-edit/assets/icon.json --samples 1024 -o ../rtc-edit/assets/icon.ppm` and `cargo run --release -- assets/icon.json --samples 1024 -o assets/icon.ppm`
Expected: two P6 files of 1024×1024 (`head -c 15` shows `P6\n1024 1024`).

- [ ] **Step 3: Look at both** (Read the PNGs from `sips -s format png` into the scratch folder): a red ball and a blue glass ball on checkers. If the blue one reads as black or washed out, adjust its material and re-render; ledger the values.

- [ ] **Step 4: Commit**

`ENHANCEMENT: 1024 px app icons for rtc-edit (red ball) and rtc-metal (blue glass) (milestone 31)`

### Task 3: xtask's pure parts

**Files:**
- Create: `RUST/xtask/Cargo.toml` (package `xtask`, 0.1.0, edition 2024, no dependencies), `RUST/xtask/src/main.rs` (only `mod bundle;` and an empty `fn main` for now), `RUST/xtask/src/bundle.rs`

**Interfaces:**
- Produces:
  - `pub struct App { pub id: String, pub name: String, pub version: String }`
  - `pub fn info_plist(app: &App) -> String`
  - `pub struct Paths { pub contents, pub exe, pub plist, pub icns, pub helper_contents, pub helper_exe, pub helper_plist, pub helper_icns: PathBuf }`; `pub fn layout(app_dir: &Path) -> Paths` (names fixed: `rtc-edit`, `rtc-metal`)
  - `pub fn icon_sizes() -> Vec<(String, u32)>`
  - `pub fn crate_version(cargo_toml: &str) -> Option<String>`

- [ ] **Step 1: Write the failing tests in `bundle.rs`**

```rust
#[test] fn the_plist_has_every_key()
// app = App{id:"no.sb21.rtc-edit", name:"rtc-edit", version:"0.1.0"}; text = info_plist(&app)
// starts_with("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<!DOCTYPE plist PUBLIC \"-//Apple//DTD PLIST 1.0//EN\" \"http://www.apple.com/DTDs/PropertyList-1.0.dtd\">")
// contains each "<key>K</key>\n\t<string>V</string>" for: CFBundleIdentifier no.sb21.rtc-edit,
// CFBundleName rtc-edit, CFBundleExecutable rtc-edit, CFBundleIconFile rtc-edit.icns,
// CFBundleShortVersionString 0.1.0, CFBundleVersion 0.1.0, CFBundlePackageType APPL,
// CFBundleInfoDictionaryVersion 6.0, LSMinimumSystemVersion 14.4;
// contains "<key>NSHighResolutionCapable</key>\n\t<true/>"; ends_with("</dict>\n</plist>\n")
#[test] fn plist_values_are_escaped()
// name "a&b<c" → contains "<string>a&amp;b&lt;c</string>" and not "a&b<c"
#[test] fn the_helper_sits_inside_the_editor()
// layout("/b/rtc-edit.app"): contents "/b/rtc-edit.app/Contents", exe ".../Contents/MacOS/rtc-edit",
// plist ".../Contents/Info.plist", icns ".../Contents/Resources/rtc-edit.icns",
// helper_contents ".../Contents/Helpers/rtc-metal.app/Contents", helper_exe ".../MacOS/rtc-metal",
// helper_plist ".../Info.plist", helper_icns ".../Resources/rtc-metal.icns"
#[test] fn icon_sizes_are_what_iconutil_wants()
// == [("icon_16x16.png",16),("icon_16x16@2x.png",32),("icon_32x32.png",32),("icon_32x32@2x.png",64),
//     ("icon_128x128.png",128),("icon_128x128@2x.png",256),("icon_256x256.png",256),
//     ("icon_256x256@2x.png",512),("icon_512x512.png",512),("icon_512x512@2x.png",1024)]
#[test] fn crate_version_reads_the_package()
// "[package]\nname = \"x\"\nversion = \"0.1.0\"\n\n[dependencies]\nfoo = { version = \"9\" }\n" → Some("0.1.0")
// "[dependencies]\nversion = \"9\"\n" → None; "" → None
```

Plist keys in the order of the spec's table, one per line, tab-indented.

- [ ] **Step 2: Run them to verify they fail**

Run: `cd RUST/xtask && cargo test`
Expected: compile errors for the missing items.

- [ ] **Step 3: Implement** — `crate_version` reads lines from `[package]` until the next `[`, takes `version = "…"`. Escape `&`, `<`, `>`.

- [ ] **Step 4: Run**

Run: `cd RUST/xtask && cargo test && cargo clippy --all-targets -- -D warnings`
Expected: 5 passed (dead-code warnings for items `main` doesn't use yet are allowed only until Task 4; if clippy fails on them, add `#![allow(dead_code)]` to `bundle.rs` and remove it in Task 4 — ledger it).

- [ ] **Step 5: Commit**

`ENHANCEMENT: xtask's bundle layout, Info.plist and icon sizes (milestone 31)` (with `Cargo.lock`)

### Task 4: `cargo xtask bundle [--install]`

**Files:**
- Modify: `RUST/xtask/src/main.rs`
- Create: `RUST/.cargo/config.toml`
- Modify: `.gitignore` (`/RUST/build/`, `/RUST/xtask/target`)

**Interfaces:**
- Consumes: Task 3's `App`, `info_plist`, `layout`, `icon_sizes`, `crate_version`; Task 2's `icon.ppm` files.

- [ ] **Step 1: Implement `main`**

- `RUST = Path::new(env!("CARGO_MANIFEST_DIR")).parent()`.
- Args: `["bundle"]` or `["bundle", "--install"]`; else `eprintln!("usage: cargo xtask bundle [--install]")`, exit 2.
- `fn run(step: &str, cmd: &mut Command) -> Result<(), String>`: status and stderr; `Err(format!("{step} failed: {stderr}"))`.
- Steps per the spec; the iconset in `std::env::temp_dir()/xtask-<pid>/<name>.iconset`, removed at the end. Versions via `crate_version(read_to_string(<crate>/Cargo.toml))`, `.ok_or("… has no version")`.
- On any `Err(e)`: remove `build/rtc-edit.app`, `eprintln!("xtask: {e}")`, exit 1. On success print the bundle path (and the install path).
- `.cargo/config.toml`: `[alias]\nxtask = "run --quiet --manifest-path xtask/Cargo.toml --"`.

- [ ] **Step 2: Run it**

Run: `cd RUST && cargo xtask bundle && codesign --verify --strict --verbose=2 build/rtc-edit.app && ls -R build/rtc-edit.app`
Expected: `valid on disk`, `satisfies its Designated Requirement`; the tree matches the layout.

- [ ] **Step 3: Twice, and a failure**

Run: `cargo xtask bundle` again → valid again. Then `mv rtc-metal/assets/icon.ppm /tmp/x.ppm; cargo xtask bundle; echo $?; ls build; mv /tmp/x.ppm rtc-metal/assets/icon.ppm`
Expected: `xtask: … failed: …`, exit 1, no `rtc-edit.app` in `build/`.

- [ ] **Step 4: Bad arguments**

Run: `cargo xtask; echo $?` and `cargo xtask bundel; echo $?`
Expected: the usage line, 2.

- [ ] **Step 5: xtask tests and clippy**

Run: `cd RUST/xtask && cargo test && cargo clippy --all-targets -- -D warnings`
Expected: 5 passed; no warnings.

- [ ] **Step 6: Commit**

`ENHANCEMENT: cargo xtask bundle builds rtc-edit.app with rtc-metal.app inside (milestone 31)`

### Task 5: docs

**Files:**
- Modify: `RUST/rtc-edit/README.md`, `RUST/rtc-metal/README.md`

- [ ] **Step 1: README text** per the spec's Docs section: "Building the app" in rtc-edit's README (command, `RUST/build/rtc-edit.app`, `--install` to `~/Applications`, quit first, ad-hoc signing for this Mac and how others build their own, run from `RUST/`, opening scenes from Finder doesn't work yet); milestone list `**31:** an app bundle with the viewer inside.` and after it `Planned: Open Recent; a scene file type, opened from Finder; the viewer's samples and traversal.`; layout table `../xtask`, `assets/icon.*`. rtc-metal's README: its icon file and re-render command.

- [ ] **Step 2: Commit**

`DOC: building rtc-edit.app in the READMEs (milestone 31)`

The by-eye list from the spec's Tests section goes to the user after the final review, with `cargo xtask bundle --install` run first.
