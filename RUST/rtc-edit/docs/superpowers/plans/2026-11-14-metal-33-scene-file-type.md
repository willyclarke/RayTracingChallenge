# Scene File Type Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `.rtc` scene files that RTC Edit owns, open from Finder (double-click, Open With, Dock drop), with `.json` still accepted.

**Architecture:** `xtask` writes the type declaration and document types into RTC Edit's `Info.plist`. `paths.rs` owns the extension rules (untitled name, Save As name, panel filters). A new `finder.rs` adds `application:openURLs:` to winit's delegate class and hands paths to iced through a queue-then-channel `Inbox`; `main.rs` turns them into `Message::OpenFiles` and opens the first like Open….

**Tech Stack:** Rust 2024, iced 0.14, objc2 0.6.5 (`class_addMethod`, `Imp`), objc2-app-kit 0.3.2, objc2-foundation 0.3.2 (`NSURL`, `NSArray`).

**Spec:** `rtc-edit/docs/superpowers/specs/2026-11-14-metal-33-scene-file-type-design.md`

## Global Constraints

- Type `no.sb21.rtc.scene`, description `RTC Scene`, conforms to `public.json`, extension `rtc`.
- Document types: `RTC Scene` / Editor / `Owner` / `no.sb21.rtc.scene`; `JSON` / Editor / `Alternate` / `public.json`. RTC Viewer's plist gets neither key.
- New scene: `untitled.rtc`. Save As suggests `<stem>.rtc`. Panels accept `["rtc", "json"]`. Plain Save keeps the file's own extension.
- Failure to install the Finder hook prints `rtc-edit: opening from Finder is off: <why>` to stderr and the editor carries on.
- Several files: the first opens; status `<name> (N more not opened)`.
- After a Finder open the editor window takes focus; the viewer is not raised.
- All commands run from `RUST/`. Commits: `ENHANCEMENT:`/`FIX:`/`DOC:` prefix, `(milestone 33)`, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. Do not push.
- rtc-metal is not changed. Repo example scenes stay `.json`.

## Deviations from the spec

- **No `file_paths(&[String])` percent decoder.** The added method asks each `NSURL` for `isFileURL` and `path` (AppKit decodes `%20` itself, as the spike did). A hand-written decoder would be code to get wrong. The spec's `file_paths` tests go; `Inbox` tests stay.
- **`more_not_opened` lives in `finder.rs`**, not `main.rs`: `main.rs` has no test module (the binary's glue is checked by eye), and the function is about Finder opens.
- **`open_path(&Path)`** is pulled out of `open()` and `open_recent()` (the same `from_path` → `switch` / `remember` / `note_recent` or error-status match) so `open_files` uses it too rather than a third copy.

## Review Focus

1. A Finder open arriving while an rfd panel or the "Save changes?" alert is up (modal run loop): the open should be queued in the channel and handled after, not lost or re-entrant.
2. A cold launch where `application:openURLs:` fires before the subscription starts: the file must still open (Inbox queue) — pinned by `inbox_hands_queued_paths_over_first`.
3. The subscription restarting (new receiver, old one dropped): later opens must reach the new one and not be dropped — pinned by `inbox_queues_again_when_the_receiver_is_gone`.
4. A `.rtc` that is not valid JSON or not a scene: status error, current scene untouched, editor still focused.
5. A file name with spaces or non-ASCII letters (`My Scene ø.rtc`): opens, and shows correctly in the status line and Open Recent (checked by eye, item 1 uses such a name).

---

### Task 1: RTC Edit's plist declares the scene type

**Files:**
- Modify: `xtask/src/bundle.rs` (`App`, `info_plist`, tests)
- Modify: `xtask/src/main.rs:87-104` (the two `App` literals)

**Interfaces:**
- Produces: `App { id, name, display, version, documents: bool }`.

- [ ] **Step 1: Write the failing tests** in `bundle.rs`'s test module. The test helper `app(display)` gains `documents: true`.

```rust
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
    let text = info_plist(&App { documents: false, ..app("RTC Viewer") });
    assert!(!text.contains("UTExportedTypeDeclarations"));
    assert!(!text.contains("CFBundleDocumentTypes"));
}
```

- [ ] **Step 2: Run them**

Run: `cargo test -p xtask`
Expected: compile error, `App` has no field `documents`.

- [ ] **Step 3: Implement.** Add `pub documents: bool` (doc: "Whether Finder opens scenes with this app: the `.rtc` type and the document types.") to `App`; `info_plist` appends the two keys exactly as in the spec's XML block, before `NSHighResolutionCapable`, when `app.documents`. Strings go through `escape`. Set `documents: true` for rtc-edit and `false` for rtc-metal in `xtask/src/main.rs`.

- [ ] **Step 4: Run the tests**

Run: `cargo test -p xtask`
Expected: all pass, including `the_plist_has_every_key` and `plist_values_are_escaped`.

- [ ] **Step 5: Check the real plist is valid**

Run: `cargo xtask bundle && plutil -lint "build/RTC Edit.app/Contents/Info.plist" "build/RTC Edit.app/Contents/Helpers/RTC Viewer.app/Contents/Info.plist"`
Expected: both `OK`.

- [ ] **Step 6: Commit**

```bash
git add xtask/src
git commit -m "ENHANCEMENT: RTC Edit's Info.plist declares the .rtc scene type and opens JSON as an alternate (milestone 33)"
```

### Task 2: `.rtc` names and panel filters

**Files:**
- Modify: `rtc-edit/src/paths.rs:169-174` and its tests
- Modify: `rtc-edit/src/main.rs:30, 450, 1027-1028`

**Interfaces:**
- Produces: `paths::SCENE_EXTENSIONS: [&str; 2]`, `paths::save_as_name(Option<&Path>) -> String`; `scene_name(None)` is now `"untitled.rtc"`.

- [ ] **Step 1: Write the failing tests** in `paths.rs` (change the existing `scene_name_is_the_file_name_or_untitled` assertion to `untitled.rtc`, add):

```rust
#[test]
fn save_as_suggests_the_stem_with_rtc() {
    assert_eq!(save_as_name(None), "untitled.rtc");
    assert_eq!(save_as_name(Some(Path::new("/a/pencil.json"))), "pencil.rtc");
    assert_eq!(save_as_name(Some(Path::new("/a/pencil.rtc"))), "pencil.rtc");
    assert_eq!(save_as_name(Some(Path::new("/a/notes"))), "notes.rtc");
    assert_eq!(save_as_name(Some(Path::new("/a/a.b.json"))), "a.b.rtc");
}

#[test]
fn panels_accept_rtc_first_then_json() {
    assert_eq!(SCENE_EXTENSIONS, ["rtc", "json"]);
}
```

- [ ] **Step 2: Run them**

Run: `cargo test -p rtc-edit --lib paths`
Expected: compile error (`save_as_name`, `SCENE_EXTENSIONS` missing).

- [ ] **Step 3: Implement** the constant and `save_as_name` (`Path::file_stem` + `.rtc`, `untitled.rtc` without a path or stem) as the spec's signatures; change `scene_name`'s default and doc comment to `untitled.rtc`.

- [ ] **Step 4: Run them**

Run: `cargo test -p rtc-edit --lib paths`
Expected: pass.

- [ ] **Step 5: Wire into `main.rs`.** Both `add_filter("scene", &["json"])` → `add_filter("scene", &SCENE_EXTENSIONS)`; Save As's `.set_file_name(self.name())` → `.set_file_name(save_as_name(self.path.as_deref()))`. Import both from `rtc_edit::paths`.

- [ ] **Step 6: Whole crate**

Run: `cargo test -p rtc-edit && cargo clippy -p rtc-edit --all-targets -- -D warnings`
Expected: all pass, no warnings.

- [ ] **Step 7: Commit**

```bash
git add rtc-edit/src/paths.rs rtc-edit/src/main.rs
git commit -m "ENHANCEMENT: new scenes are untitled.rtc, Save As suggests .rtc, panels take .rtc and .json (milestone 33)"
```

### Task 3: `finder.rs` — the hook and the inbox

**Files:**
- Create: `rtc-edit/src/finder.rs`
- Modify: `rtc-edit/src/lib.rs` (`pub mod finder;`)
- Modify: `rtc-edit/Cargo.toml` (objc2-foundation / objc2-app-kit features only if the compiler asks: `NSURL`, `NSArray`, `NSString`)

**Interfaces:**
- Produces:
  - `pub fn install(mtm: MainThreadMarker)`
  - `pub fn opens() -> UnboundedReceiver<Vec<PathBuf>>` (iced's `futures::channel::mpsc`)
  - `pub fn more_not_opened(n: usize) -> String`
  - `struct Inbox { queue: Vec<Vec<PathBuf>>, sender: Option<UnboundedSender<Vec<PathBuf>>> }` with `fn push(&mut self, paths: Vec<PathBuf>)` and `fn take_over(&mut self, sender: UnboundedSender<Vec<PathBuf>>)`, plus `static INBOX: Mutex<Inbox>`.

- [ ] **Step 1: Write the failing tests** in `finder.rs`:

```rust
fn p(s: &str) -> Vec<PathBuf> { vec![PathBuf::from(s)] }

#[test]
fn inbox_hands_queued_paths_over_first() {
    let mut inbox = Inbox::default();
    inbox.push(p("/a.rtc"));
    inbox.push(p("/b.rtc"));
    let (tx, mut rx) = unbounded();
    inbox.take_over(tx);
    inbox.push(p("/c.rtc"));
    let got: Vec<_> = std::iter::from_fn(|| rx.try_next().ok().flatten()).collect();
    assert_eq!(got, [p("/a.rtc"), p("/b.rtc"), p("/c.rtc")]);
    assert!(inbox.queue.is_empty());
}

#[test]
fn inbox_queues_again_when_the_receiver_is_gone() {
    let mut inbox = Inbox::default();
    let (tx, rx) = unbounded();
    inbox.take_over(tx);
    drop(rx);
    inbox.push(p("/a.rtc"));
    let (tx, mut rx) = unbounded();
    inbox.take_over(tx);
    assert_eq!(rx.try_next().ok().flatten(), Some(p("/a.rtc")));
}

#[test]
fn more_not_opened_counts_the_rest() {
    assert_eq!(more_not_opened(0), "");
    assert_eq!(more_not_opened(1), "");
    assert_eq!(more_not_opened(2), " (1 more not opened)");
    assert_eq!(more_not_opened(3), " (2 more not opened)");
}
```

(`more_not_opened` returns the suffix with its leading space, so the caller writes `format!("{name}{}", more_not_opened(n))`.)

- [ ] **Step 2: Run them**

Run: `cargo test -p rtc-edit --lib finder`
Expected: compile error (module items missing).

- [ ] **Step 3: Implement `Inbox`, `INBOX`, `opens()` and `more_not_opened`.** `push`: try `sender.unbounded_send`; on error (or no sender) drop the sender and queue. `take_over`: send every queued entry in order to the new sender (re-queue on failure), then keep it. `opens()`: make a channel, `INBOX.lock()…take_over(tx)`, return `rx`. A poisoned lock is taken with `into_inner` (`unwrap_or_else(PoisonError::into_inner)`).

- [ ] **Step 4: Run them**

Run: `cargo test -p rtc-edit --lib finder`
Expected: 3 pass.

- [ ] **Step 5: Implement `install(mtm)`** from the spike (memory `winit-finder-open-spike`):
  - `NSApplication::sharedApplication(mtm).delegate()`; none → warn `rtc-edit: opening from Finder is off: no application delegate`.
  - `extern "C-unwind" fn open_urls(_this: &AnyObject, _cmd: Sel, _app: &NSApplication, urls: &NSArray<NSURL>)`: collect `url.path()` for each URL with `url.isFileURL()`, as `PathBuf`; when not empty, `INBOX.lock().push(paths)`. Nothing else.
  - `transmute` to `objc2::runtime::Imp`, `objc2::ffi::class_addMethod(class, sel!(application:openURLs:), imp, c"v@:@@".as_ptr())`; `false` → warn `rtc-edit: opening from Finder is off: application:openURLs: is already defined`.
  - SAFETY comments on the transmute (signature matches `v@:@@`: void, self, _cmd, two objects), the `class_addMethod` call (class pointer from a live object, main thread, before AppKit sends opens), and the module doc saying why not `setDelegate:` or an Apple Event handler (spike results, one line each).

- [ ] **Step 6: Whole crate**

Run: `cargo test -p rtc-edit && cargo clippy -p rtc-edit --all-targets -- -D warnings`
Expected: all pass, no warnings.

- [ ] **Step 7: Commit**

```bash
git add rtc-edit/src/finder.rs rtc-edit/src/lib.rs rtc-edit/Cargo.toml
git commit -m "ENHANCEMENT: finder.rs receives application:openURLs: on winit's delegate and queues the paths (milestone 33)"
```

### Task 4: the editor opens Finder files; README

**Files:**
- Modify: `rtc-edit/src/main.rs` (`Message`, `boot`, `open`, `open_recent`, new `open_path` / `open_files`, `handle`, `subscription`)
- Modify: `rtc-edit/README.md`

**Interfaces:**
- Consumes: `finder::{install, opens, more_not_opened}` (Task 3); `SCENE_EXTENSIONS`, `save_as_name` (Task 2).
- Produces: `Message::OpenFiles(Vec<PathBuf>)`; `Editor::open_path(&mut self, path: &Path) -> bool`.

The glue needs a window and AppKit; there is no test module in `main.rs`. It is verified by the build, the existing suite, and the eye check.

- [ ] **Step 1: `open_path`.** Pull the shared `match SceneDescription::from_path(path)` out of `open()` and `open_recent()` into `fn open_path(&mut self, path: &Path) -> bool` (ok → `switch`, `remember(SCENE_FOLDER, folder_value(path))`, `note_recent`, `true`; error → status `<name>: <error>`, `false`). `open()` and `open_recent()` call it. Run `cargo test -p rtc-edit` — expected: all pass (no behaviour change).

- [ ] **Step 2: `open_files`.** `Message::OpenFiles(Vec<PathBuf>)`; `fn open_files(&mut self, paths: Vec<PathBuf>)`: empty → return; `ask_unsaved()` false → return; `open_path(&paths[0])` true and `paths.len() > 1` → `self.status = format!("{}{}", file_name(&paths[0]), finder::more_not_opened(paths.len()))`. In `handle`, the arm calls it and returns `window::gain_focus(e.main)` (whether or not it opened: the user just asked for the editor).

- [ ] **Step 3: Hook up.** In `boot()`, before the launch scene is loaded: `if let Some(mtm) = MainThreadMarker::new() { finder::install(mtm); }`. In `subscription()`, next to `menu_events`: `Subscription::run(finder_opens)` with `fn finder_opens() -> impl Stream<Item = Message>` = `finder::opens().map(Message::OpenFiles)`. Add `Message::OpenFiles(_)` wherever `Message::MenuEvent(_)` is listed among messages that don't mark the scene (line ~1226) if that list is exhaustive for such messages.

- [ ] **Step 4: Whole crate**

Run: `cargo test -p rtc-edit && cargo clippy -p rtc-edit --all-targets -- -D warnings`
Expected: all pass, no warnings.

- [ ] **Step 5: README.** Mark milestone 33 done; drop "a scene file type, opened from Finder;" from the Planned line (keep the viewer's samples and traversal and Import Scene…); under the File commands or a short "Opening scenes" paragraph: `.rtc` files open from Finder by double-click or a drop on the Dock icon, `.json` scenes through Open With ▸ RTC Edit; new scenes are saved as `.rtc`, `.json` still opens and saves; a `src/finder.rs` row in the layout table.

- [ ] **Step 6: Commit**

```bash
git add rtc-edit/src/main.rs rtc-edit/README.md
git commit -m "ENHANCEMENT: scenes opened from Finder open like Open…, the first of several, with the editor in front (milestone 33)"
```

- [ ] **Step 7: Install for the eye check**

Run: `cargo xtask bundle --install && /System/Library/Frameworks/CoreServices.framework/Frameworks/LaunchServices.framework/Support/lsregister -dump | grep -c no.sb21.rtc.scene`
Expected: the install succeeds; the count is at least 1.

Then hand the user the spec's 7 checks by eye, using `~/rtc/scenes/My Scene ø.rtc` (a copy of an example scene) for check 1.
