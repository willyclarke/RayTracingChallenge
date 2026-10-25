//! rtc-edit: `rtc-edit [scene.json]` shows the scene as a tree, starts
//! `rtc-metal --window` on a preview copy of it, and saves it
//! (Cmd-S, and Cmd-Shift-S to save elsewhere with asset paths rebased).
//! Cmd-R starts the viewer again after it was closed. Simple shapes,
//! point lights and the camera are edited in fixed fields; Cmd-L loads the
//! edits into the viewer. Cmd-N and Cmd-O start a new scene or open another
//! in the same viewer, which keeps the focus off itself, and shapes and lights are added and deleted from a button bar and keys.
//! The menu bar lists every command with its Cmd key, and owns those keys.

use iced::futures::channel::mpsc::{UnboundedReceiver, unbounded};
use iced::futures::{Stream, StreamExt};
use iced::widget::operation;
use iced::widget::{
    Column, Row, button, checkbox, column, container, row, scrollable, space, text, text_input,
};
use iced::{Element, Font, Length, Point, Rectangle, Size, Subscription, Task, keyboard, window};
use muda::MenuEvent;
use rtc_edit::edit::{self, Item, ancestors, new_scene};
use rtc_edit::fields::{
    Field, closed, fields, is_cut, lens_radius, next_box, set_closed, set_remembering,
};
use rtc_edit::follow::Follower;
use rtc_edit::history::{Group, History};
use rtc_edit::menu::{self, Command, Menus, State};
use rtc_edit::paths::{dir_of, rebase, scene_name};
use rtc_edit::preview::Preview;
use rtc_edit::tree::{Node, json, next_row, rows};
use rtc_edit::viewer::{Viewer, camera_start, exit_source, find, held_fresh, share_path};
use rtc_rust::scene::{FocalBlur, SceneDescription};
use rtc_share::{Shared, Snapshot};
use std::collections::{HashMap, HashSet};
use std::ffi::OsStr;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

/// The preview and the viewer, which must not outlive the editor. They live
/// here rather than in `Editor` because Quit in the app menu (Cmd-Q) ends
/// the process with `exit()`, which skips `Drop`; the `atexit` hook
/// `clean_up` empties this on every normal exit instead.
static CHILDREN: Mutex<Children> = Mutex::new(Children {
    preview: None,
    viewer: None,
});

struct Children {
    preview: Option<Preview>,
    viewer: Option<Viewer>,
}

unsafe extern "C" {
    fn atexit(f: extern "C" fn()) -> std::ffi::c_int;
}

/// Close the viewer and remove the preview folder and the shared file.
extern "C" fn clean_up() {
    if let Ok(mut children) = CHILDREN.lock() {
        children.viewer = None;
        children.preview = None;
    }
    let _ = std::fs::remove_file(own_share_path());
}

/// This editor's shared file, made by the first Start Viewer.
fn own_share_path() -> PathBuf {
    share_path(&std::env::temp_dir(), std::process::id())
}

struct Editor {
    /// The scene's file; none for a new scene not yet saved.
    path: Option<PathBuf>,
    desc: SceneDescription,
    folded: HashSet<Node>,
    selected: Option<Node>,
    /// Why there is no preview file, when there is none.
    preview_error: Option<String>,
    status: String,
    /// The text typed in each box of the selected item.
    texts: HashMap<Field, String>,
    /// Boxes whose text doesn't parse, with why; the scene keeps their last
    /// good value.
    errors: HashMap<Field, String>,
    /// Edits not yet saved, shown as ` •` in the title.
    unsaved: bool,
    /// Edits the viewer hasn't loaded; Cmd-L loads them.
    unloaded: bool,
    /// Cmd is held, so a box's one-character insert is a shortcut's letter.
    command: bool,
    /// The shaded box: the arrows and h/j/k/l move it, Enter, i and c edit
    /// it, x empties it. With neither it nor `editing`, keys go to the tree.
    left: Option<Field>,
    /// The box being edited, which takes ←/→ and letters itself.
    editing: Option<Field>,
    /// The focal blur that emptying the selected camera's f-number removed,
    /// so retyping one brings its focal distance and samples back.
    removed_blur: Option<FocalBlur>,
    /// A key waiting for its second: `a` for the item to add, `d` for a
    /// second `d` that deletes.
    pending: Option<char>,
    /// The window's outer top-left and inner size, in points from the main
    /// screen's top-left; each viewer starts at its top-right corner.
    frame: Rectangle,
    /// The file shared with every viewer: their camera, the scenes sent to
    /// them; made by the first Start Viewer.
    share: Option<Shared>,
    /// Whether a viewer is running: set when one starts, cleared when its
    /// exit is reported.
    live: bool,
    /// The running viewer's pid, for its exit source.
    viewer_pid: u32,
    /// The camera write the running viewer's camera stream waits from.
    camera_start: u32,
    /// The viewer's camera as last read from `share`.
    seen: Option<Snapshot>,
    /// A load held until the viewer says it takes loads: whether fresh.
    held: Option<bool>,
    /// Whether the camera's `from` follows the viewer's.
    follower: Follower,
    /// The steps Cmd-Z undoes and Cmd-Shift-Z redoes.
    history: History,
    /// The menu bar, once the window has opened.
    menus: Option<Menus>,
}

#[derive(Debug, Clone)]
enum Message {
    Select(Node),
    Toggle(Node),
    Save,
    SaveAs,
    New,
    Open,
    ShowViewer,
    Edit(Field, String),
    /// Pasted text, which arrives while Cmd is held.
    Paste(Field, String),
    Closed(bool),
    Load,
    /// Enter in a box: leave it.
    Submit(Field),
    /// Tab and Shift-Tab.
    FocusNext,
    FocusPrevious,
    /// An arrow no box took: right and down, each -1, 0 or 1.
    Move(isize, isize),
    /// Enter no box took.
    Enter,
    /// A plain character no box took: h j k l i c x u a d s p y o.
    Key(char),
    /// Add an item next to the selected one.
    Add(Item),
    /// Delete the selected item.
    Delete,
    /// Esc; `true` when a box took it and lost focus.
    Escape(bool),
    /// A click; `true` when a widget took it, such as a box it focused.
    Clicked(bool),
    /// Whether a box has focus, asked after a click.
    Focused(Field, bool),
    /// Whether Cmd is held.
    Command(bool),
    /// The window opened, moved or was resized.
    Window(window::Event),
    /// The viewer with this pid exited.
    ViewerExited(u32),
    /// The viewer wrote its camera.
    Camera(Snapshot),
    /// Waiting for the viewer's camera failed, and why.
    CameraLost(String),
    /// The follow checkbox.
    Follow(bool),
    /// Keep the followed camera (the button or Cmd-K).
    Keep,
    /// Cmd-Z and Cmd-Shift-Z.
    Undo,
    Redo,
    /// A menu bar click or key, from muda's event handler.
    MenuEvent(MenuEvent),
}

/// The editor window's size at start-up.
const WIDTH: f32 = 900.0;
const HEIGHT: f32 = 640.0;

fn main() -> iced::Result {
    iced::application(boot, update, view)
        .title(|e: &Editor| {
            let mark = if e.unsaved { " •" } else { "" };
            format!("rtc-edit — {}{mark}", e.name())
        })
        .subscription(subscription)
        .window_size((WIDTH, HEIGHT))
        .position(window::Position::Specific(Point::ORIGIN))
        .run()
}

fn file_name(path: &Path) -> String {
    path.file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default()
}

fn boot() -> Editor {
    let loaded = std::env::args().nth(1).map(|arg| {
        let path = PathBuf::from(arg);
        let desc = SceneDescription::from_path(&path).unwrap_or_else(|e| {
            eprintln!("rtc-edit: {e}");
            std::process::exit(1);
        });
        (desc, path)
    });
    // SAFETY: registers a plain function with no arguments; it runs once, at exit.
    unsafe { atexit(clean_up) };
    let mut editor = Editor {
        path: None,
        desc: new_scene(),
        folded: HashSet::new(),
        selected: None,
        preview_error: None,
        status: String::new(),
        texts: HashMap::new(),
        errors: HashMap::new(),
        unsaved: false,
        unloaded: false,
        command: false,
        left: None,
        editing: None,
        removed_blur: None,
        pending: None,
        frame: Rectangle::new(Point::ORIGIN, Size::new(WIDTH, HEIGHT)),
        share: None,
        live: false,
        seen: None,
        held: None,
        viewer_pid: 0,
        camera_start: 0,
        follower: Follower::default(),
        // `switch` below starts the history on the scene it opens.
        history: History::new(&new_scene()),
        menus: None,
    };
    match loaded {
        Some((desc, path)) => editor.switch(desc, Some(path)),
        None => editor.switch(new_scene(), None),
    }
    editor.show_viewer();
    editor
}

impl Editor {
    /// What the menu bar's greying and tick need.
    fn menu_state(&self) -> State {
        State {
            can_undo: self.history.can_undo(),
            can_redo: self.history.can_redo(),
            can_delete: self.left.is_none()
                && self.editing.is_none()
                && self
                    .selected
                    .as_ref()
                    .is_some_and(|n| edit::can_delete(&self.desc, n).is_ok()),
            viewer_running: self.live,
            following: self.follower.following(),
        }
    }

    /// The scene's file name, or `untitled.json`.
    fn name(&self) -> String {
        scene_name(self.path.as_deref())
    }

    /// The scene's file as an absolute path, `..` resolved where it exists;
    /// `None` for an untitled scene.
    fn full_path(&self) -> Option<PathBuf> {
        let path = self.path.as_deref()?;
        Some(path.canonicalize().unwrap_or_else(|_| path.to_path_buf()))
    }

    /// The status after a save: where the file went, in full.
    fn saved(&self) -> String {
        let path = self.full_path().unwrap_or_default();
        self.warn(format!("saved {}", path.display()))
    }

    /// The scene's folder, or the working directory for an untitled scene.
    fn dir(&self) -> PathBuf {
        dir_for(self.path.as_deref())
    }

    /// Replace the scene: the old preview goes, a preview of `desc` is
    /// made, and the running viewer shows it from its own camera.
    fn switch(&mut self, desc: SceneDescription, path: Option<PathBuf>) {
        let sent = {
            let mut children = CHILDREN.lock().unwrap();
            self.seen = None;
            self.follower.stop();
            // Before the new preview: both live in the same folder, which
            // the old one's drop removes.
            children.preview = None;
            let name = scene_name(path.as_deref());
            let preview = Preview::new(&std::env::temp_dir(), OsStr::new(&name))
                .and_then(|p| p.write(&desc, &dir_for(path.as_deref())).map(|()| p));
            self.preview_error = None;
            match preview {
                Ok(p) => children.preview = Some(p),
                Err(e) => self.preview_error = Some(e),
            }
            let running = children.viewer.as_mut().is_some_and(Viewer::running);
            match (&children.preview, &self.share) {
                (Some(p), Some(share)) if running => {
                    Some(send_or_hold(share, &mut self.held, p.path(), true))
                }
                _ => None,
            }
        };
        if let Some(Err(e)) = sent {
            self.status = e;
        }
        self.desc = desc;
        self.path = path;
        self.folded.clear();
        self.selected = None;
        self.texts.clear();
        self.errors.clear();
        self.left = None;
        self.editing = None;
        self.removed_blur = None;
        self.pending = None;
        self.unsaved = false;
        self.unloaded = false;
        self.history.clear(&self.desc);
    }

    /// Before New or Open: true when the scene has no unsaved edits, or they
    /// were saved or given up; false when the user cancelled or saving failed.
    fn ask_unsaved(&mut self) -> bool {
        if !self.unsaved {
            return true;
        }
        let mut question = format!("Save changes to {}?", self.name());
        if let Some(path) = self.full_path() {
            question += &format!("\n\nin {}", dir_of(&path).display());
        }
        let answer = rfd::MessageDialog::new()
            .set_level(rfd::MessageLevel::Warning)
            .set_description(question)
            .set_buttons(rfd::MessageButtons::YesNoCancelCustom(
                "Save".into(),
                "Don't Save".into(),
                "Cancel".into(),
            ))
            .show();
        // The key-up of Cmd went to the alert.
        self.command = false;
        match answer {
            rfd::MessageDialogResult::Custom(s) if s == "Save" => {}
            rfd::MessageDialogResult::Yes => {}
            rfd::MessageDialogResult::Custom(s) if s == "Don't Save" => return true,
            rfd::MessageDialogResult::No => return true,
            _ => return false,
        }
        self.save();
        !self.unsaved
    }

    /// Cmd-N.
    fn new_file(&mut self) {
        if self.ask_unsaved() {
            self.switch(new_scene(), None);
        }
    }

    /// Cmd-O.
    fn open(&mut self) {
        if !self.ask_unsaved() {
            return;
        }
        let picked = rfd::FileDialog::new()
            .set_directory(self.dir())
            .add_filter("scene", &["json"])
            .pick_file();
        self.command = false;
        let Some(path) = picked else { return };
        match SceneDescription::from_path(&path) {
            Ok(desc) => self.switch(desc, Some(path)),
            Err(e) => self.status = format!("{}: {e}", file_name(&path)),
        }
    }

    /// Start the viewer on the scene as it is now, unless it is still
    /// open; the shared file is made the first time.
    fn show_viewer(&mut self) {
        let mut children = CHILDREN.lock().unwrap();
        if children.viewer.as_mut().is_some_and(Viewer::running) {
            self.status = "viewer already running".into();
            return;
        }
        self.viewer_ended();
        if self.share.is_none() {
            match Shared::create(&own_share_path()) {
                Ok(s) => self.share = Some(s),
                Err(e) => {
                    self.status = format!("no shared file: {e}");
                    return;
                }
            }
        }
        let Some(preview) = &children.preview else {
            let e = self.preview_error.as_deref().unwrap_or("none");
            self.status = format!("no preview file: {e}");
            return;
        };
        // The new viewer shows the scene as it is now, including a view
        // followed from the last one and never loaded.
        if let Err(e) = preview.write(&self.desc, &self.dir()) {
            self.status = e;
            return;
        }
        self.unloaded = false;
        // The new viewer's camera stream skips the camera the last one left.
        self.camera_start = self.share.as_ref().map_or(0, camera_start);
        let started = find(
            std::env::var("RTC_METAL").ok().as_deref(),
            Path::new(env!("CARGO_MANIFEST_DIR")),
        )
        .and_then(|exe| {
            let at = (self.frame.x + self.frame.width, self.frame.y);
            let at = (f64::from(at.0), f64::from(at.1));
            Viewer::start(&exe, preview.path(), at, &own_share_path())
        });
        match started {
            Ok(v) => {
                self.viewer_pid = v.pid();
                self.live = true;
                children.viewer = Some(v);
                self.status = "viewer started".into();
            }
            Err(e) => self.status = e,
        }
    }

    /// The selected item's boxes, in display order.
    fn order(&self) -> Vec<Field> {
        let Some(node) = &self.selected else {
            return vec![];
        };
        fields(&self.desc, node)
            .into_iter()
            .map(|(f, _)| f)
            .collect()
    }

    /// Tab (`step` 1) or Shift-Tab (-1): edit the box next to the one being
    /// edited or shaded, wrapping; from the tree, the first or last box.
    fn tab(&mut self, step: isize) -> Task<Message> {
        let order = self.order();
        if order.is_empty() {
            return Task::none();
        }
        let n = order.len() as isize;
        let next = match self.editing.or(self.left) {
            Some(f) => order
                .iter()
                .position(|o| *o == f)
                .map_or(0, |i| (i as isize + step).rem_euclid(n) as usize),
            None if step > 0 => 0,
            None => order.len() - 1,
        };
        self.start_editing(order[next])
    }

    /// Focus `field` for typing; the cursor goes to the end.
    fn start_editing(&mut self, field: Field) -> Task<Message> {
        self.left = None;
        self.editing = Some(field);
        operation::focus(box_id(field))
    }

    /// An arrow or h/j/k/l: move the shading, or else the tree.
    fn arrow(&mut self, right: isize, down: isize) {
        if self.editing.is_some() {
            return;
        }
        if let Some(left) = self.left {
            self.left = Some(next_box(&self.order(), left, right, down));
            return;
        }
        if down != 0 {
            let rows = rows(&self.desc, &self.folded);
            if let Some(node) = next_row(&rows, self.selected.as_ref(), down) {
                self.select(node);
            }
            return;
        }
        let Some(node) = self.selected.clone() else {
            return;
        };
        let foldable = rows(&self.desc, &self.folded)
            .iter()
            .any(|r| r.node == node && r.foldable);
        if foldable && right < 0 {
            self.folded.insert(node);
        } else if foldable {
            self.folded.remove(&node);
        }
    }

    /// A plain character no box took: the second key of `a` or `d`, else a
    /// move, a letter for the shaded box, or the first key of `a` or `d`.
    fn key(&mut self, c: char) -> Task<Message> {
        match self.pending.take() {
            Some('a') => {
                match Item::from_key(c) {
                    Some(item) => self.add_item(item),
                    None => self.status.clear(),
                }
                return Task::none();
            }
            Some('d') => {
                if c == 'd' {
                    self.delete_selected();
                }
                return Task::none();
            }
            _ => {}
        }
        match c {
            'h' => self.arrow(-1, 0),
            'l' => self.arrow(1, 0),
            'k' => self.arrow(0, -1),
            'j' => self.arrow(0, 1),
            'i' | 'c' | 'x' | 'u' => return self.letter(c),
            'a' | 'd' if self.left.is_none() && self.editing.is_none() => {
                if c == 'a' {
                    self.pending = Some('a');
                    self.status =
                        "add: s sphere  c cube  p plane  y cylinder  o cone  l light".into();
                } else if self.selected.is_some() {
                    self.pending = Some('d');
                }
            }
            _ => {}
        }
        Task::none()
    }

    /// Add `item` next to the selected one and select it, unfolding what
    /// hides it.
    fn add_item(&mut self, item: Item) {
        let before = self.desc.clone();
        let node = edit::add(&mut self.desc, self.selected.as_ref(), item);
        let label = format!("add {}", self.item_name(&node));
        self.history.record(&before, Group::Once, &label);
        if let Node::Shape(p) = &node {
            rtc_edit::edit::shift(&mut self.folded, p, 1);
        }
        for a in ancestors(&node) {
            self.folded.remove(&a);
        }
        self.select(node);
        self.status.clear();
        self.unsaved = true;
        self.unloaded = true;
    }

    /// Delete the selected item and select the one that takes its place, or
    /// say why it can't go.
    fn delete_selected(&mut self) {
        let Some(node) = self.selected.clone() else {
            return;
        };
        let before = self.desc.clone();
        let label = format!("delete {}", self.item_name(&node));
        match edit::delete(&mut self.desc, &node) {
            Ok(next) => {
                self.history.record(&before, Group::Once, &label);
                if let Node::Shape(p) = &node {
                    rtc_edit::edit::shift(&mut self.folded, p, -1);
                }
                self.select(next);
                self.status.clear();
                self.unsaved = true;
                self.unloaded = true;
            }
            Err(reason) => self.status = reason,
        }
    }

    /// i, c, x or u on the shaded box.
    fn letter(&mut self, c: char) -> Task<Message> {
        let (Some(field), Some(node)) = (self.left, &self.selected) else {
            return Task::none();
        };
        match c {
            'i' => self
                .start_editing(field)
                .chain(operation::move_cursor_to_front(box_id(field))),
            'c' => {
                self.edit_as(field, String::new(), "change");
                self.start_editing(field)
            }
            'x' => {
                // Only a box that may be empty; one that needs a number stays.
                let before = self.desc.clone();
                if set_remembering(&mut self.desc, node, field, "", &mut self.removed_blur).is_ok()
                {
                    if self.desc != before {
                        let label = format!("empty {}", row_label(field));
                        self.history.record(&before, Group::Once, &label);
                    }
                    self.texts.insert(field, String::new());
                    self.errors.remove(&field);
                    self.unsaved = true;
                    self.unloaded = true;
                    self.refresh_camera(Some(field));
                }
                Task::none()
            }
            'u' => {
                self.reset(field);
                Task::none()
            }
            _ => Task::none(),
        }
    }

    /// Write the preview and send it to the viewer, which reloads the
    /// edits keeping its view.
    fn load(&mut self) {
        let children = CHILDREN.lock().unwrap();
        let Some(preview) = &children.preview else {
            let e = self.preview_error.as_deref().unwrap_or("none");
            self.status = format!("no preview file: {e}");
            return;
        };
        let sent =
            preview
                .write(&self.desc, &self.dir())
                .and_then(|()| match (&self.share, self.live) {
                    (Some(share), true) => {
                        send_or_hold(share, &mut self.held, preview.path(), false)
                    }
                    _ => Err("no viewer running".into()),
                });
        self.status = match sent {
            Ok(()) => {
                self.unloaded = false;
                self.warn("loaded into viewer".into())
            }
            Err(e) => e,
        };
    }

    fn save(&mut self) {
        let Some(path) = &self.path else {
            self.save_as();
            return;
        };
        self.status = match self.desc.save(path) {
            Ok(()) => {
                self.unsaved = false;
                self.history.saved(&self.desc);
                self.saved()
            }
            Err(e) => e.to_string(),
        };
    }

    /// `message`, with a warning when a box holds text that isn't in the scene.
    fn warn(&self, message: String) -> String {
        if self.errors.is_empty() {
            message
        } else {
            format!("{message} (a field has an error)")
        }
    }

    fn select(&mut self, node: Node) {
        self.left = None;
        self.editing = None;
        self.texts = fields(&self.desc, &node).into_iter().collect();
        self.errors.clear();
        self.removed_blur = None;
        self.selected = Some(node);
    }

    fn edit(&mut self, field: Field, value: String) {
        self.edit_as(field, value, "typing in");
    }

    /// `edit`, with the undo step named `<verb> <row label>`.
    fn edit_as(&mut self, field: Field, value: String, verb: &str) {
        self.left = None;
        self.editing = Some(field);
        let Some(node) = self.selected.clone() else {
            return;
        };
        let before = self.desc.clone();
        match set_remembering(&mut self.desc, &node, field, &value, &mut self.removed_blur) {
            Ok(()) => {
                if self.desc != before {
                    let label = format!("{verb} {}", row_label(field));
                    self.history
                        .record(&before, Group::Box(node, field), &label);
                }
                self.errors.remove(&field);
                self.unsaved = true;
                self.unloaded = true;
                self.refresh_camera(Some(field));
            }
            Err(e) => {
                self.errors.insert(field, e);
            }
        }
        self.texts.insert(field, value);
    }

    /// Some camera boxes depend on each other: an f-number creates or removes
    /// the other focal blur boxes. Show the scene's values in all of them but
    /// `edited` and the red ones.
    fn refresh_camera(&mut self, edited: Option<Field>) {
        let Some(node @ Node::Camera) = &self.selected else {
            return;
        };
        for (f, t) in fields(&self.desc, node) {
            if Some(f) != edited && !self.errors.contains_key(&f) {
                self.texts.insert(f, t);
            }
        }
    }

    /// Undo the last step (`back`) or redo the last undone one, and show
    /// the scene as it now is.
    fn undo(&mut self, back: bool) {
        self.follower.stop();
        let done = if back {
            self.history.undo(&mut self.desc)
        } else {
            self.history.redo(&mut self.desc)
        };
        let Some(label) = done else {
            self.status = if back {
                "nothing to undo"
            } else {
                "nothing to redo"
            }
            .into();
            return;
        };
        self.status = format!("{} {label}", if back { "undid:" } else { "redid:" });
        self.unsaved = self.history.differs(&self.desc);
        self.unloaded = true;
        self.editing = None;
        self.left = None;
        self.errors.clear();
        match self.selected.clone() {
            Some(node)
                if rows(&self.desc, &self.folded)
                    .iter()
                    .any(|r| r.node == node) =>
            {
                self.select(node)
            }
            _ => {
                self.selected = None;
                self.texts.clear();
            }
        }
    }

    /// An item's tree label without its comment, such as `sphere`.
    fn item_name(&self, node: &Node) -> String {
        rows(&self.desc, &HashSet::new())
            .into_iter()
            .find(|r| r.node == *node)
            .map(|r| r.label.split(" — ").next().unwrap_or("").to_string())
            .unwrap_or_else(|| "item".into())
    }

    /// The viewer with `pid` exited: if it is the running one, it ended.
    fn viewer_exited(&mut self, pid: u32) {
        if pid == self.viewer_pid {
            self.viewer_ended();
        }
    }

    /// No viewer runs: end its camera stream (once), say so, and stop
    /// following.
    fn viewer_ended(&mut self) {
        if self.live {
            if let Some(share) = &self.share {
                share.end_viewer();
            }
            self.status = "viewer closed".into();
        }
        self.live = false;
        self.seen = None;
        // The next viewer starts on the preview as it is.
        self.held = None;
        self.follower.stop();
    }

    /// The viewer wrote its camera: on a new write, show it and, following,
    /// copy it into the scene.
    fn camera(&mut self, shot: Snapshot) {
        // Until the viewer is ready the camera may be the last viewer's.
        if !self.share.as_ref().is_some_and(Shared::is_ready) {
            return;
        }
        // The new viewer says it is ready before its first camera write:
        // the load it may have missed goes now.
        if self.held.is_some()
            && let Some(share) = &self.share
            && let Some(preview) = &CHILDREN.lock().unwrap().preview
            && let Err(e) = send_or_hold(share, &mut self.held, preview.path(), false)
        {
            self.status = e;
        }
        if self.seen.is_some_and(|s| s.seq == shot.seq) {
            return;
        }
        self.seen = Some(shot);
        self.follow(shot);
    }

    /// Copy `shot` into the camera while following.
    fn follow(&mut self, shot: Snapshot) {
        let before = self.desc.clone();
        if self.follower.update(&mut self.desc.camera, &shot) {
            self.history
                .record(&before, Group::Follow, "follow the viewer");
            self.unsaved = true;
            self.show_camera();
        }
    }

    /// Show the camera's values in its boxes after following changed them,
    /// clearing the errors of the boxes it wrote.
    fn show_camera(&mut self) {
        for f in [0, 1, 2].map(Field::CameraFrom) {
            self.errors.remove(&f);
        }
        self.errors.remove(&Field::FocalDistance);
        self.refresh_camera(None);
    }

    /// Put the scene's value back in `field` and clear its error.
    fn reset(&mut self, field: Field) {
        let Some(node) = &self.selected else { return };
        if let Some((_, t)) = fields(&self.desc, node)
            .into_iter()
            .find(|(f, _)| *f == field)
        {
            self.texts.insert(field, t);
        }
        self.errors.remove(&field);
    }

    /// The status line, or why the shaded or edited box is red.
    fn shown_status(&self) -> String {
        let reason = self
            .left
            .or(self.editing)
            .and_then(|f| self.errors.get(&f).map(|why| (f, why)));
        match reason {
            Some((field, why)) => {
                let (_, label, axis) = describe(field, "");
                match axis {
                    Some(a) => format!("{label} {}: {why}", ["x", "y", "z"][a]),
                    None => format!("{label}: {why}"),
                }
            }
            None => self.status.clone(),
        }
    }

    /// Give the preview the scene's name, and have the viewer reload it
    /// keeping its view, so its title and the image name it offers follow
    /// a Save As. The preview keeps what the viewer last loaded.
    fn rename_preview(&mut self) {
        let name = scene_name(self.path.as_deref());
        let mut children = CHILDREN.lock().unwrap();
        let running = children.viewer.as_mut().is_some_and(Viewer::running);
        let Some(preview) = children.preview.as_mut() else {
            return;
        };
        if preview.path().file_name() == Some(OsStr::new(&name)) {
            return;
        }
        let sent = preview
            .rename(OsStr::new(&name))
            .and_then(|()| match (&self.share, running) {
                (Some(share), true) => send_or_hold(share, &mut self.held, preview.path(), false),
                _ => Ok(()),
            });
        if let Err(e) = sent {
            self.status = e;
        }
    }

    /// Ask where to save, rebasing asset paths; true when it saved.
    fn save_as(&mut self) -> bool {
        let dir = self.dir();
        let picked = rfd::FileDialog::new()
            .set_directory(&dir)
            .set_file_name(self.name())
            .add_filter("scene", &["json"])
            .save_file();
        // The key-up of Cmd went to the dialog.
        self.command = false;
        let Some(new) = picked else {
            return false;
        };
        let mut desc = self.desc.clone();
        let saved = rebase(&mut desc, &dir, dir_of(&new))
            .and_then(|()| desc.save(&new).map_err(|e| e.to_string()));
        match saved {
            Ok(()) => {
                // The undo and redo copies get the same new asset paths. A
                // copy whose paths can't be rewritten keeps them.
                self.history.for_each_scene(|d| {
                    let _ = rebase(d, &dir, dir_of(&new));
                });
                self.desc = desc;
                self.path = Some(new);
                self.unsaved = false;
                self.history.saved(&self.desc);
                // Rebasing touches only asset paths, which have no boxes,
                // so the boxes and any red text stay as they are.
                self.status = self.saved();
                self.rename_preview();
                true
            }
            Err(e) => {
                self.status = e;
                false
            }
        }
    }
}

/// Send the viewer the load of `preview`, or hold it in `held` while the
/// viewer isn't ready, as a load sent then may be missed.
fn send_or_hold(
    share: &Shared,
    held: &mut Option<bool>,
    preview: &Path,
    fresh: bool,
) -> Result<(), String> {
    let fresh = held_fresh(held.take(), fresh);
    if share.is_ready() {
        share.send_load(preview, fresh)
    } else {
        *held = Some(fresh);
        Ok(())
    }
}

/// `path`'s folder, or the working directory for a scene with no file.
fn dir_for(path: Option<&Path>) -> PathBuf {
    match path {
        Some(p) => dir_of(p).to_path_buf(),
        None => std::env::current_dir().unwrap_or_else(|_| PathBuf::from(".")),
    }
}

/// `Message::Key` for a key that does nothing but cancel a wait for a second one.
const OTHER: char = '\0';

fn update(e: &mut Editor, message: Message) -> Task<Message> {
    let editing = e.editing;
    let task = handle(e, message);
    // Leaving a box, for none or another, ends its undo step. Entering one
    // from none doesn't: its first keystroke both enters it and records.
    if editing.is_some() && e.editing != editing {
        e.history.close();
    }
    let state = e.menu_state();
    if let Some(menus) = &mut e.menus {
        menus.refresh(&state);
    }
    task
}

/// The message a menu command sends.
fn menu_message(command: Command, following: bool) -> Message {
    match command {
        Command::New => Message::New,
        Command::Open => Message::Open,
        Command::Save => Message::Save,
        Command::SaveAs => Message::SaveAs,
        Command::Undo => Message::Undo,
        Command::Redo => Message::Redo,
        Command::Delete => Message::Delete,
        Command::Add(item) => Message::Add(item),
        Command::Load => Message::Load,
        Command::StartViewer => Message::ShowViewer,
        Command::Follow => Message::Follow(!following),
        Command::Keep => Message::Keep,
    }
}

fn handle(e: &mut Editor, message: Message) -> Task<Message> {
    // Anything but a key, Cmd or a focus answer ends a wait for a second key;
    // Esc only does that.
    if !matches!(
        message,
        Message::Key(_)
            | Message::Command(_)
            | Message::Focused(..)
            | Message::Window(_)
            | Message::ViewerExited(_)
            | Message::Camera(_)
            | Message::CameraLost(_)
            | Message::MenuEvent(_)
    ) && let Some(waiting) = e.pending.take()
    {
        if waiting == 'a' {
            e.status.clear();
        }
        if matches!(message, Message::Escape(_)) {
            return Task::none();
        }
    }
    match message {
        Message::Select(node) => e.select(node),
        Message::Toggle(node) => {
            if !e.folded.remove(&node) {
                e.folded.insert(node);
            }
        }
        Message::Save => e.save(),
        Message::SaveAs => {
            e.save_as();
        }
        Message::New => e.new_file(),
        Message::Open => e.open(),
        Message::ShowViewer => e.show_viewer(),
        Message::Edit(field, value) => {
            // A focused box also takes Cmd-L's `l` as typing, over any
            // selection; with Cmd held, keep only a cut.
            let old = e.texts.get(&field).map_or("", String::as_str);
            if !e.command || is_cut(old, &value) {
                e.edit(field, value);
            }
        }
        Message::Closed(value) => {
            if let Some(node) = &e.selected {
                let before = e.desc.clone();
                set_closed(&mut e.desc, node, value);
                e.history.record(&before, Group::Once, "closed");
                e.unsaved = true;
                e.unloaded = true;
            }
        }
        Message::Paste(field, value) => e.edit(field, value),
        Message::Load => e.load(),
        Message::ViewerExited(pid) => e.viewer_exited(pid),
        Message::Camera(shot) => e.camera(shot),
        Message::CameraLost(why) => {
            e.follower.stop();
            e.status = format!("lost the viewer's camera: {why}");
        }
        Message::Follow(true) => {
            if e.live {
                e.follower.start(&e.desc.camera);
                if let Some(shot) = e.seen {
                    e.follow(shot);
                }
            }
        }
        Message::Follow(false) => {
            let before = e.desc.clone();
            if e.follower.untick(&mut e.desc.camera) {
                e.history.record(&before, Group::Once, "untick follow");
                e.unsaved = true;
                e.show_camera();
            }
        }
        Message::Keep => {
            e.follower.keep();
            e.history.close();
        }
        Message::Undo => e.undo(true),
        Message::Redo => e.undo(false),
        Message::MenuEvent(event) => {
            if let Some(c) = menu::command(&event.id.0) {
                let message = menu_message(c, e.follower.following());
                return handle(e, message);
            }
        }
        Message::Submit(field) => {
            e.left = Some(field);
            e.editing = None;
            // Focusing an id no widget has unfocuses them all.
            return operation::focus(iced::widget::Id::new("none"));
        }
        Message::FocusNext => return e.tab(1),
        Message::FocusPrevious => return e.tab(-1),
        Message::Move(right, down) => e.arrow(right, down),
        Message::Enter => {
            if let Some(field) = e.left {
                return e.start_editing(field);
            }
            if e.editing.is_none() {
                e.left = e.order().first().copied();
            }
        }
        Message::Key(c) => return e.key(c),
        Message::Add(item) => e.add_item(item),
        Message::Delete => e.delete_selected(),
        Message::Escape(true) => {
            if let Some(field) = e.editing.filter(|f| e.errors.contains_key(f)) {
                e.reset(field);
            }
            e.left = e.editing.take();
        }
        Message::Escape(false) => e.left = None,
        Message::Clicked(taken) => {
            if taken {
                e.left = None;
            }
            e.editing = None;
            // iced doesn't say which box a click focused; ask each one.
            return Task::batch(e.order().into_iter().map(|f| {
                operation::is_focused(box_id(f)).map(move |yes| Message::Focused(f, yes))
            }));
        }
        Message::Focused(field, yes) => {
            if yes {
                e.editing = Some(field);
            }
        }
        Message::Command(held) => e.command = held,
        Message::Window(event) => match event {
            window::Event::Opened { position, size } => {
                if e.menus.is_none() {
                    let menus = Menus::new();
                    menus.install();
                    e.menus = Some(menus);
                }
                e.frame.x = position.map_or(e.frame.x, |p| p.x);
                e.frame.y = position.map_or(e.frame.y, |p| p.y);
                (e.frame.width, e.frame.height) = (size.width, size.height);
            }
            window::Event::Moved(p) => (e.frame.x, e.frame.y) = (p.x, p.y),
            window::Event::Resized(size) => {
                (e.frame.width, e.frame.height) = (size.width, size.height)
            }
            _ => {}
        },
    }
    Task::none()
}

/// Tab and Shift-Tab, also while a box has focus and takes the key itself;
/// whether Cmd is held; Esc; the arrows, h/j/k/l, Enter and the letters
/// a c d h i j k l o p s u x y that no box took, and any other such key as one
/// that only cancels; and clicks. The Cmd keys belong to the menu bar.
fn shortcut(
    event: iced::Event,
    status: iced::event::Status,
    _: iced::window::Id,
) -> Option<Message> {
    let ignored = status == iced::event::Status::Ignored;
    let event = match event {
        iced::Event::Keyboard(event) => event,
        iced::Event::Mouse(iced::mouse::Event::ButtonPressed(_)) => {
            return Some(Message::Clicked(!ignored));
        }
        _ => return None,
    };
    let (key, modifiers) = match event {
        keyboard::Event::ModifiersChanged(m) => return Some(Message::Command(m.command())),
        keyboard::Event::KeyPressed { key, modifiers, .. } => (key, modifiers),
        _ => return None,
    };
    use keyboard::key::Named;
    let plain = ignored && modifiers.is_empty();
    match key {
        keyboard::Key::Named(Named::Tab) if modifiers.shift() => Some(Message::FocusPrevious),
        keyboard::Key::Named(Named::Tab) => Some(Message::FocusNext),
        keyboard::Key::Named(Named::Escape) => Some(Message::Escape(!ignored)),
        keyboard::Key::Named(named) if plain => match named {
            Named::Enter => Some(Message::Enter),
            Named::ArrowLeft => Some(Message::Move(-1, 0)),
            Named::ArrowRight => Some(Message::Move(1, 0)),
            Named::ArrowUp => Some(Message::Move(0, -1)),
            Named::ArrowDown => Some(Message::Move(0, 1)),
            Named::Shift | Named::Control | Named::Alt | Named::Super | Named::Meta => None,
            // Any other key only cancels a wait for a second key.
            _ => Some(Message::Key(OTHER)),
        },
        keyboard::Key::Character(c) if ignored && !modifiers.command() => match c.as_str() {
            "h" | "j" | "k" | "l" | "i" | "c" | "x" | "u" | "a" | "d" | "s" | "p" | "y" | "o"
                if modifiers.is_empty() =>
            {
                c.chars().next().map(Message::Key)
            }
            _ => Some(Message::Key(OTHER)),
        },
        _ => None,
    }
}

fn subscription(e: &Editor) -> Subscription<Message> {
    let mut subscriptions = vec![
        iced::event::listen_with(shortcut),
        window::events().map(|(_, event)| Message::Window(event)),
        Subscription::run(menu_events),
    ];
    if e.live
        && let Some(share) = &e.share
    {
        let key = (own_share_path(), share.generation(), e.camera_start);
        subscriptions.push(Subscription::run_with(key, camera_stream));
        subscriptions.push(Subscription::run_with(e.viewer_pid, viewer_exit));
    }
    Subscription::batch(subscriptions)
}

/// The viewer's camera writes, from the shared file at `key.0`, while the
/// viewer of generation `key.1` runs: a thread sleeps on the file and sends
/// each write on.
fn camera_stream(key: &(PathBuf, u32, u32)) -> UnboundedReceiver<Message> {
    let (path, generation, start) = key.clone();
    let (tx, rx) = unbounded();
    std::thread::spawn(move || {
        let shared = match Shared::open(&path) {
            Ok(s) => s,
            Err(e) => {
                let _ = tx.unbounded_send(Message::CameraLost(e));
                return;
            }
        };
        let mut last = start;
        loop {
            match shared.wait_camera(last, generation) {
                Ok(Some(shot)) => {
                    last = shot.seq;
                    if tx.unbounded_send(Message::Camera(shot)).is_err() {
                        return;
                    }
                }
                Ok(None) => return,
                Err(e) => {
                    let _ = tx.unbounded_send(Message::CameraLost(e));
                    return;
                }
            }
        }
    });
    rx
}

/// The menu bar's clicks and keys: muda calls its event handler on the
/// main thread for each, which forwards it here. A restarted subscription
/// sets the handler again for its own channel.
fn menu_events() -> UnboundedReceiver<Message> {
    let (tx, rx) = unbounded();
    MenuEvent::set_event_handler(Some(move |event| {
        let _ = tx.unbounded_send(Message::MenuEvent(event));
    }));
    rx
}

/// `ViewerExited(pid)` when the viewer `pid` exits, from a GCD proc
/// source kept alive by the stream.
fn viewer_exit(pid: &u32) -> impl Stream<Item = Message> + use<> {
    let pid = *pid;
    let (tx, rx) = unbounded();
    let notify = tx.clone();
    let source = exit_source(pid, move || {
        let _ = notify.unbounded_send(Message::ViewerExited(pid));
    });
    // An exit before the source was active goes unreported: look once.
    let running = CHILDREN
        .lock()
        .unwrap()
        .viewer
        .as_mut()
        .is_some_and(|v| v.pid() == pid && v.running());
    if !running {
        let _ = tx.unbounded_send(Message::ViewerExited(pid));
    }
    rx.map(move |m| {
        let _source = &source;
        m
    })
}

fn view(e: &Editor) -> Element<'_, Message> {
    let mut tree = Column::new().spacing(1);
    for r in rows(&e.desc, &e.folded) {
        let indent = "    ".repeat(r.depth);
        let marker: Element<Message> = if r.foldable {
            let open = !e.folded.contains(&r.node);
            button(text(format!("{indent}{}", if open { "▾" } else { "▸" })).size(14))
                .on_press(Message::Toggle(r.node.clone()))
                .style(button::text)
                .into()
        } else {
            button(text(format!("{indent}  ")).size(14))
                .style(button::text)
                .into()
        };
        let selected = e.selected.as_ref() == Some(&r.node);
        let label = button(text(r.label).size(14))
            .on_press(Message::Select(r.node))
            .width(Length::Fill)
            .style(if selected {
                button::primary
            } else {
                button::text
            });
        tree = tree.push(row![marker, label]);
    }
    let inspector: Element<Message> = match &e.selected {
        Some(node) => column![
            boxes(e, node),
            text(json(&e.desc, node)).font(Font::MONOSPACE).size(13)
        ]
        .spacing(12)
        .into(),
        None => text("select an item").into(),
    };
    let can_delete = e
        .selected
        .as_ref()
        .is_some_and(|n| edit::can_delete(&e.desc, n).is_ok());
    let add = |label: &'static str, item| {
        button(text(label).size(13))
            .on_press(Message::Add(item))
            .padding([2, 8])
    };
    let bar = row![
        add("+ sphere", Item::Sphere),
        add("+ cube", Item::Cube),
        add("+ plane", Item::Plane),
        add("+ cylinder", Item::Cylinder),
        add("+ cone", Item::Cone),
        add("+ light", Item::Light),
        space::horizontal(),
        button(text("Delete").size(13))
            .on_press_maybe(can_delete.then_some(Message::Delete))
            .padding([2, 8])
            .style(button::danger),
    ]
    .spacing(4);
    column![
        container(bar).padding([4, 8]),
        row![
            container(scrollable(tree)).width(260).padding(8),
            container(scrollable(inspector).width(Length::Fill))
                .width(Length::Fill)
                .padding(8),
        ]
        .height(Length::Fill),
        container(row![
            text(e.shown_status()).size(13),
            space::horizontal(),
            text(if e.unloaded {
                "not in viewer — Cmd-L to load"
            } else {
                ""
            })
            .size(13),
        ])
        .padding([4, 8]),
    ]
    .into()
}

/// The selected item's boxes, under a heading per group; a vector's three
/// axes share a row.
fn boxes<'a>(e: &'a Editor, node: &Node) -> Element<'a, Message> {
    let mut out = Column::new().spacing(4);
    let mut heading = "";
    let mut current: Option<Row<'a, Message>> = None;
    let is_cone = serde_json::from_str::<serde_json::Value>(&json(&e.desc, node))
        .is_ok_and(|v| v["type"] == "cone");
    let bounds = if is_cone { "Cone" } else { "Cylinder" };
    for (field, _) in fields(&e.desc, node) {
        let (group, label, axis) = describe(field, bounds);
        if group != heading {
            if let Some(r) = current.take() {
                out = out.push(r);
            }
            if heading == bounds {
                out = out.push(closed_box(e, node));
            }
            out = out.push(text(group).size(14));
            heading = group;
        }
        if axis.is_none_or(|a| a == 0) {
            if let Some(r) = current.take() {
                out = out.push(r);
            }
            if field == Field::CameraFrom(0) && e.live {
                out = out.push(follow_controls(e));
            }
            current = Some(row![text(label).size(13).width(90)].spacing(4));
        }
        let input = text_input("", e.texts.get(&field).map_or("", String::as_str))
            .on_submit(Message::Submit(field))
            .id(box_id(field))
            .size(13);
        // While following, the viewer writes `from`; its boxes are read-only.
        let input = if e.follower.following() && matches!(field, Field::CameraFrom(_)) {
            input
        } else {
            input
                .on_input(move |s| Message::Edit(field, s))
                .on_paste(move |s| Message::Paste(field, s))
        };
        // A comment is free text, wide as a row of three number boxes.
        let input = if field == Field::Comment {
            input.width(248)
        } else {
            input.font(Font::MONOSPACE).width(80)
        };
        let error = e.errors.contains_key(&field);
        let marked = e.left == Some(field);
        let input = input.style(move |theme: &iced::Theme, status| {
            let mut style = text_input::default(theme, status);
            if marked {
                style.background = theme.palette().primary.scale_alpha(0.3).into();
            }
            if error {
                style.border.color = theme.palette().danger;
                style.border.width = 2.0;
            }
            style
        });
        current = current.map(|r| r.push(input));
        if let Field::CameraFrom(2) | Field::CameraTo(2) | Field::CameraUp(2) = field {
            current = current.map(|r| r.push(viewer_column(e, field)));
        }
    }
    if let Some(r) = current.take() {
        out = out.push(r);
    }
    if let Some(radius) = lens_radius(&e.desc).filter(|_| *node == Node::Camera) {
        let radius = (radius * 1e6).round() / 1e6;
        out = out.push(text(format!("lens radius {radius}")).size(13));
    }
    out.into()
}

/// Beside a camera vector's boxes: the viewer's value, greyed. Without a
/// viewer, or before it writes, one word on the `from` row.
fn viewer_column<'a>(e: &Editor, field: Field) -> Element<'a, Message> {
    let is_from = matches!(field, Field::CameraFrom(_));
    let grey = |theme: &iced::Theme| text::Style {
        color: Some(theme.palette().text.scale_alpha(0.5)),
    };
    let shown = match (e.live, e.seen) {
        (true, Some(s)) => {
            let v = match field {
                Field::CameraFrom(_) => s.from,
                Field::CameraTo(_) => s.to,
                _ => s.up,
            };
            format!("{:>8.3} {:>8.3} {:>8.3}", v[0], v[1], v[2])
        }
        (true, None) if is_from => "waiting".into(),
        (false, _) if is_from => "no viewer".into(),
        _ => String::new(),
    };
    row![
        text("viewer").size(13).style(grey),
        text(shown).font(Font::MONOSPACE).size(13).style(grey),
    ]
    .spacing(8)
    .padding(iced::Padding::ZERO.left(12))
    .into()
}

/// Above the `from` row, with a viewer running: the follow checkbox and,
/// while following, Keep.
fn follow_controls<'a>(e: &Editor) -> Element<'a, Message> {
    let mut r = row![
        checkbox(e.follower.following())
            .label("follow the viewer")
            .on_toggle(Message::Follow)
            .size(14)
    ]
    .spacing(12)
    .align_y(iced::Alignment::Center);
    if e.follower.following() {
        r = r.push(button(text("Keep").size(13)).on_press(Message::Keep));
    }
    r.into()
}

fn closed_box<'a>(e: &Editor, node: &Node) -> Element<'a, Message> {
    checkbox(closed(&e.desc, node).unwrap_or(false))
        .label("closed")
        .on_toggle(Message::Closed)
        .size(14)
        .into()
}

/// A box's row label, such as `position`, for undo steps.
fn row_label(field: Field) -> &'static str {
    describe(field, "bounds").1
}

/// A box's heading, its row label, and its axis when it is part of a vector.
/// `bounds` is the heading for min and max: `Cylinder` or `Cone`.
fn describe(field: Field, bounds: &'static str) -> (&'static str, &'static str, Option<usize>) {
    use Field::*;
    match field {
        Comment => ("", "comment", None),
        Position(a) => ("Transform", "position", Some(a)),
        Rotation(a) => ("Transform", "rotation°", Some(a)),
        Scale(a) => ("Transform", "scale", Some(a)),
        Min => (bounds, "min", None),
        Max => (bounds, "max", None),
        Color(a) => ("Material", "color", Some(a)),
        Ambient => ("Material", "ambient", None),
        Diffuse => ("Material", "diffuse", None),
        Specular => ("Material", "specular", None),
        Shininess => ("Material", "shininess", None),
        Reflective => ("Material", "reflective", None),
        Transparency => ("Material", "transparency", None),
        RefractiveIndex => ("Material", "refractive index", None),
        LightPosition(a) => ("Light", "position", Some(a)),
        LightIntensity(a) => ("Light", "intensity", Some(a)),
        CameraFrom(a) => ("View", "from", Some(a)),
        CameraTo(a) => ("View", "to", Some(a)),
        CameraUp(a) => ("View", "up", Some(a)),
        FieldOfView => ("Lens", "field of view°", None),
        Width => ("Lens", "width", None),
        Height => ("Lens", "height", None),
        Antialias => ("Antialias", "n", None),
        FNumber => ("Focal blur", "f-number", None),
        FocalDistance => ("Focal blur", "focal distance", None),
        FocalSamples => ("Focal blur", "samples", None),
    }
}

/// The widget id of a box, so Tab can focus it.
fn box_id(field: Field) -> iced::widget::Id {
    format!("{field:?}").into()
}
