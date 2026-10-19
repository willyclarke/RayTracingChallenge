//! rtc-edit: `rtc-edit [scene.json]` shows the scene as a tree, starts
//! `rtc-metal --window --watch` on a preview copy of it, and saves it
//! (Cmd-S, and Cmd-Shift-S to save elsewhere with asset paths rebased).
//! Cmd-R starts the viewer again after it was closed. Simple shapes,
//! point lights and the camera are edited in fixed fields; Cmd-L loads the
//! edits into the viewer.

use iced::widget::operation;
use iced::widget::{
    Column, Row, button, checkbox, column, container, row, scrollable, space, text, text_input,
};
use iced::{Element, Font, Length, Subscription, Task, keyboard};
use rtc_edit::edit::{self, Item, ancestors, new_scene};
use rtc_edit::fields::{
    Field, closed, fields, is_cut, lens_radius, next_box, set_closed, set_remembering,
};
use rtc_edit::paths::{dir_of, rebase, scene_name};
use rtc_edit::preview::Preview;
use rtc_edit::tree::{Node, json, next_row, rows};
use rtc_edit::viewer::{Viewer, find};
use rtc_rust::scene::{FocalBlur, SceneDescription};
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

/// Close the viewer and remove the preview folder.
extern "C" fn clean_up() {
    if let Ok(mut children) = CHILDREN.lock() {
        children.viewer = None;
        children.preview = None;
    }
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
}

fn main() -> iced::Result {
    iced::application(boot, update, view)
        .title(|e: &Editor| {
            let mark = if e.unsaved { " •" } else { "" };
            format!("rtc-edit — {}{mark}", e.name())
        })
        .subscription(subscription)
        .window_size((900.0, 640.0))
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
    };
    match loaded {
        Some((desc, path)) => editor.switch(desc, Some(path)),
        None => editor.switch(new_scene(), None),
    }
    editor
}

impl Editor {
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

    /// Replace the scene: the old viewer ends (and offers its image), the old
    /// preview goes, and a preview and viewer start on `desc`.
    fn switch(&mut self, desc: SceneDescription, path: Option<PathBuf>) {
        {
            let mut children = CHILDREN.lock().unwrap();
            children.viewer = None;
            children.preview = None;
            let name = scene_name(path.as_deref());
            let preview = Preview::new(&std::env::temp_dir(), OsStr::new(&name))
                .and_then(|p| p.write(&desc, &dir_for(path.as_deref())).map(|()| p));
            self.preview_error = None;
            match preview {
                Ok(p) => children.preview = Some(p),
                Err(e) => self.preview_error = Some(e),
            }
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
        self.show_viewer();
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

    /// Start the viewer unless it is still open.
    fn show_viewer(&mut self) {
        let mut children = CHILDREN.lock().unwrap();
        if children.viewer.as_mut().is_some_and(Viewer::running) {
            return;
        }
        let Some(preview) = &children.preview else {
            let e = self.preview_error.as_deref().unwrap_or("none");
            self.status = format!("no preview file: {e}");
            return;
        };
        let started = find(
            std::env::var("RTC_METAL").ok().as_deref(),
            Path::new(env!("CARGO_MANIFEST_DIR")),
        )
        .and_then(|exe| Viewer::start(&exe, preview.path()));
        match started {
            Ok(v) => {
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
        let node = edit::add(&mut self.desc, self.selected.as_ref(), item);
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
        match edit::delete(&mut self.desc, &node) {
            Ok(next) => {
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
                self.edit(field, String::new());
                self.start_editing(field)
            }
            'x' => {
                // Only a box that may be empty; one that needs a number stays.
                if set_remembering(&mut self.desc, node, field, "", &mut self.removed_blur).is_ok()
                {
                    self.texts.insert(field, String::new());
                    self.errors.remove(&field);
                    self.unsaved = true;
                    self.unloaded = true;
                    self.refresh_camera(field);
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

    /// Write the preview, so the viewer reloads the edits.
    fn load(&mut self) {
        let children = CHILDREN.lock().unwrap();
        let Some(preview) = &children.preview else {
            let e = self.preview_error.as_deref().unwrap_or("none");
            self.status = format!("no preview file: {e}");
            return;
        };
        self.status = match preview.write(&self.desc, &self.dir()) {
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
        self.left = None;
        self.editing = Some(field);
        let Some(node) = &self.selected else { return };
        match set_remembering(&mut self.desc, node, field, &value, &mut self.removed_blur) {
            Ok(()) => {
                self.errors.remove(&field);
                self.unsaved = true;
                self.unloaded = true;
                self.refresh_camera(field);
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
    fn refresh_camera(&mut self, edited: Field) {
        let Some(node @ Node::Camera) = &self.selected else {
            return;
        };
        for (f, t) in fields(&self.desc, node) {
            if f != edited && !self.errors.contains_key(&f) {
                self.texts.insert(f, t);
            }
        }
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
                self.desc = desc;
                self.path = Some(new);
                self.unsaved = false;
                // Rebasing touches only asset paths, which have no boxes,
                // so the boxes and any red text stay as they are.
                self.status = self.saved();
                true
            }
            Err(e) => {
                self.status = e;
                false
            }
        }
    }
}

/// `path`'s folder, or the working directory for a scene with no file.
fn dir_for(path: Option<&Path>) -> PathBuf {
    match path {
        Some(p) => dir_of(p).to_path_buf(),
        None => std::env::current_dir().unwrap_or_else(|_| PathBuf::from(".")),
    }
}

fn update(e: &mut Editor, message: Message) -> Task<Message> {
    // Anything but a key, Cmd or a focus answer ends a wait for a second key;
    // Esc only does that.
    if !matches!(
        message,
        Message::Key(_) | Message::Command(_) | Message::Focused(..)
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
                set_closed(&mut e.desc, node, value);
                e.unsaved = true;
                e.unloaded = true;
            }
        }
        Message::Paste(field, value) => e.edit(field, value),
        Message::Load => e.load(),
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
    }
    Task::none()
}

/// Cmd-S, Cmd-Shift-S, Cmd-N, Cmd-O, Cmd-R, Cmd-L, Tab and Shift-Tab, also while a box
/// has focus and takes the key itself; whether Cmd is held; Esc; the arrows,
/// h/j/k/l, Enter and the letters a c d h i j k l o p s u x y that no box took; and clicks.
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
    let c = match key {
        keyboard::Key::Named(Named::Tab) if modifiers.shift() => {
            return Some(Message::FocusPrevious);
        }
        keyboard::Key::Named(Named::Tab) => return Some(Message::FocusNext),
        keyboard::Key::Named(Named::Escape) => return Some(Message::Escape(!ignored)),
        keyboard::Key::Named(named) if plain => {
            return match named {
                Named::Enter => Some(Message::Enter),
                Named::ArrowLeft => Some(Message::Move(-1, 0)),
                Named::ArrowRight => Some(Message::Move(1, 0)),
                Named::ArrowUp => Some(Message::Move(0, -1)),
                Named::ArrowDown => Some(Message::Move(0, 1)),
                _ => None,
            };
        }
        keyboard::Key::Character(c) if plain => {
            return match c.as_str() {
                "h" | "j" | "k" | "l" | "i" | "c" | "x" | "u" | "a" | "d" | "s" | "p" | "y"
                | "o" => c.chars().next().map(Message::Key),
                _ => None,
            };
        }
        keyboard::Key::Character(c) if modifiers.command() => c,
        _ => return None,
    };
    match c.to_lowercase().as_str() {
        "s" if modifiers.shift() => Some(Message::SaveAs),
        "s" => Some(Message::Save),
        "n" => Some(Message::New),
        "o" => Some(Message::Open),
        "r" => Some(Message::ShowViewer),
        "l" => Some(Message::Load),
        _ => None,
    }
}

fn subscription(_: &Editor) -> Subscription<Message> {
    iced::event::listen_with(shortcut)
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
            current = Some(row![text(label).size(13).width(90)].spacing(4));
        }
        let input = text_input("", e.texts.get(&field).map_or("", String::as_str))
            .on_input(move |s| Message::Edit(field, s))
            .on_paste(move |s| Message::Paste(field, s))
            .on_submit(Message::Submit(field))
            .id(box_id(field))
            .font(Font::MONOSPACE)
            .size(13)
            .width(80);
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

fn closed_box<'a>(e: &Editor, node: &Node) -> Element<'a, Message> {
    checkbox(closed(&e.desc, node).unwrap_or(false))
        .label("closed")
        .on_toggle(Message::Closed)
        .size(14)
        .into()
}

/// A box's heading, its row label, and its axis when it is part of a vector.
/// `bounds` is the heading for min and max: `Cylinder` or `Cone`.
fn describe(field: Field, bounds: &'static str) -> (&'static str, &'static str, Option<usize>) {
    use Field::*;
    match field {
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
