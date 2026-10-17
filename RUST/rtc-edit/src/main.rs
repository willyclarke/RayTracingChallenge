//! rtc-edit: `rtc-edit <scene.json>` shows the scene as a tree, starts
//! `rtc-metal --window --watch` on a preview copy of it, and saves it
//! (Cmd-S, and Cmd-Shift-S to save elsewhere with asset paths rebased).
//! Cmd-R starts the viewer again after it was closed. Milestone 17 is
//! read-only; editing comes next.

use iced::widget::{Column, button, column, container, row, scrollable, text};
use iced::{Element, Font, Length, Subscription, Task, keyboard};
use rtc_edit::paths::{dir_of, rebase};
use rtc_edit::preview::Preview;
use rtc_edit::tree::{Node, json, rows};
use rtc_edit::viewer::{Viewer, find};
use rtc_rust::scene::SceneDescription;
use std::collections::HashSet;
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
    path: PathBuf,
    desc: SceneDescription,
    folded: HashSet<Node>,
    selected: Option<Node>,
    /// Why there is no preview file, when there is none.
    preview_error: Option<String>,
    status: String,
}

#[derive(Debug, Clone)]
enum Message {
    Select(Node),
    Toggle(Node),
    Save,
    SaveAs,
    ShowViewer,
}

fn main() -> iced::Result {
    iced::application(boot, update, view)
        .title(|e: &Editor| format!("rtc-edit — {}", file_name(&e.path)))
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
    let Some(arg) = std::env::args().nth(1) else {
        eprintln!("usage: rtc-edit <scene.json>");
        std::process::exit(2);
    };
    let path = PathBuf::from(arg);
    let desc = SceneDescription::from_path(&path).unwrap_or_else(|e| {
        eprintln!("rtc-edit: {e}");
        std::process::exit(1);
    });
    // SAFETY: registers a plain function with no arguments; it runs once, at exit.
    unsafe { atexit(clean_up) };
    let preview = path
        .file_name()
        .ok_or_else(|| format!("{} names no file", path.display()))
        .and_then(|name| Preview::new(&std::env::temp_dir(), name))
        .and_then(|p| p.write(&desc, dir_of(&path)).map(|()| p));
    let mut editor = Editor {
        path,
        desc,
        folded: HashSet::new(),
        selected: None,
        preview_error: None,
        status: String::new(),
    };
    match preview {
        Ok(p) => CHILDREN.lock().unwrap().preview = Some(p),
        Err(e) => editor.preview_error = Some(e),
    }
    editor.show_viewer();
    editor
}

impl Editor {
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

    fn save(&mut self) {
        self.status = match self.desc.save(&self.path) {
            Ok(()) => format!("saved {}", file_name(&self.path)),
            Err(e) => e.to_string(),
        };
    }

    fn save_as(&mut self) {
        let Some(new) = rfd::FileDialog::new()
            .set_directory(dir_of(&self.path))
            .set_file_name(file_name(&self.path))
            .add_filter("scene", &["json"])
            .save_file()
        else {
            return;
        };
        let mut desc = self.desc.clone();
        let saved = rebase(&mut desc, dir_of(&self.path), dir_of(&new))
            .and_then(|()| desc.save(&new).map_err(|e| e.to_string()));
        self.status = match saved {
            Ok(()) => {
                self.desc = desc;
                self.path = new;
                format!("saved {}", file_name(&self.path))
            }
            Err(e) => e,
        };
    }
}

fn update(e: &mut Editor, message: Message) -> Task<Message> {
    match message {
        Message::Select(node) => e.selected = Some(node),
        Message::Toggle(node) => {
            if !e.folded.remove(&node) {
                e.folded.insert(node);
            }
        }
        Message::Save => e.save(),
        Message::SaveAs => e.save_as(),
        Message::ShowViewer => e.show_viewer(),
    }
    Task::none()
}

/// Cmd-S, Cmd-Shift-S and Cmd-R.
fn shortcut(event: keyboard::Event) -> Option<Message> {
    let keyboard::Event::KeyPressed { key, modifiers, .. } = event else {
        return None;
    };
    let keyboard::Key::Character(c) = key else {
        return None;
    };
    if !modifiers.command() {
        return None;
    }
    match c.to_lowercase().as_str() {
        "s" if modifiers.shift() => Some(Message::SaveAs),
        "s" => Some(Message::Save),
        "r" => Some(Message::ShowViewer),
        _ => None,
    }
}

fn subscription(_: &Editor) -> Subscription<Message> {
    keyboard::listen().filter_map(shortcut)
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
        Some(node) => text(json(&e.desc, node))
            .font(Font::MONOSPACE)
            .size(13)
            .into(),
        None => text("select an item").into(),
    };
    column![
        row![
            container(scrollable(tree)).width(260).padding(8),
            container(scrollable(inspector).width(Length::Fill))
                .width(Length::Fill)
                .padding(8),
        ]
        .height(Length::Fill),
        container(text(&e.status).size(13)).padding([4, 8]),
    ]
    .into()
}
