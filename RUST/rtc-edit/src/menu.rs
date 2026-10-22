//! The menu bar: its commands, the ids its items carry, and when each can
//! act.

use crate::edit::Item;
use muda::accelerator::{Accelerator, Code, Modifiers};
use muda::{AboutMetadata, CheckMenuItem, IsMenuItem, Menu, MenuItem, PredefinedMenuItem, Submenu};

/// What a menu item asks the editor to do.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Command {
    New,
    Open,
    Save,
    SaveAs,
    Undo,
    Redo,
    Delete,
    Add(Item),
    Load,
    StartViewer,
    Follow,
    Keep,
}

/// Every command, in menu order.
pub const ALL: [Command; 17] = [
    Command::New,
    Command::Open,
    Command::Save,
    Command::SaveAs,
    Command::Undo,
    Command::Redo,
    Command::Delete,
    Command::Add(Item::Sphere),
    Command::Add(Item::Cube),
    Command::Add(Item::Plane),
    Command::Add(Item::Cylinder),
    Command::Add(Item::Cone),
    Command::Add(Item::Light),
    Command::Load,
    Command::StartViewer,
    Command::Follow,
    Command::Keep,
];

/// The id the command's menu item carries.
pub fn id(command: Command) -> &'static str {
    match command {
        Command::New => "new",
        Command::Open => "open",
        Command::Save => "save",
        Command::SaveAs => "save-as",
        Command::Undo => "undo",
        Command::Redo => "redo",
        Command::Delete => "delete",
        Command::Add(Item::Sphere) => "add-sphere",
        Command::Add(Item::Cube) => "add-cube",
        Command::Add(Item::Plane) => "add-plane",
        Command::Add(Item::Cylinder) => "add-cylinder",
        Command::Add(Item::Cone) => "add-cone",
        Command::Add(Item::Light) => "add-light",
        Command::Load => "load",
        Command::StartViewer => "start-viewer",
        Command::Follow => "follow",
        Command::Keep => "keep",
    }
}

/// The command a menu item's id names.
pub fn command(id: &str) -> Option<Command> {
    ALL.into_iter().find(|&c| self::id(c) == id)
}

/// What the greying needs to know about the editor.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct State {
    pub can_undo: bool,
    pub can_redo: bool,
    /// An item is selected, it may be deleted, and no box is shaded or edited.
    pub can_delete: bool,
    pub viewer_running: bool,
    pub following: bool,
}

/// Whether `command`'s item can act, or is greyed.
pub fn enabled(command: Command, state: &State) -> bool {
    match command {
        Command::Undo => state.can_undo,
        Command::Redo => state.can_redo,
        Command::Delete => state.can_delete,
        Command::Load | Command::Follow => state.viewer_running,
        Command::StartViewer => !state.viewer_running,
        Command::Keep => state.following,
        Command::New | Command::Open | Command::Save | Command::SaveAs | Command::Add(_) => true,
    }
}

/// The menu bar, with the items whose greying or tick changes.
pub struct Menus {
    menu: Menu,
    /// Every item but Follow the Viewer, with its command.
    items: Vec<(Command, MenuItem)>,
    follow: CheckMenuItem,
    /// The state the items were last set for.
    last: Option<State>,
}

impl Default for Menus {
    fn default() -> Menus {
        Menus::new()
    }
}

impl Menus {
    /// rtc-edit, File, Edit, Scene and Viewer.
    pub fn new() -> Menus {
        let mut items = Vec::new();
        let mut item = |c: Command, label: &str, key: Option<(Modifiers, Code)>| {
            let accel = key.map(|(m, k)| Accelerator::new(m, k));
            let i = MenuItem::with_id(id(c), label, true, accel);
            items.push((c, i.clone()));
            i
        };
        let cmd = Modifiers::META;
        let shift = Modifiers::META | Modifiers::SHIFT;
        let file = [
            item(Command::New, "New", Some((cmd, Code::KeyN))),
            item(Command::Open, "Open…", Some((cmd, Code::KeyO))),
            item(Command::Save, "Save", Some((cmd, Code::KeyS))),
            item(Command::SaveAs, "Save As…", Some((shift, Code::KeyS))),
        ];
        let edit = [
            item(Command::Undo, "Undo", Some((cmd, Code::KeyZ))),
            item(Command::Redo, "Redo", Some((shift, Code::KeyZ))),
            item(Command::Delete, "Delete", Some((cmd, Code::Backspace))),
        ];
        let scene = [
            item(Command::Add(Item::Sphere), "Add Sphere", None),
            item(Command::Add(Item::Cube), "Add Cube", None),
            item(Command::Add(Item::Plane), "Add Plane", None),
            item(Command::Add(Item::Cylinder), "Add Cylinder", None),
            item(Command::Add(Item::Cone), "Add Cone", None),
            item(Command::Add(Item::Light), "Add Light", None),
        ];
        let load = item(Command::Load, "Load Edits", Some((cmd, Code::KeyL)));
        let start = item(
            Command::StartViewer,
            "Start Viewer",
            Some((cmd, Code::KeyR)),
        );
        let keep = item(Command::Keep, "Keep Followed View", Some((cmd, Code::KeyK)));
        let follow =
            CheckMenuItem::with_id(id(Command::Follow), "Follow the Viewer", true, false, None);

        let about = PredefinedMenuItem::about(
            Some("About rtc-edit"),
            Some(AboutMetadata {
                name: Some("rtc-edit".into()),
                version: Some(env!("CARGO_PKG_VERSION").into()),
                ..Default::default()
            }),
        );
        let app = Submenu::with_items(
            "rtc-edit",
            true,
            &[
                &about,
                &PredefinedMenuItem::separator(),
                &PredefinedMenuItem::quit(None),
            ],
        )
        .expect("rtc-edit menu");
        let submenu = |name: &str, items: &[MenuItem]| {
            let items: Vec<&dyn IsMenuItem> = items.iter().map(|i| i as &dyn IsMenuItem).collect();
            Submenu::with_items(name, true, &items).expect(name)
        };
        let viewer = Submenu::with_items("Viewer", true, &[&load, &start, &follow, &keep])
            .expect("Viewer menu");
        let menu = Menu::new();
        menu.append_items(&[
            &app,
            &submenu("File", &file),
            &submenu("Edit", &edit),
            &submenu("Scene", &scene),
            &viewer,
        ])
        .expect("menu bar");
        Menus {
            menu,
            items,
            follow,
            last: None,
        }
    }

    /// Make it the app's menu bar; once, when the window has opened.
    pub fn install(&self) {
        self.menu.init_for_nsapp();
    }

    /// Grey the items that can't act in `state`, and tick Follow the Viewer
    /// while following.
    pub fn refresh(&mut self, state: &State) {
        if self.last != Some(*state) {
            for (c, i) in &self.items {
                i.set_enabled(enabled(*c, state));
            }
            self.follow.set_enabled(enabled(Command::Follow, state));
            self.last = Some(*state);
        }
        // A click flips the tick itself, even when following doesn't start.
        if self.follow.is_checked() != state.following {
            self.follow.set_checked(state.following);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_command_id_maps_back() {
        let mut ids = std::collections::HashSet::new();
        for c in ALL {
            assert_eq!(command(id(c)), Some(c));
            assert!(ids.insert(id(c)), "{c:?}");
        }
        assert_eq!(command("nonsense"), None);
    }

    #[test]
    fn undo_and_redo_need_a_step() {
        let none = State::default();
        assert!(!enabled(Command::Undo, &none) && !enabled(Command::Redo, &none));
        let undo = State {
            can_undo: true,
            ..none
        };
        assert!(enabled(Command::Undo, &undo) && !enabled(Command::Redo, &undo));
        let redo = State {
            can_redo: true,
            ..none
        };
        assert!(!enabled(Command::Undo, &redo) && enabled(Command::Redo, &redo));
    }

    #[test]
    fn delete_needs_a_deletable_item() {
        assert!(!enabled(Command::Delete, &State::default()));
        let state = State {
            can_delete: true,
            ..State::default()
        };
        assert!(enabled(Command::Delete, &state));
    }

    #[test]
    fn viewer_items_follow_the_viewer() {
        for running in [false, true] {
            let state = State {
                viewer_running: running,
                ..State::default()
            };
            assert_eq!(enabled(Command::Load, &state), running);
            assert_eq!(enabled(Command::Follow, &state), running);
            assert_eq!(enabled(Command::StartViewer, &state), !running);
        }
    }

    #[test]
    fn keep_needs_following() {
        for following in [false, true] {
            let state = State {
                following,
                ..State::default()
            };
            assert_eq!(enabled(Command::Keep, &state), following);
        }
    }

    #[test]
    fn the_rest_are_always_enabled() {
        let state = State::default();
        for c in ALL {
            if matches!(
                c,
                Command::New | Command::Open | Command::Save | Command::SaveAs | Command::Add(_)
            ) {
                assert!(enabled(c, &state), "{c:?}");
            }
        }
    }
}
