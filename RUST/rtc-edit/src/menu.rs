//! The menu bar: its commands, the ids its items carry, and when each can
//! act.

use crate::edit::Item;

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
