//! The menu bar: its commands, the ids its items carry, and when each can
//! act.

use crate::about;
use crate::edit::Item;
use crate::recent;
use muda::accelerator::{Accelerator, Code, Modifiers};
use muda::{
    AboutMetadata, CheckMenuItem, Icon, IsMenuItem, Menu, MenuItem, PredefinedMenuItem, Submenu,
};

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
    Help,
    ReportIssue,
    /// File ▸ Open Recent's entry `n`, newest first.
    OpenRecent(usize),
    /// File ▸ Open Recent ▸ Clear Menu.
    ClearRecent,
}

/// The ids of Open Recent's entries, one per place on the list.
const RECENT_IDS: [&str; recent::MAX] = [
    "recent-0", "recent-1", "recent-2", "recent-3", "recent-4", "recent-5", "recent-6", "recent-7",
    "recent-8", "recent-9",
];

/// Every command with a fixed item, in menu order; Open Recent's entries and
/// Clear Menu come and go with the list.
pub const ALL: [Command; 19] = [
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
    Command::Help,
    Command::ReportIssue,
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
        Command::Help => "help",
        Command::ReportIssue => "report-issue",
        Command::OpenRecent(n) => RECENT_IDS[n],
        Command::ClearRecent => "clear-recent",
    }
}

/// The command's menu item text.
pub fn label(command: Command) -> &'static str {
    match command {
        Command::New => "New",
        Command::Open => "Open…",
        Command::Save => "Save",
        Command::SaveAs => "Save As…",
        Command::Undo => "Undo",
        Command::Redo => "Redo",
        Command::Delete => "Delete",
        Command::Add(Item::Sphere) => "Add Sphere",
        Command::Add(Item::Cube) => "Add Cube",
        Command::Add(Item::Plane) => "Add Plane",
        Command::Add(Item::Cylinder) => "Add Cylinder",
        Command::Add(Item::Cone) => "Add Cone",
        Command::Add(Item::Light) => "Add Light",
        Command::Load => "Load Edits",
        Command::StartViewer => "Start Viewer",
        Command::Follow => "Follow the Viewer",
        Command::Keep => "Keep Followed View",
        Command::Help => "RTC Edit Help",
        Command::ReportIssue => "Report an Issue…",
        // The submenu's title; its entries are labelled by `recent::labels`.
        Command::OpenRecent(_) => "Open Recent",
        Command::ClearRecent => "Clear Menu",
    }
}

/// The command's key, if it has one.
pub fn key(command: Command) -> Option<(Modifiers, Code)> {
    let cmd = Modifiers::META;
    let shift = Modifiers::META | Modifiers::SHIFT;
    match command {
        Command::New => Some((cmd, Code::KeyN)),
        Command::Open => Some((cmd, Code::KeyO)),
        Command::Save => Some((cmd, Code::KeyS)),
        Command::SaveAs => Some((shift, Code::KeyS)),
        Command::Undo => Some((cmd, Code::KeyZ)),
        Command::Redo => Some((shift, Code::KeyZ)),
        Command::Delete => Some((cmd, Code::Backspace)),
        Command::Load => Some((cmd, Code::KeyL)),
        Command::StartViewer => Some((cmd, Code::KeyR)),
        Command::Keep => Some((cmd, Code::KeyK)),
        // macOS keeps Cmd-? for its own Help-menu search, which finds this
        // item; a key of ours would never arrive.
        Command::Add(_)
        | Command::Follow
        | Command::Help
        | Command::ReportIssue
        | Command::OpenRecent(_)
        | Command::ClearRecent => None,
    }
}

/// The menu the command's item is in.
pub fn menu_of(command: Command) -> &'static str {
    match command {
        Command::New
        | Command::Open
        | Command::Save
        | Command::SaveAs
        | Command::OpenRecent(_)
        | Command::ClearRecent => "File",
        Command::Undo | Command::Redo | Command::Delete => "Edit",
        Command::Add(_) => "Scene",
        Command::Load | Command::StartViewer | Command::Follow | Command::Keep => "Viewer",
        Command::Help | Command::ReportIssue => "Help",
    }
}

/// The command a menu item's id names.
pub fn command(id: &str) -> Option<Command> {
    ALL.into_iter()
        .chain([Command::ClearRecent])
        .find(|&c| self::id(c) == id)
        .or_else(|| {
            RECENT_IDS
                .iter()
                .position(|&r| r == id)
                .map(Command::OpenRecent)
        })
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
        Command::New
        | Command::Open
        | Command::Save
        | Command::SaveAs
        | Command::Add(_)
        | Command::Help
        | Command::ReportIssue
        | Command::OpenRecent(_)
        | Command::ClearRecent => true,
    }
}

/// What About RTC Edit shows: the icon, the build and the credits.
pub fn about_metadata() -> AboutMetadata {
    AboutMetadata {
        name: Some(about::NAME.into()),
        // AppKit shows "Version <version> (<short_version>)".
        version: Some(env!("CARGO_PKG_VERSION").into()),
        short_version: Some(about::version().into()),
        copyright: Some(about::COPYRIGHT.into()),
        credits: Some(about::CREDITS.into()),
        icon: about::icon(about::ICON_PPM)
            .and_then(|(rgba, w, h)| Icon::from_rgba(rgba, w, h).ok()),
        ..Default::default()
    }
}

/// The menu bar, with the items whose greying or tick changes.
pub struct Menus {
    menu: Menu,
    /// The Help menu, which macOS gives a search box.
    help: Submenu,
    /// Every item but Follow the Viewer, with its command.
    items: Vec<(Command, MenuItem)>,
    follow: CheckMenuItem,
    /// File ▸ Open Recent, rebuilt by [`Menus::set_recent`].
    recent: Submenu,
    /// The state the items were last set for.
    last: Option<State>,
}

impl Default for Menus {
    fn default() -> Menus {
        Menus::new()
    }
}

impl Menus {
    /// rtc-edit, File, Edit, Scene, Viewer and Help.
    pub fn new() -> Menus {
        let mut items = Vec::new();
        let mut item = |c: Command| {
            let accel = key(c).map(|(m, k)| Accelerator::new(m, k));
            let i = MenuItem::with_id(id(c), label(c), true, accel);
            items.push((c, i.clone()));
            i
        };
        let new = item(Command::New);
        let open = item(Command::Open);
        let save = item(Command::Save);
        let save_as = item(Command::SaveAs);
        let edit = [
            item(Command::Undo),
            item(Command::Redo),
            item(Command::Delete),
        ];
        let scene = [
            item(Command::Add(Item::Sphere)),
            item(Command::Add(Item::Cube)),
            item(Command::Add(Item::Plane)),
            item(Command::Add(Item::Cylinder)),
            item(Command::Add(Item::Cone)),
            item(Command::Add(Item::Light)),
        ];
        let load = item(Command::Load);
        let start = item(Command::StartViewer);
        let keep = item(Command::Keep);
        let help = [item(Command::Help), item(Command::ReportIssue)];
        let follow = CheckMenuItem::with_id(
            id(Command::Follow),
            label(Command::Follow),
            true,
            false,
            None,
        );

        let about = PredefinedMenuItem::about(
            Some(&format!("About {}", about::NAME)),
            Some(about_metadata()),
        );
        let app = Submenu::with_items(
            about::NAME,
            true,
            &[
                &about,
                &PredefinedMenuItem::separator(),
                &PredefinedMenuItem::services(None),
                &PredefinedMenuItem::separator(),
                &PredefinedMenuItem::hide(Some(&format!("Hide {}", about::NAME))),
                &PredefinedMenuItem::hide_others(None),
                &PredefinedMenuItem::show_all(None),
                &PredefinedMenuItem::separator(),
                &PredefinedMenuItem::quit(Some(&format!("Quit {}", about::NAME))),
            ],
        )
        .expect("rtc-edit menu");
        let submenu = |name: &str, items: &[MenuItem]| {
            let items: Vec<&dyn IsMenuItem> = items.iter().map(|i| i as &dyn IsMenuItem).collect();
            Submenu::with_items(name, true, &items).expect(name)
        };
        let recent = Submenu::new(label(Command::OpenRecent(0)), true);
        let file = Submenu::with_items("File", true, &[&new, &open, &recent, &save, &save_as])
            .expect("File menu");
        let viewer = Submenu::with_items("Viewer", true, &[&load, &start, &follow, &keep])
            .expect("Viewer menu");
        let help = submenu("Help", &help);
        let menu = Menu::new();
        menu.append_items(&[
            &app,
            &file,
            &submenu("Edit", &edit),
            &submenu("Scene", &scene),
            &viewer,
            &help,
        ])
        .expect("menu bar");
        let mut menus = Menus {
            menu,
            help,
            items,
            follow,
            recent,
            last: None,
        };
        menus.set_recent(&[]);
        menus
    }

    /// Fill File ▸ Open Recent: one entry per label, newest first, then
    /// Clear Menu, greyed when there are none.
    pub fn set_recent(&mut self, labels: &[String]) {
        while self.recent.remove_at(0).is_some() {}
        for (id, text) in RECENT_IDS.iter().zip(labels) {
            let _ = self
                .recent
                .append(&MenuItem::with_id(*id, text, true, None));
        }
        if !labels.is_empty() {
            let _ = self.recent.append(&PredefinedMenuItem::separator());
        }
        let clear = Command::ClearRecent;
        let _ = self.recent.append(&MenuItem::with_id(
            id(clear),
            label(clear),
            !labels.is_empty(),
            None,
        ));
    }

    /// Make it the app's menu bar; once, when the window has opened.
    pub fn install(&self) {
        self.menu.init_for_nsapp();
        self.help.set_as_help_menu_for_nsapp();
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
    fn labels_and_keys_are_the_menus() {
        let cmd = Modifiers::META;
        let shift = Modifiers::META | Modifiers::SHIFT;
        assert_eq!(label(Command::New), "New");
        assert_eq!(key(Command::New), Some((cmd, Code::KeyN)));
        assert_eq!(label(Command::Open), "Open…");
        assert_eq!(key(Command::SaveAs), Some((shift, Code::KeyS)));
        assert_eq!(key(Command::Redo), Some((shift, Code::KeyZ)));
        assert_eq!(key(Command::Delete), Some((cmd, Code::Backspace)));
        assert_eq!(label(Command::Add(Item::Cone)), "Add Cone");
        assert_eq!(key(Command::Add(Item::Cone)), None);
        assert_eq!(label(Command::Load), "Load Edits");
        assert_eq!(key(Command::StartViewer), Some((cmd, Code::KeyR)));
        assert_eq!(label(Command::Follow), "Follow the Viewer");
        assert_eq!(key(Command::Follow), None);
        assert_eq!(key(Command::Keep), Some((cmd, Code::KeyK)));
        assert_eq!(label(Command::Help), "RTC Edit Help");
        // macOS keeps Cmd-? for its own Help-menu search.
        assert_eq!(key(Command::Help), None);
    }

    #[test]
    fn every_command_has_a_label_and_a_menu() {
        for c in ALL {
            assert!(!label(c).is_empty(), "{c:?}");
            assert!(
                ["File", "Edit", "Scene", "Viewer", "Help"].contains(&menu_of(c)),
                "{c:?}"
            );
        }
        assert_eq!(menu_of(Command::SaveAs), "File");
        assert_eq!(menu_of(Command::Delete), "Edit");
        assert_eq!(menu_of(Command::Add(Item::Light)), "Scene");
        assert_eq!(menu_of(Command::Keep), "Viewer");
        assert_eq!(menu_of(Command::Help), "Help");
        assert!(enabled(Command::Help, &State::default()));
    }

    #[test]
    fn report_an_issue_is_in_help_and_always_enabled() {
        assert_eq!(label(Command::ReportIssue), "Report an Issue…");
        assert_eq!(menu_of(Command::ReportIssue), "Help");
        assert_eq!(key(Command::ReportIssue), None);
        assert_eq!(id(Command::ReportIssue), "report-issue");
        assert!(enabled(Command::ReportIssue, &State::default()));
        assert_eq!(ALL.last(), Some(&Command::ReportIssue));
    }

    #[test]
    fn about_shows_the_crate_version_then_the_build() {
        // AppKit shows `version` as "Version X" and `short_version` in
        // parentheses after it.
        let about = about_metadata();
        assert_eq!(about.name.as_deref(), Some("RTC Edit"));
        assert_eq!(about.version.as_deref(), Some(env!("CARGO_PKG_VERSION")));
        assert_eq!(about.short_version.as_deref(), Some(about::version()));
        assert_eq!(about.credits.as_deref(), Some(about::CREDITS));
        assert_eq!(about.copyright.as_deref(), Some(about::COPYRIGHT));
        assert!(about.icon.is_some());
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

    #[test]
    fn recent_ids_map_back_only_in_range() {
        assert_eq!(command("recent-0"), Some(Command::OpenRecent(0)));
        assert_eq!(command("recent-9"), Some(Command::OpenRecent(9)));
        assert_eq!(command("clear-recent"), Some(Command::ClearRecent));
        for bad in ["recent-10", "recent-", "recent-x", "recent--1"] {
            assert_eq!(command(bad), None, "{bad}");
        }
        assert_eq!(id(Command::OpenRecent(3)), "recent-3");
    }

    #[test]
    fn recent_commands_are_file_items_outside_all() {
        for c in [Command::OpenRecent(0), Command::ClearRecent] {
            assert!(!ALL.contains(&c), "{c:?}");
            assert_eq!(menu_of(c), "File");
            assert_eq!(key(c), None);
            assert!(enabled(c, &State::default()), "{c:?}");
        }
        assert_eq!(label(Command::ClearRecent), "Clear Menu");
    }
}
