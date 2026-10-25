//! The help window's page, as data: the menu commands built from the menu
//! itself, and fixed tables for the editor's keys, the viewer and the
//! command line.

use crate::menu::{ALL, key, label, menu_of};
use iced::theme::palette::{Extended, Pair};
use iced::widget::{Column, Row, container, scrollable, text};
use iced::{Alignment, Element, Font, Length, Theme, border, font};
use muda::accelerator::{Code, Modifiers};

/// One piece of a key or mouse action, drawn as a chip.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Chip {
    /// ⌃ ⌥ ⇧ ⌘.
    Modifier(&'static str),
    Key(&'static str),
    Mouse(&'static str),
    /// "then", "or": drawn without a box.
    Word(&'static str),
}

/// A line of the page: its chips, then what they do.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Entry {
    pub chips: Vec<Chip>,
    pub text: &'static str,
}

/// Entries under a small title: a menu, a place, a program.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Group {
    pub title: &'static str,
    pub entries: Vec<Entry>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Section {
    pub title: &'static str,
    pub groups: Vec<Group>,
    /// Drawn in a monospace font: the command line.
    pub monospace: bool,
}

/// A menu key's chips: the modifiers in the Mac's order, then the key.
pub fn chips(modifiers: Modifiers, code: Code) -> Vec<Chip> {
    let mut chips: Vec<Chip> = [
        (Modifiers::CONTROL, "⌃"),
        (Modifiers::ALT, "⌥"),
        (Modifiers::SHIFT, "⇧"),
        (Modifiers::META, "⌘"),
    ]
    .into_iter()
    .filter(|(m, _)| modifiers.contains(*m))
    .map(|(_, s)| Chip::Modifier(s))
    .collect();
    chips.push(Chip::Key(match code {
        Code::KeyK => "K",
        Code::KeyL => "L",
        Code::KeyN => "N",
        Code::KeyO => "O",
        Code::KeyR => "R",
        Code::KeyS => "S",
        Code::KeyZ => "Z",
        Code::Backspace => "⌫",
        _ => "…",
    }));
    chips
}

fn entry(chips: &[Chip], text: &'static str) -> Entry {
    Entry {
        chips: chips.to_vec(),
        text,
    }
}

fn group(title: &'static str, entries: Vec<Entry>) -> Group {
    Group { title, entries }
}

/// The page: menu commands, editor keys, viewer, command line.
pub fn sections() -> Vec<Section> {
    use Chip::{Key as K, Modifier as M, Mouse, Word as W};
    let mut menus: Vec<Group> = Vec::new();
    for c in ALL {
        let e = entry(&key(c).map_or(vec![], |(m, k)| chips(m, k)), label(c));
        match menus.last_mut() {
            Some(g) if g.title == menu_of(c) => g.entries.push(e),
            _ => menus.push(group(menu_of(c), vec![e])),
        }
    }
    // The README's key table, in its order.
    let editor = vec![
        group(
            "Tree",
            vec![
                entry(
                    &[K("↑"), W("/"), K("↓"), W("or"), K("k"), W("/"), K("j")],
                    "select the previous or next row",
                ),
                entry(
                    &[K("←"), W("or"), K("h"), W(","), K("→"), W("or"), K("l")],
                    "fold or unfold the selected item",
                ),
                entry(&[K("Enter")], "shade the item's first box"),
                entry(
                    &[
                        K("a"),
                        W("then"),
                        K("s"),
                        K("c"),
                        K("p"),
                        K("y"),
                        K("o"),
                        K("l"),
                    ],
                    "add a sphere, cube, plane, cylinder, cone or light",
                ),
                entry(&[K("d"), K("d")], "delete the selected shape or light"),
            ],
        ),
        group(
            "Shaded box",
            vec![
                entry(
                    &[K("arrows"), W("or"), K("h"), K("j"), K("k"), K("l")],
                    "move the shading",
                ),
                entry(
                    &[K("Enter"), W(","), K("i")],
                    "edit it, with the cursor at the end or the start",
                ),
                entry(&[K("c")], "empty it and edit it"),
                entry(&[K("x")], "empty it, if the box may be empty"),
                entry(&[K("u")], "put the scene's value back, clearing a red box"),
                entry(&[K("Esc")], "back to the tree"),
            ],
        ),
        group(
            "Editing",
            vec![entry(
                &[K("Enter"), W("or"), K("Esc")],
                "leave the box, shaded; Esc on a red box also puts the scene's value back",
            )],
        ),
        group(
            "Anywhere",
            vec![
                entry(
                    &[K("Tab"), W(","), M("⇧"), K("Tab")],
                    "edit the next or previous box",
                ),
                entry(&[], "the letters act only with no modifier held"),
            ],
        ),
    ];
    // By hand from rtc-metal's src/window.rs; keep it in step.
    let viewer = vec![
        group(
            "Mouse",
            vec![
                entry(&[Mouse("drag")], "orbit the camera"),
                entry(&[Mouse("scroll")], "zoom"),
            ],
        ),
        group(
            "Keys",
            vec![
                entry(&[K("Space")], "return to the scene's view"),
                entry(&[K("i")], "show or hide the status overlay"),
            ],
        ),
        group(
            "Menus",
            vec![
                entry(&[M("⌘"), K("S")], "Save Image…"),
                entry(&[M("⌘"), K("1")], "BVH traversal"),
                entry(&[M("⌘"), K("2")], "Metal traversal"),
                entry(&[M("⌘"), K("]")], "double the samples"),
                entry(&[M("⌘"), K("[")], "halve the samples"),
                entry(&[M("⌘"), K("H")], "hide the viewer"),
            ],
        ),
    ];
    let command_line = vec![
        group(
            "rtc-edit",
            vec![
                entry(&[], "rtc-edit [scene.json]"),
                entry(&[], "RTC_METAL=/path/to/rtc-metal rtc-edit scene.json"),
            ],
        ),
        group(
            "rtc-metal",
            vec![
                entry(
                    &[],
                    "rtc-metal scene.json [-o out.ppm] [--samples N] [--window [--watch]] [--traversal bvh|metal]",
                ),
                entry(&[], "rtc-metal scene.json -o out.ppm"),
                entry(&[], "rtc-metal scene.json --window --samples 1024"),
                entry(&[], "rtc-metal scene.json --window --watch"),
            ],
        ),
    ];
    let section = |title, groups, monospace| Section {
        title,
        groups,
        monospace,
    };
    vec![
        section("Menu commands", menus, false),
        section("Editor keys", editor, false),
        section("Viewer", viewer, false),
        section("Command line", command_line, true),
    ]
}

/// The page, drawn: a scrolling column of sections, groups and entries,
/// with each chip in its kind's colour.
pub fn view<M: 'static>(sections: Vec<Section>) -> Element<'static, M> {
    let bold = Font {
        weight: font::Weight::Bold,
        ..Font::DEFAULT
    };
    let mut page = Column::new().spacing(6).padding(16);
    for section in sections {
        page = page.push(text(section.title).size(18).font(bold));
        for g in section.groups {
            if !g.title.is_empty() {
                page = page.push(text(g.title).size(14).font(bold));
            }
            for e in g.entries {
                let mut line = Row::new().spacing(6).align_y(Alignment::Center);
                for c in e.chips {
                    line = line.push(chip(c));
                }
                let words = if section.monospace {
                    text(e.text).size(13).font(Font::MONOSPACE)
                } else {
                    text(e.text).size(13)
                };
                page = page.push(line.push(words));
            }
        }
    }
    scrollable(page).width(Length::Fill).into()
}

/// A chip: a small rounded box coloured by kind; a word is plain text.
fn chip<M: 'static>(c: Chip) -> Element<'static, M> {
    let (s, pick): (_, fn(&Extended) -> Pair) = match c {
        Chip::Word(s) => return text(s).size(12).into(),
        Chip::Modifier(s) => (s, |p| p.secondary.weak),
        Chip::Key(s) => (s, |p| p.primary.weak),
        Chip::Mouse(s) => (s, |p| p.success.weak),
    };
    container(text(s).size(12))
        .padding([2, 6])
        .style(move |theme: &Theme| {
            let pair = pick(theme.extended_palette());
            container::Style {
                background: Some(pair.color.into()),
                text_color: Some(pair.text),
                border: border::rounded(4),
                ..container::Style::default()
            }
        })
        .into()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_page_draws() {
        let _: iced::Element<'static, ()> = view(sections());
    }

    #[test]
    fn chips_put_modifiers_in_mac_order() {
        use Chip::*;
        let cmd = Modifiers::META;
        let shift = Modifiers::META | Modifiers::SHIFT;
        assert_eq!(
            chips(shift, Code::KeyS),
            [Modifier("⇧"), Modifier("⌘"), Key("S")]
        );
        assert_eq!(chips(cmd, Code::Backspace), [Modifier("⌘"), Key("⌫")]);
        let all = Modifiers::CONTROL | Modifiers::ALT | Modifiers::SHIFT | Modifiers::META;
        assert_eq!(
            chips(all, Code::KeyZ),
            [
                Modifier("⌃"),
                Modifier("⌥"),
                Modifier("⇧"),
                Modifier("⌘"),
                Key("Z")
            ]
        );
    }

    #[test]
    fn every_menu_command_appears_once() {
        let menus = &sections()[0];
        let entries: Vec<&Entry> = menus.groups.iter().flat_map(|g| &g.entries).collect();
        assert_eq!(entries.len(), ALL.len());
        for c in ALL {
            let group = menus
                .groups
                .iter()
                .find(|g| g.title == menu_of(c))
                .unwrap_or_else(|| panic!("{}", menu_of(c)));
            let found: Vec<_> = group
                .entries
                .iter()
                .filter(|e| e.text == label(c))
                .collect();
            assert_eq!(found.len(), 1, "{c:?}");
            let want = key(c).map_or(vec![], |(m, k)| chips(m, k));
            assert_eq!(found[0].chips, want, "{c:?}");
        }
        let titles: Vec<_> = menus.groups.iter().map(|g| g.title).collect();
        assert_eq!(titles, ["File", "Edit", "Scene", "Viewer", "Help"]);
    }

    #[test]
    fn every_menu_key_has_a_name() {
        for c in ALL {
            if let Some((m, k)) = key(c) {
                assert!(!chips(m, k).contains(&Chip::Key("…")), "{c:?}");
            }
        }
    }

    #[test]
    fn sections_come_in_order_and_none_is_empty() {
        let s = sections();
        let titles: Vec<_> = s.iter().map(|s| s.title).collect();
        assert_eq!(
            titles,
            ["Menu commands", "Editor keys", "Viewer", "Command line"]
        );
        for section in &s {
            assert!(!section.groups.is_empty(), "{}", section.title);
            for g in &section.groups {
                assert!(!g.entries.is_empty(), "{} / {}", section.title, g.title);
            }
        }
        assert_eq!(
            s.iter().map(|s| s.monospace).collect::<Vec<_>>(),
            [false, false, false, true]
        );
    }
}
