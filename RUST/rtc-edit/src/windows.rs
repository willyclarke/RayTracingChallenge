//! Which window a key or a close belongs to: the editor's window, or the
//! help window.

use iced::window::Id;
use iced::{Point, Rectangle, Size};

/// The help window's size when it opens.
pub const HELP_SIZE: Size = Size::new(560.0, 720.0);

/// Where help opens: at the editor's top-right corner.
pub fn help_at(editor: Rectangle) -> Point {
    Point::new(editor.x + editor.width, editor.y)
}

/// What kind of input a key or click is, for routing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Input {
    Escape,
    /// Cmd pressed or released: the editor tracks it from any window.
    Modifiers,
    Other,
}

/// Where an input goes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Route {
    Editor,
    CloseHelp,
    Drop,
}

/// Inputs in the editor's window act there; in help, Esc closes it, Cmd is
/// tracked, and the rest are ignored.
pub fn route(window: Id, main: Id, help: Option<Id>, input: Input) -> Route {
    if window == main || input == Input::Modifiers && help == Some(window) {
        Route::Editor
    } else if help == Some(window) && input == Input::Escape {
        Route::CloseHelp
    } else {
        Route::Drop
    }
}

/// What a window's closing means.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Closed {
    /// The editor closed: the app quits.
    Quit,
    /// Help closed: forget it.
    Help,
    Other,
}

pub fn closed(window: Id, main: Id, help: Option<Id>) -> Closed {
    if window == main {
        Closed::Quit
    } else if help == Some(window) {
        Closed::Help
    } else {
        Closed::Other
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn help_opens_at_the_editors_top_right() {
        let frame = Rectangle::new(Point::new(10.0, 20.0), Size::new(900.0, 640.0));
        assert_eq!(help_at(frame), Point::new(910.0, 20.0));
    }

    #[test]
    fn keys_act_in_the_editor_and_esc_closes_help() {
        let (main, help, other) = (Id::unique(), Id::unique(), Id::unique());
        for input in [Input::Escape, Input::Modifiers, Input::Other] {
            assert_eq!(route(main, main, Some(help), input), Route::Editor);
            assert_eq!(route(main, main, None, input), Route::Editor);
        }
        assert_eq!(
            route(help, main, Some(help), Input::Escape),
            Route::CloseHelp
        );
        assert_eq!(route(help, main, Some(help), Input::Other), Route::Drop);
        assert_eq!(
            route(help, main, Some(help), Input::Modifiers),
            Route::Editor
        );
        assert_eq!(route(other, main, Some(help), Input::Escape), Route::Drop);
    }

    #[test]
    fn closing_the_editor_quits_and_closing_help_forgets_it() {
        let (main, help, other) = (Id::unique(), Id::unique(), Id::unique());
        assert_eq!(closed(main, main, Some(help)), Closed::Quit);
        assert_eq!(closed(help, main, Some(help)), Closed::Help);
        assert_eq!(closed(other, main, Some(help)), Closed::Other);
        assert_eq!(closed(other, main, None), Closed::Other);
    }
}
