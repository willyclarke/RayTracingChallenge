//! Which window a key or a close belongs to: the editor's window, or the
//! help window.

use iced::window::{Id, Position};
use iced::{Point, Rectangle, Size};
use rtc_config::{Map, Value};

/// The help window's size when it opens.
pub const HELP_SIZE: Size = Size::new(560.0, 720.0);

/// Where help opens: at the editor's top-right corner.
pub fn help_at(editor: Rectangle) -> Point {
    Point::new(editor.x + editor.width, editor.y)
}

/// Where help opens: beside the editor when it fits on `screen` (the
/// editor's monitor), else centred, so it never opens out of sight.
pub fn help_position(editor: Rectangle, screen: Option<Size>) -> Position {
    match screen {
        Some(screen) if editor.x + editor.width + HELP_SIZE.width <= screen.width => {
            Position::Specific(help_at(editor))
        }
        _ => Position::Centered,
    }
}

/// The smallest and largest saved size taken, in points: a smaller window
/// is hard to find, and a larger one is bigger than any screen and may be
/// more than a Metal surface can hold.
const FRAME_MIN: Size = Size::new(200.0, 150.0);
const FRAME_MAX: f32 = 8192.0;

/// The editor's frame as the settings file keeps it: four finite numbers,
/// the size within `FRAME_MIN`..`FRAME_MAX`; anything else is `None`.
pub fn frame_from(config: &Map<String, Value>) -> Option<Rectangle> {
    let frame = config.get(rtc_config::EDITOR_FRAME)?;
    let n = |k| {
        frame
            .get(k)
            .and_then(Value::as_f64)
            .map(|v| v as f32)
            .filter(|v| v.is_finite())
    };
    let (x, y, width, height) = (n("x")?, n("y")?, n("width")?, n("height")?);
    ((FRAME_MIN.width..=FRAME_MAX).contains(&width)
        && (FRAME_MIN.height..=FRAME_MAX).contains(&height))
    .then(|| Rectangle::new(Point::new(x, y), Size::new(width, height)))
}

/// `frame` as the settings file keeps it.
pub fn frame_json(frame: Rectangle) -> Value {
    rtc_config::json!({
        "x": frame.x, "y": frame.y, "width": frame.width, "height": frame.height
    })
}

/// Whether `frame` lies wholly on `screen` (the window's monitor, its
/// top-left at the origin); no screen, no fit.
pub fn fits(frame: Rectangle, screen: Option<Size>) -> bool {
    screen.is_some_and(|s| {
        frame.x >= 0.0
            && frame.y >= 0.0
            && frame.x + frame.width <= s.width
            && frame.y + frame.height <= s.height
    })
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
    fn help_is_centred_when_there_is_no_room_beside_the_editor() {
        let frame = Rectangle::new(Point::new(0.0, 0.0), Size::new(900.0, 640.0));
        let wide = Some(Size::new(1920.0, 1080.0));
        assert!(matches!(
            help_position(frame, wide),
            Position::Specific(p) if p == Point::new(900.0, 0.0)
        ));
        // 900 + 560 > 1440: help would poke off the right edge.
        let narrow = Some(Size::new(1440.0, 900.0));
        assert!(matches!(help_position(frame, narrow), Position::Centered));
        let full = Rectangle::new(Point::new(0.0, 0.0), Size::new(1920.0, 1050.0));
        assert!(matches!(help_position(full, wide), Position::Centered));
        assert!(matches!(help_position(frame, None), Position::Centered));
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

    #[test]
    fn a_saved_frame_needs_four_numbers_and_a_size() {
        use rtc_config::{Map, Value, json};
        let map = |v: Value| Map::from_iter([("editor_frame".to_string(), v)]);
        let full = json!({"x": 10, "y": 25, "width": 900, "height": 640});
        let frame = Rectangle::new(Point::new(10.0, 25.0), Size::new(900.0, 640.0));
        assert_eq!(frame_from(&map(full)), Some(frame));
        assert_eq!(frame_from(&map(frame_json(frame))), Some(frame));
        assert_eq!(frame_from(&Map::new()), None);
        for bad in [
            json!({"x": 0, "y": 0, "width": 900}),
            json!({"x": "0", "y": 0, "width": 900, "height": 640}),
            json!({"x": 0, "y": 0, "width": 0, "height": 640}),
            json!({"x": 0, "y": 0, "width": -5, "height": 640}),
            json!([1, 2, 3, 4]),
            // Too big for any screen or surface, or too small to find.
            json!({"x": 0, "y": 0, "width": 1e39, "height": 640}),
            json!({"x": 0, "y": 0, "width": 20000, "height": 640}),
            json!({"x": 0, "y": 0, "width": 1, "height": 1}),
            json!({"x": 1e39, "y": 0, "width": 900, "height": 640}),
        ] {
            assert_eq!(frame_from(&map(bad.clone())), None, "{bad}");
        }
    }

    #[test]
    fn a_frame_fits_only_inside_its_screen() {
        let screen = Some(Size::new(1440.0, 900.0));
        let at = |x, y, w, h| Rectangle::new(Point::new(x, y), Size::new(w, h));
        assert!(fits(at(0.0, 25.0, 900.0, 640.0), screen));
        // Touching the right and bottom edges.
        assert!(fits(at(540.0, 260.0, 900.0, 640.0), screen));
        assert!(!fits(at(600.0, 25.0, 900.0, 640.0), screen));
        assert!(!fits(at(0.0, 300.0, 900.0, 640.0), screen));
        assert!(!fits(at(-1.0, 25.0, 900.0, 640.0), screen));
        assert!(!fits(at(0.0, 25.0, 900.0, 640.0), None));
    }
}
