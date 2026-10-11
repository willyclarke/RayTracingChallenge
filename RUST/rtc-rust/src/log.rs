use std::fmt;
use std::time::{SystemTime, UNIX_EPOCH};

pub enum Color {
    Reset,
    Red,
    Green,
    Yellow,
    Blue,
    Violet,
    Cyan,
    Gray,
    Gray1,
    White1,
    White2,
}

impl fmt::Display for Color {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let code = match self {
            Color::Reset => "\x1b[0m",
            Color::Red => "\x1b[31m",
            Color::Green => "\x1b[32m",
            Color::Yellow => "\x1b[33m",
            Color::Blue => "\x1b[34m",
            Color::Violet => "\x1b[35m",
            Color::Cyan => "\x1b[36m",
            Color::Gray => "\x1b[90m",
            Color::Gray1 => "\x1b[37m",
            Color::White1 => "\x1b[38m",
            Color::White2 => "\x1b[39m",
        };

        write!(f, "{}", code)
    }
}

pub fn now_ms() -> u128 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_millis()
}

#[macro_export]
macro_rules! __log_internal {
    ($level:expr, $color:expr, $what:expr, $($arg:tt)*) => {
        println!(
            "{}[{}]{} {}[{}]{} {}:{} | {}{}{} | WHAT: {} | WHY: {}",
            $crate::log::Color::Yellow,
            $crate::log::now_ms(),
            $crate::log::Color::Reset,
            $crate::log::Color::Cyan,
            module_path!(),
            $crate::log::Color::Reset,
            file!(),
            line!(),
            $color,
            $level,
            $crate::log::Color::Reset,
            $what,
            format_args!($($arg)*)
        );
    };
}

// Debug log (compiled out in release)
#[macro_export]
macro_rules! logd {
    ($what:expr, $($arg:tt)*) => {
        #[cfg(debug_assertions)]
        {
            $crate::__log_internal!(
                "DEBUG",
                $crate::log::Color::Violet,
                $what,
                $($arg)*
            );
        }
    };
}

// Info log (always enabled)
#[macro_export]
macro_rules! logi {
    ($what:expr, $($arg:tt)*) => {
        $crate::__log_internal!(
            " INFO",
            $crate::log::Color::Green,
            $what,
            $($arg)*
        );
    };
}

// Error log (always enabled)
#[macro_export]
macro_rules! loge {
    ($what:expr, $($arg:tt)*) => {
        $crate::__log_internal!(
            "ERROR",
            $crate::log::Color::Red,
            $what,
            $($arg)*
        );
    };
}
