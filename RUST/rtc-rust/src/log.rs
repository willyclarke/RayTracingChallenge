use std::time::{SystemTime, UNIX_EPOCH};

#[allow(dead_code)]
pub(crate) const RESET: &str = "\x1b[0m";
#[allow(dead_code)]
pub(crate) const RED: &str = "\x1b[31m";
#[allow(dead_code)]
pub(crate) const GREEN: &str = "\x1b[32m";
#[allow(dead_code)]
pub(crate) const YELLOW: &str = "\x1b[33m";
#[allow(dead_code)]
pub(crate) const BLUE: &str = "\x1b[34m";
#[allow(dead_code)]
pub(crate) const VIOLET: &str = "\x1b[35m";
#[allow(dead_code)]
pub(crate) const GRAY1: &str = "\x1b[37m";
#[allow(dead_code)]
pub(crate) const WHITE1: &str = "\x1b[38m";
#[allow(dead_code)]
pub(crate) const WHITE2: &str = "\x1b[39m";
#[allow(dead_code)]
pub(crate) const CYAN: &str = "\x1b[36m";
#[allow(dead_code)]
pub(crate) const GRAY: &str = "\x1b[90m";

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
            $crate::log::YELLOW,
            $crate::log::now_ms(),
            $crate::log::RESET,
            $crate::log::CYAN,
            module_path!(),
            $crate::log::RESET,
            file!(),
            line!(),
            $color,
            $level,
            $crate::log::RESET,
            $what,
            format!($($arg)*)
        );
    };
}

// Debug log (compiled out in release)
#[macro_export]
macro_rules! logd {
    ($what:expr, $($arg:tt)*) => {
        #[cfg(debug_assertions)]
        {
            $crate::__log_internal!("DEBUG", $crate::log::VIOLET, $what, $($arg)*);
        }
    };
}

// Info log (always enabled)
#[macro_export]
macro_rules! logi {
    ($what:expr, $($arg:tt)*) => {
            $crate::__log_internal!(" INFO", $crate::log::GREEN, $what, $($arg)*);

    };
}

// Error log (always enabled)
#[macro_export]
macro_rules! loge {
    ($what:expr, $($arg:tt)*) => {
        $crate::__log_internal!("ERROR", $crate::log::RED, $what, $($arg)*);
    };
}

#[allow(unused_imports)]
pub(crate) use __log_internal;
#[allow(unused_imports)]
pub(crate) use logd;
#[allow(unused_imports)]
pub(crate) use loge;
#[allow(unused_imports)]
pub(crate) use logi;
