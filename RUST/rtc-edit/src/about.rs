//! What About rtc-edit and Help ▸ Report an Issue… show: the build, the
//! credits, the icon and the GitHub issue to open.

/// The commit rtc-edit was built from and its date, such as
/// `7f1756a 2026-10-25`, with ` dirty` for uncommitted changes; `unknown`
/// without git.
pub fn version() -> &'static str {
    env!("RTC_EDIT_COMMIT")
}

pub const COPYRIGHT: &str = "© 2026 Willy Clarke";

pub const CREDITS: &str = "Built while working through The Ray Tracer Challenge by Jamis Buck.\n\
Made with Rust, iced, muda and Metal.\n\
github.com/willyclarke/RayTracingChallenge";

/// `s` for a URL query: letters, digits and `-_.~` stay, every other UTF-8
/// byte becomes `%XX`.
pub fn percent_encode(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for b in s.bytes() {
        if b.is_ascii_alphanumeric() || b"-_.~".contains(&b) {
            out.push(b as char);
        } else {
            out.push_str(&format!("%{b:02X}"));
        }
    }
    out
}

/// A new GitHub issue for rtc-edit `version` (`0.1.0 (7f1756a …)`) on macOS
/// `macos`, its body an outline to fill in.
pub fn issue_url(version: &str, macos: &str) -> String {
    let body = format!(
        "**What happened**\n\n**What I expected**\n\n**Steps to reproduce**\n\n---\nrtc-edit {version}, macOS {macos}\n"
    );
    format!(
        "https://github.com/willyclarke/RayTracingChallenge/issues/new?title={}&body={}",
        percent_encode("rtc-edit: "),
        percent_encode(&body)
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn percent_encode_keeps_unreserved_and_encodes_the_rest() {
        assert_eq!(percent_encode("a b*\n:?"), "a%20b%2A%0A%3A%3F");
        assert_eq!(percent_encode("AZaz09-_.~"), "AZaz09-_.~");
        assert_eq!(percent_encode("é"), "%C3%A9");
    }

    #[test]
    fn the_issue_url_has_the_title_and_the_versions() {
        let url = issue_url("0.1.0 (abc 2026-10-25)", "26.0");
        let prefix = "https://github.com/willyclarke/RayTracingChallenge/issues/new?title=rtc-edit%3A%20&body=";
        assert!(url.starts_with(prefix), "{url}");
        let body = &url[prefix.len()..];
        assert!(
            body.starts_with(&percent_encode("**What happened**\n\n")),
            "{body}"
        );
        assert!(
            body.ends_with(&percent_encode(
                "---\nrtc-edit 0.1.0 (abc 2026-10-25), macOS 26.0\n"
            )),
            "{body}"
        );
    }

    #[test]
    fn the_version_is_one_line() {
        assert!(!version().is_empty());
        assert!(!version().contains('\n'));
    }
}
