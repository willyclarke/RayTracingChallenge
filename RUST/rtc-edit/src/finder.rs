//! Scenes opened from Finder: a double-click, Open With or a drop on the
//! Dock icon.
//!
//! AppKit tells the application delegate with `application:openURLs:`,
//! and winit's delegate doesn't answer it. [`install`] adds the method to
//! winit's delegate class; the paths then reach iced through [`opens`].
//! A spike ruled out the other ways: replacing winit's delegate with
//! `setDelegate:` makes winit panic when it looks its own delegate up, and
//! an `odoc` Apple Event handler set in `boot()` is replaced by AppKit's.

use iced::futures::channel::mpsc::{UnboundedReceiver, UnboundedSender, unbounded};
use objc2::runtime::{AnyObject, Imp, Sel};
use objc2::sel;
use objc2_app_kit::NSApplication;
use objc2_foundation::{MainThreadMarker, NSArray, NSURL};
use std::path::PathBuf;
use std::sync::{Mutex, PoisonError};

/// Opens waiting for the editor: AppKit can hand them over before iced's
/// subscriptions run (a cold launch from Finder), so they queue until a
/// receiver takes over.
#[derive(Default)]
struct Inbox {
    queue: Vec<Vec<PathBuf>>,
    sender: Option<UnboundedSender<Vec<PathBuf>>>,
}

impl Inbox {
    /// Send `paths` on, or queue them while nobody is listening.
    fn push(&mut self, paths: Vec<PathBuf>) {
        let paths = match &self.sender {
            Some(tx) => match tx.unbounded_send(paths) {
                Ok(()) => return,
                Err(e) => e.into_inner(),
            },
            None => paths,
        };
        self.sender = None;
        self.queue.push(paths);
    }

    /// Hand the queue to `sender`, in order, and send there from now on.
    fn take_over(&mut self, sender: UnboundedSender<Vec<PathBuf>>) {
        self.sender = Some(sender);
        for paths in std::mem::take(&mut self.queue) {
            self.push(paths);
        }
    }
}

static INBOX: Mutex<Inbox> = Mutex::new(Inbox {
    queue: Vec::new(),
    sender: None,
});

fn inbox() -> std::sync::MutexGuard<'static, Inbox> {
    INBOX.lock().unwrap_or_else(PoisonError::into_inner)
}

/// The files opened from Finder, for a subscription: any that came before
/// it started first, then each open as it comes. A restarted subscription
/// takes over from the old one.
pub fn opens() -> UnboundedReceiver<Vec<PathBuf>> {
    let (tx, rx) = unbounded();
    inbox().take_over(tx);
    rx
}

/// Add `application:openURLs:` to winit's application delegate class, so
/// files opened from Finder reach [`opens`]. When it can't, say so on
/// stderr; the editor works as before, only without Finder opens.
pub fn install(mtm: MainThreadMarker) {
    if let Err(why) = add_open_urls(mtm) {
        eprintln!("rtc-edit: opening from Finder is off: {why}");
    }
}

fn add_open_urls(mtm: MainThreadMarker) -> Result<(), &'static str> {
    let app = NSApplication::sharedApplication(mtm);
    let delegate = app.delegate().ok_or("no application delegate")?;
    let object: &AnyObject = delegate.as_ref();
    let class = object.class();
    type OpenUrls = extern "C-unwind" fn(&AnyObject, Sel, &NSApplication, &NSArray<NSURL>);
    // SAFETY: `open_urls` has the signature the type encoding `v@:@@`
    // below names (void; self, _cmd, then two objects), which is how the
    // runtime calls an `Imp` registered with it.
    let imp = unsafe { std::mem::transmute::<OpenUrls, Imp>(open_urls) };
    // SAFETY: `class` is the class of a live object (winit's delegate),
    // the selector and encoding match `open_urls`, and this runs on the
    // main thread in `boot()`, before AppKit sends any open.
    let added = unsafe {
        objc2::ffi::class_addMethod(
            class as *const _ as *mut _,
            sel!(application:openURLs:),
            imp,
            c"v@:@@".as_ptr(),
        )
    };
    if added.as_bool() {
        Ok(())
    } else {
        Err("application:openURLs: is already defined")
    }
}

/// AppKit's call with the opened files: their paths go to the inbox. It
/// touches nothing of winit's.
extern "C-unwind" fn open_urls(
    _this: &AnyObject,
    _cmd: Sel,
    _app: &NSApplication,
    urls: &NSArray<NSURL>,
) {
    let paths: Vec<PathBuf> = urls
        .iter()
        .filter(|url| url.isFileURL())
        .filter_map(|url| url.path())
        .map(|path| PathBuf::from(path.to_string()))
        .collect();
    if !paths.is_empty() {
        inbox().push(paths);
    }
}

/// What the status line adds when `n` files came at once and only the
/// first opened: ` (2 more not opened)`, or nothing for one file.
pub fn more_not_opened(n: usize) -> String {
    match n {
        0 | 1 => String::new(),
        n => format!(" ({} more not opened)", n - 1),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(s: &str) -> Vec<PathBuf> {
        vec![PathBuf::from(s)]
    }

    #[test]
    fn inbox_hands_queued_paths_over_first() {
        let mut inbox = Inbox::default();
        inbox.push(p("/a.rtc"));
        inbox.push(p("/b.rtc"));
        let (tx, mut rx) = unbounded();
        inbox.take_over(tx);
        inbox.push(p("/c.rtc"));
        let got: Vec<_> = std::iter::from_fn(|| rx.try_recv().ok()).collect();
        assert_eq!(got, [p("/a.rtc"), p("/b.rtc"), p("/c.rtc")]);
        assert!(inbox.queue.is_empty());
    }

    #[test]
    fn inbox_queues_again_when_the_receiver_is_gone() {
        let mut inbox = Inbox::default();
        let (tx, rx) = unbounded();
        inbox.take_over(tx);
        drop(rx);
        inbox.push(p("/a.rtc"));
        let (tx, mut rx) = unbounded();
        inbox.take_over(tx);
        assert_eq!(rx.try_recv().ok(), Some(p("/a.rtc")));
    }

    #[test]
    fn more_not_opened_counts_the_rest() {
        assert_eq!(more_not_opened(0), "");
        assert_eq!(more_not_opened(1), "");
        assert_eq!(more_not_opened(2), " (1 more not opened)");
        assert_eq!(more_not_opened(3), " (2 more not opened)");
    }
}
