//! The scene editor's load commands: a thread sleeps on the shared file
//! until the editor sends a scene, then hands it to the main thread.

use crate::window;
use dispatch2::DispatchQueue;
use rtc_share::Shared;
use std::path::Path;

/// Take load commands from the shared file at `share`, from the next one
/// sent: each is loaded on the main queue by `window::deliver`. The editor
/// holds its loads until this says the viewer is ready.
pub fn listen(share: &Path) -> Result<(), String> {
    let shared = Shared::open(share)?;
    let mut last = shared.load_seq();
    // From here on no load is missed, so the editor may send them.
    shared.viewer_ready();
    std::thread::spawn(move || {
        loop {
            match shared.wait_load(last) {
                Ok(load) => {
                    last = load.seq;
                    DispatchQueue::main().exec_async(move || window::deliver(load));
                }
                Err(e) => {
                    eprintln!("rtc-metal: loads stopped: {e}");
                    return;
                }
            }
        }
    });
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn listening_tells_the_editor_the_viewer_takes_loads() -> Result<(), String> {
        let p = std::env::temp_dir().join(format!("rtc-metal-loads-{}.share", std::process::id()));
        let editor = Shared::create(&p)?;
        assert!(!editor.is_ready());
        listen(&p)?;
        assert!(editor.is_ready());
        std::fs::remove_file(&p).unwrap();
        Ok(())
    }
}
