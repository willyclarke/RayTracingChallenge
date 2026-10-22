//! The scene editor's load commands: a thread sleeps on the shared file
//! until the editor sends a scene, then hands it to the main thread.

use crate::window;
use dispatch2::DispatchQueue;
use rtc_share::Shared;
use std::path::Path;

/// Take load commands from the shared file at `share`, from the next one
/// sent: each is loaded on the main queue by `window::deliver`.
pub fn listen(share: &Path) -> Result<(), String> {
    let shared = Shared::open(share)?;
    let mut last = shared.load_seq();
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
