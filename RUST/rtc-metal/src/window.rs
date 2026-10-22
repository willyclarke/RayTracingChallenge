//! Live window (`--window`): an `MTKView` that adds a batch of sampled
//! frames per display refresh until `target` frames, showing the running
//! mean; the batch grows while the GPU keeps up with the display and halves
//! when it lags, so the GPU, not the refresh rate, sets the pace. Dragging
//! orbits the camera and scrolling zooms; both restart the accumulation
//! (events with no movement leave it alone). Space returns to the scene's
//! own view. The window can be resized or made full screen: the image is
//! rendered at the view's pixel size, with the scene's field of view, and
//! starts again at each new size. An overlay in the top-left corner shows
//! the progress, a stopwatch (seconds spent drawing the current image,
//! frozen once the target is reached), frame time and batch, size and traversal
//! (`i` hides and shows it);
//! the title bar shows the same, outside full screen. The Traversal menu
//! (Cmd-1 BVH, Cmd-2 Metal) switches how the bounded shapes are searched
//! and restarts the accumulation. The Samples menu sets how many frames to
//! refine up to (presets, Unlimited, Cmd-] doubles, Cmd-[ halves) without
//! restarting: a higher target carries on, a lower one stops.
//! With `--watch` a timer checks the scene file four times a second; when it
//! changed, the scene reloads at the current size and traversal and refines
//! from one sample again, keeping the view unless the camera changed in more
//! than width and height; a new size only resizes the window.
//! A file that doesn't load keeps the old image, and the status shows why.
//! Quitting (Cmd-Q, closing the window, or SIGTERM) hands the final image
//! and its size to `write`, with the `-o` path or, without one, a path
//! the user picks in a save panel after an alert asks whether to save.
//! The alert's Cancel (Esc) keeps the viewer running, except after a
//! SIGTERM; cancelling the save panel shows the alert again.
//!
//! Everything here runs on the main thread: the view and its delegate share
//! one `Rc<RefCell<Renderer>>`.

use crate::gpu::{Accumulator, Gpu};
use crate::orbit::Orbit;
use crate::scene::GpuScene;
use crate::types::{TRAVERSAL_BVH, TRAVERSAL_METAL};
use crate::watch::Watch;
use objc2::rc::Retained;
use objc2::runtime::ProtocolObject;
use objc2::{DefinedClass, MainThreadMarker, MainThreadOnly, Message, define_class, msg_send, sel};
use objc2_app_kit::{
    NSAlert, NSAlertFirstButtonReturn, NSAlertSecondButtonReturn, NSApplication,
    NSApplicationActivationPolicy, NSApplicationDelegate, NSApplicationTerminateReply,
    NSBackingStoreType, NSColor, NSControlStateValueOff, NSControlStateValueOn, NSEvent,
    NSEventModifierFlags, NSFont, NSMenu, NSMenuItem, NSModalResponseOK, NSResponder, NSSavePanel,
    NSScreen, NSTextField, NSView, NSWindow, NSWindowCollectionBehavior, NSWindowDelegate,
    NSWindowStyleMask,
};
use objc2_foundation::{
    NSArray, NSObject, NSObjectProtocol, NSPoint, NSRect, NSSize, NSString, NSTimer, NSURL,
};
use objc2_metal::{MTLCommandBuffer, MTLCommandEncoder, MTLCommandQueue, MTLPixelFormat};
use objc2_metal_kit::{MTKView, MTKViewDelegate};
use objc2_quartz_core::CAMetalDrawable;
use rtc_rust::camera::Camera;
use rtc_rust::scene::{CameraDescription, SceneDescription};
use rtc_share::Shared;
use std::cell::{Cell, RefCell};
use std::path::{Path, PathBuf};
use std::rc::Rc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::{Duration, Instant};

/// The sample target that never stops refining.
const UNLIMITED: u32 = u32::MAX;
/// The Samples menu's fixed targets.
const PRESETS: [u32; 5] = [16, 64, 256, 1024, 4096];
/// Screen refreshes per second the view asks for.
const FPS: isize = 60;
/// Seconds between checks for SIGTERM and of a `--watch`ed scene file.
const WATCH_INTERVAL: f64 = 0.25;
/// SIGTERM's number on macOS.
const SIGTERM: i32 = 15;
/// Time between screen refreshes, in ms.
const REFRESH_MS: f64 = 1000.0 / FPS as f64;

/// Set by the SIGTERM handler; the timer quits the app when it sees it,
/// and the quit then offers no Cancel.
static QUIT: AtomicBool = AtomicBool::new(false);

unsafe extern "C" {
    /// libc's `signal(3)`; returns the previous handler (unused here).
    fn signal(sig: i32, handler: extern "C" fn(i32)) -> usize;
}

/// SIGTERM handler: only flags the quit, as a handler may do nothing that
/// isn't async-signal-safe; the main thread's timer does the quitting.
extern "C" fn on_sigterm(_sig: i32) {
    QUIT.store(true, Ordering::SeqCst);
}

/// Writes an image, its width and height, to a path: on quit and for
/// Save Image…; an `Err` is shown to the user.
pub type Writer = Box<dyn Fn(Vec<[f32; 4]>, usize, usize, &Path) -> Result<(), String>>;

/// Where the image goes on quit.
pub enum Save {
    /// Write it to this path (`-o`) without asking.
    To(String),
    /// Ask, offering this file name (`image_name` of the scene).
    Ask(String),
}

/// The default image file name for `scene`: its stem with `.ppm`.
pub fn image_name(scene: &Path) -> String {
    format!(
        "{}.ppm",
        scene.file_stem().unwrap_or_default().to_string_lossy()
    )
}

/// The file name a save panel offers: the asked-for name, or the `-o`
/// file's.
fn offered_name(save: &Save) -> String {
    match save {
        Save::Ask(name) => name.clone(),
        Save::To(path) => Path::new(path)
            .file_name()
            .map_or_else(|| path.clone(), |n| n.to_string_lossy().into_owned()),
    }
}

/// The camera `orbit` shows, as from, to and up (x, y, z), with
/// `focal_distance`: what `--share` writes.
pub fn shot(orbit: &Orbit, focal_distance: f64) -> ([f64; 3], [f64; 3], [f64; 3], f64) {
    let xyz = |t: rtc_rust::tuple::Tuple| [t.x, t.y, t.z];
    (
        xyz(orbit.from()),
        xyz(orbit.to()),
        xyz(orbit.up()),
        focal_distance,
    )
}

/// What a quit does with the image.
#[derive(Debug, PartialEq)]
enum Quit {
    /// Write it to this path, then quit.
    Save(PathBuf),
    /// Quit without writing it.
    Discard,
    /// Don't quit.
    Cancel,
}

/// The quit alert's buttons.
#[derive(Clone, Copy)]
enum Button {
    Save,
    Discard,
    Cancel,
}

/// What a quit does: nothing to write before the first frame is drawn,
/// `-o`'s path without asking, otherwise whatever `ask` answers when
/// offered the default name.
fn choose(save: &Save, frames: u32, ask: impl FnOnce(&str) -> Quit) -> Quit {
    if frames == 0 {
        return Quit::Discard;
    }
    match save {
        Save::To(path) => Quit::Save(PathBuf::from(path)),
        Save::Ask(name) => ask(name),
    }
}

/// The answer from the quit alert's `button` and, after Save…, the save
/// panel's `path`; cancelling the panel shows the alert again.
fn answer(mut button: impl FnMut() -> Button, mut path: impl FnMut() -> Option<PathBuf>) -> Quit {
    loop {
        match button() {
            Button::Discard => return Quit::Discard,
            Button::Cancel => return Quit::Cancel,
            Button::Save => {
                if let Some(path) = path() {
                    return Quit::Save(path);
                }
            }
        }
    }
}

/// Whether the quit alert offers Cancel: not when the scene editor asked
/// the viewer to quit (SIGTERM), as it has already moved on.
fn offers_cancel(sigterm: bool) -> bool {
    !sigterm
}

/// The state the view and its delegate share.
struct Renderer {
    gpu: Gpu,
    scene: GpuScene,
    acc: Accumulator,
    camera: Camera,
    orbit: Orbit,
    /// The scene's own camera and orbit, for the space key.
    home: (Camera, Orbit),
    /// Whether the camera has left `home` since the last reset.
    moved: bool,
    /// Frames to refine up to; `UNLIMITED` never stops.
    target: u32,
    last_frame: Instant,
    /// The last frame's time, for the status shown while paused.
    ms: f64,
    /// Stopwatch: the time spent drawing the current image's frames.
    elapsed: Duration,
    /// Frames to draw per screen refresh; adapts so the GPU stays busy
    /// without the window lagging.
    batch: u32,
    /// The scene file's camera the current scene was loaded with; a reload
    /// keeps the view unless this changes.
    camera_desc: CameraDescription,
    /// Why the last reload failed, shown in the status until one succeeds.
    reload_error: Option<String>,
    /// The scene file, with `--watch`: reloaded whenever it changes.
    watch: Option<Watch>,
    /// The file's new `(width, height)` after a reload changed them; the
    /// view takes it to resize the window.
    resize: Option<(usize, usize)>,
    /// With `--share`: where the camera is written whenever it changes.
    share: Option<Shared>,
}

/// Whether two cameras differ at most in image size, so a reload keeps the
/// user's orbit.
fn same_view(a: &CameraDescription, b: &CameraDescription) -> bool {
    let sized = |c: &CameraDescription| CameraDescription {
        width: 0,
        height: 0,
        ..c.clone()
    };
    sized(a) == sized(b)
}

impl Renderer {
    /// Point the camera at the orbit's new position and start accumulating
    /// again.
    fn restart(&mut self) {
        self.camera.set_transform(&self.orbit.transform());
        if let Err(e) = self.scene.set_camera(&self.camera) {
            eprintln!("rtc-metal: {e}");
        }
        self.acc.reset();
        self.last_frame = Instant::now();
        self.elapsed = Duration::ZERO;
        self.publish();
    }

    /// Write the camera to the `--share` file, if any.
    fn publish(&self) {
        if let Some(share) = &self.share {
            let (from, to, up, focus) = shot(&self.orbit, self.camera.focal_distance);
            share.write_camera(from, to, up, focus);
        }
    }

    /// Return to the scene's own view, at the current size; false (and the
    /// image kept) when the camera is already there.
    fn reset_view(&mut self) -> bool {
        if !self.moved {
            return false;
        }
        let (home, orbit) = self.home.clone();
        self.camera = resized(&home, self.camera.hsize, self.camera.vsize);
        self.orbit = orbit;
        self.moved = false;
        self.restart();
        true
    }

    /// Render at `width x height` pixels from now on; false when that is
    /// already the size (or empty, as while the window is minimised).
    fn resize(&mut self, width: usize, height: usize) -> bool {
        if width == 0 || height == 0 || (width, height) == (self.camera.hsize, self.camera.vsize) {
            return false;
        }
        self.camera = resized(&self.camera, width, height);
        if let Err(e) = self.scene.set_camera(&self.camera) {
            eprintln!("rtc-metal: {e}");
            return false;
        }
        self.rebuild()
    }

    /// Refine up to `target` frames from now on, keeping the frames already
    /// summed; true when that means drawing more.
    fn set_target(&mut self, target: u32) -> bool {
        self.target = target;
        self.last_frame = Instant::now();
        self.acc.frames() < target
    }

    /// Search the bounded shapes with `traversal` (`TRAVERSAL_BVH` or
    /// `TRAVERSAL_METAL`) from now on; false when it already does.
    fn set_traversal(&mut self, traversal: u32) -> bool {
        if traversal == self.scene.info.counts.traversal {
            return false;
        }
        self.scene.info.counts.traversal = traversal;
        self.rebuild()
    }

    /// A new, empty accumulator for the scene's current size and traversal.
    fn rebuild(&mut self) -> bool {
        match self.gpu.accumulator(&self.scene) {
            Ok(acc) => self.acc = acc,
            Err(e) => {
                eprintln!("rtc-metal: {e}");
                return false;
            }
        }
        self.last_frame = Instant::now();
        self.elapsed = Duration::ZERO;
        true
    }

    /// Swap in the scene `desc` describes (asset paths relative to
    /// `base_dir`), at the current size and traversal, and start refining
    /// again. The view is kept unless the camera changed in more than width
    /// and height, which then becomes the new home view. On error nothing
    /// changes but `reload_error`.
    fn reload(&mut self, desc: &SceneDescription, base_dir: &Path) -> Result<(), String> {
        let result = self.try_reload(desc, base_dir);
        self.reload_error = result.as_ref().err().cloned();
        result
    }

    fn try_reload(&mut self, desc: &SceneDescription, base_dir: &Path) -> Result<(), String> {
        let (world, file_camera) = desc.build(base_dir).map_err(|e| e.0)?;
        let new_view = (!same_view(&desc.camera, &self.camera_desc)).then(|| {
            let camera = resized(&file_camera, self.camera.hsize, self.camera.vsize);
            (camera, Orbit::of(&desc.camera))
        });
        let camera = new_view.as_ref().map_or(self.camera, |(c, _)| *c);
        let mut scene = GpuScene::from_world(&world, &camera)?;
        scene.info.counts.traversal = self.scene.info.counts.traversal;
        let acc = self.gpu.accumulator(&scene)?;
        self.scene = scene;
        self.acc = acc;
        let size = |c: &CameraDescription| (c.width, c.height);
        self.resize = (size(&desc.camera) != size(&self.camera_desc)).then(|| size(&desc.camera));
        self.camera_desc = desc.camera.clone();
        if let Some((camera, orbit)) = new_view {
            self.camera = camera;
            self.home = (camera, orbit.clone());
            self.orbit = orbit;
            self.moved = false;
            self.publish();
        }
        self.elapsed = Duration::ZERO;
        self.last_frame = Instant::now();
        Ok(())
    }

    /// Frames the next draw adds: the batch, but no more than the target
    /// leaves.
    fn frames_to_draw(&self) -> u32 {
        self.batch.min(self.target - self.acc.frames())
    }

    /// A draw of `frames` frames finished at `now`: the time per frame in
    /// ms, also kept in `ms`; the draw's time goes on the stopwatch and
    /// sets the next batch.
    fn tick(&mut self, now: Instant, frames: u32) -> f64 {
        let draw = now - self.last_frame;
        self.last_frame = now;
        self.elapsed += draw;
        let draw_ms = draw.as_secs_f64() * 1e3;
        self.batch = next_batch(self.batch, draw_ms);
        self.ms = draw_ms / frames as f64;
        self.ms
    }
}

/// The batch after a draw that took `draw_ms`: one more while draws keep
/// up with the screen (up to a quarter over its refresh interval), half
/// (at least 1) once they lag.
fn next_batch(batch: u32, draw_ms: f64) -> u32 {
    if draw_ms <= 1.25 * REFRESH_MS {
        batch + 1
    } else {
        (batch / 2).max(1)
    }
}

/// Twice `target`, up to 2^31; Unlimited stays Unlimited.
fn doubled(target: u32) -> u32 {
    if target == UNLIMITED {
        target
    } else {
        target.saturating_mul(2).min(1 << 31)
    }
}

/// Half `target`, at least 1; from Unlimited, the largest preset.
fn halved(target: u32) -> u32 {
    if target == UNLIMITED {
        PRESETS[PRESETS.len() - 1]
    } else {
        (target / 2).max(1)
    }
}

/// The content size for a `width` x `height` image on a screen with backing
/// `scale`: one drawable pixel per image pixel, shrunk uniformly to fit
/// `room` (the content area available, if there is a screen), never grown.
fn window_size(width: usize, height: usize, scale: f64, room: Option<NSSize>) -> NSSize {
    let (w, h) = (width as f64 / scale, height as f64 / scale);
    let fit = room.map_or(1.0, |room| (room.width / w).min(room.height / h).min(1.0));
    NSSize::new(w * fit, h * fit)
}

/// The narrowest room beside the editor worth placing the window in;
/// with less, the window is centred.
const MIN_ROOM: f64 = 200.0;

/// `--at X,Y`: two finite numbers.
pub fn parse_at(s: &str) -> Option<(f64, f64)> {
    let (x, y) = s.split_once(',')?;
    let (x, y) = (x.parse::<f64>().ok()?, y.parse::<f64>().ok()?);
    (x.is_finite() && y.is_finite()).then_some((x, y))
}

/// Where a window whose frame's top-left is `at` goes: `at` is in points,
/// right and down from the top-left of the primary screen, the one at
/// AppKit's origin. `screens` are each screen's frame and visible frame.
/// Gives the top-left in AppKit's coordinates, moved down below the menu
/// bar if above it, and the room the frame has to the visible right and
/// bottom edges; `None` when no screen holds the point or the room is
/// narrower or lower than `MIN_ROOM`.
pub fn place(at: (f64, f64), screens: &[(NSRect, NSRect)]) -> Option<(NSPoint, NSSize)> {
    let primary = screens
        .iter()
        .find(|(frame, _)| frame.origin == NSPoint::new(0.0, 0.0))
        .or(screens.first())?;
    let (x, y) = (at.0, primary.0.size.height - at.1);
    let (_, visible) = screens.iter().find(|(frame, _)| {
        let (min, max) = (frame.origin, frame.origin.x + frame.size.width);
        (min.x..=max).contains(&x) && (min.y..=min.y + frame.size.height).contains(&y)
    })?;
    let top = y.min(visible.origin.y + visible.size.height);
    let top_left = NSPoint::new(x, top);
    let room = room_from(top_left, *visible);
    (room.width >= MIN_ROOM && room.height >= MIN_ROOM).then_some((top_left, room))
}

/// `top_left` moved left and up just enough that a window frame of `size`
/// fits inside `visible`; a frame larger than `visible` keeps its top-left
/// corner on the screen.
pub fn on_screen(top_left: NSPoint, size: NSSize, visible: NSRect) -> NSPoint {
    let (left, bottom) = (visible.origin.x, visible.origin.y);
    let (right, top) = (left + visible.size.width, bottom + visible.size.height);
    NSPoint::new(
        top_left.x.min(right - size.width).max(left),
        top_left.y.max(bottom + size.height).min(top),
    )
}

/// The room for a window resized on a reload, keeping its top-left: to the
/// visible right and bottom edges, or the whole visible frame when that
/// leaves less than `MIN_ROOM` either way.
pub fn reload_room(top_left: NSPoint, visible: NSRect) -> NSSize {
    let room = room_from(top_left, visible);
    if room.width >= MIN_ROOM && room.height >= MIN_ROOM {
        room
    } else {
        visible.size
    }
}

/// The room from a window frame's top-left (AppKit coordinates) to the
/// right and bottom edges of `visible`.
pub fn room_from(top_left: NSPoint, visible: NSRect) -> NSSize {
    NSSize::new(
        visible.origin.x + visible.size.width - top_left.x,
        top_left.y - visible.origin.y,
    )
}

/// The overlay's and title's text.
fn status(
    name: &str,
    frames: u32,
    target: u32,
    secs: f64,
    (ms, batch): (f64, u32),
    (width, height): (usize, usize),
    traversal: u32,
) -> String {
    let traversal = if traversal == TRAVERSAL_METAL {
        "metal"
    } else {
        "bvh"
    };
    let target = if target == UNLIMITED {
        "∞".to_string()
    } else {
        target.to_string()
    };
    format!(
        "{name} · {frames}/{target} samples · {secs:.1} s · {ms:.1} ms/frame ×{batch} · {width}×{height} px · {traversal}"
    )
}

/// `camera` at `width x height` pixels: same field of view, view transform
/// and lens. The field of view spans the longer side, so a wider window
/// shows more to the sides.
fn resized(camera: &Camera, width: usize, height: usize) -> Camera {
    let mut c = Camera::new(width, height, camera.field_of_view);
    c.set_transform(&camera.transform());
    c.aa_samples = camera.aa_samples;
    c.aa_threshold = camera.aa_threshold;
    c.aperture = camera.aperture;
    c.focal_distance = camera.focal_distance;
    c.dof_samples = camera.dof_samples;
    c.motion_samples = camera.motion_samples;
    c
}

define_class!(
    // SAFETY: MTKView has no subclassing requirements; `RtcView` does not
    // implement `Drop`.
    #[unsafe(super(MTKView, NSView, NSResponder, NSObject))]
    #[thread_kind = MainThreadOnly]
    #[ivars = Rc<RefCell<Renderer>>]
    /// The image view. Dragging orbits the camera "grabbing the scene"
    /// (dragging right orbits it left, dragging down raises it; the signs
    /// live in `Orbit::drag`); scrolling zooms; space resets the view; `i`
    /// toggles the status overlay.
    struct RtcView;

    impl RtcView {
        #[unsafe(method(acceptsFirstResponder))]
        fn accepts_first_responder(&self) -> bool {
            true
        }

        #[unsafe(method(mouseDragged:))]
        fn mouse_dragged(&self, event: &NSEvent) {
            // AppKit deltas are in points, +y = mouse down; `Orbit::drag`
            // owns the sign convention ("grab the scene": dragging right
            // orbits the camera left, dragging down raises it).
            let (dx, dy) = (event.deltaX(), event.deltaY());
            if dx == 0.0 && dy == 0.0 {
                return; // keep the converged image
            }
            self.interact(|orbit| orbit.drag(dx, dy));
        }

        #[unsafe(method(scrollWheel:))]
        fn scroll_wheel(&self, event: &NSEvent) {
            // Trackpad gesture begin/end, horizontal scrolls and tilt
            // wheels send dy = 0: keep the converged image.
            let dy = event.scrollingDeltaY();
            if dy == 0.0 {
                return;
            }
            // Trackpads report precise deltas in points, many per gesture;
            // wheels report whole lines.
            let dy = if event.hasPreciseScrollingDeltas() {
                dy * 0.1
            } else {
                dy
            };
            self.interact(|orbit| orbit.zoom(dy));
        }

        #[unsafe(method(useBvhTraversal:))]
        fn use_bvh_traversal(&self, _sender: Option<&NSObject>) {
            self.switch_traversal(TRAVERSAL_BVH);
        }

        #[unsafe(method(useMetalTraversal:))]
        fn use_metal_traversal(&self, _sender: Option<&NSObject>) {
            self.switch_traversal(TRAVERSAL_METAL);
        }

        /// A Samples preset; the item's tag is the target (0 = Unlimited).
        #[unsafe(method(useSampleTarget:))]
        fn use_sample_target(&self, sender: &NSMenuItem) {
            self.change_target(|_| tag_target(sender.tag()));
        }

        #[unsafe(method(doubleSamples:))]
        fn double_samples(&self, _sender: Option<&NSObject>) {
            self.change_target(doubled);
        }

        #[unsafe(method(halveSamples:))]
        fn halve_samples(&self, _sender: Option<&NSObject>) {
            self.change_target(halved);
        }

        /// The timer: quit after a SIGTERM, otherwise reload a `--watch`ed
        /// scene when its file changed.
        #[unsafe(method(checkQuitAndWatchedFile:))]
        fn check_quit_and_watched_file(&self, _timer: &NSTimer) {
            // QUIT stays set: the quit it starts offers no Cancel, and a
            // repeat while the alert is up is refused.
            if QUIT.load(Ordering::SeqCst) {
                NSApplication::sharedApplication(self.mtm()).terminate(None);
                return;
            }
            self.reload_if_changed();
        }

        /// AppKit asks before showing a menu: tick the current traversal and
        /// sample target.
        #[unsafe(method(validateMenuItem:))]
        fn validate_menu_item(&self, item: &NSMenuItem) -> bool {
            let r = self.ivars().borrow();
            let ticked = match item.action() {
                Some(a) if a == sel!(useBvhTraversal:) => {
                    Some(r.scene.info.counts.traversal == TRAVERSAL_BVH)
                }
                Some(a) if a == sel!(useMetalTraversal:) => {
                    Some(r.scene.info.counts.traversal == TRAVERSAL_METAL)
                }
                Some(a) if a == sel!(useSampleTarget:) => Some(tag_target(item.tag()) == r.target),
                _ => None,
            };
            if let Some(ticked) = ticked {
                item.setState(if ticked {
                    NSControlStateValueOn
                } else {
                    NSControlStateValueOff
                });
            }
            true
        }

        #[unsafe(method(keyDown:))]
        fn key_down(&self, event: &NSEvent) {
            if event.charactersIgnoringModifiers().is_some_and(|c| c.to_string() == "i") {
                // The status overlay is the view's only subview.
                for overlay in self.subviews() {
                    overlay.setHidden(!overlay.isHidden());
                }
                return;
            }
            // 49 is the space bar's virtual key code on every layout.
            if event.keyCode() != 49 {
                // SAFETY: NSResponder's keyDown: takes one NSEvent.
                unsafe { msg_send![super(self), keyDown: event] }
                return;
            }
            if self.ivars().borrow_mut().reset_view() {
                self.resume();
            }
        }
    }
);

impl RtcView {
    fn new(
        mtm: MainThreadMarker,
        renderer: Rc<RefCell<Renderer>>,
        frame: NSRect,
    ) -> Retained<Self> {
        let device = renderer.borrow().gpu.device().retain();
        let this = Self::alloc(mtm).set_ivars(renderer);
        // SAFETY: `initWithFrame:device:` is MTKView's designated initialiser,
        // with the argument types its header declares.
        unsafe { msg_send![super(this), initWithFrame: frame, device: Some(&*device)] }
    }

    /// Move the camera with `change`, restart accumulation and resume drawing.
    fn interact(&self, change: impl FnOnce(&mut Orbit)) {
        {
            let mut r = self.ivars().borrow_mut();
            let before = r.orbit.radius();
            change(&mut r.orbit);
            // Zooming keeps the lens focused on what it was focused on.
            r.camera.focal_distance *= r.orbit.radius() / before;
            r.moved = true;
            r.restart();
        }
        self.resume();
    }

    /// Set the sample target to `change(target)`; draw more if that is
    /// above the frames done, otherwise just show the new target.
    fn change_target(&self, change: impl FnOnce(u32) -> u32) {
        let more = {
            let mut r = self.ivars().borrow_mut();
            let target = change(r.target);
            r.set_target(target)
        };
        if more {
            self.resume();
        } else {
            // Paused: one draw call refreshes the status text.
            self.setNeedsDisplay(true);
        }
    }

    /// Switch traversal from the menu and restart drawing.
    fn switch_traversal(&self, traversal: u32) {
        if self.ivars().borrow_mut().set_traversal(traversal) {
            self.resume();
        }
    }

    /// Start drawing frames again after the accumulation restarted.
    fn resume(&self) {
        self.setPaused(false);
        self.setEnableSetNeedsDisplay(false);
    }

    /// Reload the watched scene file if it changed: refine the new scene, or
    /// keep the old one and show why the reload failed.
    fn reload_if_changed(&self) {
        let reloaded = {
            let mut r = self.ivars().borrow_mut();
            let Some(watch) = r.watch.as_mut() else {
                return;
            };
            if !watch.changed() {
                return;
            }
            let path = watch.path().to_path_buf();
            let base_dir = path.parent().unwrap_or(Path::new("."));
            let result = SceneDescription::from_path(&path)
                .map_err(|e| e.0)
                .and_then(|desc| r.reload(&desc, base_dir));
            if let Err(e) = &result {
                eprintln!("rtc-metal: reload failed: {e}");
                // A file that doesn't parse never reaches `reload`.
                r.reload_error = Some(e.clone());
            }
            let resize = if result.is_ok() {
                r.resize.take()
            } else {
                None
            };
            (result.is_ok(), resize)
        };
        let (reloaded, resize) = reloaded;
        if let (Some((w, h)), Some(window)) = (resize, self.window()) {
            let scale = window.backingScaleFactor();
            // The window keeps its top-left; AppKit's resize keeps the
            // bottom-left.
            let frame = window.frame();
            let top_left = NSPoint::new(frame.origin.x, frame.origin.y + frame.size.height);
            let visible = window.screen().map(|screen| screen.visibleFrame());
            let room = visible.map(|visible| {
                NSWindow::contentRectForFrameRect_styleMask(
                    NSRect::new(NSPoint::new(0.0, 0.0), reload_room(top_left, visible)),
                    window.styleMask(),
                    MainThreadMarker::from(self),
                )
                .size
            });
            window.setContentSize(window_size(w, h, scale, room));
            // Grown back to the whole screen near an edge, it moves back
            // onto the screen.
            let size = window.frame().size;
            window.setFrameTopLeftPoint(
                visible.map_or(top_left, |visible| on_screen(top_left, size, visible)),
            );
        }
        if reloaded {
            self.resume();
        } else {
            // Paused: one draw call shows the message.
            self.setNeedsDisplay(true);
        }
    }
}

define_class!(
    // SAFETY: NSObject has no subclassing requirements; `RtcDelegate` does
    // not implement `Drop`.
    #[unsafe(super(NSObject))]
    #[thread_kind = MainThreadOnly]
    #[ivars = Delegate]
    struct RtcDelegate;

    impl RtcDelegate {
        /// File ▸ Save Image…: write the image as it stands where the user
        /// picks; rendering carries on.
        #[unsafe(method(saveImage:))]
        fn save_image(&self, _sender: Option<&NSObject>) {
            let name = offered_name(&self.ivars().save);
            // No borrow is held while the panel runs: the view may draw.
            if let Some(path) = save_path(self.mtm(), &name) {
                self.write_image(&path);
            }
        }
    }

    unsafe impl NSObjectProtocol for RtcDelegate {}

    unsafe impl MTKViewDelegate for RtcDelegate {
        #[unsafe(method(drawInMTKView:))]
        fn draw_in_mtk_view(&self, view: &MTKView) {
            self.draw(view);
        }

        #[unsafe(method(mtkView:drawableSizeWillChange:))]
        fn drawable_size_will_change(&self, view: &MTKView, size: NSSize) {
            // The window was resized or went full screen; `size` is in pixels.
            let resized = self
                .ivars()
                .renderer
                .borrow_mut()
                .resize(size.width as usize, size.height as usize);
            if resized {
                view.setPaused(false);
                view.setEnableSetNeedsDisplay(false);
            }
        }
    }

    unsafe impl NSWindowDelegate for RtcDelegate {
        /// Closing the window quits, asking first, so Cancel keeps it open.
        #[unsafe(method(windowShouldClose:))]
        fn window_should_close(&self, _sender: &NSWindow) -> bool {
            NSApplication::sharedApplication(self.mtm()).terminate(None);
            false
        }
    }

    unsafe impl NSApplicationDelegate for RtcDelegate {
        #[unsafe(method(applicationShouldTerminateAfterLastWindowClosed:))]
        fn should_terminate_after_last_window_closed(&self, _sender: &NSApplication) -> bool {
            true
        }

        /// Write the image where `save` says (asking first without `-o`),
        /// then quit, unless the user cancels; a second quit while asking
        /// is cancelled.
        #[unsafe(method(applicationShouldTerminate:))]
        fn should_terminate(&self, _sender: &NSApplication) -> NSApplicationTerminateReply {
            let ivars = self.ivars();
            if ivars.asking.replace(true) {
                return NSApplicationTerminateReply::TerminateCancel;
            }
            // No borrow is held while the alert and panel run: the view may
            // draw meanwhile.
            let frames = ivars.renderer.borrow().acc.frames();
            let mtm = self.mtm();
            let cancel = offers_cancel(QUIT.load(Ordering::SeqCst));
            let quit = choose(&ivars.save, frames, |name| ask(mtm, name, cancel));
            if quit == Quit::Cancel {
                ivars.asking.set(false);
                return NSApplicationTerminateReply::TerminateCancel;
            }
            if let Quit::Save(path) = quit {
                self.write_image(&path);
            }
            NSApplicationTerminateReply::TerminateNow
        }
    }
);

struct Delegate {
    renderer: Rc<RefCell<Renderer>>,
    /// The status label over the image.
    overlay: Retained<NSTextField>,
    /// Where the image goes on quit.
    save: Save,
    /// Set while quitting asks where to save.
    asking: Cell<bool>,
    /// Writes the image to a file.
    write: Writer,
    title: String,
}

impl RtcDelegate {
    fn new(mtm: MainThreadMarker, ivars: Delegate) -> Retained<Self> {
        let this = Self::alloc(mtm).set_ivars(ivars);
        // SAFETY: NSObject's `init`.
        unsafe { msg_send![super(this), init] }
    }

    /// Write the image as it stands to `path`; a failure is shown.
    fn write_image(&self, path: &Path) {
        let ivars = self.ivars();
        let image = {
            let r = ivars.renderer.borrow();
            let size = (r.camera.hsize, r.camera.vsize);
            r.acc.resolve(&r.gpu).map(|pixels| (pixels, size))
        };
        let written =
            image.and_then(|(pixels, (width, height))| (ivars.write)(pixels, width, height, path));
        if let Err(e) = written {
            eprintln!("rtc-metal: {e}");
            show_error(self.mtm(), &e);
        }
    }

    /// Add a batch of frames and show the running mean; pause once
    /// `target` frames are in.
    fn draw(&self, view: &MTKView) {
        let ivars = self.ivars();
        let mut r = ivars.renderer.borrow_mut();
        if r.acc.frames() >= r.target {
            let ms = r.ms;
            drop(r);
            self.show_status(view, ms);
            view.setPaused(true);
            view.setEnableSetNeedsDisplay(true);
            return;
        }
        let Some(drawable) = view.currentDrawable() else {
            return;
        };
        let Some(commands) = r.gpu.queue().commandBuffer() else {
            eprintln!("rtc-metal: no command buffer");
            return;
        };
        let Some(encoder) = commands.computeCommandEncoder() else {
            eprintln!("rtc-metal: no compute encoder");
            return;
        };
        let info = r.scene.info;
        let frames = r.frames_to_draw();
        for _ in 0..frames {
            r.acc.encode_frame(&encoder, &info);
        }
        r.acc.encode_resolve_texture(&encoder, &drawable.texture());
        encoder.endEncoding();
        commands.presentDrawable(ProtocolObject::from_ref(&*drawable));
        commands.commit();

        let ms = r.tick(Instant::now(), frames);
        drop(r);
        self.show_status(view, ms);
    }

    /// Put the status text in the title and the overlay.
    fn show_status(&self, view: &MTKView, ms: f64) {
        let ivars = self.ivars();
        let r = ivars.renderer.borrow();
        let mut text = status(
            &ivars.title,
            r.acc.frames(),
            r.target,
            r.elapsed.as_secs_f64(),
            (ms, r.batch),
            (r.camera.hsize, r.camera.vsize),
            r.scene.info.counts.traversal,
        );
        if let Some(e) = &r.reload_error {
            text.push_str(&format!(" · reload failed: {e}"));
        }
        drop(r);
        let text = NSString::from_str(&text);
        if let Some(window) = view.window() {
            window.setTitle(&text);
        }
        let overlay = &ivars.overlay;
        overlay.setStringValue(&text);
        overlay.sizeToFit();
        // AppKit's y axis points up: pin the label 8 points from the top left.
        let top = view.bounds().size.height;
        let size = overlay.frame().size;
        overlay.setFrameOrigin(NSPoint::new(8.0, top - size.height - 8.0));
    }
}

/// Open a window showing `scene` refining up to `samples` frames, viewed
/// from `camera` (the scene file's `camera_desc`, at the scene's size), run
/// the app until it quits, then exit. With `watch`, the scene reloads when
/// its file changes. On quit, `write` gets the final image (row-major
/// RGBA, as `Accumulator::resolve`) and the path `save` gives (see
/// `choose`), if at least one frame was drawn and a path was chosen. A
/// SIGTERM quits as Cmd-Q does. With `share`, the camera is written there
/// at start-up and whenever it changes.
#[allow(clippy::too_many_arguments)]
pub fn run(
    gpu: Gpu,
    scene: GpuScene,
    camera: Camera,
    camera_desc: CameraDescription,
    samples: u32,
    title: String,
    watch: Option<Watch>,
    at: Option<(f64, f64)>,
    share: Option<Shared>,
    save: Save,
    write: Writer,
) -> ! {
    let mtm = MainThreadMarker::new().expect("the window must run on the main thread");
    let (width, height) = (scene.width() as f64, scene.height() as f64);
    let acc = gpu.accumulator(&scene).unwrap_or_else(|e| {
        eprintln!("rtc-metal: {e}");
        std::process::exit(1);
    });
    let orbit = Orbit::of(&camera_desc);
    let renderer = Rc::new(RefCell::new(Renderer {
        gpu,
        scene,
        acc,
        camera,
        home: (camera, orbit.clone()),
        orbit,
        moved: false,
        target: samples,
        last_frame: Instant::now(),
        ms: 0.0,
        elapsed: Duration::ZERO,
        batch: 1,
        camera_desc,
        reload_error: None,
        watch,
        resize: None,
        share,
    }));
    renderer.borrow().publish();

    let app = NSApplication::sharedApplication(mtm);
    app.setActivationPolicy(NSApplicationActivationPolicy::Regular);
    let overlay = status_label(mtm, &title);
    let delegate = RtcDelegate::new(
        mtm,
        Delegate {
            renderer: renderer.clone(),
            overlay: overlay.clone(),
            save,
            asking: Cell::new(false),
            write,
            title: title.clone(),
        },
    );
    app.setDelegate(Some(ProtocolObject::from_ref(&*delegate)));
    app.setMainMenu(Some(&main_menu(mtm, &delegate)));

    // One drawable pixel per image pixel, shrunk uniformly so the window,
    // title bar included, fits the screen.
    let style = NSWindowStyleMask::Titled
        | NSWindowStyleMask::Closable
        | NSWindowStyleMask::Miniaturizable
        | NSWindowStyleMask::Resizable;
    let screen = NSScreen::mainScreen(mtm);
    let scale = screen.as_ref().map_or(1.0, |s| s.backingScaleFactor());
    // With `--at`, the room right of that point on its screen; else the
    // main screen's.
    let placed = at.and_then(|at| {
        let screens: Vec<_> = NSScreen::screens(mtm)
            .iter()
            .map(|s| (s.frame(), s.visibleFrame()))
            .collect();
        place(at, &screens)
    });
    let room = match placed {
        Some((_, room)) => Some(
            NSWindow::contentRectForFrameRect_styleMask(
                NSRect::new(NSPoint::new(0.0, 0.0), room),
                style,
                mtm,
            )
            .size,
        ),
        None => screen.map(|screen| {
            NSWindow::contentRectForFrameRect_styleMask(screen.visibleFrame(), style, mtm).size
        }),
    };
    let frame = NSRect::new(
        NSPoint::new(0.0, 0.0),
        window_size(width as usize, height as usize, scale, room),
    );

    let view = RtcView::new(mtm, renderer, frame);
    view.setDelegate(Some(ProtocolObject::from_ref(&*delegate)));
    view.setColorPixelFormat(MTLPixelFormat::BGRA8Unorm);
    view.setFramebufferOnly(false);
    // The drawable follows the view's size in pixels; each change reaches
    // `drawableSizeWillChange`, which re-sizes the camera.
    view.setAutoResizeDrawable(true);
    view.setPreferredFramesPerSecond(FPS);
    view.addSubview(&overlay);
    // SAFETY: `checkQuitAndWatchedFile:` takes the timer, as the selector
    // expects; the view lives until the app exits, and so does the
    // repeating timer, which the run loop retains.
    let _ = unsafe {
        NSTimer::scheduledTimerWithTimeInterval_target_selector_userInfo_repeats(
            WATCH_INTERVAL,
            &view,
            sel!(checkQuitAndWatchedFile:),
            None,
            true,
        )
    };
    // SAFETY: `on_sigterm` has the signature `signal` expects and only
    // stores to an atomic, which is async-signal-safe.
    unsafe { signal(SIGTERM, on_sigterm) };

    // SAFETY: `releasedWhenClosed` is turned off right below, before the
    // window can be closed, so the `Retained` stays the only owner.
    let window = unsafe {
        NSWindow::initWithContentRect_styleMask_backing_defer(
            NSWindow::alloc(mtm),
            frame,
            style,
            NSBackingStoreType::Buffered,
            false,
        )
    };
    // SAFETY: see above.
    unsafe { window.setReleasedWhenClosed(false) };
    window.setTitle(&NSString::from_str(&title));
    window.setContentView(Some(&view));
    window.setCollectionBehavior(NSWindowCollectionBehavior::FullScreenPrimary);
    window.setDelegate(Some(ProtocolObject::from_ref(&*delegate)));
    match placed {
        Some((top_left, _)) => window.setFrameTopLeftPoint(top_left),
        None => window.center(),
    }
    window.makeKeyAndOrderFront(None);
    window.makeFirstResponder(Some(&view));
    app.activate(); // needs macOS 14+ (older: activateIgnoringOtherApps)
    app.run();
    std::process::exit(0);
}

/// Ask whether to save the image and, if so, where, offering `name` in the
/// current directory, with a Cancel button (Esc) when `cancel`.
fn ask(mtm: MainThreadMarker, name: &str, cancel: bool) -> Quit {
    let app = NSApplication::sharedApplication(mtm);
    // A SIGTERM may come while another app (the scene editor) is active:
    // the plain `activate` only asks, so force the alert to the front.
    #[allow(deprecated)]
    app.activateIgnoringOtherApps(true);
    let alert = NSAlert::new(mtm);
    alert.setMessageText(&NSString::from_str("Save the rendered image?"));
    alert.addButtonWithTitle(&NSString::from_str("Save…"));
    alert.addButtonWithTitle(&NSString::from_str("Quit Without Saving"));
    if cancel {
        // AppKit gives a button titled Cancel the Esc key.
        alert.addButtonWithTitle(&NSString::from_str("Cancel"));
    }
    let button = || match alert.runModal() {
        r if r == NSAlertFirstButtonReturn => Button::Save,
        r if r == NSAlertSecondButtonReturn => Button::Discard,
        _ => Button::Cancel,
    };
    answer(button, || save_path(mtm, name))
}

/// A save panel offering `name` in the current directory; the chosen path,
/// or `None` when cancelled.
fn save_path(mtm: MainThreadMarker, name: &str) -> Option<PathBuf> {
    let panel = NSSavePanel::savePanel(mtm);
    if let Ok(dir) = std::env::current_dir() {
        let dir = NSString::from_str(&dir.to_string_lossy());
        panel.setDirectoryURL(Some(&NSURL::fileURLWithPath_isDirectory(&dir, true)));
    }
    panel.setNameFieldStringValue(&NSString::from_str(name));
    // `allowedContentTypes` would need the UniformTypeIdentifiers crate.
    #[allow(deprecated)]
    panel.setAllowedFileTypes(Some(&NSArray::from_retained_slice(&[NSString::from_str(
        "ppm",
    )])));
    if panel.runModal() != NSModalResponseOK {
        return None;
    }
    panel
        .URL()
        .and_then(|url| url.path())
        .map(|path| PathBuf::from(path.to_string()))
}

/// Tell the user the image could not be written, and why.
fn show_error(mtm: MainThreadMarker, error: &str) {
    let alert = NSAlert::new(mtm);
    alert.setMessageText(&NSString::from_str("The image could not be saved."));
    alert.setInformativeText(&NSString::from_str(error));
    alert.runModal();
}

/// White text on a translucent black box, for the top-left overlay.
fn status_label(mtm: MainThreadMarker, text: &str) -> Retained<NSTextField> {
    let label = NSTextField::labelWithString(&NSString::from_str(text), mtm);
    label.setTextColor(Some(&NSColor::whiteColor()));
    label.setBackgroundColor(Some(&NSColor::colorWithWhite_alpha(0.0, 0.55)));
    label.setDrawsBackground(true);
    // Fixed-width digits so the counters don't jitter as they change.
    label.setFont(Some(&NSFont::monospacedDigitSystemFontOfSize_weight(
        12.0, 0.0,
    )));
    label
}

/// The main menu: the app menu with Hide, Hide Others, Show All and Quit
/// (Cmd-Q; a bare `NSApplication` has no menu, so Cmd-Q would do nothing),
/// File with Save Image… (Cmd-S, sent to `delegate`), Traversal and
/// Samples.
fn main_menu(mtm: MainThreadMarker, delegate: &RtcDelegate) -> Retained<NSMenu> {
    let item = |name: &str, action, key: &str| {
        // SAFETY: each action is NSApplication's (reached through the
        // responder chain) or `RtcDelegate`'s, taking one sender.
        unsafe {
            NSMenuItem::initWithTitle_action_keyEquivalent(
                NSMenuItem::alloc(mtm),
                &NSString::from_str(name),
                Some(action),
                &NSString::from_str(key),
            )
        }
    };
    let app_menu = NSMenu::new(mtm);
    app_menu.addItem(&item("Hide rtc-metal", sel!(hide:), "h"));
    let others = item("Hide Others", sel!(hideOtherApplications:), "h");
    others
        .setKeyEquivalentModifierMask(NSEventModifierFlags::Option | NSEventModifierFlags::Command);
    app_menu.addItem(&others);
    app_menu.addItem(&item("Show All", sel!(unhideAllApplications:), ""));
    app_menu.addItem(&NSMenuItem::separatorItem(mtm));
    app_menu.addItem(&item("Quit", sel!(terminate:), "q"));
    let app_item = NSMenuItem::new(mtm);
    app_item.setSubmenu(Some(&app_menu));
    let file_menu = NSMenu::initWithTitle(NSMenu::alloc(mtm), &NSString::from_str("File"));
    let save = item("Save Image…", sel!(saveImage:), "s");
    // SAFETY: the delegate lives until the app exits.
    unsafe { save.setTarget(Some(delegate)) };
    file_menu.addItem(&save);
    let file_item = NSMenuItem::new(mtm);
    file_item.setSubmenu(Some(&file_menu));
    let menu = NSMenu::new(mtm);
    menu.addItem(&app_item);
    menu.addItem(&file_item);
    menu.addItem(&traversal_menu(mtm));
    menu.addItem(&samples_menu(mtm));
    menu
}

/// BVH (Cmd-1) and Metal (Cmd-2). The actions have no target, so they go
/// to the first responder, the image view.
fn traversal_menu(mtm: MainThreadMarker) -> Retained<NSMenuItem> {
    let title = NSString::from_str("Traversal");
    let submenu = NSMenu::initWithTitle(NSMenu::alloc(mtm), &title);
    for (name, action, key) in [
        ("BVH", sel!(useBvhTraversal:), "1"),
        ("Metal", sel!(useMetalTraversal:), "2"),
    ] {
        // SAFETY: both actions are `RtcView` methods taking one sender.
        let item = unsafe {
            NSMenuItem::initWithTitle_action_keyEquivalent(
                NSMenuItem::alloc(mtm),
                &NSString::from_str(name),
                Some(action),
                &NSString::from_str(key),
            )
        };
        submenu.addItem(&item);
    }
    let item = NSMenuItem::new(mtm);
    item.setSubmenu(Some(&submenu));
    item
}

/// A Samples preset item's tag as a target: 0 is Unlimited.
fn tag_target(tag: isize) -> u32 {
    if tag == 0 { UNLIMITED } else { tag as u32 }
}

/// The presets and Unlimited (ticked when current), then Double (Cmd-])
/// and Halve (Cmd-[).
fn samples_menu(mtm: MainThreadMarker) -> Retained<NSMenuItem> {
    let submenu = NSMenu::initWithTitle(NSMenu::alloc(mtm), &NSString::from_str("Samples"));
    let item = |name: &str, action, key: &str| {
        // SAFETY: every action is an `RtcView` method taking one sender.
        unsafe {
            NSMenuItem::initWithTitle_action_keyEquivalent(
                NSMenuItem::alloc(mtm),
                &NSString::from_str(name),
                Some(action),
                &NSString::from_str(key),
            )
        }
    };
    for target in PRESETS {
        let preset = item(&target.to_string(), sel!(useSampleTarget:), "");
        preset.setTag(target as isize);
        submenu.addItem(&preset);
    }
    let unlimited = item("Unlimited", sel!(useSampleTarget:), "");
    unlimited.setTag(0);
    submenu.addItem(&unlimited);
    submenu.addItem(&NSMenuItem::separatorItem(mtm));
    submenu.addItem(&item("Double", sel!(doubleSamples:), "]"));
    submenu.addItem(&item("Halve", sel!(halveSamples:), "["));
    let menu_item = NSMenuItem::new(mtm);
    menu_item.setSubmenu(Some(&submenu));
    menu_item
}

#[cfg(test)]
mod tests {
    use super::*;
    use rtc_rust::matrix::Matrix4;
    use rtc_rust::tuple::Tuple;

    #[test]
    fn shot_gives_the_orbits_camera() {
        let mut o = Orbit::new(
            Tuple::point(0.0, 1.0, -5.0),
            Tuple::point(0.0, 1.0, 0.0),
            Tuple::vector(0.0, 1.0, 0.0),
        );
        let close = |a: [f64; 3], b: [f64; 3]| a.iter().zip(b).all(|(a, b)| (a - b).abs() < 1e-9);
        let (from, to, up, focus) = shot(&o, 4.5);
        assert!(close(from, [0.0, 1.0, -5.0]), "{from:?}");
        assert!(close(to, [0.0, 1.0, 0.0]) && close(up, [0.0, 1.0, 0.0]));
        assert_eq!(focus, 4.5);
        o.zoom(3.0);
        let f = o.from();
        assert!(close(shot(&o, 4.5).0, [f.x, f.y, f.z]));
    }
    use rtc_rust::scene::LightDescription;

    fn renderer() -> Result<Renderer, String> {
        let world = rtc_rust::world::World::default_world();
        let camera = Camera::new(8, 6, 1.0);
        let gpu = Gpu::new()?;
        let scene = GpuScene::from_world(&world, &camera)?;
        let acc = gpu.accumulator(&scene)?;
        let orbit = Orbit::new(
            rtc_rust::tuple::Tuple::point(0.0, 0.0, -5.0),
            rtc_rust::tuple::Tuple::point(0.0, 0.0, 0.0),
            rtc_rust::tuple::Tuple::vector(0.0, 1.0, 0.0),
        );
        Ok(Renderer {
            gpu,
            scene,
            acc,
            camera,
            home: (camera, orbit.clone()),
            orbit,
            moved: false,
            target: 4,
            last_frame: Instant::now(),
            ms: 0.0,
            elapsed: Duration::ZERO,
            batch: 1,
            camera_desc: CameraDescription {
                width: 8,
                height: 6,
                field_of_view: 1.0,
                from: [0.0, 0.0, -5.0],
                to: [0.0, 0.0, 0.0],
                up: [0.0, 1.0, 0.0],
                antialias: None,
                focal_blur: None,
                motion_blur: None,
            },
            reload_error: None,
            watch: None,
            resize: None,
            share: None,
        })
    }

    #[test]
    fn batch_grows_while_draws_keep_up_and_halves_when_they_lag() {
        // Up to a quarter over the refresh interval still counts as keeping up.
        assert_eq!(next_batch(1, REFRESH_MS), 2);
        assert_eq!(next_batch(5, 1.25 * REFRESH_MS), 6);
        assert_eq!(next_batch(6, 1.3 * REFRESH_MS), 3);
        assert_eq!(next_batch(7, 3.0 * REFRESH_MS), 3);
        assert_eq!(next_batch(1, 3.0 * REFRESH_MS), 1);
    }

    #[test]
    fn a_draw_takes_the_batch_but_never_passes_the_target() -> Result<(), String> {
        let mut r = renderer()?;
        r.batch = 3;
        assert_eq!(r.frames_to_draw(), 3);
        r.acc.add_frame(&r.gpu, &r.scene.info)?;
        r.acc.add_frame(&r.gpu, &r.scene.info)?;
        // Target 4, 2 frames in: only 2 more.
        assert_eq!(r.frames_to_draw(), 2);
        r.set_target(UNLIMITED);
        assert_eq!(r.frames_to_draw(), 3);
        Ok(())
    }

    #[test]
    fn tick_shares_the_draw_time_among_its_frames_and_adapts_the_batch() -> Result<(), String> {
        let mut r = renderer()?;
        r.batch = 4;
        let start = r.last_frame;
        assert_eq!(r.tick(start + Duration::from_millis(16), 4), 4.0);
        assert_eq!((r.batch, r.elapsed), (5, Duration::from_millis(16)));
        r.tick(start + Duration::from_millis(56), 5);
        assert_eq!((r.ms, r.batch), (8.0, 2));
        Ok(())
    }

    #[test]
    fn status_shows_the_batch() {
        let text = status("a.json", 8, 256, 1.0, (7.14, 2), (3, 4), TRAVERSAL_BVH);
        assert!(text.contains("· 7.1 ms/frame ×2 ·"));
    }

    #[test]
    fn stopwatch_sums_frame_times_and_restarts_with_the_image() -> Result<(), String> {
        let mut r = renderer()?;
        let start = r.last_frame;
        assert_eq!(r.tick(start + Duration::from_millis(30), 1), 30.0);
        r.tick(start + Duration::from_millis(50), 1);
        assert_eq!((r.ms, r.elapsed), (20.0, Duration::from_millis(50)));
        // Raising the target carries the stopwatch on, without the pause.
        r.set_target(8);
        let resumed = r.last_frame;
        r.tick(resumed + Duration::from_millis(10), 1);
        assert_eq!(r.elapsed, Duration::from_millis(60));
        // A new image starts from zero.
        r.restart();
        assert_eq!(r.elapsed, Duration::ZERO);
        r.tick(r.last_frame + Duration::from_millis(10), 1);
        r.set_traversal(TRAVERSAL_METAL);
        assert_eq!(r.elapsed, Duration::ZERO);
        Ok(())
    }

    #[test]
    fn status_shows_the_stopwatch() {
        let text = status("a.json", 256, 256, 72.24, (8.3, 1), (3, 4), TRAVERSAL_BVH);
        assert!(text.contains("· 256/256 samples · 72.2 s · 8.3 ms/frame ×1 ·"));
    }

    /// A renderer for `scenes/chapter7.json` at 40x20, as `--window` builds it.
    fn chapter7() -> Result<(Renderer, SceneDescription), String> {
        let desc = SceneDescription::from_path("scenes/chapter7.json").map_err(|e| e.0)?;
        let (world, camera) = desc.build(Path::new("scenes")).map_err(|e| e.0)?;
        let camera = resized(&camera, 40, 20);
        let gpu = Gpu::new()?;
        let scene = GpuScene::from_world(&world, &camera)?;
        let acc = gpu.accumulator(&scene)?;
        let orbit = Orbit::of(&desc.camera);
        let mut r = Renderer {
            gpu,
            scene,
            acc,
            camera,
            home: (camera, orbit.clone()),
            orbit,
            moved: false,
            target: 4,
            last_frame: Instant::now(),
            ms: 0.0,
            elapsed: Duration::ZERO,
            batch: 1,
            camera_desc: desc.camera.clone(),
            reload_error: None,
            watch: None,
            resize: None,
            share: None,
        };
        r.acc.add_frame(&r.gpu, &r.scene.info)?;
        r.acc.add_frame(&r.gpu, &r.scene.info)?;
        Ok((r, desc))
    }

    #[test]
    fn reload_with_the_same_camera_keeps_the_view() -> Result<(), String> {
        let (mut r, mut desc) = chapter7()?;
        r.orbit.drag(0.3, 0.0);
        r.moved = true;
        r.restart();
        let view = r.camera.transform();
        desc.shapes[3].material.get_or_insert_default().color = Some([1.0, 0.0, 0.0]);
        r.reload(&desc, Path::new("scenes"))?;
        assert_eq!(r.acc.frames(), 0);
        assert_eq!(r.elapsed, Duration::ZERO);
        assert_eq!((r.camera.hsize, r.camera.vsize), (40, 20));
        assert!(r.camera.transform().approx_eq(view));
        assert!(r.moved);
        Ok(())
    }

    #[test]
    fn same_view_ignores_the_size() -> Result<(), String> {
        let (_, desc) = chapter7()?;
        let base = desc.camera;
        let mut resized = base.clone();
        resized.width = 999;
        resized.height = 3;
        assert!(same_view(&base, &resized));
        let mut moved = base.clone();
        moved.from = [0.0, 3.0, -8.0];
        assert!(!same_view(&base, &moved));
        let mut zoomed = base.clone();
        zoomed.field_of_view += 0.1;
        assert!(!same_view(&base, &zoomed));
        let mut blurred = base.clone();
        blurred.focal_blur = Some(rtc_rust::scene::FocalBlur {
            aperture: 0.1,
            focal_distance: 5.0,
            samples: 4,
        });
        assert!(!same_view(&base, &blurred));
        Ok(())
    }

    #[test]
    fn reload_with_only_a_new_size_keeps_the_view() -> Result<(), String> {
        let (mut r, mut desc) = chapter7()?;
        r.orbit.drag(0.3, 0.0);
        r.moved = true;
        r.restart();
        let view = r.camera.transform();
        desc.camera.width = 77;
        r.reload(&desc, Path::new("scenes"))?;
        assert_eq!((r.camera.hsize, r.camera.vsize), (40, 20));
        assert!(r.camera.transform().approx_eq(view));
        assert!(r.moved);
        assert_eq!(r.camera_desc, desc.camera);
        Ok(())
    }

    #[test]
    fn reload_with_a_new_size_asks_for_that_window_size() -> Result<(), String> {
        let (mut r, mut desc) = chapter7()?;
        r.orbit.drag(0.3, 0.0);
        r.moved = true;
        r.restart();
        let view = r.camera.transform();
        desc.camera.width = 77;
        r.reload(&desc, Path::new("scenes"))?;
        assert_eq!(r.resize, Some((77, desc.camera.height)));
        assert!(r.camera.transform().approx_eq(view));
        assert!(r.moved);
        Ok(())
    }

    #[test]
    fn reload_with_the_same_size_asks_for_no_resize() -> Result<(), String> {
        let (mut r, mut desc) = chapter7()?;
        desc.camera.focal_blur = None;
        desc.camera.from = [0.0, 3.0, -8.0];
        r.reload(&desc, Path::new("scenes"))?;
        assert_eq!(r.resize, None);
        Ok(())
    }

    #[test]
    fn parse_at_reads_two_numbers() {
        assert_eq!(parse_at("900,0"), Some((900.0, 0.0)));
        assert_eq!(parse_at("-1440,12.5"), Some((-1440.0, 12.5)));
    }

    #[test]
    fn parse_at_rejects_bad_values() {
        for bad in ["5", "a,b", "1,2,3", "inf,0", "NaN,0", ""] {
            assert_eq!(parse_at(bad), None, "{bad}");
        }
    }

    /// A 1440x900 main screen with a 25-point menu bar.
    fn main_screen() -> (NSRect, NSRect) {
        (
            NSRect::new(NSPoint::new(0.0, 0.0), NSSize::new(1440.0, 900.0)),
            NSRect::new(NSPoint::new(0.0, 0.0), NSSize::new(1440.0, 875.0)),
        )
    }

    #[test]
    fn place_puts_the_window_right_of_the_point() {
        assert_eq!(
            place((900.0, 30.0), &[main_screen()]),
            Some((NSPoint::new(900.0, 870.0), NSSize::new(540.0, 870.0)))
        );
    }

    #[test]
    fn place_moves_a_point_under_the_menu_bar_down() {
        assert_eq!(
            place((900.0, 0.0), &[main_screen()]),
            Some((NSPoint::new(900.0, 875.0), NSSize::new(540.0, 875.0)))
        );
    }

    #[test]
    fn place_finds_the_screen_holding_the_point() {
        let left = NSRect::new(NSPoint::new(-1280.0, 0.0), NSSize::new(1280.0, 1024.0));
        assert_eq!(
            place((-380.0, 0.0), &[main_screen(), (left, left)]),
            Some((NSPoint::new(-380.0, 900.0), NSSize::new(380.0, 900.0)))
        );
    }

    #[test]
    fn place_gives_none_without_200_points_of_room() {
        assert_eq!(place((1300.0, 30.0), &[main_screen()]), None);
        assert!(place((1240.0, 30.0), &[main_screen()]).is_some());
    }

    #[test]
    fn place_gives_none_off_every_screen() {
        assert_eq!(place((5000.0, 30.0), &[main_screen()]), None);
    }

    #[test]
    fn place_measures_from_the_screen_at_the_origin() {
        // The key screen comes first in no particular order; the y flip
        // uses the screen at (0, 0), whatever its place in the list.
        let left = NSRect::new(NSPoint::new(-1280.0, 0.0), NSSize::new(1280.0, 1024.0));
        assert_eq!(
            place((900.0, 30.0), &[(left, left), main_screen()]),
            Some((NSPoint::new(900.0, 870.0), NSSize::new(540.0, 870.0)))
        );
    }

    #[test]
    fn place_gives_none_without_200_points_of_height() {
        assert_eq!(place((900.0, 800.0), &[main_screen()]), None);
        assert!(place((900.0, 700.0), &[main_screen()]).is_some());
    }

    #[test]
    fn reload_room_falls_back_to_the_whole_screen_near_an_edge() {
        let visible = main_screen().1;
        assert_eq!(
            reload_room(NSPoint::new(100.0, 800.0), visible),
            NSSize::new(1340.0, 800.0)
        );
        for near in [
            NSPoint::new(1300.0, 800.0),
            NSPoint::new(100.0, 150.0),
            NSPoint::new(100.0, -20.0),
        ] {
            assert_eq!(reload_room(near, visible), visible.size, "{near:?}");
        }
    }

    #[test]
    fn on_screen_keeps_a_top_left_that_fits() {
        let visible = main_screen().1;
        let size = NSSize::new(400.0, 300.0);
        assert_eq!(
            on_screen(NSPoint::new(100.0, 800.0), size, visible),
            NSPoint::new(100.0, 800.0)
        );
    }

    #[test]
    fn on_screen_moves_a_window_back_from_the_right_and_bottom() {
        let visible = main_screen().1;
        let size = NSSize::new(400.0, 300.0);
        assert_eq!(
            on_screen(NSPoint::new(1300.0, 150.0), size, visible),
            NSPoint::new(1040.0, 300.0)
        );
        // Larger than the screen: its top-left corner stays visible.
        assert_eq!(
            on_screen(
                NSPoint::new(1300.0, 150.0),
                NSSize::new(2000.0, 1000.0),
                visible
            ),
            NSPoint::new(0.0, 875.0)
        );
    }

    #[test]
    fn room_from_measures_to_the_visible_right_and_bottom() {
        assert_eq!(
            room_from(NSPoint::new(100.0, 800.0), main_screen().1),
            NSSize::new(1340.0, 800.0)
        );
    }

    #[test]
    fn window_size_fits_the_room_and_never_grows() {
        let size = |room| window_size(800, 600, 2.0, room);
        assert_eq!(size(None), NSSize::new(400.0, 300.0));
        assert_eq!(
            size(Some(NSSize::new(200.0, 300.0))),
            NSSize::new(200.0, 150.0)
        );
        assert_eq!(
            size(Some(NSSize::new(4000.0, 4000.0))),
            NSSize::new(400.0, 300.0)
        );
    }

    #[test]
    fn reload_with_a_new_camera_takes_it_as_home() -> Result<(), String> {
        let (mut r, mut desc) = chapter7()?;
        r.moved = true;
        desc.camera.from = [0.0, 3.0, -8.0];
        r.reload(&desc, Path::new("scenes"))?;
        let view = rtc_rust::world::view_transform(
            rtc_rust::tuple::Tuple::point(0.0, 3.0, -8.0),
            rtc_rust::tuple::Tuple::point(0.0, 1.0, 0.0),
            rtc_rust::tuple::Tuple::vector(0.0, 1.0, 0.0),
        );
        assert_eq!((r.camera.hsize, r.camera.vsize), (40, 20));
        assert!(r.camera.transform().approx_eq(view));
        assert!(r.home.0.transform().approx_eq(view));
        assert!(!r.moved);
        assert_eq!(r.camera_desc, desc.camera);
        Ok(())
    }

    #[test]
    fn a_broken_reload_keeps_the_old_scene() -> Result<(), String> {
        let (mut r, desc) = chapter7()?;
        let mut broken = desc.clone();
        broken.light = None;
        broken.lights = Some(vec![]);
        broken.camera.from = [0.0, 3.0, -8.0];
        let err = r
            .reload(&broken, Path::new("scenes"))
            .err()
            .ok_or("expected an error")?;
        assert!(err.contains("lights"), "{err}");
        assert_eq!(r.reload_error.as_deref(), Some(err.as_str()));
        assert_eq!(r.acc.frames(), 2);
        assert_eq!(r.camera_desc, desc.camera);
        r.reload(&desc, Path::new("scenes"))?;
        assert_eq!(r.reload_error, None);
        Ok(())
    }

    #[test]
    fn an_unsupported_reload_keeps_the_old_scene() -> Result<(), String> {
        let (mut r, mut desc) = chapter7()?;
        desc.light = Some(LightDescription::Spot {
            position: [-10.0, 10.0, -10.0],
            target: [0.0, 0.0, 0.0],
            intensity: [1.0, 1.0, 1.0],
            inner_angle: 0.2,
            outer_angle: 0.4,
        });
        let err = r
            .reload(&desc, Path::new("scenes"))
            .err()
            .ok_or("expected an error")?;
        assert!(err.contains("spot"), "{err}");
        assert_eq!(r.acc.frames(), 2);
        Ok(())
    }

    #[test]
    fn reload_keeps_the_traversal() -> Result<(), String> {
        let (mut r, desc) = chapter7()?;
        assert!(r.set_traversal(TRAVERSAL_METAL));
        r.reload(&desc, Path::new("scenes"))?;
        assert_eq!(r.scene.info.counts.traversal, TRAVERSAL_METAL);
        r.acc.add_frame(&r.gpu, &r.scene.info)?;
        Ok(())
    }

    #[test]
    fn doubling_and_halving_the_target() {
        assert_eq!(doubled(256), 512);
        assert_eq!(doubled(1 << 31), 1 << 31);
        assert_eq!(doubled(UNLIMITED), UNLIMITED);
        assert_eq!(halved(256), 128);
        assert_eq!(halved(1), 1);
        // From Unlimited, halving lands on the largest preset.
        assert_eq!(halved(UNLIMITED), 4096);
    }

    #[test]
    fn set_target_resumes_only_below_it() -> Result<(), String> {
        let mut r = renderer()?;
        for _ in 0..3 {
            r.acc.add_frame(&r.gpu, &r.scene.info)?;
        }
        // Raising the target carries on from the frames already summed.
        assert!(r.set_target(8));
        assert_eq!((r.target, r.acc.frames()), (8, 3));
        // Lowering it below them stops at the current image.
        assert!(!r.set_target(2));
        assert_eq!((r.target, r.acc.frames()), (2, 3));
        assert!(r.set_target(UNLIMITED));
        Ok(())
    }

    #[test]
    fn status_shows_unlimited_as_infinity() {
        let text = status(
            "a.json",
            37,
            UNLIMITED,
            0.0,
            (1.0, 1),
            (3, 4),
            TRAVERSAL_BVH,
        );
        assert!(text.contains("· 37/∞ samples ·"));
    }

    #[test]
    fn set_traversal_rebuilds_the_accumulator() -> Result<(), String> {
        let mut r = renderer()?;
        r.acc.add_frame(&r.gpu, &r.scene.info)?;
        assert!(r.set_traversal(TRAVERSAL_METAL));
        assert_eq!(r.scene.info.counts.traversal, TRAVERSAL_METAL);
        assert_eq!(r.acc.frames(), 0);
        // The new accumulator takes Metal frames; the old one refused them.
        r.acc.add_frame(&r.gpu, &r.scene.info)?;
        assert!(!r.set_traversal(TRAVERSAL_METAL));
        assert!(r.set_traversal(TRAVERSAL_BVH));
        r.acc.add_frame(&r.gpu, &r.scene.info)?;
        Ok(())
    }

    #[test]
    fn status_shows_progress_size_and_traversal() {
        let text = status(
            "cornell.json",
            37,
            256,
            4.56,
            (12.34, 3),
            (3456, 2234),
            TRAVERSAL_BVH,
        );
        assert_eq!(
            text,
            "cornell.json · 37/256 samples · 4.6 s · 12.3 ms/frame ×3 · 3456×2234 px · bvh"
        );
        assert!(
            status("a.json", 1, 2, 0.0, (1.0, 1), (3, 4), TRAVERSAL_METAL).ends_with("· metal")
        );
    }

    #[test]
    fn resized_camera_keeps_everything_but_the_size() {
        let mut camera = Camera::new(800, 600, 0.68).with_focal_blur(0.1, 3.5, 8);
        camera.set_transform(&Matrix4::translation(1.0, 2.0, 3.0));
        let wide = resized(&camera, 1600, 600);
        assert_eq!((wide.hsize, wide.vsize), (1600, 600));
        assert_eq!(wide.field_of_view, camera.field_of_view);
        assert!(wide.transform().approx_eq(camera.transform()));
        assert_eq!(
            (wide.aperture, wide.focal_distance, wide.dof_samples),
            (0.1, 3.5, 8)
        );
        // The field of view spans the longer side: a wider window shows more
        // to the sides, at the same height.
        assert_eq!(wide.half_width, camera.half_width);
        assert!(wide.half_height < camera.half_height);
        let fresh = Camera::new(1600, 600, 0.68);
        assert_eq!(wide.pixel_size, fresh.pixel_size);
    }

    #[test]
    fn image_name_is_the_scene_stem_with_ppm() {
        assert_eq!(
            image_name(Path::new("scenes/chapter7.json")),
            "chapter7.ppm"
        );
        assert_eq!(image_name(Path::new("/tmp/rtc-edit-1/a.b.json")), "a.b.ppm");
    }

    #[test]
    fn offered_name_is_the_ask_name_or_the_out_file_name() {
        assert_eq!(offered_name(&Save::Ask("a.ppm".into())), "a.ppm");
        assert_eq!(offered_name(&Save::To("dir/b.ppm".into())), "b.ppm");
    }

    #[test]
    fn a_given_path_is_used_without_asking() {
        let save = Save::To("x.ppm".to_string());
        assert_eq!(
            choose(&save, 3, |_| panic!("asked")),
            Quit::Save(PathBuf::from("x.ppm"))
        );
    }

    #[test]
    fn asking_returns_the_answer() {
        let save = Save::Ask("a.ppm".to_string());
        let answer = choose(&save, 3, |name| {
            assert_eq!(name, "a.ppm");
            Quit::Save(PathBuf::from("/d/b.ppm"))
        });
        assert_eq!(answer, Quit::Save(PathBuf::from("/d/b.ppm")));
        assert_eq!(choose(&save, 3, |_| Quit::Discard), Quit::Discard);
        assert_eq!(choose(&save, 3, |_| Quit::Cancel), Quit::Cancel);
    }

    #[test]
    fn nothing_is_saved_before_a_frame_is_drawn() {
        for save in [
            Save::To("x.ppm".to_string()),
            Save::Ask("a.ppm".to_string()),
        ] {
            assert_eq!(choose(&save, 0, |_| panic!("asked")), Quit::Discard);
        }
    }

    /// Runs `answer` with the alert's buttons in turn and the panel's paths
    /// in turn, panicking if either is asked once too often.
    fn answer_with(buttons: &[Button], paths: &[Option<&str>]) -> Quit {
        let (mut buttons, mut paths) = (buttons.iter(), paths.iter());
        answer(
            || *buttons.next().expect("alert asked too often"),
            || {
                paths
                    .next()
                    .expect("panel asked too often")
                    .map(PathBuf::from)
            },
        )
    }

    #[test]
    fn the_alert_buttons_give_their_answers() {
        assert_eq!(answer_with(&[Button::Discard], &[]), Quit::Discard);
        assert_eq!(answer_with(&[Button::Cancel], &[]), Quit::Cancel);
        assert_eq!(
            answer_with(&[Button::Save], &[Some("/d/a.ppm")]),
            Quit::Save(PathBuf::from("/d/a.ppm"))
        );
    }

    #[test]
    fn cancelling_the_save_panel_asks_again() {
        assert_eq!(
            answer_with(&[Button::Save, Button::Discard], &[None]),
            Quit::Discard
        );
        assert_eq!(
            answer_with(&[Button::Save, Button::Save], &[None, Some("/d/b.ppm")]),
            Quit::Save(PathBuf::from("/d/b.ppm"))
        );
    }

    #[test]
    fn cancel_is_offered_unless_the_editor_asked_to_quit() {
        assert!(offers_cancel(false));
        assert!(!offers_cancel(true));
    }
}
