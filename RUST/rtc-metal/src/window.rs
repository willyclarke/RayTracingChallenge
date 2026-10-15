//! Live window (`--window`): an `MTKView` that adds one sampled frame per
//! display refresh until `target` frames, showing the running mean. Dragging
//! orbits the camera and scrolling zooms; both restart the accumulation
//! (events with no movement leave it alone). Space returns to the scene's
//! own view. The window can be resized or made full screen: the image is
//! rendered at the view's pixel size, with the scene's field of view, and
//! starts again at each new size.
//! Quitting hands the final image and its size to `on_close`.
//!
//! Everything here runs on the main thread: the view and its delegate share
//! one `Rc<RefCell<Renderer>>`.

use crate::gpu::{Accumulator, Gpu};
use crate::orbit::Orbit;
use crate::scene::GpuScene;
use objc2::rc::Retained;
use objc2::runtime::ProtocolObject;
use objc2::{DefinedClass, MainThreadMarker, MainThreadOnly, Message, define_class, msg_send, sel};
use objc2_app_kit::{
    NSApplication, NSApplicationActivationPolicy, NSApplicationDelegate, NSBackingStoreType,
    NSEvent, NSMenu, NSMenuItem, NSResponder, NSScreen, NSView, NSWindow,
    NSWindowCollectionBehavior, NSWindowStyleMask,
};
use objc2_foundation::{
    NSNotification, NSObject, NSObjectProtocol, NSPoint, NSRect, NSSize, NSString,
};
use objc2_metal::{MTLCommandBuffer, MTLCommandEncoder, MTLCommandQueue, MTLPixelFormat};
use objc2_metal_kit::{MTKView, MTKViewDelegate};
use objc2_quartz_core::CAMetalDrawable;
use rtc_rust::camera::Camera;
use std::cell::{Cell, RefCell};
use std::rc::Rc;
use std::time::Instant;

/// Called on quit with the final image and its width and height.
type OnClose = Box<dyn FnOnce(Vec<[f32; 4]>, usize, usize)>;

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
    target: u32,
    last_frame: Instant,
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
        match self.gpu.accumulator(&self.scene) {
            Ok(acc) => self.acc = acc,
            Err(e) => {
                eprintln!("rtc-metal: {e}");
                return false;
            }
        }
        self.last_frame = Instant::now();
        true
    }
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
    /// live in `Orbit::drag`); scrolling zooms; space resets the view.
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

        #[unsafe(method(keyDown:))]
        fn key_down(&self, event: &NSEvent) {
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

    /// Start drawing frames again after the accumulation restarted.
    fn resume(&self) {
        self.setPaused(false);
        self.setEnableSetNeedsDisplay(false);
    }
}

define_class!(
    // SAFETY: NSObject has no subclassing requirements; `RtcDelegate` does
    // not implement `Drop`.
    #[unsafe(super(NSObject))]
    #[thread_kind = MainThreadOnly]
    #[ivars = Delegate]
    struct RtcDelegate;

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

    unsafe impl NSApplicationDelegate for RtcDelegate {
        #[unsafe(method(applicationShouldTerminateAfterLastWindowClosed:))]
        fn should_terminate_after_last_window_closed(&self, _sender: &NSApplication) -> bool {
            true
        }

        #[unsafe(method(applicationWillTerminate:))]
        fn will_terminate(&self, _notification: &NSNotification) {
            let ivars = self.ivars();
            let r = ivars.renderer.borrow();
            match r.acc.resolve(&r.gpu) {
                Ok(pixels) => {
                    if let Some(on_close) = ivars.on_close.take() {
                        on_close(pixels, r.camera.hsize, r.camera.vsize);
                    }
                }
                Err(e) => eprintln!("rtc-metal: {e}"),
            }
        }
    }
);

struct Delegate {
    renderer: Rc<RefCell<Renderer>>,
    on_close: Cell<Option<OnClose>>,
    title: String,
}

impl RtcDelegate {
    fn new(mtm: MainThreadMarker, ivars: Delegate) -> Retained<Self> {
        let this = Self::alloc(mtm).set_ivars(ivars);
        // SAFETY: NSObject's `init`.
        unsafe { msg_send![super(this), init] }
    }

    /// Add one frame and show the running mean; pause once `target` frames
    /// are in.
    fn draw(&self, view: &MTKView) {
        let ivars = self.ivars();
        let mut r = ivars.renderer.borrow_mut();
        if r.acc.frames() >= r.target {
            drop(r);
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
        r.acc.encode_frame(&encoder, &info);
        r.acc.encode_resolve_texture(&encoder, &drawable.texture());
        encoder.endEncoding();
        commands.presentDrawable(ProtocolObject::from_ref(&*drawable));
        commands.commit();

        let now = Instant::now();
        let ms = (now - r.last_frame).as_secs_f64() * 1e3;
        r.last_frame = now;
        let title = format!(
            "{} — {}/{} samples, {ms:.1} ms/frame",
            ivars.title,
            r.acc.frames(),
            r.target
        );
        drop(r);
        if let Some(window) = view.window() {
            window.setTitle(&NSString::from_str(&title));
        }
    }
}

/// Open a window showing `scene` refining up to `samples` frames, run the
/// app until it quits, then exit. `on_close` gets the final image (row-major
/// RGBA, as `Accumulator::resolve`) if at least one frame was drawn.
pub fn run(
    gpu: Gpu,
    scene: GpuScene,
    camera: Camera,
    orbit: Orbit,
    samples: u32,
    title: String,
    on_close: OnClose,
) -> ! {
    let mtm = MainThreadMarker::new().expect("the window must run on the main thread");
    let (width, height) = (scene.width() as f64, scene.height() as f64);
    let acc = gpu.accumulator(&scene).unwrap_or_else(|e| {
        eprintln!("rtc-metal: {e}");
        std::process::exit(1);
    });
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
    }));

    let app = NSApplication::sharedApplication(mtm);
    app.setActivationPolicy(NSApplicationActivationPolicy::Regular);
    app.setMainMenu(Some(&quit_menu(mtm)));
    let delegate = RtcDelegate::new(
        mtm,
        Delegate {
            renderer: renderer.clone(),
            on_close: Cell::new(Some(on_close)),
            title: title.clone(),
        },
    );
    app.setDelegate(Some(ProtocolObject::from_ref(&*delegate)));

    // One drawable pixel per image pixel, shrunk uniformly so the window,
    // title bar included, fits the screen.
    let style = NSWindowStyleMask::Titled
        | NSWindowStyleMask::Closable
        | NSWindowStyleMask::Miniaturizable
        | NSWindowStyleMask::Resizable;
    let screen = NSScreen::mainScreen(mtm);
    let scale = screen.as_ref().map_or(1.0, |s| s.backingScaleFactor());
    let (mut w, mut h) = (width / scale, height / scale);
    if let Some(screen) = screen {
        let room = NSWindow::contentRectForFrameRect_styleMask(screen.visibleFrame(), style, mtm);
        let fit = (room.size.width / w).min(room.size.height / h).min(1.0);
        (w, h) = (w * fit, h * fit);
    }
    let frame = NSRect::new(NSPoint::new(0.0, 0.0), NSSize::new(w, h));

    let view = RtcView::new(mtm, renderer, frame);
    view.setDelegate(Some(ProtocolObject::from_ref(&*delegate)));
    view.setColorPixelFormat(MTLPixelFormat::BGRA8Unorm);
    view.setFramebufferOnly(false);
    // The drawable follows the view's size in pixels; each change reaches
    // `drawableSizeWillChange`, which re-sizes the camera.
    view.setAutoResizeDrawable(true);
    view.setPreferredFramesPerSecond(60);

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
    window.center();
    window.makeKeyAndOrderFront(None);
    window.makeFirstResponder(Some(&view));
    app.activate(); // needs macOS 14+ (older: activateIgnoringOtherApps)
    app.run();
    std::process::exit(0);
}

/// A main menu holding just the app menu with Quit (Cmd-Q): a bare
/// `NSApplication` has no menu, so Cmd-Q would do nothing.
fn quit_menu(mtm: MainThreadMarker) -> Retained<NSMenu> {
    let quit = NSString::from_str("Quit");
    let q = NSString::from_str("q");
    // SAFETY: `terminate:` is NSApplication's action method, reached through
    // the responder chain.
    let quit_item = unsafe {
        NSMenuItem::initWithTitle_action_keyEquivalent(
            NSMenuItem::alloc(mtm),
            &quit,
            Some(sel!(terminate:)),
            &q,
        )
    };
    let app_menu = NSMenu::new(mtm);
    app_menu.addItem(&quit_item);
    let app_item = NSMenuItem::new(mtm);
    app_item.setSubmenu(Some(&app_menu));
    let menu = NSMenu::new(mtm);
    menu.addItem(&app_item);
    menu
}

#[cfg(test)]
mod tests {
    use super::*;
    use rtc_rust::matrix::Matrix4;

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
}
