//! Live window (`--window`): an `MTKView` that adds one sampled frame per
//! display refresh until `target` frames, showing the running mean. Dragging
//! orbits the camera and scrolling zooms; both restart the accumulation
//! (events with no movement leave it alone).
//! Quitting hands the final image to `on_close`.
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
    NSEvent, NSMenu, NSMenuItem, NSResponder, NSScreen, NSView, NSWindow, NSWindowStyleMask,
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

type OnClose = Box<dyn FnOnce(Vec<[f32; 4]>)>;

/// The state the view and its delegate share.
struct Renderer {
    gpu: Gpu,
    scene: GpuScene,
    acc: Accumulator,
    camera: Camera,
    orbit: Orbit,
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
}

define_class!(
    // SAFETY: MTKView has no subclassing requirements; `RtcView` does not
    // implement `Drop`.
    #[unsafe(super(MTKView, NSView, NSResponder, NSObject))]
    #[thread_kind = MainThreadOnly]
    #[ivars = Rc<RefCell<Renderer>>]
    /// The image view. Dragging orbits the camera "grabbing the scene"
    /// (dragging right orbits it left, dragging down raises it; the signs
    /// live in `Orbit::drag`); scrolling zooms.
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
            r.restart();
        }
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
        fn drawable_size_will_change(&self, _view: &MTKView, _size: NSSize) {
            // The drawable size is fixed to the camera's.
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
                        on_close(pixels);
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
        orbit,
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
    let style =
        NSWindowStyleMask::Titled | NSWindowStyleMask::Closable | NSWindowStyleMask::Miniaturizable;
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
    view.setAutoResizeDrawable(false);
    view.setDrawableSize(NSSize::new(width, height));
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
