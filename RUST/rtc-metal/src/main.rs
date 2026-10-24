//! rtc-metal: render on the GPU via Metal, write PPM through rtc-rust's Canvas.
//!
//! Usage: rtc-metal [scene.json] [-o out.ppm] [--samples N] [--window [--watch] [--at X,Y]
//!                  [--share PATH]] [--traversal bvh|metal]
//!
//! With a scene, renders it (scene format: rtc-rust's README); without one,
//! writes the milestone 1 test gradient. `--samples N` averages N sampled
//! frames (random light samples, jittered pixels) instead of rendering
//! rtc-rust's deterministic image. `--window` shows the sampled frames
//! refining live (default 256, or `--samples N`); drag to orbit the camera,
//! scroll to zoom, space to return to the scene's view, `i` to toggle the
//! status overlay (which times each image until it reaches the target); resizing the window or going full screen renders at the
//! new size. On quit (or SIGTERM) the image, at the size last shown, is
//! written to `-o` or, without it, wherever the user picks after an alert
//! asks whether to save it. With
//! `--watch` the window reloads the scene whenever its file changes (a scene
//! editor, or a text editor, saves it) and refines it again; a file that
//! doesn't load keeps the old image and shows why. `--share PATH` writes the
//! window's camera to a file made by the scene editor (rtc-share) whenever
//! it changes, loads the scenes the editor sends through it, and leaves the
//! focus with the editor at start. `--traversal metal` searches the bounded
//! shapes with Metal's ray tracing API instead of the BVH walk
//! (milestone 10).

use rtc_metal::gpu::Gpu;
use rtc_metal::scene::GpuScene;
use rtc_metal::types::{TRAVERSAL_BVH, TRAVERSAL_METAL};
use rtc_metal::watch::Watch;
use rtc_metal::window::{Save, image_name};
use rtc_rust::canvas::Canvas;
use rtc_rust::scene::SceneDescription;
use rtc_rust::tuple::Tuple;
use std::path::Path;
use std::time::Instant;

const GRADIENT_WIDTH: usize = 400;
const GRADIENT_HEIGHT: usize = 300;
/// Frames the window accumulates without `--samples`.
const WINDOW_SAMPLES: u32 = 256;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut out = None;
    let mut scene_path = None;
    let mut samples = None;
    let mut window = false;
    let mut watch = false;
    let mut at = None;
    let mut share = None;
    let mut traversal = TRAVERSAL_BVH;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "-o" => {
                out = Some(args.get(i + 1).cloned().unwrap_or_else(|| usage()));
                i += 1;
            }
            "--samples" => {
                let n = args.get(i + 1).and_then(|n| n.parse::<u32>().ok());
                samples = Some(n.filter(|&n| n >= 1).unwrap_or_else(|| usage()));
                i += 1;
            }
            "--window" => window = true,
            "--watch" => watch = true,
            "--at" => {
                let xy = args.get(i + 1).and_then(|a| rtc_metal::window::parse_at(a));
                at = Some(xy.unwrap_or_else(|| usage()));
                i += 1;
            }
            "--share" => {
                share = Some(args.get(i + 1).cloned().unwrap_or_else(|| usage()));
                i += 1;
            }
            "--traversal" => {
                traversal = match args.get(i + 1).map(String::as_str) {
                    Some("bvh") => TRAVERSAL_BVH,
                    Some("metal") => TRAVERSAL_METAL,
                    _ => usage(),
                };
                i += 1;
            }
            a if a.starts_with('-') => usage(),
            a => scene_path = Some(a.to_string()),
        }
        i += 1;
    }

    if (watch || at.is_some() || share.is_some()) && !window {
        usage();
    }

    let gpu = Gpu::new().unwrap_or_else(|e| fail(&e));
    if window {
        let path = scene_path.unwrap_or_else(|| usage());
        let save = match out {
            Some(out) => Save::To(out),
            None => Save::Ask(image_name(Path::new(&path))),
        };
        open_window(
            gpu,
            &path,
            samples.unwrap_or(WINDOW_SAMPLES),
            traversal,
            watch,
            at,
            share.as_deref(),
            save,
        );
    }
    let out = out.unwrap_or_else(|| match &scene_path {
        Some(path) => image_name(Path::new(path)),
        None => "gradient.ppm".to_string(),
    });
    let (pixels, width, height) = match &scene_path {
        Some(path) => {
            let (world, camera) =
                rtc_rust::scene::load(path).unwrap_or_else(|e| fail(&e.to_string()));
            let mut scene = GpuScene::from_world(&world, &camera).unwrap_or_else(|e| fail(&e));
            scene.info.counts.traversal = traversal;
            let start = Instant::now();
            let pixels = match samples {
                Some(n) => accumulate(&gpu, &scene, n).unwrap_or_else(|e| fail(&e)),
                None => gpu.render(&scene).unwrap_or_else(|e| fail(&e)),
            };
            eprintln!(
                "GPU render: {:.1} ms ({})",
                start.elapsed().as_secs_f64() * 1e3,
                if traversal == TRAVERSAL_METAL {
                    "metal"
                } else {
                    "bvh"
                }
            );
            (pixels, scene.width(), scene.height())
        }
        None if samples.is_some() => usage(),
        None => {
            let pixels = gpu
                .gradient(GRADIENT_WIDTH, GRADIENT_HEIGHT)
                .unwrap_or_else(|e| fail(&e));
            (pixels, GRADIENT_WIDTH, GRADIENT_HEIGHT)
        }
    };
    write_ppm(pixels, width, height, Path::new(&out)).unwrap_or_else(|e| fail(&e));
}

/// Show the scene at `path` refining live until `samples` frames, reloading
/// it when the file changes if `watch`; on quit, write the image where
/// `save` says.
#[allow(clippy::too_many_arguments)]
fn open_window(
    gpu: Gpu,
    path: &str,
    samples: u32,
    traversal: u32,
    watch: bool,
    at: Option<(f64, f64)>,
    share: Option<&str>,
    save: Save,
) -> ! {
    let path = Path::new(path);
    // Watch before reading, so a save while the scene loads still reloads.
    let watch = watch.then(|| Watch::new(path));
    let desc = SceneDescription::from_path(path).unwrap_or_else(|e| fail(&e.to_string()));
    let (world, camera) = desc
        .build(path.parent().unwrap_or(Path::new(".")))
        .unwrap_or_else(|e| fail(&e.to_string()));
    let mut scene = GpuScene::from_world(&world, &camera).unwrap_or_else(|e| fail(&e));
    scene.info.counts.traversal = traversal;
    let title = path.file_name().unwrap().to_string_lossy().into_owned();
    // Without the file the editor's loads would never arrive: exit, and
    // the editor sees the viewer end.
    let shared = share.map(|p| {
        rtc_metal::loads::listen(Path::new(p)).unwrap_or_else(|e| fail(&e));
        rtc_share::Shared::open(Path::new(p)).unwrap_or_else(|e| fail(&e))
    });
    rtc_metal::window::run(
        gpu,
        scene,
        camera,
        desc.camera.clone(),
        samples,
        title,
        watch,
        at,
        shared,
        save,
        Box::new(write_ppm),
        rtc_metal::window::activates_at_start(share.is_some()),
    )
}

/// `n` sampled frames of `scene`, averaged; progress on stderr at most once
/// a second.
fn accumulate(gpu: &Gpu, scene: &GpuScene, n: u32) -> Result<Vec<[f32; 4]>, String> {
    let mut acc = gpu.accumulator(scene)?;
    let mut shown = Instant::now();
    for done in 1..=n {
        acc.add_frame(gpu, &scene.info)?;
        if shown.elapsed().as_secs_f64() >= 1.0 {
            eprint!("\r{done}/{n} samples");
            shown = Instant::now();
        }
    }
    eprintln!("\r{n}/{n} samples");
    acc.resolve(gpu)
}

/// Write row-major RGBA `pixels` as a PPM through rtc-rust's Canvas.
fn write_ppm(pixels: Vec<[f32; 4]>, width: usize, height: usize, out: &Path) -> Result<(), String> {
    let mut canvas = Canvas::new(width, height);
    for (i, [r, g, b, _]) in pixels.into_iter().enumerate() {
        canvas.write_pixel(
            i % width,
            i / width,
            Tuple::color(r as f64, g as f64, b as f64),
        );
    }
    canvas
        .write_ppm(out)
        .map_err(|e| format!("{}: {e}", out.display()))?;
    println!("wrote {} ({width}x{height})", out.display());
    Ok(())
}

fn usage() -> ! {
    eprintln!(
        "usage: rtc-metal [scene.json] [-o out.ppm] [--samples N] [--window [--watch] [--at X,Y] [--share PATH]] [--traversal bvh|metal]"
    );
    eprintln!("  --samples N  average N >= 1 sampled frames (soft shadows, anti-aliasing)");
    eprintln!(
        "  --window     refine live (default {WINDOW_SAMPLES} samples); drag orbits, scroll zooms,\n               space resets the view,\n               i toggles the status overlay (with seconds to reach the target),\n               resizing re-renders,\n               quit asks whether and where to save the image at the window's
               size (with -o, writes there without asking);
               Samples menu: presets, Unlimited, Cmd-] doubles, Cmd-[ halves"
    );
    eprintln!("  --watch      with --window: reload the scene whenever its file changes");
    eprintln!(
        "  --at X,Y     with --window: the window's top-left, in points right and down from\n               the main screen's top-left, fitted to the room right of it\n               (centred when under 200 points are left)"
    );
    eprintln!(
        "  --share PATH with --window: write the camera to PATH, a file the scene editor made,\n               whenever it changes, and load the scenes it sends; don't take\n               the focus at start; exits when it can't use PATH"
    );
    eprintln!("  --traversal  bvh (default): our BVH walk; metal: Metal's ray tracing API");
    eprintln!("               (in --window, also the Traversal menu: Cmd-1 bvh, Cmd-2 metal)");
    eprintln!("  the scene's antialias.n is ignored on the GPU");
    std::process::exit(2);
}

fn fail(msg: &str) -> ! {
    eprintln!("rtc-metal: {msg}");
    std::process::exit(1);
}
