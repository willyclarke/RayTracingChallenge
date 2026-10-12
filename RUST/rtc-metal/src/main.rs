//! rtc-metal: render on the GPU via Metal, write PPM through rtc-rust's Canvas.
//!
//! Usage: rtc-metal [scene.json] [-o out.ppm] [--samples N] [--window]
//!
//! With a scene, renders it (scene format: rtc-rust's README); without one,
//! writes the milestone 1 test gradient. `--samples N` averages N sampled
//! frames (random light samples, jittered pixels) instead of rendering
//! rtc-rust's deterministic image.

use rtc_metal::gpu::Gpu;
use rtc_metal::scene::GpuScene;
use rtc_rust::canvas::Canvas;
use rtc_rust::tuple::Tuple;
use std::time::Instant;

const GRADIENT_WIDTH: usize = 400;
const GRADIENT_HEIGHT: usize = 300;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut out = None;
    let mut scene_path = None;
    let mut samples = None;
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
            a if a.starts_with('-') => usage(),
            a => scene_path = Some(a.to_string()),
        }
        i += 1;
    }

    let gpu = Gpu::new().unwrap_or_else(|e| fail(&e));
    let (pixels, width, height) = match &scene_path {
        Some(path) => {
            let (world, camera) =
                rtc_rust::scene::load(path).unwrap_or_else(|e| fail(&e.to_string()));
            let scene = GpuScene::from_world(&world, &camera).unwrap_or_else(|e| fail(&e));
            let start = Instant::now();
            let pixels = match samples {
                Some(n) => accumulate(&gpu, &scene, n).unwrap_or_else(|e| fail(&e)),
                None => gpu.render(&scene).unwrap_or_else(|e| fail(&e)),
            };
            eprintln!("GPU render: {:.1} ms", start.elapsed().as_secs_f64() * 1e3);
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

    let out = out.unwrap_or_else(|| match &scene_path {
        Some(path) => format!(
            "{}.ppm",
            std::path::Path::new(path)
                .file_stem()
                .unwrap()
                .to_string_lossy()
        ),
        None => "gradient.ppm".to_string(),
    });
    write_ppm(pixels, width, height, &out);
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
fn write_ppm(pixels: Vec<[f32; 4]>, width: usize, height: usize, out: &str) {
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
        .unwrap_or_else(|e| fail(&e.to_string()));
    println!("wrote {out} ({width}x{height})");
}

fn usage() -> ! {
    eprintln!("usage: rtc-metal [scene.json] [-o out.ppm] [--samples N] [--window]");
    eprintln!("  --samples N  average N >= 1 sampled frames (soft shadows, anti-aliasing)");
    eprintln!("  the scene's antialias.n is ignored on the GPU");
    std::process::exit(2);
}

fn fail(msg: &str) -> ! {
    eprintln!("rtc-metal: {msg}");
    std::process::exit(1);
}
