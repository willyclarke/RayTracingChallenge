//! `rtc` — render a JSON scene description to a PPM image.
//!
//! ```text
//! rtc <scene.json> [-o out.ppm]
//! ```
//!
//! The scene format is documented in `README.md`; `scenes/cover.json` is a
//! complete example. Without `-o`, the output lands next to the current
//! working directory as `<scene-stem>.ppm`.

use std::path::PathBuf;
use std::process::ExitCode;
use std::time::Instant;

use rtc_rust::scene;

fn usage() {
    eprintln!("usage: rtc <scene.json> [-o out.ppm]");
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut scene_path: Option<PathBuf> = None;
    let mut out_path: Option<PathBuf> = None;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "--help" => {
                usage();
                return ExitCode::SUCCESS;
            }
            "-o" | "--output" => {
                i += 1;
                match args.get(i) {
                    Some(p) => out_path = Some(PathBuf::from(p)),
                    None => {
                        usage();
                        return ExitCode::FAILURE;
                    }
                }
            }
            p if scene_path.is_none() && !p.starts_with('-') => {
                scene_path = Some(PathBuf::from(p));
            }
            _ => {
                usage();
                return ExitCode::FAILURE;
            }
        }
        i += 1;
    }

    let Some(scene_path) = scene_path else {
        usage();
        return ExitCode::FAILURE;
    };
    let out_path = out_path.unwrap_or_else(|| {
        let stem = scene_path
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "render".into());
        PathBuf::from(format!("{stem}.ppm"))
    });

    let (world, camera) = match scene::load(&scene_path) {
        Ok(built) => built,
        Err(e) => {
            eprintln!("rtc: {e}");
            return ExitCode::FAILURE;
        }
    };

    let start = Instant::now();
    let canvas = world.render_parallel(camera);
    let elapsed = start.elapsed();

    if let Err(e) = canvas.write_ppm(&out_path) {
        eprintln!("rtc: writing {}: {e}", out_path.display());
        return ExitCode::FAILURE;
    }
    println!(
        "rendered {}x{} in {:.3} s -> {}",
        camera.hsize,
        camera.vsize,
        elapsed.as_secs_f64(),
        out_path.display()
    );
    ExitCode::SUCCESS
}
