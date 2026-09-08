//! `rtc` — render a JSON scene description to a PPM image.
//!
//! ```text
//! rtc <scene.json> [-o out.ppm] [--orbit <frames>]
//! ```
//!
//! The scene format is documented in `README.md`; `scenes/cover.json` is a
//! complete example. Without `-o`, the output lands next to the current
//! working directory as `<scene-stem>.ppm`.
//!
//! `--orbit <frames>` renders a film instead of a still: the camera circles
//! the scene's look-at point once (rotating about the camera's `up` axis),
//! writing `<stem>_0000.ppm` … `<stem>_{frames-1}.ppm`. Frame `frames` would
//! coincide with frame 0, so the sequence loops seamlessly. The frames are
//! merged into a video with ffmpeg; `rtc` prints the command afterwards.

use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Instant;

use rtc_rust::scene::SceneDescription;
use rtc_rust::tuple::Tuple;
use rtc_rust::world::{orbit_from, view_transform};

mod help;

fn usage() {
    eprintln!("usage: rtc <scene.json> [-o out.ppm] [--orbit <frames>]");
    eprintln!("       rtc help [topic]      scene-format reference ('rtc help' lists topics)");
    eprintln!();
    eprintln!(
        "  -o, --output <file>   output PPM (default: <scene-stem>.ppm in the working directory)"
    );
    eprintln!("      --orbit <frames>  circle the look-at point once, writing <stem>_0000.ppm …");
    eprintln!(
        "                        (10 s at 24 fps = 240 frames); prints the ffmpeg merge command"
    );
}

/// Numbered frame path next to `out`: `frames/dice.ppm` + 7 -> `frames/dice_0007.ppm`.
fn frame_path(out: &Path, index: usize) -> PathBuf {
    let stem = out
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "render".into());
    out.with_file_name(format!("{stem}_{index:04}.ppm"))
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();

    if args.first().map(String::as_str) == Some("help") {
        return match help::print(args.get(1).map(String::as_str)) {
            true => ExitCode::SUCCESS,
            false => ExitCode::FAILURE,
        };
    }
    let mut scene_path: Option<PathBuf> = None;
    let mut out_path: Option<PathBuf> = None;
    let mut orbit_frames: Option<usize> = None;

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
            "--orbit" => {
                i += 1;
                match args.get(i).and_then(|n| n.parse::<usize>().ok()) {
                    Some(n) if n >= 1 => orbit_frames = Some(n),
                    _ => {
                        eprintln!("rtc: --orbit needs a frame count of 1 or more");
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

    let base_dir = scene_path.parent().unwrap_or(Path::new("."));
    let desc = match SceneDescription::from_path(&scene_path) {
        Ok(desc) => desc,
        Err(e) => {
            eprintln!("rtc: {e}");
            return ExitCode::FAILURE;
        }
    };
    let (world, mut camera) = match desc.build(base_dir) {
        Ok(built) => built,
        Err(e) => {
            eprintln!("rtc: {e}");
            return ExitCode::FAILURE;
        }
    };

    let Some(frames) = orbit_frames else {
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
        return ExitCode::SUCCESS;
    };

    if let Some(dir) = out_path.parent().filter(|d| !d.as_os_str().is_empty())
        && let Err(e) = std::fs::create_dir_all(dir)
    {
        eprintln!("rtc: creating {}: {e}", dir.display());
        return ExitCode::FAILURE;
    }

    let [fx, fy, fz] = desc.camera.from;
    let [tx, ty, tz] = desc.camera.to;
    let [ux, uy, uz] = desc.camera.up;
    let (from, to, up) = (
        Tuple::point(fx, fy, fz),
        Tuple::point(tx, ty, tz),
        Tuple::vector(ux, uy, uz),
    );

    let total = Instant::now();
    for k in 0..frames {
        let angle = 2.0 * std::f64::consts::PI * k as f64 / frames as f64;
        camera.set_transform(&view_transform(orbit_from(from, to, up, angle), to, up));

        let start = Instant::now();
        let canvas = world.render_parallel(camera);
        let elapsed = start.elapsed();

        let path = frame_path(&out_path, k);
        if let Err(e) = canvas.write_ppm(&path) {
            eprintln!("rtc: writing {}: {e}", path.display());
            return ExitCode::FAILURE;
        }
        println!(
            "frame {}/{} rendered in {:.3} s -> {}",
            k + 1,
            frames,
            elapsed.as_secs_f64(),
            path.display()
        );
    }

    let pattern = frame_path(&out_path, 0)
        .display()
        .to_string()
        .replace("_0000.ppm", "_%04d.ppm");
    let video = out_path.with_extension("mp4");
    println!(
        "rendered {} frames of {}x{} in {:.3} s; merge them with:",
        frames,
        camera.hsize,
        camera.vsize,
        total.elapsed().as_secs_f64()
    );
    println!(
        "ffmpeg -framerate 24 -i {pattern} -vf \"pad=ceil(iw/2)*2:ceil(ih/2)*2\" -c:v libx264 -pix_fmt yuv420p {}",
        video.display()
    );
    ExitCode::SUCCESS
}
