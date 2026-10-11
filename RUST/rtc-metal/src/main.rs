//! rtc-metal: render on the GPU via Metal, write PPM through rtc-rust's Canvas.
//!
//! Usage: rtc-metal [-o out.ppm]

use rtc_metal::gpu::Gpu;
use rtc_rust::canvas::Canvas;
use rtc_rust::tuple::Tuple;

const WIDTH: usize = 400;
const HEIGHT: usize = 300;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let out = match args.iter().position(|a| a == "-o") {
        Some(i) => args.get(i + 1).cloned().unwrap_or_else(|| usage()),
        None => "gradient.ppm".to_string(),
    };

    let gpu = Gpu::new().unwrap_or_else(|e| fail(&e));
    let pixels = gpu.gradient(WIDTH, HEIGHT).unwrap_or_else(|e| fail(&e));

    let mut canvas = Canvas::new(WIDTH, HEIGHT);
    for (i, [r, g, b, _]) in pixels.into_iter().enumerate() {
        let color = Tuple::color(r as f64, g as f64, b as f64);
        canvas.write_pixel(i % WIDTH, i / WIDTH, color);
    }
    canvas
        .write_ppm(&out)
        .unwrap_or_else(|e| fail(&e.to_string()));
    println!("wrote {out} ({WIDTH}x{HEIGHT})");
}

fn usage() -> ! {
    eprintln!("usage: rtc-metal [-o out.ppm]");
    std::process::exit(2);
}

fn fail(msg: &str) -> ! {
    eprintln!("rtc-metal: {msg}");
    std::process::exit(1);
}
