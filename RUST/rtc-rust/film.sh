#!/bin/sh
# Render an orbit film: the camera circles the scene's look-at point once.
#
#   ./film.sh <scene.json> [frames=240] [fps=24]     # or positional: film.sh scene.json 240 24
#
# Frames land in frames/<stem>/<stem>_NNNN.ppm and are merged into <stem>.mp4
# in the working directory with ffmpeg (240 frames at 24 fps = 10 s).
set -eu

usage="usage: film.sh <scene.json> [frames=240] [fps=24]"
scene=${1:?$usage}
shift
frames=240
fps=24
positional=0
for arg in "$@"; do
    case $arg in
        frames=*) frames=${arg#frames=} ;;
        fps=*) fps=${arg#fps=} ;;
        *=*) echo "film.sh: unknown option '$arg'" >&2; echo "$usage" >&2; exit 1 ;;
        *) positional=$((positional + 1))
           case $positional in
               1) frames=$arg ;;
               2) fps=$arg ;;
               *) echo "film.sh: too many arguments" >&2; echo "$usage" >&2; exit 1 ;;
           esac ;;
    esac
done
for v in "$frames" "$fps"; do
    case $v in
        ''|*[!0-9]*|0) echo "film.sh: frames and fps must be whole numbers of 1 or more, got '$v'" >&2
                       echo "$usage" >&2; exit 1 ;;
    esac
done
stem=$(basename "$scene" .json)

cargo build --release
rm -rf "frames/$stem"   # stale frames from a longer run would otherwise end up in the video
target/release/rtc "$scene" --orbit "$frames" -o "frames/$stem/$stem.ppm"
# HEVC through Apple's hardware encoder: much lighter to decode than H.264 at
# Retina-sized 48 fps frames, which QuickTime and the screen saver stalled on.
# -tag:v hvc1 is the four-char code Apple players require; -g <fps> puts a
# keyframe every second so a looping player restarts and recovers quickly
# (x264's default of 250 left a 10 s film with two); +faststart moves the
# index to the front of the file.
ffmpeg -y -framerate "$fps" -i "frames/$stem/${stem}_%04d.ppm" \
    -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -c:v hevc_videotoolbox -q:v 65 -tag:v hvc1 \
    -g "$fps" -pix_fmt yuv420p -movflags +faststart "$stem.mp4"
echo "wrote $stem.mp4 ($frames frames at $fps fps)"
