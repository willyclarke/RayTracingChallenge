#!/bin/sh
# Build RtcSaver.saver, a macOS screen saver that loops an rtc orbit film.
#
#   screensaver/build.sh [film.mp4] [install]
#
# Defaults to cover.mp4 in the crate root (make one with ./film.sh).
# Writes screensaver/build/RtcSaver.saver; `install` copies it to
# ~/Library/Screen Savers/ and quits System Settings so it reloads the bundle.
# Needs Xcode's command-line tools (swiftc, codesign).
set -eu

here=$(cd "$(dirname "$0")" && pwd)
film=${1:-$here/../cover.mp4}
install=${2:-}
[ -f "$film" ] || { echo "build.sh: film '$film' not found (render one with ./film.sh)" >&2; exit 1; }

saver="$here/build/RtcSaver.saver"
rm -rf "$saver"
mkdir -p "$saver/Contents/MacOS" "$saver/Contents/Resources"

swiftc -O -module-name RtcSaver \
    -target "$(uname -m)-apple-macos14.0" \
    -sdk "$(xcrun --show-sdk-path)" \
    -emit-library -o "$saver/Contents/MacOS/RtcSaver" \
    -framework ScreenSaver -framework AVFoundation \
    "$here/RtcSaverView.swift"
cp "$here/Info.plist" "$saver/Contents/Info.plist"
cp "$film" "$saver/Contents/Resources/film.mp4"
# Thumbnail for the settings tile: the film's first frame, aspect-filled to
# Apple's 90x58 shape (doubled so it stays crisp in the wallpaper pane).
if command -v ffmpeg >/dev/null; then
    for spec in thumbnail.png:180:116 thumbnail@2x.png:360:232; do
        name=${spec%%:*}; size=${spec#*:}; w=${size%:*}; h=${size#*:}
        ffmpeg -v error -y -i "$film" -frames:v 1 \
            -vf "scale=$w:$h:force_original_aspect_ratio=increase,crop=$w:$h" \
            "$saver/Contents/Resources/$name"
    done
else
    echo "build.sh: ffmpeg not found, skipping the thumbnail" >&2
fi
codesign --force --sign - "$saver"
echo "built $saver"

if [ "$install" = install ]; then
    dest="$HOME/Library/Screen Savers"
    mkdir -p "$dest"
    rm -rf "$dest/RtcSaver.saver"
    cp -R "$saver" "$dest/"
    osascript -e 'tell application "System Settings" to quit' >/dev/null 2>&1 || true
    echo "installed to $dest/RtcSaver.saver; pick RtcSaver in System Settings > Wallpaper (screen saver section, under Other)"
else
    echo "install with: $0 $film install"
fi
