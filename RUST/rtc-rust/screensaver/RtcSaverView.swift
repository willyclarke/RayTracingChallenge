// A macOS screen saver that loops a film rendered by `rtc --orbit`.
//
// The principal class of the `.saver` bundle built by `build.sh`. The film
// lives inside the bundle as `film.mp4` (the screen saver host is sandboxed,
// so it cannot read from the home directory). One instance runs per display.

import AVFoundation
import ScreenSaver

@objc(RtcSaverView)
final class RtcSaverView: ScreenSaverView {
    private var player: AVQueuePlayer?
    private var looper: AVPlayerLooper?
    private var playerLayer: AVPlayerLayer?

    override init?(frame: NSRect, isPreview: Bool) {
        super.init(frame: frame, isPreview: isPreview)
        setUp()
    }

    required init?(coder: NSCoder) {
        super.init(coder: coder)
        setUp()
    }

    private func setUp() {
        wantsLayer = true
        layer?.backgroundColor = NSColor.black.cgColor

        let bundle = Bundle(for: RtcSaverView.self)
        guard let url = bundle.url(forResource: "film", withExtension: "mp4") else {
            NSLog("RtcSaver: film.mp4 is missing from %@", bundle.bundlePath)
            return
        }

        let queue = AVQueuePlayer()
        queue.isMuted = true
        // AVPlayerLooper re-queues the item so the orbit repeats without a gap.
        looper = AVPlayerLooper(player: queue, templateItem: AVPlayerItem(url: url))

        let videoLayer = AVPlayerLayer(player: queue)
        videoLayer.videoGravity = .resizeAspectFill
        videoLayer.frame = bounds
        layer?.addSublayer(videoLayer)

        player = queue
        playerLayer = videoLayer
    }

    override func layout() {
        super.layout()
        playerLayer?.frame = bounds
    }

    override func startAnimation() {
        super.startAnimation()
        player?.play()
    }

    override func stopAnimation() {
        super.stopAnimation()
        player?.pause()
    }

    // The player layer draws the frames; nothing to do per tick.
    override func animateOneFrame() {}

    override var hasConfigureSheet: Bool { false }
}
