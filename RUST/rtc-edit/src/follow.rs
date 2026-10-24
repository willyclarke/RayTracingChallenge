//! Following the viewer's camera into the scene's: ticking **follow**
//! remembers the camera's `from` (and focal distance), each new view from
//! the viewer is copied in, Keep keeps it, and unticking puts the
//! remembered values back.

use rtc_rust::scene::CameraDescription;
use rtc_share::Snapshot;

/// The camera's `from` and, with focal blur, its focal distance and aperture.
type Cached = ([f64; 3], Option<(f64, f64)>);

/// Whether the scene's camera follows the viewer's, and what to put back.
#[derive(Default)]
pub struct Follower {
    /// While following: the `from` and, with focal blur, the focal distance
    /// and aperture the camera had when following started.
    cache: Option<Cached>,
}

/// To 4 decimals, so the boxes and the file stay readable.
fn rounded(v: f64) -> f64 {
    (v * 1e4).round() / 1e4
}

impl Follower {
    pub fn following(&self) -> bool {
        self.cache.is_some()
    }

    /// Follow from now on, remembering `camera`'s values.
    pub fn start(&mut self, camera: &CameraDescription) {
        let blur = camera
            .focal_blur
            .as_ref()
            .map(|b| (b.focal_distance, b.aperture));
        self.cache = Some((camera.from, blur));
    }

    /// While following, copy `shot`'s `from`, and its focal distance when
    /// `camera` has focal blur, rounded; whether `camera` changed. The
    /// aperture stays, as the viewer keeps its lens when it zooms, so the
    /// scene renders as the viewer shows it.
    pub fn update(&mut self, camera: &mut CameraDescription, shot: &Snapshot) -> bool {
        if !self.following() {
            return false;
        }
        let from = shot.from.map(rounded);
        let mut changed = camera.from != from;
        camera.from = from;
        let distance = rounded(shot.focal_distance);
        if let Some(blur) = &mut camera.focal_blur
            && distance > 0.0
            && blur.focal_distance != distance
        {
            blur.focal_distance = distance;
            changed = true;
        }
        changed
    }

    /// Stop following, keeping the values followed.
    pub fn keep(&mut self) {
        self.cache = None;
    }

    /// Stop following and put back the values from when it started; whether
    /// `camera` changed.
    pub fn untick(&mut self, camera: &mut CameraDescription) -> bool {
        let Some((from, cached)) = self.cache.take() else {
            return false;
        };
        let mut changed = camera.from != from;
        camera.from = from;
        if let (Some(blur), Some((distance, aperture))) = (&mut camera.focal_blur, cached)
            && (blur.focal_distance, blur.aperture) != (distance, aperture)
        {
            (blur.focal_distance, blur.aperture) = (distance, aperture);
            changed = true;
        }
        changed
    }

    /// Stop following without putting anything back: the scene was
    /// replaced or the viewer ended.
    pub fn stop(&mut self) {
        self.cache = None;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rtc_rust::scene::FocalBlur;

    fn camera() -> CameraDescription {
        CameraDescription {
            width: 100,
            height: 50,
            field_of_view: 1.0,
            from: [0.0, 1.5, -5.0],
            to: [0.0, 1.0, 0.0],
            up: [0.0, 1.0, 0.0],
            antialias: None,
            focal_blur: None,
            motion_blur: None,
        }
    }

    fn blurred() -> CameraDescription {
        CameraDescription {
            focal_blur: Some(FocalBlur {
                focal_distance: 5.0,
                aperture: 0.25,
                samples: 32,
            }),
            ..camera()
        }
    }

    fn shot(seq: u32) -> Snapshot {
        Snapshot {
            seq,
            from: [0.81234, 2.1, -4.87],
            to: [0.0, 1.0, 0.0],
            up: [0.0, 1.0, 0.0],
            focal_distance: 6.0,
        }
    }

    #[test]
    fn nothing_follows_until_started() {
        let mut c = camera();
        assert!(!Follower::default().update(&mut c, &shot(1)));
        assert_eq!(c, camera());
    }

    #[test]
    fn following_copies_from_rounded() {
        let (mut c, mut f) = (camera(), Follower::default());
        f.start(&c);
        assert!(f.following());
        assert!(f.update(&mut c, &shot(1)));
        assert_eq!(c.from, [0.8123, 2.1, -4.87]);
    }

    #[test]
    fn an_unchanged_from_is_no_edit() {
        let (mut c, mut f) = (camera(), Follower::default());
        f.start(&c);
        f.update(&mut c, &shot(1));
        assert!(!f.update(&mut c, &shot(2)));
    }

    #[test]
    fn focal_distance_follows_only_with_focal_blur() {
        let (mut c, mut f) = (camera(), Follower::default());
        f.start(&c);
        f.update(&mut c, &shot(1));
        assert_eq!(c.focal_blur, None);

        let (mut c, mut f) = (blurred(), Follower::default());
        f.start(&c);
        f.update(&mut c, &shot(1));
        let blur = c.focal_blur.unwrap();
        assert_eq!(blur.focal_distance, 6.0);
        // The viewer keeps its lens radius as it zooms, so the scene does too.
        assert_eq!(blur.aperture, 0.25);
    }

    #[test]
    fn unticking_restores_from_and_focal_distance() {
        let (mut c, mut f) = (blurred(), Follower::default());
        f.start(&c);
        f.update(&mut c, &shot(1));
        assert!(f.untick(&mut c));
        assert!(!f.following());
        assert_eq!(c, blurred());
    }

    #[test]
    fn unticking_puts_back_the_exact_blur() {
        let (mut c, mut f) = (blurred(), Follower::default());
        f.start(&c);
        for d in [6.0, 7.3, 2.9] {
            f.update(
                &mut c,
                &Snapshot {
                    focal_distance: d,
                    ..shot(1)
                },
            );
        }
        f.untick(&mut c);
        assert_eq!(c.focal_blur, blurred().focal_blur);
    }

    #[test]
    fn keeping_keeps_the_values() {
        let (mut c, mut f) = (camera(), Follower::default());
        f.start(&c);
        f.update(&mut c, &shot(1));
        f.keep();
        assert!(!f.following());
        assert_eq!(c.from, [0.8123, 2.1, -4.87]);
        let mut moved = shot(2);
        moved.from = [9.0, 9.0, 9.0];
        assert!(!f.update(&mut c, &moved));
    }

    #[test]
    fn stopping_keeps_the_values() {
        let (mut c, mut f) = (camera(), Follower::default());
        f.start(&c);
        f.update(&mut c, &shot(1));
        f.stop();
        assert!(!f.following());
        assert_eq!(c.from, [0.8123, 2.1, -4.87]);
    }
}
