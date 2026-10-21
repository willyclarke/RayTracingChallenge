//! Undo and redo: copies of the whole scene, one per step, each with a
//! label for the status line. Changes of one open group (typing in one box,
//! a run of follow updates) are one step.

use crate::fields::Field;
use crate::tree::Node;
use rtc_rust::scene::SceneDescription;

/// Steps kept; older ones are dropped.
const DEPTH: usize = 100;

/// What a change belongs to: changes of the open group merge into one step.
#[derive(Clone, PartialEq, Debug)]
pub enum Group {
    /// Typing in one box, from entering it to leaving it.
    Box(Node, Field),
    /// A run of follow updates.
    Follow,
    /// A step of its own.
    Once,
}

pub struct History {
    /// The scenes before each step, oldest first, with the step's label.
    undo: Vec<(SceneDescription, String)>,
    /// The scenes undone, most recently undone last.
    redo: Vec<(SceneDescription, String)>,
    /// The group the last step is still open for.
    open: Option<Group>,
    /// The scene as last saved or opened.
    saved: SceneDescription,
}

impl History {
    /// No steps; `scene` is the saved one.
    pub fn new(scene: &SceneDescription) -> History {
        History {
            undo: Vec::new(),
            redo: Vec::new(),
            open: None,
            saved: scene.clone(),
        }
    }

    /// Note a change, given the scene `before` it: part of the open step
    /// when `group` is the open group, else a new step labelled `label`.
    pub fn record(&mut self, before: &SceneDescription, group: Group, label: &str) {
        if group != Group::Once && self.open.as_ref() == Some(&group) {
            return;
        }
        self.undo.push((before.clone(), label.to_string()));
        if self.undo.len() > DEPTH {
            self.undo.remove(0);
        }
        self.redo.clear();
        self.open = (group != Group::Once).then_some(group);
    }

    /// End the open group: the next change is a new step.
    pub fn close(&mut self) {
        self.open = None;
    }

    /// Put `scene` back as it was before the last step; its label.
    pub fn undo(&mut self, scene: &mut SceneDescription) -> Option<String> {
        Self::swap(&mut self.undo, &mut self.redo, scene, &mut self.open)
    }

    /// Re-apply the last step undone; its label.
    pub fn redo(&mut self, scene: &mut SceneDescription) -> Option<String> {
        Self::swap(&mut self.redo, &mut self.undo, scene, &mut self.open)
    }

    fn swap(
        from: &mut Vec<(SceneDescription, String)>,
        to: &mut Vec<(SceneDescription, String)>,
        scene: &mut SceneDescription,
        open: &mut Option<Group>,
    ) -> Option<String> {
        *open = None;
        let (other, label) = from.pop()?;
        to.push((std::mem::replace(scene, other), label.clone()));
        Some(label)
    }

    /// Apply `f` to every scene kept (the undo and redo copies and the
    /// saved one), as Save As rewrites asset paths for a new folder.
    pub fn for_each_scene(&mut self, mut f: impl FnMut(&mut SceneDescription)) {
        for (scene, _) in self.undo.iter_mut().chain(self.redo.iter_mut()) {
            f(scene);
        }
        f(&mut self.saved);
    }

    /// `scene` was saved.
    pub fn saved(&mut self, scene: &SceneDescription) {
        self.saved = scene.clone();
    }

    /// Whether `scene` differs from the one last saved or opened.
    pub fn differs(&self, scene: &SceneDescription) -> bool {
        *scene != self.saved
    }

    /// Forget every step; `scene` is the saved one.
    pub fn clear(&mut self, scene: &SceneDescription) {
        *self = History::new(scene);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::edit::new_scene;

    /// `s` with its field of view set to `fov`.
    fn with_fov(s: &SceneDescription, fov: f64) -> SceneDescription {
        let mut s = s.clone();
        s.camera.field_of_view = fov;
        s
    }

    /// Record a change of `scene`'s field of view to `fov` in `group`.
    fn change(h: &mut History, scene: &mut SceneDescription, fov: f64, group: Group) {
        h.record(scene, group, "a");
        *scene = with_fov(scene, fov);
    }

    fn fov_box() -> Group {
        Group::Box(Node::Camera, Field::FieldOfView)
    }

    #[test]
    fn undo_restores_and_redo_reapplies() {
        let mut s = with_fov(&new_scene(), 1.0);
        let mut h = History::new(&s);
        change(&mut h, &mut s, 2.0, Group::Once);
        assert_eq!(h.undo(&mut s).as_deref(), Some("a"));
        assert_eq!(s.camera.field_of_view, 1.0);
        assert_eq!(h.redo(&mut s).as_deref(), Some("a"));
        assert_eq!(s.camera.field_of_view, 2.0);
    }

    #[test]
    fn one_box_group_is_one_step() {
        let mut s = with_fov(&new_scene(), 1.0);
        let mut h = History::new(&s);
        change(&mut h, &mut s, 2.0, fov_box());
        change(&mut h, &mut s, 3.0, fov_box());
        h.undo(&mut s);
        assert_eq!(s.camera.field_of_view, 1.0);
        assert_eq!(h.undo(&mut s), None);
    }

    #[test]
    fn another_group_or_close_splits() {
        let width = Group::Box(Node::Camera, Field::Width);
        for (first, second, close) in [
            (fov_box(), width, false),
            (fov_box(), fov_box(), true),
            (Group::Once, Group::Once, false),
        ] {
            let mut s = with_fov(&new_scene(), 1.0);
            let mut h = History::new(&s);
            change(&mut h, &mut s, 2.0, first);
            if close {
                h.close();
            }
            change(&mut h, &mut s, 3.0, second);
            h.undo(&mut s);
            assert_eq!(s.camera.field_of_view, 2.0);
            h.undo(&mut s);
            assert_eq!(s.camera.field_of_view, 1.0);
        }
    }

    #[test]
    fn a_new_record_clears_redo() {
        let mut s = with_fov(&new_scene(), 1.0);
        let mut h = History::new(&s);
        change(&mut h, &mut s, 2.0, Group::Once);
        h.undo(&mut s);
        change(&mut h, &mut s, 3.0, Group::Once);
        assert_eq!(h.redo(&mut s), None);
        assert_eq!(s.camera.field_of_view, 3.0);
    }

    #[test]
    fn the_oldest_step_goes_after_100() {
        let mut s = with_fov(&new_scene(), 0.0);
        let mut h = History::new(&s);
        for fov in 1..=101 {
            change(&mut h, &mut s, fov as f64, Group::Once);
        }
        for _ in 0..100 {
            assert!(h.undo(&mut s).is_some());
        }
        assert_eq!(h.undo(&mut s), None);
        assert_eq!(s.camera.field_of_view, 1.0);
    }

    #[test]
    fn undo_back_to_the_saved_scene_is_not_unsaved() {
        let mut s = with_fov(&new_scene(), 1.0);
        let mut h = History::new(&s);
        change(&mut h, &mut s, 2.0, Group::Once);
        assert!(h.differs(&s));
        h.undo(&mut s);
        assert!(!h.differs(&s));
        h.redo(&mut s);
        assert!(h.differs(&s));
        h.saved(&s);
        assert!(!h.differs(&s));
    }

    #[test]
    fn clear_forgets_everything() {
        let mut s = with_fov(&new_scene(), 1.0);
        let mut h = History::new(&s);
        change(&mut h, &mut s, 2.0, Group::Once);
        let other = with_fov(&new_scene(), 5.0);
        h.clear(&other);
        let mut t = other.clone();
        assert_eq!(h.undo(&mut t), None);
        assert!(!h.differs(&other));
    }

    #[test]
    fn every_copy_can_be_rewritten() {
        let mut s = with_fov(&new_scene(), 1.0);
        let mut h = History::new(&s);
        change(&mut h, &mut s, 2.0, Group::Once);
        change(&mut h, &mut s, 3.0, Group::Once);
        h.undo(&mut s);
        // As Save As rewrites asset paths: here, every copy's width.
        h.for_each_scene(|d| d.camera.width = 7);
        s.camera.width = 7;
        h.undo(&mut s);
        assert_eq!((s.camera.field_of_view, s.camera.width), (1.0, 7));
        // The saved scene was rewritten too.
        assert!(!h.differs(&s));
        h.redo(&mut s);
        h.redo(&mut s);
        assert_eq!((s.camera.field_of_view, s.camera.width), (3.0, 7));
    }

    #[test]
    fn nothing_to_undo_or_redo_leaves_the_scene() {
        let mut s = with_fov(&new_scene(), 1.0);
        let mut h = History::new(&s);
        assert_eq!(h.undo(&mut s), None);
        assert_eq!(h.redo(&mut s), None);
        assert_eq!(s.camera.field_of_view, 1.0);
    }
}
