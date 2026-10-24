//! The preview file the viewer watches: the scene with every asset path
//! absolute, in a folder of its own that is removed when the editor ends.

use crate::paths::absolutize;
use rtc_rust::scene::SceneDescription;
use std::ffi::OsStr;
use std::path::{Path, PathBuf};

pub struct Preview {
    dir: PathBuf,
    path: PathBuf,
}

impl Preview {
    /// Make `root/rtc-edit-<pid>/`; the preview is `scene_name` in it.
    pub fn new(root: &Path, scene_name: &OsStr) -> Result<Preview, String> {
        let dir = root.join(format!("rtc-edit-{}", std::process::id()));
        std::fs::create_dir_all(&dir).map_err(|e| format!("{}: {e}", dir.display()))?;
        Ok(Preview {
            path: dir.join(scene_name),
            dir,
        })
    }

    /// The preview's folder, removed with it.
    pub fn dir(&self) -> &Path {
        &self.dir
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    /// Save `desc`, asset paths made absolute against `scene_dir`.
    pub fn write(&self, desc: &SceneDescription, scene_dir: &Path) -> Result<(), String> {
        let mut desc = desc.clone();
        absolutize(&mut desc, scene_dir)?;
        desc.save(&self.path).map_err(|e| e.to_string())
    }

    /// Move the preview to `scene_name` in its folder, keeping what it holds.
    pub fn rename(&mut self, scene_name: &OsStr) -> Result<(), String> {
        let path = self.dir.join(scene_name);
        std::fs::rename(&self.path, &path).map_err(|e| format!("{}: {e}", path.display()))?;
        self.path = path;
        Ok(())
    }
}

impl Drop for Preview {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    fn fresh(name: &str) -> PathBuf {
        let root =
            std::env::temp_dir().join(format!("rtc_edit_preview_{name}_{}", std::process::id()));
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(&root).unwrap();
        root
    }

    #[test]
    fn write_puts_a_loadable_scene_in_the_folder() {
        let root = fresh("write");
        let dir = Path::new("../rtc-metal/scenes");
        let desc = SceneDescription::from_path(dir.join("chapter15.json")).unwrap();
        let p = Preview::new(&root, OsStr::new("chapter15.json")).unwrap();
        p.write(&desc, dir).unwrap();
        assert!(p.path().ends_with("chapter15.json"));
        let folder = p
            .path()
            .parent()
            .unwrap()
            .file_name()
            .unwrap()
            .to_string_lossy();
        assert!(folder.starts_with("rtc-edit-"), "{folder}");
        let (world, _) = rtc_rust::scene::load(p.path()).unwrap();
        assert!(world.shapes.len() > 2, "teapots missing");
    }

    #[test]
    fn dropping_the_preview_removes_its_folder() {
        let root = fresh("drop");
        let p = Preview::new(&root, OsStr::new("x.json")).unwrap();
        let folder = p.path().parent().unwrap().to_path_buf();
        assert!(folder.is_dir());
        drop(p);
        assert!(!folder.exists());
    }

    #[test]
    fn renaming_the_preview_moves_its_file_and_keeps_its_contents() {
        let root = fresh("rename");
        let mut p = Preview::new(&root, OsStr::new("untitled.json")).unwrap();
        fs::write(p.path(), "the scene").unwrap();
        let old = p.path().to_path_buf();
        p.rename(OsStr::new("boxes.json")).unwrap();
        assert!(p.path().ends_with("boxes.json"));
        assert_eq!(fs::read_to_string(p.path()).unwrap(), "the scene");
        assert!(!old.exists());
    }
}
