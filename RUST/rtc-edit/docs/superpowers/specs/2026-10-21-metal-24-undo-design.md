# Milestone 24: undo and redo in the scene editor

Branch `metal/24-undo`, merged into develop with `--no-ff` when done.
Milestones 17–21 and 23 built the scene editor; 22 gave the viewer's quit
alert a Cancel button.

## Goal

The editor has no undo: a wrong delete or a bad edit means reopening the
file and losing the other edits. Cmd-Z should take back the last edit and
Cmd-Shift-Z put it back, as in any Mac app.

Success:
- Delete a sphere, Cmd-Z: it is back, with its boxes.
- Type `1.25` into a box, leave it, Cmd-Z: the box shows its old value.
- Follow the viewer, orbit, Keep, Cmd-Z: the camera's `from` is back where
  it was before following.
- Edit, save, edit again, Cmd-Z: the ` •` disappears.

## Out of scope

- Undoing the viewer's orbit, or anything in a file on disk.
- Undo of folding, selection or the shaded box: only changes to the scene
  are steps.
- A history that survives New, Open or quitting.
- An undo menu or list of steps; the status line names each one.

## Decisions

### What one step is

| Action | One step |
|---|---|
| Typing in a box | everything from entering the box to leaving it |
| `x`, `c`, `u` on a shaded box | each one |
| The `closed` checkbox | each click |
| Add, delete (button, `a`, `d d`, Cmd-Delete) | each item |
| Following the viewer | the whole run of follow updates, until following stops or another edit is made |
| Untick (puts back the old `from`) | one step |
| Keep | no step; it ends the follow step |

- **Leaving a box:** the box being edited changes, to none or another box
  (Enter, Esc, Tab, a click, a selection).
- **Cmd-Z while following** stops following first (without restoring),
  then undoes.
- **Redo** (Cmd-Shift-Z) re-applies undone steps; any new edit clears
  them.
- **Depth:** the last 100 steps; older ones are dropped.
- **New and Open** clear the history; Save doesn't.

### The ` •` and the viewer

- The editor keeps the scene as last saved or opened. After an undo or
  redo the ` •` shows only when the scene differs from it
  (`SceneDescription: PartialEq`).
- After an undo or redo the viewer is behind (Cmd-L loads it), as after
  any edit.

### After an undo or redo

- The selection stays when that item is still a row of the tree;
  otherwise nothing is selected (the tree has the keys). The boxes show
  the scene's values; red boxes are cleared.
- The status line names the step: `undid: <label>` or `redid: <label>`;
  with nothing to undo or redo, `nothing to undo` / `nothing to redo`.
- Labels: `typing in <row label>` (e.g. `typing in position`),
  `closed`, `add <item>` (e.g. `add sphere`), `delete <item>`,
  `follow the viewer`, `untick follow`, and for `x`, `c`, `u`
  `empty <row label>`, `change <row label>`, `reset <row label>`.

## Code

- `src/history.rs`, no iced:
  - `pub enum Group { Box(Node, Field), Follow, Once }`;
  - `History` holds the undo and redo lists of `(SceneDescription,
    String)`, the open group, and the saved scene;
  - `record(&mut self, scene: &SceneDescription, group: Group, label:
    &str)`, called before a change: the same open group (not `Once`)
    merges, anything else pushes a copy of `scene`, clears redo, and
    opens `group`;
  - `close(&mut self)`; `undo` / `redo(&mut self, scene: &mut
    SceneDescription) -> Option<String>` (close first, swap with the top
    copy, return its label); `saved(&mut self, scene)`; `differs(&self,
    scene) -> bool`; `clear(&mut self, scene)` (empty, and `scene` is the
    saved one).
- `src/main.rs`:
  - one `record` before each change: `edit` (`Box` for typing, `Once` for
    `x`/`c`/`u`), `Closed`, `add_item`, `delete_selected`, `follow`
    (`Follow`), untick (`Once`);
  - after every message, a change of the box being edited calls
    `close()`;
  - Cmd-Z → undo, Cmd-Shift-Z → redo; undo and redo stop following;
  - Keep calls `close()`, ending the follow step;
  - Save and Open call `saved`; New and Open call `clear`.
- `README.md`: the keys and what one step is.

## Error handling

Nothing can fail: undo and redo with empty lists only say so.

## Testing

`history.rs` (unit tests):
- undo restores the scene before the change, redo re-applies it;
- the same `Box` group merges; another group, `close`, or `Once` splits;
- a new record clears redo;
- the 101st step drops the oldest;
- undoing back to the saved scene makes `differs` false;
- `clear` forgets everything;
- undo and redo with nothing to do return `None` and leave the scene.

By eye: the four success checks above, Cmd-Shift-Z after each, and `x`
then Cmd-Z on a shaded box.
