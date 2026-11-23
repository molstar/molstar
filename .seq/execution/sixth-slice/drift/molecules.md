# Drift: molecules.md

## `block.form` in the context menu

**What**: the design's interaction matrix and its Consequences put `block.form` "in the context menu (as `copies`)".
The viewer has it as a command only: the `editing` behavior's `block.form`, in no built-in menu, with no default key;
a host's `on: 'element'` header entry may name it (viewer spec, Commands).

**Where**: `packages/seq/viewer/src/` (the editing behavior's commands); the viewer spec says what is built.

**Options**: (a) the design follows the code (a host decides which menus offer it); (b) add it to the default
element menu.

**Recommendation**: (a), matching the library-not-viewer line (composition is the application's); the stories show it
through a host's header entry.
