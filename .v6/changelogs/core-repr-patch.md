## Core: patches for the repr tree

- `@molstar/core/repr/patch`: git's model for repr trees. A patch of six ops (set, unset, insert, remove, move with
  rename, replace) is applied to an immutable root by `applyPatch`, which normalizes and validates what changed, keeps
  refs following a move, and refuses a dangling or silently retargeted ref unless the patch allows it. It returns the
  new root, its net `Changes` and an inverse. `diff` derives a patch between two roots; `patchLanguage(language)`
  describes patches as repr trees. `@molstar/core/repr/patch/conformance` is the law suite (apply, inverse, diff) a
  language can run.
