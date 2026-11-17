# Drift: lane G's choices in the viewer's models (to fold into viewer-models.md at the stage 4 gate)

- `SpecModel.apply` takes a whole new root (or a function of the old one): no patch format in the tree core yet; switch
  to `repr` patches when they land.
- A viewport is offsets along the axis in display order (not positions): handles a section's range and a wrapped cyclic
  range. Edit the design.
- A scroll group shares the range only; each section keeps its own vertical scroll and width.
- `Follow.point` is valid within one run only (an intron breaks it); followers use the exact `followRange` hull adjusted
  for fractions. Note it in the core `_spec`.
- A `hover` gesture in the binding table; `key` gestures go to `@molstar/ui`'s `Keymap` (reported if put in the table).
- The hit type carries `extent`, `section` and `target`; lane F's hit type must match it (or one replaces the other at
  the stage 4 merge).
- Config as a plain default object with added items (`hover.closeCodons`, `viewport.maxPxPerPosition`,
  `viewport.groupEqualAxes`, `viewport.onResize` `keep-zoom`, `follow.followerMoves` `ignore`, `navigate.*`); the
  `ParamDefinition`s come later.
- After the fix round: `follow.cyclicArc` (`wrap`: a cyclic follower's viewport runs past the axis total and `pieces`
  splits it; `hull` the alternative); follow is continuous via `point` within a run, the `followRange` hull only where
  an end maps nowhere or the ends cross; `modifiers: '*'` bindings tried after an exact match, hover bound with `*`,
  Shift+wheel pans, Alt unbound; `SpecModel.attachHistory` (changes outside a run invalidate the history); selection
  bindings in `SelectionBindings`, in no preset until the selection lane; a saved `session()` stores positions, not
  offsets.
- The orchestrator added `./models` and `./behaviors` entries to `seqx-viewer` (packages/seq AGENTS.md: across folders
  name the entry).

