## Performance: many structures in one scene

Measured with up to 500 separate entries of 1.2k–11.5k atoms in one plugin (`tests/src/pipeline/many.ts`, loads in one
data transaction as `OpenFiles` and Prism do). Loading 500 to a full scene: 74 s → 16 s; the representation presets
alone 60 s → 10 s, allocating 72 → 4.2 GB.

- State updates with more than one edit (a preset's component and representation, an empty component deleted) no
  longer walk the whole state tree twice: the reconciler takes deletes and update roots from what the builder's tree
  recorded (`StateBuilder.EditInfo.edits`), the same as the whole-tree diff, which still runs after an update that
  threw (core)
- The structure info behavior (model and structure indices) keeps its root models and structures instead of walking
  the state for each created model or structure (on top of #2008's skip for other objects): model and structure
  creation for 500 entries 7.1 s → 0.36 s
- `updateFocusRepr` (every representation preset) skips themes the focus representation already has; its behavior
  update walks the state once and commits nothing without focus representations. Such an update no longer clears the
  data state's undo history as a side effect
- SelectLoci no longer marks an empty selection on every representation for each representation created
- The volume and particle hierarchies skip their walk while the state has no volume or particle list: with #2008's
  measurement sync, a single theme change, overpaint or removal with 500 structures 122–435 ms → 45 ms; overpaint of
  all 6.9 s → 0.19 s
