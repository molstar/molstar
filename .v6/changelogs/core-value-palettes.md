## Core: value palettes shared by MolViewSpec and Seq\*

- MolViewSpec's palettes (categorical, discrete, continuous; their params, split colors and named lists) are one
  implementation in `@molstar/core/util/color/value-palette` and `value-palette-params`, which MVS's runtime and Seq\*
  both use (Seq\* had its own copy). It takes the values to color rather than an MVS annotation, so any per-residue
  values can be colored through a palette; its continuous and discrete functions allocate nothing per call
  (`palette.at`). MVS and SVS colors are unchanged.
- `@molstar/mvs/components/annotation-color-theme`'s `makeContinuousPaletteCheckpoints` is removed (see the breaking
  changes).
- Fixes: a categorical palette's value named like an `Object.prototype` member (`constructor`) no longer breaks MVS's
  color lookup; color names are only the table's own.
