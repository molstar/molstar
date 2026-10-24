## Fixes from the perf audit

- Fix disk-style links in the impostor and line builders being drawn from another edge's midpoint
- Fix complex and shape visuals reading the previous geometry's positions after a geometry change
- Fix `computeGaussianDensityTexture3d` computing a 2d texture
- Fix visuals recreated with a layer set (GPU ↔ CPU or backend switch) leaving a merged surface unmerged, every member
  smoothed
- Performance: sphere impostors evaluate the size theme only for a new bounding sphere (rebuild of 1.1M atoms: 44 →
  25 ms); new spheres fill their position-group texture once (1M: 24 → 21 ms)
- Performance: an unchanged merged group skips the merge check at each update (10k members: 39 → 9 ms); a visual's
  first update after a build skips the foreign-resource walk (first cartoon recolor of 1.1M atoms: 65 → 44 ms)
- Performance: inter-unit bond marking walks the smaller side of a loci and a unit pair; the segment mask tests each
  voxel of its box once (124³ box: 100 → 60 ms)
- Performance: the cylinder primitive cache is bounded and checks its keys (20k gizmo radii: 153 → 1.4 MB retained)
- Performance: the illustrative theme lightens each style color once (186k atoms: 26 → 13 ms); cartoon params are
  made once per structure (500 structures, one update: 45 → 0.1 ms)
- `PluginContext.dataTransaction` documents that scripts committing per structure should run inside it
