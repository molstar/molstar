## Particles: nested lists

- Nested particle lists: a target of kind `particles` (a template) places a whole particle list at each of its
  particles, two levels (`model/particles/particle-nesting.ts`); a repeated part stored once, its placements a copy
  each (the yeast cell on a 2x2 grid: a pack of 714 MB, not 2.8 GB written flat). Every expanded instance a particle of
  its own: picked, labelled, highlighted, selected and layered alone, in one copy or in every copy (loci by path:
  `Loci.partial`, `ParticleNesting.ofCopy`, `ofInstance`, `ofTemplateParticles`).
- Particle packs with templates: a template a document of its own in the same file, written and read once;
  `scripts/particle-grid-pack.mjs` writes a cell's grid nested. A template tile-ordered (tiles about a store page's
  particles, by target within), with its coarse levels and boundary stored (read, not sampled: the yeast nested 2x2's
  representations 1.33 -> 0.82 s); a template an older writer wrote is read with a warning (write it again).
- The target representation draws a list with templates by its copies, not flattened: on WebGPU culled and composed on
  the GPU (the page store's cull: copies near first, groups culled whole, a copy the occlusion pyramid hides left out),
  on WebGL2 expanded in the vertex stage, on WebGL1 flattened. Far copies drawn at their template's coarse levels; a
  copy whose leaves are not loaded drawn at its coarse level meanwhile. The isg cell on a 16x16x16 grid, near and
  turning (Chrome): 45 -> 23 ms a frame.
- A template's leaf pages loaded where the GPU's cull wants them (a copy at full detail, its page in view), let go when
  no longer wanted; their GPU memory within the device's storage buffer limit (the yeast nested 2x2 whole: 1,728 ->
  773 MB of GPU storage).
- A nested pack opened from a file or URL by range where the GPU's cull draws it: its templates' particles read by
  page as the view wants them, kept in a 256 MB cache once on the GPU (the yeast nested 2x2 opened: 304 of 784 MB read,
  371 MB held instead of 855, 173 ms instead of 275). What needs all of a template's particles (flattening, spacefill,
  orientation, fibers, writing a pack, WebGL's copies) refuses such a template; a loci of part of a copy is bounded by
  its tiles' boxes (a single particle's sphere a tile's).
