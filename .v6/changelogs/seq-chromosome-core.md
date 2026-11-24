## Seq* (next): the chromosome's core

- A ruler's `units: bases` labels a genome's positions in bp, kb or Mb by the tick step (p53's chromosome uses it).
- A numbering's ruler ticks its step's round values and labels insertion codes (52A, 100A–E, 111.1) where they have
  room.
- The ruler rail draws just below the ruler on the canvas, never over its ticks or labels; exports are unchanged.
- Inline annotations over a long system read one tile when zoomed out instead of one per 1,024 positions.
- A template's track shows only while its node's region is in view, so browsing no longer leaves tracks behind.
- The exported SVG names a pairs layer's arcs.
- ESPript's `gaps` setting is now wired to its identity colours (EC-135), and its default `count` changes them: every
  gapped column of an existing ESPript view now colours differently.

## Core: `formatGenomic`

- `formatGenomic(n, step)` in `@molstar/core`'s `util/number`: a genomic position in bases, kb or Mb, chosen by a
  ruler's tick step.
