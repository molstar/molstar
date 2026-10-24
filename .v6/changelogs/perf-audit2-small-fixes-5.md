## Fixes from the perf audit

- Fix: long parses (text CIF, CSV, Cube, DX, MOL2, PSF) yield and report progress at about the update rate; the gaps between yields
  reached ~2× the rate, and several times it when the parse rate changed
- Fix: a kinemage that re-opens a `@group` gives that group's lists its own masters, not the previous group's
- Fix: on a device without `fragDepth` (WebGL1 without EXT_frag_depth) the mesh fallback of spheres, bonds, backbone,
  nucleotide, interactions, particle and volume-dot visuals is kept; it was destroyed and rebuilt at every update (a
  recolor of 186k-atom interactions: 0.4–1 s)
- Fix: water bridges with `includeParent` follow the child by `parentDisplay`, as the other interactions do; every
  bridge of the whole structure was drawn
- Fix: overlapping attaches of a structure property with the same params compute once (a selection query during an
  ASA-colored build computed ASA twice: 0.9 s at 186k atoms)
- Performance: interaction feature grids get an exact box without the 98-direction boundary, and a tester's subset is
  made once per type set. 186k atoms: interactions 389–401 → 351–357 ms; the grids of a computation 89–92 → 21 ms
- Performance: water and ion units read their struct_conn partners once for the inter-unit bonds. 1.1M atoms:
  bonds-inter 173–180 → 151–158 ms
- Performance: polymer and other chain-level component queries take whole units without copying their elements.
  1.1M atoms: 4.8 → 1.7 ms, 31 → 2 MB allocated
- Performance: a component of whole units shares the parent's inter-unit bond lists. 7776 units, 95k edges: 25 →
  14 ms, 76 → 43 MB allocated
- Performance: model building buckets atoms reading each key once (sortAtomSite). 1.1M atoms: 45–49 → 32–33 ms
- Performance: clipping and wiggle layers skip the visuals a change does not reach (Prism cuts). 1344 visuals: one
  texture update per reached visual instead of 1344 per change
