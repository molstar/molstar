## Fixes from the perf audit

- Fix molecular surfaces of large structures at coarse resolution computing for seconds without a progress update or
  an abort check (res 8: 1 chunk per pass → 1 per 903 atoms)
- Fix `PeriodicBox.minimumImage` giving a farther image in skewed boxes (rhombic dodecahedron, L = 10, d = (0, 0, 4.5):
  7.52 → 4.5 Å), which put a chain of a centred complex at a wrong image
- Pipeline bench: coarse-only entries run through; a stage that throws fails the run (exit 1, no numbers in its row)
- Performance: MessagePack `bin` values are copied natively (`decodeMsgPack` of a 32 MB BinaryCIF: 14 → 2 ms)
- Performance: builders keep their optimized code across major GCs (`anchorMap` in `core/util/map-anchor`; 7qpd
  cartoon alloc 72 → 35 MB, gaussian surface 84 → 73 MB)
- Performance: a recolor or alpha change of a merged representation skips the merge check (7qpd: 7/7 → 0 members)
- Performance: the `.cif` variant sniff reads a plain mmCIF 3 times instead of 5 (124 MB: 237 → 143 ms)
- Performance: a normalized value palette finds its range without strings (1.1M values: 360 → 8 ms, 310 → 60 MB)
- Performance: `getAtomicPolymerElements` drops an always-true test per residue (flat-32: ~1.5 → ~0.8 ms)
- Periodic frames made whole are cached per topology and let go with it (100k atoms: 82 → 9 MB kept after the
  trajectory is dropped)
