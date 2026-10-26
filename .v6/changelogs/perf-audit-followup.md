## Fixes from the perf audit

- Performance: a hover or click marks only the representations of the hovered structure (and its equivalent copies),
  not every representation in the scene (200 structures: hover 1.7–2.1 → 0.40 ms, select 1.8–3.7 → 0.25 ms per
  event); representations can report the structure roots they mark (`Representation.getStructureRootHashes`)
- Performance: deselecting all (a click on empty space) visits only the representations that were highlighted or
  selected (200 structures: 12.6 → 0.6 ms)
