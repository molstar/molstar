## Fixes from the perf audit

- Fix a builder's `updateState` evaluating its node, and a rebased builder losing its tag, `dependsOn` and flag edits
- Fix `setDependsOn` not updating the tree's dependencies: deleting the dependency left the dependent in error
- Fix `OrderedSet.forEachSegment` emitting an interval's last segment twice, `IntervalIterator` for intervals not at
  0, `Segmentation.projectValue` inside a segment, `ChunkedArray.is`, Column `toArray` from the last row, `deepClone`
  of typed arrays, Map and Set, `LRUCache.set` of an existing key, the UTF-8 fallback decoder (2- and 4-byte
  sequences), `Sphere3D.copy`/`expandBySphere` keeping stale extrema, a reused `MinimizeRmsd` state writing into a
  caller's centre, and Tensor `getCoords` of rank 4+ spaces writing a stray key
- Performance: builder `updateState` no longer scans the builder's actions. 10k nodes: 257 → 4.6 ms
- Performance: a builder commit carries the flags written after its `build()` without visiting every flagged node.
  10k flagged nodes: 1.1 → 0.01 ms per commit
- Performance: a state update walks the evaluated subtrees once for their status and the scheduled set; queued updates
  wake only the next one (1000 queued: 5.6 → 0.5 ms)
- Performance: `OrderedSet.subtract` of a sorted array from an interval walks both once (200k − 100k: 4.9 → 0.3 ms);
  `SortedRanges` finds a unit's first polymer range by binary search (2500 units × 2500 ranges: 24 → 0.7 ms)
- Performance: `PD.isParamEqual` compares array values by index (100k: 3.3 → 0.04 ms); `distinctColors` keeps the
  lists of its last 16 calls ('generate' palettes rebuilt per theme: 3.3 → 0.01 ms for 20 colors)
- Performance: `floodfill3` of a plain space walks by strides from the faces (220×227×233: 451 → 200 ms)
