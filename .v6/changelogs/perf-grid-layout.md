## Core: faster grid lookups

- `GridLookup3D` queries are faster with the same hits in the same order and the same squared distances, and no
  extra memory: lookup indices are resolved without the per-candidate `OrderedSet.getAt` dispatch, and a contiguous
  `SortedArray` (the elements of a typical unit) is treated as a range. `find` 1.25–1.35× (one find per atom, 186k
  atoms: r 4 Å 108 → 82 ms, r 6.4 Å 199 → 153 ms; 7qpd r 4 Å 8.4 → 6.5 ms), `approxNearest` ~1.2×, `nearest` ~6 %
  (#1992)
