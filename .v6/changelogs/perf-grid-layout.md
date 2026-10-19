## Core: faster grid lookups

- `GridLookup3D` queries (`find`, `check`, `approxNearest`, `nearest`) are 1.25–1.3× faster with the same hits in the
  same order and the same squared distances, and no extra memory: the lookup indices are resolved without the
  per-candidate `OrderedSet.getAt` dispatch (one find per atom, 186k atoms: r 4 Å 107 → 85 ms, r 6.4 Å 197 → 158 ms;
  7qpd r 4 Å 8.2 → 6.4 ms) (#1992)
