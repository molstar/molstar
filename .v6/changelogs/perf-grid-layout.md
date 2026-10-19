## Core: faster grid lookups

- `GridLookup3D` queries (`find`, `check`, `approxNearest`, `nearest`) are 1.5–1.7× faster with the same hits in the
  same order and the same squared distances: the build keeps bucket-ordered copies of the positions (and radii), so a
  query scans them linearly (one find per atom, 186k atoms: r 4 Å 109 → 65 ms, r 6.4 Å 199 → 125 ms; 7qpd r 4 Å 8.4
  → 5.0 ms). Costs 12 B per element retained (4 B more with radii); coordinates edited in place after the build are
  not seen
