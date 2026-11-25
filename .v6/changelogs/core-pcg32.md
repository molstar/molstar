## Core: `PCG` is a real PCG32

- `PCG` (`@molstar/core`'s `data/util/hash-functions`) is now a real PCG32 (64-bit state, XSH RR, the reference
  `pcg32_srandom(seed, 54)` sequence). The old one kept 32 bits of state and was biased (floats' eighths from 8.6% to
  21.1%). Seeded sequences change: the gray shades of distinct group colors around a gray base, SSAO's sample noise, the
  volume dot representation's and the streamlines' seeded placements.
- Seq* (next): the rows' majority summaries count blocks under 64 rows exactly.
