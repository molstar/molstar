# Performance: fixtures, goldens, bench

## 1. Fixtures

`src/_test/fixtures/{features,motm1,exosome}.mvsj`: the `index.mvsj` of real MVSX files, without their assets (tests
never download). motm1 and exosome are long narrated animations (2,700–2,900 distinct frames); features covers most node
kinds and one long rotation. Excluded from Biome (`biome.json`): kept byte for byte.

## 2. Goldens

`src/_test/load-animation.golden.test.ts`, in the unit run: per snapshot, a hash of the canonical JSON of every MVS
frame and of every state snapshot frame (the root transform's version left out: it is the plugin's), with the durations;
`Math.random` seeded. A change of the load's output — a value, a key, a ref, a version, which steps merge — fails at the
frame where it starts. Regenerate only for an intended change of the output (`MVS_UPDATE_GOLDENS=1`), and say why in the
commit.

The seeded `Math.random` stream is shared by everything that calls it during the load: a change in who else consumes it
(as the old `UUID.create22` did) moves the noise values of motm1 without the output being wrong. Regenerate then, and
check the expanded frames (computed before any other consumer runs) did not change.

## 3. Bench

`src/_test/load-animation.bench.test.ts`, skipped unless `MVS_BENCH=1` (`--silent=false` to see the table;
`MVS_BENCH_PROFILE=<dir>` writes a `.cpuprofile` per fixture; `MVS_BENCH_FIXTURES` takes `.mvsj` paths too, for
documents too large to keep; `MVS_BENCH_LOAD_ONLY=1` measures `loadMVS` alone, for those). Node, no canvas,
`keepCamera`: the camera and canvas work of a browser load is not in it.

| Fixture  | Stage                     | Baseline (`0050f1464`) |             Now |
| -------- | ------------------------- | ---------------------: | --------------: |
| exosome  | `generateStateTransition` |         396 ms, 795 MB |    20 ms, 24 MB |
| exosome  | `loadMVS`                 |      1,803 ms, 3.27 GB | ~350 ms, 868 MB |
| exosome  | retained                  |                  77 MB |         14.7 MB |
| motm1    | `generateStateTransition` |         140 ms, 312 MB |    11 ms, 15 MB |
| motm1    | `loadMVS`                 |        736 ms, 1.40 GB | ~150 ms, 447 MB |
| motm1    | retained                  |                  25 MB |          4.2 MB |
| features | `loadMVS`                 |        587 ms, 1.29 GB | ~285 ms, 733 MB |
| story\*  | `loadMVS`                 |        39 s\*\*, 43 GB | ~2.45 s, 5.4 GB |
| story\*  | retained                  |                1.56 GB |           84 MB |

\* "Terms of entrapment", an MVSX of 16 snapshots and a 118 MB `index.mvsj` (22,546 frames; a primitive of 10,596 lines
whose colors flash in 360 fades): not a fixture, measured from a local copy (`MVS_BENCH_FIXTURES=<path>`,
`MVS_BENCH_LOAD_ONLY=1`, `NODE_OPTIONS=--max-old-space-size=16000`). \*\* after steps 1–4 and the missing-start-color
fix; it did not load before the fix.

Times vary ±5 % between runs. Update the "Now" column with every change that moves it.
