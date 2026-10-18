# Testing: fixtures, goldens, bench

## 1. Fixtures

`src/_test/fixtures/{features,motm1,exosome}.mvsj`: the `index.mvsj` of real MVSX files, without assets (tests never
download). motm1 and exosome are long animations; features covers most node kinds and frames with `include_camera`.
Excluded from Biome: kept byte for byte.

## 2. Goldens

`src/_test/load-animation.golden.test.ts` (unit run): per snapshot, a fingerprint of every MVS frame and of every state
snapshot frame (refs, params, versions) with its camera and canvas (a stub canvas in the plugin), and the durations;
`Math.random` seeded. Any change of the load's output fails at the frame where it starts.

- **Platform-independent**: math functions (`sin`, `pow`) can differ in the last bit between platforms (CI is x64
  Linux). A fingerprint is a hash of the canonical JSON with numbers left out, plus their count and sums, compared with
  a relative tolerance of 1e-9 of the sum of `|x|`. Versions hash exact numbers: each is replaced by its index among the
  distinct versions of the snapshot and its frames; the root transform's version (the plugin's) is left out.
- Regenerate (`MVS_UPDATE_GOLDENS=1`) only for an intended change of the output, and say why in the commit.

The seeded `Math.random` stream is shared by everything that calls it during a load: a new consumer moves the noise
values without the output being wrong. Then regenerate, and check the expanded frames did not change.

## 3. Bench

`src/_test/load-animation.bench.test.ts`, skipped unless `MVS_BENCH=1` (with `--silent=false`). Prints time, allocation
and retained heap per stage, and replay counts.

| Variable                  | Effect                                                       |
| ------------------------- | ------------------------------------------------------------ |
| `MVS_BENCH_FIXTURES`      | fixture names or `.mvsj` paths (documents too large to keep) |
| `MVS_BENCH_LOAD_ONLY=1`   | `loadMVS` only (the other stages keep every frame's tree)    |
| `MVS_BENCH_RUNS`          | runs per stage (best of)                                     |
| `MVS_BENCH_PROFILE=<dir>` | a `.cpuprofile` of one `loadMVS` per fixture                 |

Node, no canvas, `keepCamera`: the camera and canvas work of a browser load is not measured.
