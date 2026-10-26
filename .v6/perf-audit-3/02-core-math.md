# Slice 02: core math (v3)

Scope: [`packages/core/src/math`](../../packages/core/src/math) (molecular surface, gaussian density, `GridLookup3D`,
primitives, `Tensor`, RMSD, graphs, the new `PeriodicBox` and `SymmetryOperator.rigidMotionTest`), plus the graphics
molsurf visual's position data ([util/molecular-surface.ts](../../packages/graphics/src/repr/structure/visual/util/molecular-surface.ts))
and model's `areBackboneConnected`. HEAD `16277d822`, v2 base `0944823d4`. Node 25.6, Apple M4.

Measurement conditions: the machine was heavily loaded during this run (load average 8–32, mostly 20–30, with other
auditors). Wall times are unusable at that load: the same molsurf stage took 2.6–15 s from one run to the next. So every
time below is **CPU time (`process.cpuUsage`) of variants interleaved in one process**, plus deterministic counters,
elements-kind checks and bit-exact field comparisons. Absolute CPU times are about 2–3× v2's quiet wall times (cache and
frequency contention). Ratios within one process are what counts. Benches: `…/scratchpad/pa3/02/` (see the end).

PRs since v2 that touched the slice: #2025 (02/P10 + 02/B5), #2033 (02/P11, 02/B1), #2034 (02/B2, 07/P6 = 02/P7), #2037
(`PeriodicBox`, new), #2046 (02/B3, 02/B4, 02/12/B2, 02/12/F1, 02/P9 partly; `rigidMotionTest` is from #2025/#2047's
05/P14 work, in `symmetry-operator.ts`).

## Hot paths now

| Path | Fixture | Cost (CPU, loaded machine) | Where the time goes |
| --- | --- | --- | --- |
| `calcMolecularSurface` torus pass | flat-32, 185k atoms, unit surfaces | 76 % of the surface at res 8, 66 % at res 2, 57 % at res 1.3 (with and without torii, same process conditions) | 3.54M tori at every resolution; per torus a 38-neighbour filter (135M iterations), 12–24 probes of which 98.6 % are obscured (119–147M occlusion iterations), 116–613M voxel visits by the 1.4 % unobscured probes |
| `calcMolecularSurface` points pass | flat-32 | 1.4 s CPU at res 8, 5.8 s at res 1.3 | `find` + `loadNeighbours` per atom (7.08M neighbours), occlusion 40–415M iterations |
| `GridLookup3D.query` | flat-32, one find per atom | r 2: 63 ms, r 4: 112–115 ms, r 6.4: 200 ms (run at load ~5) | scattered `px/py/pz[idx]` per candidate (P2), 32 elements per cell (P3) |
| `accumulateRange` (CPU gaussian) | 7qpd res 0.5 | median 296 ms CPU vs 265 ms for the split kernel | 157 OSR-exit deopts in 3 reps (P13), `fasterExp` per voxel (P8) |

## v2 → now

| v2 id | Title | Status | PR | Note |
| --- | --- | --- | --- | --- |
| P10 | Vec3 elements-kind contagion | fixed | #2025 | `kinds-head.ts`: after load, boundaries, lookups, intra/inter bonds and the structure lookup, `Vec3()`, `Vec3.create`, `Vec3.add` and `Mat4()` all have double elements (flat-32 and 7qpd). The heavy-lane guard v2 proposed was not added (no `%HasDoubleElements` check in the repo) |
| B5 | `areBackboneConnected` reads index -1 | fixed | #2025 | [atomic-ranges.ts:41](../../packages/model/src/model/structure/model/properties/utils/atomic-ranges.ts#L41) returns early; scalar distance, no Vec3 |
| P2 | Grid query: scattered x/y/z per candidate | still open (partly fixed in v1→v2) | | rerun (load ~5): r 2 63.1 → 46.5–49.0 ms, r 4 112–115 → 83–91 ms, r 6.4 200 → 152–154 ms with bucket-ordered copies; same hits |
| P3 | 32 elements per cell | still open | | count 8: r 2 63 → 43–44 ms, r 4 112–115 → 87–94 ms, r 6.4 200 → 177–180 ms; with P2: 35–39 / 71–74 / 145–146 ms |
| P11 | `projectTorus` Vec3 scratch, mispredicting filter | fixed | #2033 | bit-identical to v2's kernel (field diff 0, 7qpd res 0.5/36 and 8/12, flat-32 8/12). CPU, interleaved: 7qpd res 8 median 409 → 353 ms (−14 %), flat-32 res 8 6740 → 5730 ms (−15 %), 7qpd res 0.5 1896 → 1826 ms (−4 %, noise level) |
| P12 | Torus pass ignores resolution; second `find` per atom | still open | | counters above; the second `find` + `loadNeighbours` per atom is still there ([molecular-surface.ts:426](../../packages/core/src/math/geometry/molecular-surface.ts#L426)). With the compact unit input (#2034) `OrderedSet.getAt` in `loadNeighbours` is an `Interval` add, so that part of v2's note is gone |
| P13 | `accumulateRange` OSR exit per call | still open | | 157 "exit from OSR'd inner loop" deopts in 3 reps on 7qpd (v2: 212; #2034's bigger chunks mean fewer calls). v2's split kernel: median 296 → 265 ms CPU (−10 %), interleaved |
| 12/F1 | `floodfill3` through closures | fixed | #2046 | `floodfill3Strided` ([tensor.ts:498](../../packages/core/src/math/linear-algebra/tensor.ts#L498)): 220×227×233, 9.6M visited: 189–194 ms vs v2's prototype 342–420 ms in the same process, visited diff 0 |
| P8 | `fasterExp` per voxel; separable kernel | still open | | [cpu.ts:115](../../packages/core/src/math/geometry/gaussian-density/cpu.ts#L115) unchanged |
| P9 | `getEdgeIndices` `[]` per miss; edges as objects | partly fixed | #2046 | the miss returns `emptyArray` ([inter-unit-graph.ts:59](../../packages/core/src/math/graph/inter-unit-graph.ts#L59)); `Map.get` per vertex and the object edges remain (per-unit CSR is an API change) |
| P7 | Per-unit molsurf radii sized `OrderedSet.end` | fixed | #2034 | unit path is compact (x/y/z Float64, radii, `Interval` indices). The structure path ([util/molecular-surface.ts:95](../../packages/graphics/src/repr/structure/visual/util/molecular-surface.ts#L95)) still sizes by `end`, but its indices are structure serials, so `end` ≈ n: no waste |
| B1 | `approxNearest` returns the nearest last | fixed | #2033 | grid ([grid.ts:463](../../packages/core/src/math/geometry/lookup3d/grid.ts#L463)) and `StructureLookup3D._approxNearest` (one result over units) |
| B2 | `updateChunk = 1e5/(maxRadius⁹·scale)` | partly fixed | #2034 | small radii now yield; **molsurf at typical radii and coarse resolution yields 8× less often than before: new B6** |
| B3 | `Sphere3D.copy`/`expandBySphere` stale extrema | fixed | #2046 | `else out.extrema = void 0` in both |
| B4 | RMSD state aliases a caller's centre | fixed | #2046 | the state owns `centerA/B`; `setCenters` copies in, and gives the input a clone |
| 12/B2 | `getCoords` stray key | fixed | #2046 | line removed |
| 12/layouts | Surface fields lay out `[0,1,2]` by hand | still open | | [molecular-surface.ts:476](../../packages/core/src/math/geometry/molecular-surface.ts#L476), [cpu.ts:46](../../packages/core/src/math/geometry/gaussian-density/cpu.ts#L46); ownership only |
| P6 | `getBoundary` 98 directions | dropped (v2) | | v2 showed the cost was P10. Not re-measured |

Counts (18 v2 items): fixed 9 (P7, P10, P11, 12/F1, B1, B3, B4, B5, 12/B2), partly fixed 2 (P9, B2), still open 6 (P2, P3,
P8, P12, P13, 12/layouts), regressed 0, dropped 1 (P6).

New bugs or costs from fixes since v2:

- **#2034 (02/B2 fix) → B6**: the voxel-based chunk ignores the per-atom fixed cost (`find`, neighbour load, ~19 tori) and
  under-counts the points loop (ceil, `ngPoints`). Molsurf yields far less often than in v2 for normal radii (1.7 Å +
  1.4 probe): longest synchronous stretch, 186k-atom single input, 66 → 195 ms (res 1), 82 → 372 ms (res 2), 183 →
  1400–1500 ms (res 8). The gaussian is not affected (25–33 ms longest).
- **#2037 (`PeriodicBox`) → B7**: `minimumImage` is not the minimum image for triclinic boxes, contrary to its doc.
- Reviewed with no finding: #2025 `atomic-ranges` (same boolean, no Vec3); #2033 `projectTorus` scalars (bit-identical,
  `tX[tc]` stays in bounds since `tc ≤ j < nbCount ≤ nbCapacity`; NaN `rInt` keeps none, as before) and `approxNearest`;
  #2034 compact unit input (`id` is ordinal-indexed as `calcMolecularSurface` reads it, `id[i]`; boundary unchanged since
  coordinates are copied exactly as Float64); #2046 `floodfill3Strided` (stack grows before the ≤ 6 pushes; 1-voxel
  dimensions handled; `!(v < t)` matches NaN handling), `setCenters` (tm-align builds a fresh input per call, so no
  stale cached centre), Sphere3D extrema, `emptyArray`; `rigidMotionTest` (a moved first pair fixes `M`; a later equal
  pair under a moved `M` correctly fails).
- Pipeline A/B v2 base vs HEAD (`ab-pipe.sh`, 7qpd, molsurf + gaussian stages, alternated): HEAD was not slower in any
  pair (molsurf 2.6–6.5 s vs base 5.8–15 s; gaussian 1.4–3.5 s vs 2.3–8.2 s), but at load 25–30 that only rules out a
  large regression.

## Open performance findings (ranked)

| # | Id | Status | Title | File | Kind | Estimate | Effort | Frozen |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | P2 | still open | Grid query: bucket-ordered coordinate copies | [grid.ts:360](../../packages/core/src/math/geometry/lookup3d/grid.ts#L360) | time | measured: query −22–27 % (r 4 112–115 → 83–91 ms, flat-32); +12 B per element | S | no |
| 2 | P12 | still open | Molsurf torus pass ignores resolution; second `find` per atom | [molecular-surface.ts:423](../../packages/core/src/math/geometry/molecular-surface.ts#L423) | time | measured: torus pass 57–76 % of the surface CPU (flat-32, res 1.3–8); coarse skip v2: `lowest` −75 % | M | no |
| 3 | P3 | still open | 32 elements per grid cell | [grid.ts:253](../../packages/core/src/math/geometry/lookup3d/grid.ts#L253) | time | measured: query −10–30 % (count 8) | S | no |
| 4 | P13 | still open | `accumulateRange` OSR exit per call | [cpu.ts:70](../../packages/core/src/math/geometry/gaussian-density/cpu.ts#L70) | time | measured: CPU gaussian −10 % (7qpd res 0.5, 296 → 265 ms CPU median) | S | no |
| 5 | P14 | new | Torus fully inside one neighbour: exact early exit | [molecular-surface.ts:354](../../packages/core/src/math/geometry/molecular-surface.ts#L354) | time | measured: 7qpd res 8 −6–10 % (283 → 264 ms CPU median), res 0.5 within noise; bit-identical; 37 % of 7qpd's tori | S | no |
| 6 | P8 | still open | Gaussian `fasterExp` per voxel; separable kernel | [cpu.ts:115](../../packages/core/src/math/geometry/gaussian-density/cpu.ts#L115) | time | v1: `accumulateRange` −20–25 % (not re-measured) | S | no |
| 7 | P9 | partly fixed | `getEdgeIndices`: `Map.get` per vertex, edges as objects | [inter-unit-graph.ts:59](../../packages/core/src/math/graph/inter-unit-graph.ts#L59) | time+memory | v2 estimate: valence model −10–20 % with per-unit CSR | M | no |

### P2 + P3: grid query (re-measured)

`bench-grid.ts` (flat-32, one `find` per atom, 256 unit lookups, hits identical in all variants). First run at load ~5:

| r | current (32) | count 8 | copies (32) | copies + 8 |
| --- | --- | --- | --- | --- |
| 2 | 63.1–63.3 ms | 43.3–43.7 | 46.5–49.0 | 34.9–38.5 |
| 4 | 112.4–114.9 | 87.2–93.6 | 82.7–90.9 | 70.6–74.1 |
| 6.4 | 200.1–200.3 | 177.4–179.6 | 151.9–154.1 | 145.4–146.0 |

A second run at load ~13 kept the ratios (r 4: 212 → 116 ms with copies). Same as v2: the decision is the maintainers'
(+12 B per element, +4 with radii; #1992 dropped the copies on purpose; count 8 changes result order, so the pinned
lookup and bond snapshots change). These are still the biggest exact wins in the slice: the grid query is 32–55 % self of
bonds-intra/inter, interactions and ASA (v2 profile).

### P12: torus pass (counters at HEAD)

`bench-cnt.ts` with an instrumented copy (`ms-cnt.ts`), flat-32, compact unit inputs as the visual builds them:

| res / probes | CPU total | CPU without torii | tori | filter iter. | probes | unobscured | occlusion iter. (torus) | torus voxel visits |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 8 / 12 | 5850 ms | 1400 ms | 3.54M | ~135M (38/torus, 23.3 kept) | 42.5M | 0.60M | 119M | 116M |
| 2 / 24 | 9129 | 3147 | 3.54M | ~135M | 85.0M | 1.21M | 147M | 259M |
| 1.3 / 24 | 13472 | 5771 | 3.54M | ~135M | 85.0M | 1.21M | 147M | 613M |

The tori and the filter do not depend on resolution; 98.6 % of probe positions are obscured. At res 8 one unobscured
probe visits a 6³ window (`ngTorus = 2 + ⌊1.4/8⌋`) of 8 Å voxels. Shrinking that window is not exact: voxels still at
the "visited, never projected" value 1001 are filled by whatever probe reaches them, so the window is part of the output.
The fixes are v2's: (1) exact: keep each atom's neighbour ordinals from the points pass (CSR, freed at return) instead of
the second `find` + `loadNeighbours` (7.08M neighbours loaded twice); (2) output-changing: skip or thin the torus pass
when resolution > ~2× probe radius.

### P14 (new): torus probe circle inside one neighbour

- Mechanism: when the whole probe circle (centre `mid`, radius `rInt`) lies inside one neighbour `c ≠ b`, every probe is
  obscured and nothing is written, but the kernel still finishes the filter and tests 12–36 probes.
- Prototype (`ms-bury.ts`): in the filter loop, `inner = 2·nbR[j] − reach` (= `nbR − rInt(1+1e-6) − 1e-4`); `if (inner > 0
  && d² < inner² && nbIndex[j] !== b) return;`. The 1e-4 Å margin covers the probe-table rounding v2's filter comment
  already bounds, so every skipped probe would have been obscured by `c` (strict `<` on `r²`).
- Evidence (`bench-ab2.ts`, interleaved CPU, 7qpd): 82k tori exit early; res 8/12 median 282 → 264 ms (−6 %; per round
  −6–18 %), res 0.5/36 1154 → 1191 ms (noise). Fields bit-identical with v2's kernel in both.
- Verdict: small. Worth folding into any P12 rework, not on its own.

### P13: gaussian OSR exits (re-measured)

`bench-gauss-ab.ts`, 7qpd res 0.5, 8 interleaved rounds: HEAD median 296 ms CPU, v2's split kernel (`g-split.ts`:
per-atom kernel split out, `min[0..2]` and `id ? id[i] : i` hoisted) 265 ms (−10 %). `--trace-deopt`: 157 OSR-exit deopts
in 3 reps. CPU gaussian runs headless or on GPU fallback only.

## Bugs

| Id | Severity | Status | Title | File | Failure scenario |
| --- | --- | --- | --- | --- | --- |
| B6 | medium | new (regression, #2034); **fixed in #2050** | Molsurf chunks by voxels only: no yield for seconds at coarse resolution | [molecular-surface.ts:489](../../packages/core/src/math/geometry/molecular-surface.ts#L489) | a 186k-atom structure-level molecular surface at res 8: the torus pass runs as 2 synchronous chunks of 1.4–1.5 s each (v2: 56 updates, longest 183 ms); progress and abort wait for each |
| B7 | medium | new (#2037); **fixed in #2050** (shifts bounded by the inverse rows, exact for any box) | `PeriodicBox.minimumImage` is not the minimum image in triclinic boxes | [periodic-box.ts:47](../../packages/core/src/math/geometry/periodic-box.ts#L47) | GROMACS rhombic dodecahedron, L = 10 Å: d = (0, 0, 4.5) gives (−5, −5, −2.57), \|m\| 7.52 instead of 4.5; 7.7 % of random \|d\| < L/2 get a farther image, by up to 6.1 Å. In `periodic.ts:200` (centre pinning) a chain of the centred complex is put at a wrong image, so the complex stays split |

### B6: molsurf yield regression

- Cause: #2034 set `updateChunk = max(1, ceil(2e6 / (2·(maxRadius + probe)·scale + 2)³))`. At res 8 (scale 0.125) the
  denominator is 21.4, so a chunk is 93.5k atoms (old formula: 6.7k). But the points loop visits `(2·(⌈rad·scale⌉ +
  ngPoints) + 2)³` = 216 voxels per atom there, not 21, and every atom also pays a `find`, a neighbour load and ~19 tori
  (38-neighbour filter each), whatever the resolution. The torus pass is the work at coarse resolution, and its chunk
  uses the same atom count.
- Evidence (`probe-chunk.ts` / `probe-chunk-v2.ts`: flat-32's atoms merged into one input, `ctx.update` gaps, load 8–18):

  | res / probes | v2 updates, longest | HEAD updates, longest |
  | --- | --- | --- |
  | 1 / 24 | 442, 66 ms | 104, 195 ms (points) |
  | 2 / 24 | 222, 82 ms | 26, 372 ms (torii) |
  | 8 / 12 | 56, 183 ms | 4, 1511 ms (torii); points 408 ms |

  Unit surfaces (~720 atoms per unit on flat-32) are one chunk either way; the structure visuals
  (`structure-molecular-surface-mesh`/`-wireframe`) and very large units are affected, and large structures get coarse
  resolutions from quality `auto`. The CPU gaussian is fine (HEAD longest 25–33 ms vs v2 6–26 ms).
- Fix: count the real loop and a per-atom fixed cost, e.g. `ng = ⌈(maxRadius + probe)·scale⌉ + ngPoints`, `perAtom = (2·ng
  + 2)³ + 2000` (the constant ≈ `find` + tori in voxel-visit units), `updateChunk = max(1, ceil(2e6 / perAtom))`; or cap
  at a few thousand atoms. Same formula for both passes is fine once the constant is in.

### B7: `minimumImage` in triclinic boxes

- Cause: rounding each fractional coordinate gives the lattice point of the box-shaped cell around `d`. That is the
  nearest image only for a rectangular box. For a skewed box a fractional coordinate can exceed 0.5 while `|d|` is far
  under half the shortest box vector (the c-height of the dodecahedron is 0.707 L), so the doc's "exact … while |d| is
  under half the shortest box vector" is false.
- Probe (`probe-mi.ts`): box a = (L,0,0), b = (0,L,0), c = (L/2, L/2, L√2/2) (GROMACS `-bt dodecahedron`), brute force over
  5³ lattice shifts: 4026 of 52443 random `|d| < L/2` not nearest, worst excess 6.09 Å at L = 10.
- Impact: only caller is centre pinning in `Model.Periodic` ([periodic.ts:200](../../packages/model/src/model/structure/model/periodic.ts#L200)):
  molecule centres tens of Å apart are mapped to the image nearest the first molecule. Dodecahedron and truncated
  octahedron boxes are common in MD. `applyWhole` (the inlined copy, bond-length displacements) and the no-jump rounding
  (frame-to-frame centre motion) work on small displacements and are fine.
- Fix: after rounding, try the 27 neighbouring lattice shifts (−1..1 per vector) and keep the shortest. That is exact for
  GROMACS-reduced boxes and costs nothing at per-molecule rate. Add the dodecahedron case to `periodic-box.test.ts`. Owner of
  the caller: slice 05.

## Already good, re-checked

- P10 stays fixed through load, boundaries, lookups and bonds (`kinds-head.ts`, flat-32 and 7qpd).
- `floodfill3Strided` beats v2's own prototype (190 vs 345 ms) with identical visited sets.
- The molsurf compact unit input (#2034) is a little faster as well as smaller: 7qpd res 0.5 1826 → 1654 ms CPU median
  (−9 %), flat-32 res 8 5730 → 5489 ms (−4 %), bit-identical fields.

## Not audited / not measured

- Wall-clock stage times: load 20–30 for most of the run; only CPU-time ratios, counters and kinds are reported.
- P8 and P9 not re-measured (code unchanged or partly changed as stated).
- `rigidMotionTest` and `PeriodicBox` callers in model (slice 05) only as far as B7 needs.

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/31ee827f-4715-4967-bec3-085f5013d11f/scratchpad/pa3/02/`,
bundled with `node build.mjs <in.ts> <out.mjs>` (esbuild against the worktree, `molstar-src`). Fixtures in `load.ts`.

| Script | What |
| --- | --- |
| `bench-grid.ts` + `grid-opt.ts` | P2/P3 (v2's, rebuilt against HEAD) |
| `ms-cur.ts` (`NOTORII=1`), `ms-v2.ts`, `bench-ab.ts`, `bench-ab2.ts` + `ms-bury.ts` (`mk-bury.mjs`) | P11 status, compact input, P14; interleaved CPU, bit-exact field diff vs v2's kernel |
| `ms-cnt.ts` (`mk-cnt.mjs`), `bench-cnt.ts` | P12 work counters per resolution |
| `probe-chunk.ts`, `probe-chunk-v2.ts`, `probe-chunk-g.ts` | B6: longest stretch between `ctx.update` calls |
| `probe-mi.ts` | B7: minimum image vs brute force |
| `kinds-head.ts` (`--allow-natives-syntax`) | P10 elements kinds at HEAD |
| `bench-gauss-ab.ts`, `bench-gauss-cur.ts` (`--trace-deopt`) | P13 |
| `bench-flood.ts` | 12/F1 |
| `mkbase.sh`, `ab-pipe.sh` | v2-base tree and alternated pipeline runs |
