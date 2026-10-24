# Slice 02: core math (audit v2)

Scope: [`packages/core/src/math`](../../packages/core/src/math): geometry (molecular surface, gaussian density,
`GridLookup3D`, boundary helpers, primitives), linear algebra (`Vec3`, `Tensor`, RMSD), graphs. It also covers two
callers where the cost lands in core functions: the graphics molsurf visual (P7) and model's `areBackboneConnected`
(P10/B5, owned by slice 05). Run at `0944823d4` with Node 25.6 on an Apple M4. The machine was shared (load 2-5), so
times carry ±10-20 % noise; allocation numbers reproduce exactly. Fixtures: [00-profile.md](00-profile.md). MVS is out
of scope (see `../mvs-optimizations/research/`). No finding touches Canvas3D or `gl`.

## Hot paths now

| Path | Where it is paid | Cost |
| --- | --- | --- |
| `calcMolecularSurface` | `repr:molecular-surface` | `projectTorus` plus `obscured` are ~65 % of the work. The torus pass is resolution-independent: 3.5M tori on flat-32 even at `lowest`, which is 2.5 of 3.0 s at res 8 |
| `GridLookup3D.query` | bonds-intra/inter, interactions, ASA, the molsurf neighbour search | 32-55 % self of those stages; flat-32 one find per atom at r 4: 108 ms |
| `getBoundary` (98-direction EPOS) | every unit, the lookup+boundary stage | flat-32: 256 units in 62-72 ms with Vec3 contagion vs 26-31 ms without (P10) |
| `accumulateRange` | CPU gaussian surface (headless, GPU fallback) | 36 % of `repr:gaussian-surface` on 7qpd; 212 OSR-exit deopts in 3 runs |
| `InterUnitGraph.getEdgeIndices` | valence model, bond iteration | 23 % of valence-model (asm-192) |

The main new result is P10. On 1tup entries one `Vec3.create(undefined, …)` at model creation makes every `Vec3` that
later goes through the shared helpers generic (`PACKED_ELEMENTS`). Most of the 00-profile's projectTorus allocation
(its P3) and of `Vec3.dot`'s share of lookup+boundary (its P6) are this one call. In v1 the big-bcif fixture
(1tup-2x2x2) had the same trigger, so v1's P4 and P6 numbers were partly contagion too.

Changed since v1: #1988 (molsurf compact neighbour lists, P1), #1987 (`Tensor.strides`, marching cubes, P5), #1992
(grid queries without `getAt`; it dropped the bucket-ordered copies, P2), #1990 + #2004 (mesh builder with inline and
`fromFloat32Array` reads, interactions as impostors, P4), #1994 (`Box3D` helpers, `Mat4.isRigid`,
`Vec3.fromFloat32Array`, 14-direction structure spheres above 10k extrema).

## v1 → now

Fixed 3, partly fixed 1, still open 11, regressed 0, dropped 1 (16 v1 items: P1-P9, B1-B4, and 12/F1, 12/B2 and
12/layouts in `tensor.ts`).

| v1 id | Title | Status | PR | Note |
| --- | --- | --- | --- | --- |
| P1 | Molsurf occlusion rescans the full neighbour list via `getAt` | fixed | #1988 | `loadNeighbours` (:120), per-torus pruning (:331), occlusion tested only if the voxel would change (:248/:256); bit-identical (pinned test) |
| P2 | Grid query: `OrderedSet` dispatch and scattered x/y/z | partly-fixed | #1992 | indices resolved once (:155-159), one loop per index kind. The PR added bucket-ordered copies and then dropped them ("no extra memory"); they would give a further −27 % |
| P3 | Default 32 elements per cell is too coarse | still-open | | `grid.ts:253`; count 8: r 4 108 → 84 ms |
| P4 | `Vec3.fromArray` megamorphic in mesh-builder loops | fixed | #1990, #2004 | `mesh-builder.ts:26-34`; the helpers keep loads monomorphic but do not protect the out arrays' elements kind (P10) |
| P5 | Marching cubes reads through a `Space.get` closure | fixed | #1987 | `algorithm.ts:291` `Tensor.strides`; 220³ 151 → 73 ms (12-tensor) |
| P6 | `getBoundary` 98-direction EPOS | dropped | | Code unchanged, but the cost was P10's contagion: with it removed, the flattened pass is 20-22 ms against the current 20-22 ms (`v/bb-san.ts`). #1994 moved structure spheres above 10k extrema to 14 directions |
| P7 | Per-unit molsurf radius array sized `OrderedSet.end` | still-open | | `graphics/.../util/molecular-surface.ts:48` and :85 |
| P8 | Gaussian `fasterExp` per voxel; separable kernel | still-open | | `cpu.ts:112` unchanged |
| P9 | `getEdgeIndices` `[]` per miss; edges as objects in nested Maps | still-open | | `inter-unit-graph.ts:58` unchanged |
| B1 | `approxNearest` puts the nearest last; callers read [0] | still-open | | `grid.ts:469`, :481; `lookup3d.ts:312`; consumers `blob-surface.ts:182`, `external-structure.ts:150` |
| B2 | `updateChunk = 1e5/(maxRadius⁹·scale)` | still-open | | `molecular-surface.ts:462`, `cpu.ts:65` |
| B3 | `Sphere3D.copy` / `expandBySphere` keep stale extrema | still-open | | `sphere3d.ts:59`, :240 |
| B4 | RMSD state aliases and overwrites a caller's centre | still-open | | `minimize-rmsd.ts:92-96`; tm-align and superposition-ligand reuse state without centres |
| 12/F1 | `floodfill3` through accessor closures | still-open | | `tensor.ts:494`; 220³ 312-352 vs 220-224 ms strided |
| 12/B2 | `getCoords` (rank ≥ 4) writes a stray key | still-open | | `tensor.ts:379` |
| 12/layouts | Surface fields lay out `[0,1,2]` by hand | still-open | | `molecular-surface.ts:445`, `cpu.ts:53`, `blob-surface.ts`; ownership only, no time |

## Open findings (confirmed and plausible, ranked)

| # | Id | vs v1 | Title | File | Kind | Impact | Estimate | Effort | Verdict | Frozen |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | P10 | new; **fixed** #2025 | Vec3 elements-kind contagion from one `Vec3.create(undefined)` | [atomic-ranges.ts:40](../../packages/model/src/model/structure/model/properties/utils/atomic-ranges.ts#L40) | time+memory | high | flat-32: unit boundaries −40-60 %, interactions −35 % alloc and −5-15 % time, molsurf −1.4 GB alloc (−7.4 GB at flat-192) and −4 % time; **measured** (flat-32, A/B): molsurf alloc 1519–1566 → 214–251 MB per build, unit boundaries alloc 50–56 → 16–33 MB and 174–184 → 96–116 ms (round 1; round 2 within noise), `interactions` alloc ~−35 %; times on a loaded machine | S | confirmed | no |
| 2 | P2 | partly-fixed | Grid query: scattered x/y/z per candidate; copies dropped in #1992 | [grid.ts:358](../../packages/core/src/math/geometry/lookup3d/grid.ts#L358) | time | high | query −27 % (r 4 109 → 78 ms); bonds-intra/interactions/ASA −10-15 %; +12 B per element | S | confirmed | no |
| 3 | P11 | new; **fixed** #2033 | `projectTorus`: Vec3 scratch through shared helpers, mispredicting filter | [molecular-surface.ts:282](../../packages/core/src/math/geometry/molecular-surface.ts#L282) | time | medium | after P10: `lowest` −6-12 % (2.88-2.92 → 2.53-2.75 s); with P10 −9-16 %; bit-identical | S | confirmed | no |
| 4 | P3 | still-open | Default 32 elements per grid cell | [grid.ts:253](../../packages/core/src/math/geometry/lookup3d/grid.ts#L253) | time | medium | query −10-30 %; with P2 r 2 61 → 34 ms | S | confirmed | no |
| 5 | P12 | new | Molsurf torus pass ignores resolution; second `find` per atom | [molecular-surface.ts:398](../../packages/core/src/math/geometry/molecular-surface.ts#L398) | time | medium | exact: −5-7 %; coarse-quality skip: `lowest` 2.5 → 0.58 s on flat-32 (0.6 % of voxels change side) | M | confirmed | no |
| 6 | P13 | new | `accumulateRange` / `projectPointsRange` re-OSR and deopt per call | [cpu.ts:67](../../packages/core/src/math/geometry/gaussian-density/cpu.ts#L67) | time | low | CPU gaussian −8-10 % (7qpd 157-163 → 143-145 ms) | S | confirmed | no |
| 7 | 12/F1 | still-open; **fixed** #2046 | `floodfill3` through closures, whole-grid seed scan | [tensor.ts:494](../../packages/core/src/math/linear-algebra/tensor.ts#L494) | time | low | floodfill −28-40 % (220³ 312-352 → 220-224 ms) | S | confirmed | no |
| 8 | P8 | still-open | Gaussian `fasterExp` per voxel; separable kernel | [cpu.ts:112](../../packages/core/src/math/geometry/gaussian-density/cpu.ts#L112) | time | low | `accumulateRange` −20-25 % (v1) | S | confirmed | no |
| 9 | P9 | still-open; partly fixed #2046 (the miss returns `emptyArray`; per-unit CSR open) | `getEdgeIndices` `[]` per miss; edges as objects | [inter-unit-graph.ts:58](../../packages/core/src/math/graph/inter-unit-graph.ts#L58) | time+memory | low | valence model −10-20 % (CSR); emptyArray alone small | S | confirmed | no |
| 10 | P7 | still-open; **fixed** #2034 (= 07/P6) | Per-unit molsurf radius array sized `OrderedSet.end` | [molecular-surface.ts:48](../../packages/graphics/src/repr/structure/visual/util/molecular-surface.ts#L48) | memory | low | hundreds of MB of external memory on flat entries | S | confirmed | no |

### P10: Vec3 elements-kind contagion

- Trigger: `areBackboneConnected` (atomic-ranges.ts:40-41) runs even when the coarse fallback also returns -1. It then
  builds `Vec3.create(x[-1], y[-1], z[-1])` = `[undefined×3]`, and `Vec3.distance(...)` is NaN, so the result is false.
  On 1tup it runs at model creation (getAtomicRanges :98 cyclic check, :123 per residue); 7qpd has no such residue.
- Spread: the stores in `Vec3.create`/`set`/`add`/… now have feedback for both the double and the generic map.
  Optimized code transitions every double receiver to `PACKED_ELEMENTS`. Checked right after loading flat-32,
  `Vec3.create(0.5, 1, 2)` already returns a generic array (`v/real.ts`; double when sanitized, double on 7qpd). In a
  minimal bundle even a fresh `Vec3()` is generic after the trigger (`v/mech2.ts`). From then on each double stored
  into a Vec3 is a new 16 B HeapNumber, for the rest of the session.
- Evidence (`bench-contagion*.ts`, `SANITIZE=1` maps undefined to NaN in `Vec3.create`, flat-32):

  | Stage | Contagion | Sanitized |
  | --- | --- | --- |
  | molsurf res 1.3/24, sampled alloc | 1420 MB (projectTorus 1358, normalToLine 57) | 5 MB |
  | molsurf time (auditor, quieter run) | 6.32-6.74 s | 6.12-6.19 s |
  | unit boundaries ×256 | 62-72 ms, 32 MB | 26-31 ms, 4 MB |
  | computeInteractions | 573-731 ms, 429-432 MB | 515-669 ms, 278-280 MB |
  | intra bonds | 121-133 ms | 121-132 ms (no Vec3 on that path) |

- Fix: (1) root: `if (eiStart === -1 || eiEnd === -1) return false;` after the coarse fallback (same boolean). (2)
  Hardening: in the hottest kernels use scalar locals or `Float64Array` scratch instead of Vec3 helpers (P11 for
  projectTorus); typed arrays never change their elements kind. (3) Guard: a heavy-lane child process run with
  `--allow-natives-syntax` that loads 1tup-flat-32 and asserts `%HasDoubleElements(Vec3.create(0.5, 1, 2))`. Keep
  Vec3 as `number[]`.
- Risk: none for (1). Any other undefined or object that reaches a fresh Vec3 starts the contagion again (an
  out-of-range `fromArray`, `Vec3.create(...json)`), hence the guard. App-side sources other than this one were not
  surveyed; with the sanitizer, load, boundary, lookup, bonds and molsurf feed no other generic array (kinds8.ts).

### P2: grid query coordinate layout (rest)

- Mechanism: #1992 resolved the indices once, but each candidate still loads `px[idx]`, `py[idx]`, `pz[idx]` (and
  `radius[idx]`) in element order, a scattered gather per bucket.
- Evidence (`bench-grid.ts`, rerun, flat-32, one find per atom, same hits): r 2 60.9 → 44.5 ms, r 4 108.9 → 78.4 ms, r
  6.4 193.7 → 144.2 ms.
- Decision to revisit: the first commit of #1992 had exactly these copies (1.5-1.7×). The PR dropped them for "no
  extra memory". The remaining gain costs +12 B per element (+4 with radii), about 2-4 MB at 186k per cached lookup.
  Raise it with the maintainers rather than re-landing it silently.
- Fix: `bx/by/bz` (and `br`) in bucket order in `_build`, as `Float32Array` for a Float32 source, which is exact.
  `query`, `approxQueryNearest` and `queryNearest` scan them linearly.
- Risk: coordinate edits after the build are not seen (`Unit.remapModel` drops `lookup3d`).

### P11: projectTorus scratch and filter

- Mechanism: per torus, `Vec3.normalize`/`normalToLine`/`cross`/`scale` on closure-held `atob`, `n1`, `n2`. Under P10
  that is ~24 HeapNumbers per torus. The candidate filter (:331-344) keeps about half of ~44 neighbours behind a
  data-dependent branch and writes `tCount`, a closure context slot, per kept neighbour.
- Evidence (`v/bench-p11.ts` = bench-msf3 with the sanitizer, flat-32 res 8/12): orig 2.88-2.92 s, `ms-fix3.ts`
  (scalar a→b, n1, n2, mid; branchless filter, local counter) 2.53-2.75 s. Unsanitized orig 3.00-3.03 s. Bit-identical:
  `maxFieldDiff 0 idDiff 0` on 7qpd res 0.5/36 and flat-32 res 1.3/24.
- Fix: scalar locals computed exactly as the helpers compute them. Write `tX/tY/tZ/tRSq[tc]` unconditionally and add
  `(d² < reach² && nbIndex[j] !== b) ? 1 : 0` to a local counter. Keep molecular-surface.test.ts as the oracle.
- Risk: none measured. This also makes the kernel immune to P10.

### P3: default grid cell count

- Evidence (rerun, flat-32, hits identical): count 8, current query: r 2 60.9 → 41.8, r 4 108.9 → 83.9, r 6.4 193.7 →
  173.0 ms. On top of P2's copies: 44.5 → 33.6, 78.4 → 65.7, 144.2 → 138.7 ms.
- Fix: default 8, or ~4 Å cells for `unit.lookup3d`.
- Risk: result order changes, so bond and contact orders change (equivalent graphs); pinned lookup and bond snapshots
  need regenerating.

### P12: torus pass vs resolution

- Mechanism: the torus count is set by the atom pairs (~19 per atom on flat-32), not by the grid. `projectToriiRange`
  repeats `projectPointsRange`'s `find` + `loadNeighbours` for the same atom and radius (:401 vs :204).
  `loadNeighbours` still calls `OrderedSet.getAt` per neighbour (:126).
- Evidence (`bench-nt.ts`, rerun, flat-32 res 8/12, with P11): 2.50-2.51 s with the torus pass, 0.57-0.60 s
  without; 208 of 35,904 voxels change side of the iso value (0.6 %).
- Fix: (1) exact: keep each atom's neighbour ordinals from projectPoints (CSR, freed at return) and reuse them, and
  resolve them with the grid's resolved indices. Store ordinals only (4 B each) and re-gather coordinates: a
  structure-level surface at 186k atoms × 44 neighbours with coordinates would be ~230 MB. (2) Output-changing: skip
  or thin the torus pass when resolution > ~2× probe radius (`lower`, `lowest`). It needs a visual check and an entry
  in the repr quality table.
- Risk: (1) none. (2) Coarse molsurf snapshots change.

### Low

- **P13**: `accumulateRange` OSR-enters at an inner loop and deopts on leaving it ("exit from OSR'd inner loop"), once
  per call, i.e. per unit. Rerun (`bench-gauss`, 7qpd res 0.5): `--trace-deopt` OSR lines 212 → 3. Time 157-163 →
  143-145 ms with the per-atom kernel split out, `min[0..2]` and `id ? id[i] : i` hoisted (`g-split.ts`). Same in
  `projectPointsRange`. CPU gaussian runs headless or on GPU fallback only.
- **12/F1**: `floodfill3` makes one `o(x,y,z)` closure call per voxel for seeding, and `getCoords` plus 6 `o(...)` per
  visited voxel. Rerun (`bench-flood`): 220×227×233, 9.6M visited: 312-352 → 220-224 ms strided (seed the 6 faces,
  `o ± s`, Int32Array stack), visited diff 0. Optionally materialize the filled field so marching cubes keeps strides.
  Opt-in floodfill only.
- **P8**: per-atom 1D exp tables and row skip. v1: 155-173 → 122-137 ms, max |Δ| 0.057 (the `fasterExp` error).
  Combine with P13.
- **P9**: return the module `emptyArray`. The real cost is `Map.get` per vertex; per-unit CSR is the fix, and it is a
  public API change (`breaking-v6-changes.md`).
- **P7**: build compact x/y/z/r with `indices = Interval(0, n)` and `id = ordinal`, as `GaussianDensityCPU` does for
  radii.

## Bugs

| Id | vs v1 | Severity | Title | File | Verdict |
| --- | --- | --- | --- | --- | --- |
| B5 | new; **fixed** #2025 | medium | `areBackboneConnected` reads coordinates at index -1 (root of P10) | [atomic-ranges.ts:40](../../packages/model/src/model/structure/model/properties/utils/atomic-ranges.ts#L40) | confirmed |
| B1 | still-open; **fixed** #2033 | low | `approxNearest` returns the nearest last; callers read [0] | [grid.ts:469](../../packages/core/src/math/geometry/lookup3d/grid.ts#L469) | confirmed |
| B2 | still-open; **fixed** #2034 | low | `updateChunk = 1e5/(maxRadius⁹·scale)`: no yield for small radii | [molecular-surface.ts:462](../../packages/core/src/math/geometry/molecular-surface.ts#L462) | confirmed |
| B3 | still-open; **fixed** #2046 | low (latent) | `Sphere3D.copy` / `expandBySphere` keep stale extrema | [sphere3d.ts:59](../../packages/core/src/math/geometry/primitives/sphere3d.ts#L59) | confirmed (code) |
| B4 | still-open; **fixed** #2046 | low (latent) | RMSD state reuse aliases and overwrites a caller's centre | [minimize-rmsd.ts:92](../../packages/core/src/math/linear-algebra/3d/minimize-rmsd.ts#L92) | confirmed (code) |
| 12/B2 | still-open; **fixed** #2046 | low | `getCoords` (rank ≥ 4) writes `coords['undefined']` | [tensor.ts:379](../../packages/core/src/math/linear-algebra/tensor.ts#L379) | confirmed (code) |

- **B5**: a polymer residue with a trace atom but neither backbone N/C nor the coarse backbone atom gives
  `eiStart`/`eiEnd = -1`. The function then reads `x[-1]` (undefined) and compares NaN, which gives the right false by
  accident. Its consequence is P10: GB-scale extra allocation and slower Vec3 code for the session. Fix: return false
  after the fallback. It belongs to slice 05 too.
- **B1**: in the first non-empty cell `approxQueryNearest` adds each improving candidate, so `indices` runs in
  decreasing distance with the nearest last. `StructureLookup3D._approxNearest` appends across units the same way.
  Failure: atoms at 2.9 and 0.3 Å in one cell give `[far, near]`. The external-structure theme colours from
  `rf.indices[0]`, the farther atom, and blob k-means assigns to the farther centroid. Fix: track the best and add it
  once.
- **B2**: maxRadius 0.8 at res 0.5 gives a chunk of 373k atoms: a structure-level surface runs in one synchronous
  chunk, with no progress and no abort. Fix: chunk by work, `max(1, ceil(2e6 / (2·maxRadius·scale + 2)³))`, in both
  files.
- **B3**: `copy(out, a)` with `a` lacking extrema keeps `out`'s old extrema, and `includeSphere`, `radiusSphere`,
  `Box3D.fromSphere3D` and `includes` prefer them. No exercised path: `clip-volumes.ts:142` passes a fresh `out`. Fix:
  `else out.extrema = undefined`, and clear in `expandBySphere`.
- **B4**: create a state with `data1.centerA = v`, then reuse it with `data2` without centres: `computeCenter` writes
  `data2`'s centroid into the caller's `v`. tm-align and superposition-ligand never pass centres. Fix: the state owns
  its Vec3s and copies into them; no write-back into `data`.
- **12/B2**: `coords[ao[d + 1]]` is past the axis order, so it writes key `undefined` (a map change on the caller's
  array). Values are right. Fix: delete the line.

## Checked and dropped

- P6 (v1 and this run's "flatten the 98-direction pass"): with P10 sanitized, current 20-22 ms = flattened 20-22 ms
  over 256 units (`v/bb-san.ts`, box and sphere maxDiff 0). The v1 gain and the 00-profile's `Vec3.dot` share were
  contagion. The flattened prototype is still slow under contagion (28-30 ms: it builds Vec3s in the radius pass), so
  "makes it immune" was not shown either.
- P10 "the mechanism is only polymorphic store transitions": that holds in the real load. In a minimal bundle the
  trigger also flips the allocation site of `Vec3.zero`'s literal, so fresh Vec3s are born generic (`v/mech2.ts`; not
  reproduced in plain JS, `v/site*.mjs`). Either way the fix is the same.
- P11 "allocation 1.36 GB → 5 MB": that is P10's gain. P11 alone is the time item (impact high → medium).
- P10 "StructureLookup3D's pivot turns generic": not seen in this run (`pivot D` after boundaries). It is not needed
  for the finding.

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/1f4e2fbb-da6c-4948-aacb-385a4dab465c/scratchpad/perf-audit-2/core-math/`,
bundled with `node build.mjs <in.ts> <out.mjs>` (esbuild, `molstar-src`). Fixtures in `load.ts` (`7qpd`, `flat-32`,
`flat-192`); natives checks need `--allow-natives-syntax`.

| Script | What |
| --- | --- |
| `patch3.ts` (`SANITIZE=1`), `patch2.ts` | P10/B5: sanitize `Vec3.create`; report the first generic array passed to a Vec3 helper |
| `kinds*.ts`, `v/real.ts`, `v/mech.ts`, `v/mech2.ts`, `v/site*.mjs` | P10: elements kinds after load / trigger |
| `bench-contagion.ts` (`molsurf`/`boundary`/`asa`), `bench-contagion2.ts` | P10: time and sampled allocation with and without the sanitizer |
| `bench-msf3.ts` + `ms-fix3.ts`, `v/bench-p11.ts` | P11: orig vs scalar/branchless projectTorus; `CMP=1` field diff |
| `bench-nt.ts` + `ms-nt.ts` (`NOTORII=1`) | P12: molsurf without the torus pass, voxels that change side |
| `bench-grid.ts` + `grid-opt.ts` | P2/P3: current vs bucket-ordered copies, count 32/8 |
| `bench-boundary.ts`, `v/bb-san.ts`, `bench-boundary2.ts` | P6: current vs flattened, with and without the sanitizer |
| `bench-gauss.ts` + `g-orig.ts`/`g-split.ts` | P13: OSR deopts and time |
| `bench-flood.ts` | 12/F1: floodfill vs strided |
| `heapsample.ts` | sampled allocation per function |
