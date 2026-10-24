# Slice 08: graphics geo (v3)

Scope: `packages/graphics/src/geo`, which covers marching cubes, `MeshBuilder` and its adders, the
Cylinders/Lines/Spheres/Text builders, mesh color smoothing, group mapping, the font atlas and the text layout (+
`core/util/text-measure.ts`, which the text builder calls). Out of scope: GPU gaussian density and `gl/` (frozen). HEAD
`16277d822`. v2 is [../perf-audit-2/08-graphics-geo.md](../perf-audit-2/08-graphics-geo.md).

Measurement conditions: Node 26.12.1 (v2 used 25.6), Apple M4. The machine ran at **load 23–31 on 10 cores** for the
whole session, with the other auditors. Wall times were inflated 3–10× and swung by 2× between back-to-back runs, so
they carry no weight on their own. The evidence is therefore allocation (`total_allocated_bytes`, sampled heap),
counters, `--trace-deopt`, and A/B ratios measured in the same process or back to back. The base for the A/B is
`0944823d4`, built with `git archive` and bundled with the same script and an alias plugin. Benches are in
`/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/31ee827f-4715-4967-bec3-085f5013d11f/scratchpad/pa3/08/`
(see the end of this file).

PRs since v2 that touched the slice:

- #2028: one font atlas per font, LRU of 4 (B6).
- #2030: font table with lazy web fonts, kerned layout (`TextMeasure`), ink-sized glyph cells, atlas texture version.
- #2033: P16 cylinder fast path, B7.
- #2034: B3.
- #2044: B2, B4.
- #2048: P9, P11, B5.

## Hot paths now

| Path                                          | Cost on fixture                                                                                                                                                            | Where the time goes                                                                                                                                                                 |
| --------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| MC `interpolate`/`stridedRow`/`processCell`   | 120³ blob field (145k vertices): 17–20 MB alloc per run after the first; 4.4–6.5 MB with the state map kept alive. 7qpd `repr:gaussian-surface` median alloc 77 MB (fix: 45) | optimized code dropped at every major GC (P15), unchanged since v2; boxing in `interpolate`                                                                                         |
| `addTube`/`addSheet` (cartoon)                | 3000 tubes: alternating 29–46 ms / 35–38 MB and 11–14 ms / 0 MB; with an anchor 6.5–8.5 ms / 0 MB every run. 7qpd `repr:cartoon` median alloc 70.7 MB (fix: 33.1)          | P12, unchanged since v2                                                                                                                                                             |
| `calcMeshColorSmoothing`/`splatVertex`        | 145k vertices, 2000 groups: 81–116 ms against a cached-CSR recolor of 47–81 ms (under load; v2: 35 against 11 ms)                                                          | P1, unchanged                                                                                                                                                                       |
| `FontAtlas.draw` (glyph SDF)                  | 7–10 ms per glyph at the default `fontQuality` 3, 12–18 ms at 4; 94 ASCII glyphs take 0.66–0.97 s (base ≈ HEAD)                                                            | `edt1d` 65 %, `edt` 19 %, the alpha/`Math.pow` loops in `draw` 11 % (P18)                                                                                                           |
| `TextBuilder.add`                             | 50k labels of ~14 chars: 58 MB allocated, against 18 MB at base; time +10–30 %                                                                                             | `TextMeasure.layout` per label: a `split('')` array, a `Float64Array`, the result object, and a fresh pair string per grapheme for the width map (P20)                              |
| `CylindersBuilder.add`                        | 300k cylinders: 135–211 ms against base 208–323 ms (same process kinds, under load)                                                                                        | fixed (#2033); the remaining time is the caller's boxed doubles in the bench                                                                                                        |

## v2 → now

| v2 id | Title                                                          | Status            | PR                       | Note                                                                                                                                                                                                                                                    |
| ----- | -------------------------------------------------------------- | ----------------- | ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| B6    | font atlas cache keyed by every Text prop                      | fixed             | #2028, #2030             | keyed by the CSS font string; LRU of 4 (`font-atlas.ts:24`). Re-read after #2030: the key now includes the family's font list from `FontFamilyFonts`, so `setFontFamilyFonts` makes a new atlas (by design)                                            |
| P15   | MC re-tiers after every major GC                               | still open        |                          | `algorithm.ts` unchanged. Re-run on Node 26: "embedded weak objects cleared" for `interpolate` ×8, `processCell` ×8, `stridedRow` ×6, `get3dOffsetFromEdgeInfo` ×6, `slice` ×5, `processRow` ×4 per 8 runs. Alloc 17–20 → 4.4–6.5 MB with the last state kept |
| P12   | cartoon tube/sheet code dropped at every major GC              | still open        |                          | `mesh-builder.ts:120` unchanged. Same deopt (`addTube` ×2, plus `add3`, `addRepeat` of ChunkedArray). Anchor: 35–38 MB → 0 MB per cold build                                                                                                          |
| P1    | recolor re-splats color-independent weights                    | partly fixed (v1) | –                        | unchanged since v2; CSR cache still ~1.5–3× (ratio under load 1.4–1.7×; v2 idle 3×)                                                                                                                                                                    |
| B8    | smoothed `volumeInstance` layers not recomputed on a transform | still open        |                          | `Visual.setTransform` (`visual.ts:794`) still rewrites `aTransform` and bounds only; layers are still set before the transform (`units-representation.ts:372` before `:402–409`)                                                                                                               |
| P16   | `CylindersBuilder.add` 36 adder calls                          | fixed             | #2033                    | one capacity check per array, direct writes (`cylinders-builder.ts:121`). Reviewed: correct, with a fallback when any chunk ends. A/B 300k: ~1.4–1.5× faster under load                                                                               |
| P7    | cylinder impostor attributes ×6                                | still open        |                          | frozen (gl schema)                                                                                                                                                                                                                                      |
| P6    | constant mapping/index buffers per geometry                    | still open        |                          | `cylinders-builder.ts:258`, `lines-builder.ts:112`                                                                                                                                                                                                      |
| P17   | one polymorphic ChunkedArray adder site                        | still open        |                          | P16 took cylinders off it; Lines/Text/Mesh/MC still use it (see P19)                                                                                                                                                                                    |
| P9    | cylinder primitive cache unbounded, 32-bit hash                | fixed             | #2048                    | buckets compared on all 9 values, 256 entries, FIFO by hash (`builder/cylinder.ts:23–24`). Reviewed: correct; FIFO can evict the hot b&s primitive during a gizmo drag (a few µs to rebuild, not a finding)                                              |
| B7    | dashed cylinder mutates `props.topCap`                         | fixed             | #2033                    | copy for the last dash                                                                                                                                                                                                                                  |
| P11   | Spheres `positionGroup` filled twice                           | fixed             | #2048                    | `fromArrays` no longer fills. Reviewed: `update` still fills (needed: `updateValues` refills only on an LOD/size change), and `createValues` fills, so no reader sees an empty texture                                                                |
| P10   | `smoothEdges` rebuilds neighbour maps                          | still open        |                          | `mesh.ts:683–697`                                                                                                                                                                                                                                       |
| B1    | `fillEdges` never adds a triangle                              | still open        |                          | `mesh.ts:566` `< 2`                                                                                                                                                                                                                                     |
| B2    | `createGroupMapping` sized by `dataCount`                      | fixed             | #2044                    | `util.ts:204`. Reviewed: `arrayMax` scans all of `groups`, including the stale tail of a reused buffer (`compactTrimmed` keeps a ≥80 % full chunk whole), so it can only oversize; no wrong mapping                                                  |
| B3    | `getTrilinearlyInterpolated` 0/0                               | fixed             | #2034                    |                                                                                                                                                                                                                                                         |
| B4    | MC lines ids in a Float32Array                                 | fixed             | #2044                    |                                                                                                                                                                                                                                                         |
| B5    | `computeGaussianDensityTexture3d` dispatches '2d'              | fixed             | #2048                    |                                                                                                                                                                                                                                                         |
| 13/P3 | unshared 5M allocation, builder part                           | partly fixed (v2) |                          | the geo part is done; P12/P15 boxing remains                                                                                                                                                                                                            |
| 14/N4 | cartoon meshes reserve ~1.4×                                   | still open        |                          | caller's estimate (repr, slice 07)                                                                                                                                                                                                                      |

Counts: fixed 9 (B2, B3, B4, B5, B6, B7, P9, P11, P16); partly fixed 2 (P1, 13/P3); still open 9 (P6, P7, P10, P12,
P15, P17, B1, B8, 14/N4); regressed 0; dropped 0.

New bugs or costs from fixes since v2:

- **#2030 (kerned layout)**: costs allocation (P20). `TextBuilder.add` allocates 3.2× what it did (18 → 58 MB for
  50k labels) and runs 10–30 % slower. The layout itself is correct. The cell `left` offset matches the old `buffer`
  shift for zero-bearing glyphs, and the width is the sum of the measured advances.
- **#2030 (texture version)**: `FontAtlas.get` bumps `texture.version` for each new glyph. With shared atlas
  textures, any label build that draws a new glyph re-uploads the whole atlas (43 MB R8 at quality 3), at most once per
  frame. That is in gl (frozen), recorded only.
- No new correctness bug was found in #2028, #2030, #2033, #2034, #2044 or #2048.

## Open performance findings (ranked)

| Id  | Status     | Title                                                                                       | Where                                                                                                                                                                                            | Kind        | Estimate                                                                                                    | Effort | Frozen |
| --- | ---------- | ------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ----------- | ----------------------------------------------------------------------------------------------------------- | ------ | ------ |
| P15 | still open | MC optimized code dropped at every major GC (per-run `MarchingCubesState` map)              | [algorithm.ts:117](../../packages/graphics/src/geo/util/marching-cubes/algorithm.ts#L117)                                                                                                       | time+memory | measured: 7qpd gaussian alloc 77 → 45 MB (median of 7); 120³ MC alloc 17–20 → 4.4–6.5 MB; time v2 idle −8 % to 2× | S      | no     |
| P12 | still open | cartoon `addTube`/`addSheet` dropped at every major GC (`MeshBuilder.State` map)            | [mesh-builder.ts:120](../../packages/graphics/src/geo/geometry/mesh/mesh-builder.ts#L120)                                                                                                       | time+memory | measured: 7qpd `repr:cartoon` alloc 70.7 → 33.1 MB; 3000 tubes cold 29–46 ms / 36 MB → 6.5–8.5 ms / 0 MB      | S      | no     |
| P1  | partly     | recolor re-splats the color-independent weights                                             | [color-smoothing.ts:414](../../packages/graphics/src/geo/geometry/mesh/color-smoothing.ts#L414)                                                                                                 | time        | measured v2 (idle) ~3×; re-run under load 1.4–1.7×                                                          | M      | no     |
| P18 | new        | glyph SDF draw 7–10 ms per glyph at the default quality: first labels block 0.5–1 s         | [font-atlas.ts:338](../../packages/graphics/src/geo/geometry/text/font-atlas.ts#L338), [distance-transform.ts:97](../../packages/core/src/math/geometry/distance-transform.ts#L97)              | time        | measured under load: 94 glyphs 0.66–0.97 s (q3); idle likely 2–3× less (estimated)                          | M      | no     |
| P6  | still open | constant mapping and index buffers per geometry                                             | [cylinders-builder.ts:258](../../packages/graphics/src/geo/geometry/cylinders/cylinders-builder.ts#L258)                                                                                        | memory      | v2: −120 B JS heap per cylinder                                                                             | S      | no     |
| P20 | new        | `TextMeasure.layout` allocates per label (array, Float64Array, pair strings)                | [text-measure.ts:95](../../packages/core/src/util/text-measure.ts#L95), [text-builder.ts:64](../../packages/graphics/src/geo/geometry/text/text-builder.ts#L64)                                 | memory+time | measured: 50k labels 58 → ~20 MB allocated (base 18), time −10–25 % (estimated from the A/B)                 | S      | no     |
| P19 | new        | `LinesBuilder.add` 12 ChunkedArray adder calls per line (P16's pattern)                     | [lines-builder.ts:41](../../packages/graphics/src/geo/geometry/lines/lines-builder.ts#L41)                                                                                                      | time        | estimated ~1.5–3× on the add phase (P16 analogy); wireframe isosurfaces, `line` repr                        | S      | no     |
| P7  | still open | cylinder impostor attributes ×6                                                             | [cylinders-builder.ts:121](../../packages/graphics/src/geo/geometry/cylinders/cylinders-builder.ts#L121)                                                                                        | memory      | v2: buffers ~5× smaller                                                                                     | L      | yes    |
| P17 | still open | ChunkedArray adders: one polymorphic site for all builders                                  | [chunked-array.ts:59](../../packages/core/src/data/util/chunked-array.ts#L59)                                                                                                                   | time        | v2: ~1.5× variance                                                                                          | M      | no     |
| P10 | still open | `smoothEdges` rebuilds neighbour maps per iteration                                         | [mesh.ts:683](../../packages/graphics/src/geo/geometry/mesh/mesh.ts#L683)                                                                                                                       | time+memory | v2: several × with CSR                                                                                      | M      | no     |

### P15 and P12: optimized code dropped at major GC (re-checked)

- **Status**: still open; neither file changed since v2. Re-run on Node 26.12.1 with `--trace-deopt`:
  - `mc.mjs`, 8 MC runs with a gc between them: "embedded weak objects cleared" for `interpolate` ×8, `processCell` ×8,
    `stridedRow` ×6, `get3dOffsetFromEdgeInfo` ×6, `slice` ×5 and `processRow` ×4, plus `add` and `add3` of
    ChunkedArray.
  - `tube.mjs`: the same for `addTube` ×2, `add3` and `addRepeat`.
  - The ChunkedArray adders show it too. Their code embeds the `ChunkedArray.create` literal map, and no ChunkedArray
    outlives a build: `compactTrimmed` returns typed arrays. So every builder in the app re-tiers after a major GC, not
    only MC and cartoon.
- **Evidence**:
  - MC 120³, 145k vertices: 17–20 MB per run after the first (`interpolate` 15–17 MB). With the last state kept
    (`mc-last.mjs`): 4.4–6.5 MB.
  - Tubes: cold builds 35–38 MB. With `Anchor = MeshBuilder.createState(1, 1)` (`tube-a.mjs`): 0 MB from run 1 on,
    6.5–8.5 ms against 11–46 ms.
  - Pipeline 7qpd, 7 runs, base bundle against a bundle with both patches (`bench-base.mjs`, `bench-fix.mjs`), median
    allocation:
    - `repr:cartoon`: 70.7 → 33.1 MB.
    - `repr:gaussian-surface`: 77 → 44.8 MB.
    - `repr:molecular-surface`: 252 → 264 MB (noise, no gain).
    - Times unusable at load 25–31.
- **Fix** (as v2): module-level anchors.
  - `const StateMapAnchor = createState(1, 1)` in mesh-builder.
  - A pooled or anchored `MarchingCubesState`: drop the field references after `run()`, and give every field its final
    representation.
  - Also a `ChunkedArray.create` anchor in core (slice 01). Pin all three with an allocation check for repeated small
    builds with a gc between them.

### P1: recolor re-splats (re-checked)

Unchanged since v2 apart from the B3 hunk. `cs.mjs` under load: current 81–116 ms against a cached-CSR recolor of
47–81 ms. The CSR build costs 137 ms and 13.4 MB (34.5 nnz per sample). v2's idle numbers (35 against 11 ms, about 3×)
remain the reference. The fix and its risk are as in v2. It needs the same invalidation hook as B8.

### P18 (new): glyph SDF drawing cost

- **Evidence** (`glyphs.mjs`): a stub canvas with 0.6 em advances and a filled ink box, 94 printable ASCII glyphs into a
  fresh atlas.

  | `fontQuality` | Per glyph, base | Per glyph, HEAD |
  | ------------- | --------------- | --------------- |
  | 0             | 0.40–0.44 ms    | 0.46–0.49 ms    |
  | 2             | 3.8–4.2 ms      | 3.9–4.1 ms      |
  | 3 (default)   | 7.0–8.1 ms      | 8.0–10.4 ms     |
  | 4             | 12–19 ms        | 17–18 ms        |

  HEAD ≈ base, so #2030 did not regress it. CPU profile: `edt1d` 65 %, `edt` 19 %, `draw` 11 % (the alpha → grid loop
  with `Math.pow`).
- **Cause**: the cell is about 218 × 405 px at quality 3. Each glyph runs 2 full EDTs, so about 0.7M `edt1d` cell
  visits. This runs on the main thread when the first labels are built, and again for every atlas the LRU rebuilds. The
  first label build with a web font draws the atlas twice: once with the fallback, once when the font arrives
  (`loadFont` drops the cached atlas, labels.md §3). A typical label set has 40–60 distinct glyphs: 0.3–0.6 s at load,
  roughly 0.15–0.3 s idle (estimated).
- **Fix options**:
  - EDT only over the ink bounding box plus `radius`, and fill the rest of the cell from the outside distance.
  - Skip the outer EDT for rows or columns with no ink.
  - Draw glyphs in a worker with OffscreenCanvas.
  - Prefetch the ASCII set in idle time after the first atlas is created.
  - `edt1d` takes `NumberArray`: callers pass Float64/Int16 (atlas) and Float64/Int32 (`squaredDistanceTransform3D`),
    so the keyed loads are polymorphic once both run. A typed copy per caller is a small extra win (slice 02).

### P20 (new): per-label layout allocation (#2030)

- **Evidence** (`text.mjs`, 50k labels like `ALA 123 | A`, stub canvas, base against HEAD in back-to-back processes):

  | Run          | Allocated per run | `add` phase    |
  | ------------ | ----------------- | -------------- |
  | base         | 18.0–18.4 MB      | 1.00–1.39 s    |
  | HEAD (#2030) | 57.7–58.0 MB      | 1.13–1.75 s    |

  The extra ~0.8 KB per label comes from `TextMeasure.layout`: `split('')`, a `Float64Array` (with its ArrayBuffer),
  the result object, and one `previous + grapheme` concatenation per grapheme to look up the pair width. `measureText`
  calls stay bounded (249 for the whole run), so the cache works.
- **Fix**:
  - A layout into reusable scratch (a grow-only `Float64Array` of offsets plus a count) for the builder.
  - Pair widths keyed by code-unit pair for ASCII (`(a << 8) | b` in a `Map<number, number>`, or a 128² `Float64Array`
    per atlas), with no string concatenation.
  - Keep `splitGraphemes` for the non-ASCII path only.
- **Scale**: atom labels on flat-32 (186k) allocate about 150 MB more per build. That is low until labels go
  per-atom on large structures.

### P19 (new): `LinesBuilder.add` per-vertex adders

`add`/`addVec` (:41, :57) make 12 ChunkedArray adder calls per line, the pattern P16 removed from cylinders. The A/B in
`cyl.mjs` shows lines at the same per-item cost as the old cylinders path (74–350 ms for 300k lines under load).
Callers: MC lines (wireframe isosurface: `MarchingCubesLinesBuilder.get` re-adds every line, builder.ts:110),
`line`-type structure visuals and cages. The fix is P16's: one capacity check per array, then direct writes. The gain is
estimated, not measured idle.

## Bugs

| Id  | Severity | Status     | Title                                                                         | Where                                                                                                   | Failure scenario                                                                                                                                              |
| --- | -------- | ---------- | ----------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| B8  | medium   | still open | smoothed `volumeInstance` layers/colors not recomputed on a transform change | [visual.ts:794](../../packages/graphics/src/repr/visual.ts#L794)                                        | overpaint a chain of a gaussian surface, then move the structure (Prism transform, `setState({ transform })`, unitTransforms): the paint stays at the old place |
| B10 | low      | new        | a full font atlas redraws a glyph that does not fit on every `get`            | [font-atlas.ts:302](../../packages/graphics/src/geo/geometry/text/font-atlas.ts#L302)                   | labels with more than ~320 distinct CJK glyphs: each later occurrence of an overflowing glyph costs a full SDF draw (~20 ms under load) and a console warning  |
| B1  | low      | still open | `fillEdges` never adds a triangle                                             | [mesh.ts:566](../../packages/graphics/src/geo/geometry/mesh/mesh.ts#L566)                               | `includeParent` surfaces keep edge holes                                                                                                                      |

### B10 (new): full atlas, redrawn on every get

- **Cause**: `get` draws the glyph into the scratch canvas first. It then finds no room
  (`currentY + scratchH > height`), warns "canvas to small" and returns the placeholder, without recording
  `mapped[char]`. So the next `get` of that char draws again: `measureText`, `getImageData`, 2 EDTs. The texture holds
  about 350 cells of 0.75 lineheight. One-em glyphs (CJK, 320 px cells at quality 3) fit 320.
- **Probe** (`atlasfull.mjs`, 500 distinct CJK-width glyphs, 3 passes):
  - Pass 0: 500 draws, 320 mapped, 181 warnings.
  - Passes 1 and 2: 181 draws and 181 warnings each, 3.3–3.6 s per pass under load.
  - Every later label build pays this again for the session.
- **Not a regression**: the defect predates v2. #2030 widened cells to the ink (`maxWidth` 0.75 → 1.5 lineheight), so
  emoji and wide glyphs fill the texture sooner.
- **Fix**:
  - `this.mapped[char] = this.placeholder` on overflow, and warn once per atlas.
  - Or grow the texture by a row block: a new array, with the version bump that is already there.

### B8, B1

Unchanged from v2 (see there for the evidence). B8 now also blocks the P1 cache, which needs the same "transform
changed → re-smooth" hook.

## Already good, re-checked

- `CylindersBuilder.add` fast path (#2033): the alignment holds. Chunk sizes are ×6 elements, reused buffers are
  `cylinderCount·6`, and any chunk end falls back to the generic adders. Output is identical (test
  `cylinders-builder.test.ts`).
- Cylinder cache (#2048): exact 9-value compare, so no hash collisions. The eviction count stays consistent when the
  evicted key is the one being added.
- Font loading (#2030):
  - `loadFont` returns nothing in Node and for loaded or system fonts.
  - `settle` drops only the atlas of the arrived font.
  - Failed fonts are remembered.
  - `clearFontAtlasCache` empties the LRU list, which is the only store.
  - `TextMeasure` caps at 65 536 widths per atlas (≤ ~4 MB).
- `createGroupMapping` (#2044), `getTrilinearlyInterpolated` (#2034), MC lines `Uint32Array` (#2044),
  `computeGaussianDensityTexture3d` (#2048): correct.

## Not audited / not measured

- Idle wall times: the machine stayed at load 23–31, so every time above is a range under load, and the ratios rest on
  allocation or same-process A/B.
- Browser `measureText` and `getImageData` cost (P18 and P20 use a stub canvas), and GPU re-upload of the atlas texture
  per new glyph (gl, frozen).
- flat-32/asm-192 pipeline stages: skipped under this load.
- texture-mesh color smoothing (unchanged since v2), `image/`, `direct-volume/`.

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/31ee827f-4715-4967-bec3-085f5013d11f/scratchpad/pa3/08/`.
Bundle with `node buildp.mjs X.ts X.mjs [mc] [anchor] [base]`:

- `mc` keeps the last MC state alive.
- `anchor` adds the `MeshBuilder` state anchor.
- `base` resolves `@molstar/*` to `base/packages` (the `0944823d4` archive).

Run with `node --expose-gc X.mjs`.

| Script                              | What                                                                                                    |
| ----------------------------------- | ------------------------------------------------------------------------------------------------------- |
| `mc.ts`, `tube.ts`, `tube-a.ts`     | P15/P12: repeated MC 120³ and 3000-tube builds; `mc-last.mjs` and `tube-a.mjs` are the fixes; add `--trace-deopt` |
| `bench-base.mjs`, `bench-fix.mjs`   | pipeline `bench.ts` bundled without and with both patches; run with `--fixture <abs path> --reprs cartoon,gaussian-surface,molecular-surface` |
| `cs.ts`                             | P1: `calcMeshColorSmoothing` against a cached-CSR recolor                                               |
| `glyphs.ts`                         | P18: glyph drawing per `fontQuality`, base against HEAD                                                 |
| `text.ts`                           | P20: `TextBuilder.add` for 50k labels, base against HEAD                                                |
| `atlasfull.ts`                      | B10: redraws once the atlas is full                                                                     |
| `cyl.ts`                            | P16 (fixed), P19: Cylinders and Lines builder add phase, base against HEAD                              |
