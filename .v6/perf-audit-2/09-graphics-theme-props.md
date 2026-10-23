# Slice 09 v2: `graphics/src/theme`, `graphics/src/props`

Scope: the color and size themes (about 40 color, 6 size), palettes (`ColorTheme.Palette`, core `getPalette`), the layer
types (overpaint, transparency, clipping, substance, emissive, wiggle) in
[`theme/`](../../packages/graphics/src/theme); the interactions representation (intra-unit, inter-unit and bridge
visuals), the ASA, interaction-type, SIFTS and cross-link themes and volume streamlines in
[`props/`](../../packages/graphics/src/props). HEAD `0944823d4`, 2026-10-22, Node 25.6, Apple M4, moderate shared load
(the verifier's re-runs were 1.5–2× slower than the auditor's; ratios held). v1:
[../perf-audit/09-graphics-theme-props.md](../perf-audit/09-graphics-theme-props.md). Fixtures and method:
[00-profile.md](00-profile.md). No finding touches Canvas3D or `gl`. The cylinder and mesh builders are slice 08.

## Hot paths now

| Path | Cost | Where it goes |
| --- | --- | --- |
| `repr:interactions` fresh build (impostors) | flat-32 155–168 ms, 82.8 MB buffers (v1 mesh 1.1 s, 422 MB); flat-192 1.2–1.4 s, 500 MB; asm-864 1.2 GB | `createLinkCylinderImpostors` 38–54 %, `addFixedCountDashes` 20–23 % (slice 08), frozen gl merge ~11 %, GC 11–16 %, `ChunkedArray._compact` 7 % (P15) |
| same, typed-array allocation | flat-32 436 MB, asm-192 1.28 GB, flat-192 2.63 GB for 83 / 265 / 500 MB kept | the first buffer chunk is 6× the exact size, then copied (P15) |
| interactions `sizeFactor` update | flat-32 64–66 ms, asm-192 171–173, flat-192 368–409 (alpha: 0.2–2.4 ms) | full geometry rebuild for a uniform (P16) |
| inter-unit interactions on assemblies | asm-192 254k contact edges, ~265 MB | 95 % are 1,258 within-copy edges repeated per operator, not instanced (P17) |
| `Theme.create` per representation update | ≤ 2 ms at 186k atoms | except `particle-hierarchy` 356–390 ms at 10M particles, `particle-entity` 40–44 ms (P6) |
| theme `color(location)` per 186k elements | 1–8 ms most themes | illustrative 21, operator-hkl 17 (per unit in practice), ASA 10 (P4) |
| `theme-apply:ball-and-stick` | – | chain-id closure ~20 %; the rest is `createGroupColor`/`LocationIterator` (slice 08) |
| layer merge/filter | 0.1 ms for 10 layers; 17–30 ms for 1000–3000 one-residue layers (flat-32) | fine |

Changed since v1: #2004 landed impostor interactions (P1) and fixed B1 and B4, moved H-bond anchoring to
`eachBondedAtom` (across units) and added interaction-location sizing to the physical and uncertainty size themes;
#2005 and #2012 fixed 13/P1 (layer smoothing over reached instances only). #2004 together with #2007 (08/P8) introduced
P15. No color-theme file changed since v1 (`git diff 60fe4c2f1..HEAD -- packages/graphics/src/theme/color` is empty).

## v1 → now

11 entries (09/P1–P6, 09/B1–B4, 13/P1): **fixed 4, partly fixed 0, still open 7, regressed 0, dropped 0**. P15 is a new
defect inside the P1 fix, recorded as new.

| v1 id | Title | Status | PR | Note |
| --- | --- | --- | --- | --- |
| P1 | Interactions have no cylinder-impostor visual | fixed | #2004 | impostor variants of all three visuals (`intra:252`, `inter:261`, `bridge:299`); flat-32 155–168 ms, 82.8 MB (v1 1.1 s, 422 MB); buffer sizing wrong, see P15 |
| P2 | Mesh path: dash caps 2/3 of vertices, estimate ~5× low | still-open | | `link.ts:144` estimate and `dashCap: true` (`link.ts:38`) unchanged; mesh fallback only now |
| P3 | 'set' legend built by O(count) string concatenation | still-open | | core `palette.ts:95-100` unchanged; 1M 46 ms / 53 MB, 10M 617 ms / 549 MB |
| P4 | Per-element `Color.lighten` and string keys in closures | still-open | | files unchanged; illustrative 20.8/48.8 ms, entity-id 6.8/14.7 (elem/bond, flat-32) |
| P5 | H-bond endpoint search per directed edge, closure per call | still-open | | `intra:115/:126`, `inter:103/:115`; now also walks inter-unit bonds; ~5 % of the build |
| P6 | `particle-hierarchy` rebuilds its O(N) index per update | still-open | | `particle-hierarchy.ts:134,139` unchanged; re-run 386–390 ms at 10M (v1 340–357) |
| B1 | intra `includeParent` filter reads `members[offsets[i]]` | fixed | #2004 | `intra:157,166,172` read `elements[members[i]]`; donor flag `types[a[edge]]` (:112); tests `interactions-visual.test.ts:356,372` |
| B2 | entity-source indexes polymer-only `sequences` by entity index | still-open | | `entity-source.ts:81`; re-run on `ent.cif`: throws "Cannot read properties of undefined (reading 'sequence')" |
| B3 | `ColorTheme.Palette` 'set' with offsets writes `undefined` | still-open | | `color.ts:128,140`; re-run: `[[blue,0],[red,1]]` → `[undefined, undefined]` |
| B4 | bridge visual rebuilds on every update (size theme by ref) | fixed | #2004 | `bridge:319-323` no longer compares `newTheme.size`; mesh rebuilds by `SizeTheme.areEqual` |
| 13/P1 | Instanced surfaces: layer smoothing spans every instance | fixed | #2005, #2012 | mesh `color-smoothing.ts:231` `getSplatInstances`; GPU path likewise; `repr/visual.ts:378` skips all-clear; asm-192-low 4.7 s → 20 ms |

## Open findings, ranked

| Id | vs v1 | Title | File | Kind | Impact | Estimate | Effort | Verdict | Frozen |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| P15 | new; **fixed** #2025 (= 07/P11) | Impostor interactions: exact cylinder count multiplied by 6 twice; 6× first chunk, then copied | [link.ts:391](../../packages/graphics/src/repr/structure/visual/util/link.ts#L391) | time+memory | medium | fresh build: typed-array alloc ÷4–6 (flat-192 −2 GB transient), build −15–30 %; **measured**: see 07/P11 | S | confirmed | no |
| P16 | new; **fixed** #2032 | Impostor interactions rebuild geometry on `sizeFactor`, a uniform there | [interactions-intra-unit-cylinder.ts:269](../../packages/graphics/src/props/computed/representations/interactions-intra-unit-cylinder.ts#L269) | time | medium | size slider 64–410 ms per step → ~1–3 ms | S | confirmed | no |
| P17 | new | Inter-unit interactions on assemblies: one non-instanced visual, within-copy contacts repeated per operator | [interactions-inter-unit-cylinder.ts:47](../../packages/graphics/src/props/computed/representations/interactions-inter-unit-cylinder.ts#L47) | architecture | medium | 1tup asm fixtures: inter-unit buffers ÷~17 (asm-192 265 → ~20 MB, asm-864 1.2 GB → ~70 MB); real capsids with between-copy interfaces ÷1–3 | L | confirmed | no |
| P6 | still-open; **fixed** #2034 | `particle-hierarchy` (and `particle-entity`) rebuild their O(N) index per `Theme.create` | [particle-hierarchy.ts:134](../../packages/graphics/src/theme/color/particle-hierarchy.ts#L134) | time+memory | medium | ~360–390 ms + 40–80 MB per update at 10M particles → ~0 | S | confirmed | no |
| P4 | still-open | Per-element `Color.lighten` and string keys (illustrative, entity-id, chain-id multi-model, entity-source) | [illustrative.ts:80](../../packages/graphics/src/theme/color/illustrative.ts#L80) | time | low | illustrative 20 → ~5 ms per 186k elements; entity/chain-id 2–3× | S | confirmed | no |
| P5 | still-open | H-bond anchor search per directed edge, closure per call, now across inter-unit bonds | [interactions-intra-unit-cylinder.ts:115](../../packages/graphics/src/props/computed/representations/interactions-intra-unit-cylinder.ts#L115) | time | low | ~5 % of the interactions build (5–60 ms flat-32 → flat-192) | S | confirmed | no |
| P2 | still-open | Mesh fallback: dash caps 2/3 of vertices; estimate ~5× low | [link.ts:144](../../packages/graphics/src/repr/structure/visual/util/link.ts#L144) | time+memory | low | mesh fallback ~2× faster and smaller | S | confirmed | no |
| P3 | still-open | `getPalette` 'set' legend by O(count) string concatenation | [palette.ts:100](../../packages/core/src/util/color/palette.ts#L100) | time+memory | low | 1M: 46 ms / 53 MB, 10M: 0.6 s / 549 MB → O(colors) | S | confirmed | no |
| P18 | new | volume-segment theme: `Array.indexOf` per location, O(S²) per recolor | [volume-segment.ts:46](../../packages/graphics/src/theme/color/volume-segment.ts#L46) | time | low | 10k segments 19 ms, 50k 476 ms → < 5 ms | S | confirmed | no |

### P15: impostor buffers sized 6× and copied

- **Mechanism**: #2004 added `LinkBuilderProps.cylinderCount` (`withDashedCylinderCount`,
  [shared.ts:40](../../packages/graphics/src/props/computed/representations/shared.ts#L40)) so the impostor builder
  allocates once, and calls `CylindersBuilder.create(cylinderCount() * 6, linkCount, …)` ("sizes its arrays in vertices").
  Since #2007 `create` counts cylinders and multiplies both arguments by 6 itself
  ([cylinders-builder.ts:66-69](../../packages/graphics/src/geo/geometry/cylinders/cylinders-builder.ts#L66)). The six
  `ChunkedArray`s (10 floats per vertex) get a first chunk of 36 vertices per cylinder, fill 1/6, and
  `compactTrimmed` (80 % rule, `chunked-array.ts:118`) copies the used part. The growth chunk is `linkCount` cylinders
  (×6), also oversized. Updates pass the old buffers and are unaffected.
- **Evidence**: `interbench.ts` with `build2.mjs` patch `cyl6`, Float32Array bytes counted:

  | Fixture | f32 allocated base → fix | Build base → fix | Kept buffers |
  | --- | --- | --- | --- |
  | flat-32 (verifier re-run) | 435.8 → 106.2 MB | 400–522 → 187–236 ms (loaded machine) | 82.8 MB both |
  | flat-32 (auditor) | 436 → 106 MB | 155–176 → 108–132 ms | 82.8 MB |
  | asm-192 | 1,282 → 229 MB | ~250 ms, within noise | 265 MB |
  | flat-192 | 2.63 → 0.64 GB | 1.20–1.42 → 1.02–1.14 s | 500 MB |

- **Fix**: `CylindersBuilder.create(linkBuilder.cylinderCount(), Math.max(1024, linkCount / 4), cylinders)`; fix the
  comment at `link.ts:389`. An exact count fills the single chunk, so `_compact` returns it without a copy. Builder test:
  exact count → one chunk, no copy. The no-count branch `create(linkCount * 2, …)` is 2 cylinders per link and fills
  ~50 % for solid bonds (copied too): slice 07/08.
- **Risk**: low; geometry identical. An undercount only grows.

### P16: `sizeFactor` rebuilds impostor geometry

- **Mechanism**: the impostor visuals share `setUpdateState`, which sets `createGeometry` on `sizeFactor` (intra :269,
  [inter :282](../../packages/graphics/src/props/computed/representations/interactions-inter-unit-cylinder.ts#L282),
  [bridge :320](../../packages/graphics/src/props/computed/representations/interactions-bridge-cylinder.ts#L320)). For
  impostors the size is the size theme × `uSizeFactor`, which `Cylinders.updateValues` already updates
  ([cylinders.ts:376](../../packages/graphics/src/geo/geometry/cylinders/cylinders.ts#L376)). In the impostor builder
  `radius()` (the only geometry use of `sizeFactor`) runs only for multi-cylinder styles
  ([link.ts:429](../../packages/graphics/src/repr/structure/visual/util/link.ts#L429)); interactions are all `Dashed`.
  Left: the fallback bounding sphere expansion `1 * sizeFactor`.
- **Evidence**: `interbench.ts`, `sizeFactor` 0.2 ↔ 0.25 vs alpha: auditor flat-32 64–66 vs 2.2–2.4 ms, asm-192 171–173
  vs 0.2–0.3, flat-192 368–409; verifier flat-32 111–150 vs 3–9 ms (loaded), 0.1 MB allocated.
- **Fix**: drop `sizeFactor` from the impostor `setUpdateState`, keep it in `setMeshUpdateState` (the mesh bakes it).
  Expand the fallback sphere by the parameter max (`InteractionsParams.sizeFactor` max 1,
  [interactions.ts:47](../../packages/graphics/src/props/computed/representations/interactions.ts#L47); the shared
  params allow 10, so use the representation's) or the size theme max. Bond impostors
  (`bond-intra-unit-cylinder.ts:319`) need `radius` for multi-bond offsets and keep it.
- **Risk**: bounding sphere ≤ 1 Å larger; mesh unchanged.

### P17: inter-unit interactions not instanced on assemblies

- **Mechanism**: `InteractionsInterUnitVisual` is a `ComplexCylindersVisual` over every `interactions.contacts` edge in
  world coordinates (`getInterUnitInteractionsBuilderProps`, :47–200). Contacts between chains of one copy recur
  identically under every operator; each is built, stored, colored and uploaded separately. Intra-unit contacts are
  already instanced per symmetry group.
- **Evidence**: `interop.ts`. Verifier, asm-32: 41,864 edges, 40,256 (96 %) within one operator = 1,258 distinct × 32,
  1,608 between copies. Auditor, asm-192: 254,528 edges, 241,536 (95 %) = 1,258 × 192, 12,992 between. asm-192
  interactions buffers 265 MB (~336 B per edge × 3 dashes); asm-864 1.2 GB.
- **Correction**: the 1tup asm fixtures place 9-chain copies on a grid, so copies barely touch. In a real capsid with
  one or a few chains per protomer, much of the inter-unit set is between copies (op_i, op_j) and stays in the complex
  visual: there the gain is ÷1–3, not ÷17. Impact lowered to medium.
- **Fix**: split the inter-unit edges. Same-operator edges, grouped by (symmetry group A, symmetry group B) when both
  groups have aligned operator lists: geometry for the first copy in its invariant frame, one instance per operator
  (group = edge within the pattern, instance = operator). Between-copy edges stay complex. Needs a unit-pair group
  (for operator k: unit ids of A and B) and an instance → unit-pair map for `getInteractionLoci`/`eachInteraction`.
  Fallback to the current visual when operator lists differ. Same pattern applies to inter-unit bonds (slice 07).
- **Risk**: picking and marking need the new map; overpaint/transparency are group×instance and follow.

### P6: particle color index per update

- **Mechanism**: `buildHierarchyColors` (two passes, Map lookups, `Int32Array(count)` ×1–2) and
  `buildEntityColorIndex` ([particle-entity.ts:62](../../packages/graphics/src/theme/color/particle-entity.ts#L62)) run
  in every factory call; the particle transforms call `Theme.create` on every update
  ([representation.ts:126](../../packages/plugin/core/src/state/transforms/particles/representation.ts#L126)).
- **Evidence**: `verify/ph.ts`, 10M synthetic particles: hierarchy 386–390 ms (verifier), 356–371 (auditor), v1 340–357;
  entity 41–44 ms.
- **Fix**: memoize the index part per `ParticleList` (WeakMap or custom property, as `getParticleSpacingRadius`), keyed
  by `by`; build only the O(entities) color tables per call. Same for `particle-entity`.
- **Risk**: assumes immutable particle lists; key on a version otherwise.

### Low

- **P4** ([illustrative.ts:80](../../packages/graphics/src/theme/color/illustrative.ts#L80),
  [entity-id.ts:124](../../packages/graphics/src/theme/color/entity-id.ts#L124),
  [chain-id.ts:53-58](../../packages/graphics/src/theme/color/chain-id.ts#L53),
  [entity-source.ts:173](../../packages/graphics/src/theme/color/entity-source.ts#L173)): `themes.ts` flat-32
  elem/bond ms: illustrative 20.8/48.8 vs molecule-type 2.0/3.9; entity-id 6.8/14.7 vs chain-id 3.3/7.4.
  `keymb.mjs` per 1M lookups: multi-model template key + Map 21.8 ms, single-model Map 3.7, typed serial array 1.0. Fix:
  `Map<Color, Color>` memo for `lighten`; per-model `Int32Array` serial by chain index for chain-id/entity-id;
  `Map<Model, …>` for entity-source (with B2). Pin colors with a snapshot test.
- **P5** (intra :115/:126, inter :103/:115): per directed edge, a closure and a walk over the donor's bonded atoms
  (intra- and, since #2004, inter-unit bonds). Profile self: `position` inter 2.8 %, intra 0.9 %, `getEdgeIndices` 1.1 %,
  `eachIntraBondedAtom` 0.7 %. Fix: donor-H anchors once per undirected contact into a `Float32Array(6 × contacts)`,
  module-level callback.
- **P2** ([link.ts:144](../../packages/graphics/src/repr/structure/visual/util/link.ts#L144)): only the mesh fallback
  (`tryUseImpostor` off, no `fragDepth`, export). v1: `dashCap` true 5.95M vertices, 201 MB, 360–378 ms; false 2.06M,
  97 MB, 171–181 ms. Fix: `dashCap: PD.Boolean(false)` in `InteractionsSharedParams`; style-aware estimate
  `ceil(dashCount/2) * (radialSegments+1) * (2 + (dashCap ? 2 : 0))` per link.
- **P3** ([palette.ts:100](../../packages/core/src/util/color/palette.ts#L100)): reached by `model-index` (count =
  models), or a user-picked set list for particle-index/element-index, or volume-segment (count = segments).
  `verify/pal10m`: 1M 46 ms / 53.4 MB retained, 10M 617 ms / 549 MB. Fix: cap labels per row (`… (+N)`) or lazy legend.
- **P18** ([volume-segment.ts:46](../../packages/graphics/src/theme/color/volume-segment.ts#L46)): `segments.indexOf`
  per color call, one call per segment visual, array rebuilt per `Theme.create`. `v2/segidx.mjs`: 1k 0.2 ms, 10k
  19.3 ms, 50k 476 ms; Map 0.1 / 0.6 / 4.3 ms. Small next to one visual (mesh, draw call) per segment. Fix:
  `Map<SegmentIndex, number>` built once in the factory.

## Bugs

| Id | vs v1 | Severity | Title | File | Verdict |
| --- | --- | --- | --- | --- | --- |
| B2 | still-open; **fixed** #2032 | medium | entity-source indexes polymer-only `sequence.sequences` by entity index | [entity-source.ts:81](../../packages/graphics/src/theme/color/entity-source.ts#L81) | confirmed |
| B7 | new; **fixed** #2032 | medium | volume-value `isRelative`: sigma-unit domain against absolute values | [volume-value.ts:87](../../packages/graphics/src/theme/color/volume-value.ts#L87) | confirmed |
| B6 | new | low | Bridge visual ignores the child with `includeParent` | [interactions-bridge-cylinder.ts:141](../../packages/graphics/src/props/computed/representations/interactions-bridge-cylinder.ts#L141) | confirmed |
| B3 | still-open; **fixed** #2033 | low | `ColorTheme.Palette` 'set' with offsets writes `undefined` (black) | [color.ts:128](../../packages/graphics/src/theme/color.ts#L128) | confirmed |
| B5 | new | low | Impostor-capable visuals' mesh fallback recreated on every update without `fragDepth` | [interactions-intra-unit-cylinder.ts:325](../../packages/graphics/src/props/computed/representations/interactions-intra-unit-cylinder.ts#L325) | confirmed |

- **B2**: `model.sequence.sequences` holds polymer entities only; it is indexed by `entities.getEntityIndex(entityId)`
  (:80-81); the entity-index map is `byEntityKey`. `getParams` also calls `getMaps`, so the params UI and the
  carbon/illustrative 'entity-source' style fail too. Failure (verifier re-run, `entsrc.ts` on `data/ent.cif`, entity 1
  ZN, entity 2 polymer with `entity_src_gen`): `sequences = ['2']`, `byEntityKey` keys `['1']`, factory throws "Cannot
  read properties of undefined (reading 'sequence')". Also :179: `label_seq_id` 0 or past the array → `NaN` →
  `palette.color(NaN)` → black. Fix: `byEntityKey[entityIndex]`, skip rows without one, missing src or out-of-range id →
  `DefaultColor`; `Map<Model, …>` instead of `models.indexOf` per element.
- **B7**: with domain 'auto' and `isRelative`, the domain becomes `(min-mean)/sigma, (max-mean)/sigma` (:60-63), but the
  vertex path interpolates absolute values (`makeGetTrilinearlyInterpolated(grid, 'none')`, :87) and maps
  `(value - d0)/(d1 - d0)` (:97); the direct-volume path normalizes the sigma domain against absolute `min, max` (:76).
  A custom domain with `isRelative` is read as absolute. `external-volume.ts:94-111` does it right (`'relative'`
  interpolator). Failure: EM map mean 0.01, sigma 0.05, range [-0.2, 1]: domain [-4.2, 19.8]; voxel 0.5 → 0.20 of the
  palette instead of 0.58, all colors in the first fifth; direct volume gets a transfer domain normalized as
  `(-4.2 + 0.2)/1.2 = -3.3 … 16.7`. Off by default; no test covers it. Fix: `'relative'` interpolator when `isRelative`
  (any domain); in the direct path convert back (`rel * sigma + mean`) before `normalize`; test with mean ≠ 0, sigma ≠ 1.
- **B6**: with `includeParent`, `getData` passes `structure.asParent()`
  ([interactions.ts:85](../../packages/graphics/src/props/computed/representations/interactions.ts#L85)), so bridges are
  computed for the whole parent. Intra and inter visuals filter by `structure.child` and `parentDisplay`;
  `getBridgeBuilderProps` (:141-240) and `eachBridgeInteraction` never read `child` (`ignore` only drops shared legs), and
  the bounding sphere uses the parent's boundary (:251). Failure: water bridges on, interactions of a ligand component
  with `includeParent`: every water bridge of the protein is drawn. The includeParent tests (`interactions-visual.test.ts:356,
  372`) compare impostor vs mesh or intra-unit only. Fix: when `structure.child` is set, ignore a bridge by the
  `parentDisplay` rule applied to features A and B ('stub': A or B in the child; 'full'/'between': as inter-unit); child
  boundary for the sphere; add a bridge includeParent test.
- **B3**: `count = min(10^maxPrecision, 4096)` is 1 for integer offsets → `normalize(i, 0, 0)` NaN; for `t` past the
  last stop `findPredecessorIndex` returns `off.length`, `src[j]` undefined (:139-140; the interpolate branch clamps).
  Failure (verifier re-run, `palette.ts`): set `[[blue,0],[red,1]]` → `[undefined, undefined]`; `[[blue,0.2],[red,0.6]]`
  → last 4 of 10 undefined. Reached through volume-value with a user set list with offsets. Fix: `count = max(2, …)`,
  clamp `j`; extend `theme/_test/palette.test.ts`.
- **B5**: the factory picks the mesh variant when `!(tryUseImpostor && checkCylinderImpostorSupport(graphics))`, which
  includes a context without `fragDepth` (`visual/util/common.ts:55-70`). The mesh variant's
  `mustRecreate = tryUseImpostor && !!graphics` (intra :325, inter :338, bridge :374) is then true on every update;
  `UnitsRepresentation` (:151, :201, :226) and `ComplexRepresentation` (:75) destroy and rebuild it on every color, alpha
  or size change. The same pattern is in the bond visuals, `nucleotide-atomic-{element,bond}`,
  `polymer-backbone-{cylinder,sphere}` and particle spacefill. Only WebGL1 without `EXT_frag_depth` (WebGPU reports
  `fragDepth`). Failure: on such a device a recolor of flat-32 interactions rebuilds the mesh (~0.4–1 s). Fix: mesh
  `mustRecreate = tryUseImpostor && checkCylinderImpostorSupport(graphics)`, impostor `!tryUseImpostor ||
  !checkCylinderImpostorSupport(graphics)` (sphere support check for sphere visuals); all copies at once.

## Checked and dropped

- Nothing refuted. Corrections: P17 lowered high → medium (the ÷17 is the grid fixture; real capsids keep their
  between-copy interfaces); P16's "caps sizeFactor at 1" holds for `InteractionsParams` only (shared params allow 10);
  P18 measured (was by complexity); B5 widened to every impostor/mesh pair using the same `mustRecreate` pattern;
  `particle-size.ts` and the new `interaction-location.ts` size helper read clean.

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/1f4e2fbb-da6c-4948-aacb-385a4dab465c/scratchpad/perf-audit-2/graphics-theme-props/`,
bundled with `node build.mjs in.ts out.mjs` (or `node build2.mjs in.ts out.mjs [patch]`; esbuild, `molstar-src`
condition); verifier re-builds in `v/`, run `node --expose-gc --max-old-space-size=8000 v/X.mjs`.

| Script | What |
| --- | --- |
| `interbench.ts` + `build2.mjs` | interactions repr build/update/marking, `COUNT=1` Float32Array byte counter, `PROF=` cpuprofile; patch `cyl6` = P15 fix (P15, P16) |
| `interop.ts` | inter-unit edges within vs between operators on an assembly (P17) |
| `themes.ts`, `keymb.mjs` | per-element theme closures, key-string vs typed lookups (P4) |
| `entsrc.ts` (+ `data/ent.cif`), `palette.ts` | B2, B3 reproductions |
| `verify/ph.ts`, `verify/pal10m.ts` | particle-hierarchy at 10M (P6), large set legend (P3) |
| `merge.ts` | layer merge/filter costs |
| `v2/segidx.mjs` | volume-segment `indexOf` vs Map (P18) |
