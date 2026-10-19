# Slice 09: `graphics/src/theme`, `graphics/src/props`

The color and size theme registry, about 40 color themes, palettes (`ColorTheme.Palette`, `getPalette` in core) and the
overpaint, transparency, clipping, substance, emissive and wiggle layer types (`theme/`); the interactions
representation (intra-unit, inter-unit and bridge visuals), the ASA, interaction-type, SIFTS and cross-link themes, the
cross-link representation and volume streamlines (`props/`). No `_spec`. Measured on `v6-opts-261019` (`90d717f49`),
Node, Apple M4; fixtures and method in [00-profile.md](00-profile.md). A second pass re-read each finding with its
callers and re-ran the benches. The mesh primitives these visuals build are in
[08-graphics-geo.md](08-graphics-geo.md) (P3, P4, P5); this report does not repeat them.

## Hot paths

| Path                                                    | Callers                                                                     | Share / cost                                                                              |
| ------------------------------------------------------- | --------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------- |
| `createLinkCylinderMesh` with `LinkStyle.Dashed`        | all three interactions visuals (mesh only)                                  | `repr:interactions` big-bcif 1.1 s, 1.85 GB alloc, 422 MB buffers, 623 MB retained; inter 45 %, intra 42 % |
| `Theme.create` → `provider.factory`                     | structure and particle representation transforms, on every update           | structure themes ≤ 1.9 ms at 186k atoms; `particle-hierarchy` 340 ms at 10M particles (P6) |
| `getPalette` 'set' legend loop (core)                   | every palette theme's factory                                               | O(count); only large when count is atoms, particles or models (P3)                        |
| theme `color(location)` per element/bond                | color-data fill                                                             | 2–8 ms per 186k elements for most themes; illustrative 20, operator-hkl 16 (per unit), ASA 12 |

Already fine: `particle-entity` (40 ms at 10M), `unit-index` and `polymer-index` (unit-id maps), `secondary-structure`;
layer `merge`/`remap`/`areEqual` are O(layers × units + elements) after the plugin helpers merge layers per value.
Recreating structure themes on each update costs ≤ 2 ms; what a recolor costs is surface color smoothing (slice 08 P1).

## Performance findings

| Id  | Title                                                                  | Where                                                                                                                      | Kind        | Impact | Estimate                                                                   | Effort | Verdict   | Frozen |
| --- | ---------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------- | ----------- | ------ | -------------------------------------------------------------------------- | ------ | --------- | ------ |
| P1  | Interactions have no cylinder-impostor visual (dashed mesh only)       | [interactions-intra-unit-cylinder.ts:164](../../packages/graphics/src/props/computed/representations/interactions-intra-unit-cylinder.ts#L164) | time+memory | high   | big-bcif intra builder 370 → 45 ms, 201 → 39 MB; repr ~1.1 s → ~0.15 s, buffers 422 → ~80 MB | M      | confirmed | no     |
| P2  | Mesh path: dash caps are 2/3 of the vertices; vertex estimate ~5× low  | [link.ts:139](../../packages/graphics/src/repr/structure/visual/util/link.ts#L139)                                          | time+memory | medium | mesh 2.2× faster, 2× smaller (201 → 97 MB intra, big-bcif)                  | S      | confirmed | no     |
| P6  | `particle-hierarchy` rebuilds its O(N) color index on every update     | [particle-hierarchy.ts:134](../../packages/graphics/src/theme/color/particle-hierarchy.ts#L134)                             | time+memory | medium | 340 ms + 40–80 MB per update at 10M particles → ~0                          | S      | confirmed | no     |
| P3  | 'set' legend built eagerly by O(count) string concatenation            | [palette.ts:100](../../packages/core/src/util/color/palette.ts#L100)                                                        | time+memory | low    | 1M: 41 ms, 53 MB retained; 10M: 0.6 s, 549 MB → O(colors)                   | S      | confirmed | no     |
| P4  | Per-element `Color.lighten` and string keys in theme closures          | [illustrative.ts:80](../../packages/graphics/src/theme/color/illustrative.ts#L80)                                          | time        | low    | illustrative 20 → ~5 ms, entity-id 7 → ~3.5 ms per 186k elements            | S      | confirmed | no     |
| P5  | H-bond endpoint search runs for both directed copies, closure per edge | [interactions-intra-unit-cylinder.ts:93](../../packages/graphics/src/props/computed/representations/interactions-intra-unit-cylinder.ts#L93) | time        | low    | ~2× less endpoint work on H-bond edges; a few ms (not measured)            | S      | plausible | no     |

### P1: no impostor path for interactions

- **Mechanism**: `createIntraUnitInteractionsCylinderMesh` (:164), `createInterUnitInteractionCylinderMesh`
  ([inter:194](../../packages/graphics/src/props/computed/representations/interactions-inter-unit-cylinder.ts#L194)) and
  `createBridgeCylinderMesh` ([bridge:233](../../packages/graphics/src/props/computed/representations/interactions-bridge-cylinder.ts#L233))
  all call `createLinkCylinderMesh` with `LinkStyle.Dashed`: per directed edge `ceil(6/2) = 3` capped mesh cylinders
  through `addPrimitive`. Bonds already choose `UnitsCylindersVisual` + `createLinkCylinderImpostors`, which handles
  `Dashed` through `CylindersBuilder.addFixedCountDashes`.
- **Evidence**: profile `repr:interactions` big-bcif 1101 ms / 1847 MB alloc / 643 MB retained; buffers `aPosition` and
  `aNormal` 143 MB each, index 88, group 48. `inter.ts` re-run, intra-unit only (38,336 directed edges, 416 units,
  `radialSegments` 8 = the repr's 'low' quality): mesh 360–378 ms, 5.95M vertices (155 per edge), 201 MB; impostors
  43–45 ms, 39 MB. Inter-unit adds 41,864 directed edges at the same ratio. Small structures get 'high' quality
  (`radialSegments` 20), so the mesh is heavier per edge there.
- **Fix**: mirror `IntraUnitBondCylinderVisual`: impostor variants of the three visuals (`UnitsCylindersVisual`,
  `ComplexCylindersVisual`), `tryUseImpostor: PD.Boolean(true)`, choose with
  `props.tryUseImpostor && checkCylinderImpostorSupport(graphics)`, `mustRecreate` on `tryUseImpostor`/graphics; keep
  the mesh visual as the fallback. The existing `builderProps` (position, style, ignore) feed both builders.
- **Risk**: the impostor builder ignores `radius()`: the size comes from the size theme at the group's location times
  the geometry `sizeFactor`. The mesh path bakes `min(sizeA, sizeB) * sizeFactor` from element locations, so with the
  default `uniform` size theme the two agree; with `physical` they do not unless the size iterator maps a group to that
  minimum. Check dash spacing and caps against `addFixedCountDashedCylinder` (both step `d / (count + 0.5)`), the
  `includeParent` modes, picking/marking (group = edge index, unchanged), and the no-`fragDepth` fallback.

### P2: mesh fallback sizing and dash caps

- **Mechanism**: interactions inherit `dashCap: true` ([link.ts:38](../../packages/graphics/src/repr/structure/visual/util/link.ts#L38));
  every dash gets two caps although its radius is 0.08 Å (uniform 1 × `sizeFactor` 0.2 × `dashScale` 0.4).
  `vertexCountEstimate = radialSegments * 2 * linkCount * 2` assumes one uncapped cylinder per directed link: 32 against
  155 real vertices at 8 segments, so `ChunkedArray` grows by ~15 chunks of `estimate / 4` and `compact` copies
  everything once. (Slice 08 P3 found the same estimate 2× too high for solid links; the fix must be style-aware.)
- **Evidence**: `inter.ts` intra-unit, 8 segments: `dashCap` true 5.95M vertices, 201 MB, 360–378 ms; false 2.06M,
  97 MB, 171–181 ms. Slice 08 P5 (dedupe cylinder primitive vertices) shrinks caps from 2 × 33 to 2 × 17 vertices per
  dash and reduces this gain to ~1.6×.
- **Fix**: `dashCap: PD.Boolean(false)` in `InteractionsSharedParams` (or caps only at `radialSegments >= 16`); estimate
  per directed link from the style: dashed `ceil(dashCount / 2) * (radialSegments + 1) * (2 + (dashCap ? 2 : 0))`.
- **Risk**: open dash ends are visible end-on at high zoom (no caps, single-sided mesh). Matters only for the mesh
  fallback once P1 lands (no `fragDepth`, `tryUseImpostor: false`, export).

### P6: particle color index rebuilt per update

- **Mechanism**: `buildHierarchyColors` makes two passes over all particles with 1–3 `Map` lookups each and allocates
  `Int32Array(count)` (plus a second one for `by: 'function'`) on every factory call. Both particle transforms call
  `Theme.create` on each update ([representation.ts:126](../../packages/plugin/core/src/state/transforms/particles/representation.ts#L126)),
  so alpha, size or quality changes redo it. `particle-hierarchy` is the default theme for lists with entities and
  compartments ([provider.ts:62](../../packages/plugin/core/src/state/formats/particles/provider.ts#L62)).
- **Evidence**: `verify/ph.ts`, synthetic 10M particles, 20 compartments, 500 entities: `particle-hierarchy` factory
  338–357 ms; `particle-entity` 39–44 ms.
- **Fix**: memoize the index part in a `WeakMap<ParticleList, Map<by, {colorIndex, compartment/entity maps}>>`; build
  only the `colors` table (O(entities)) from palette/variability/shift per call.
- **Risk**: assumes particle lists are immutable; if one can change in place, key on a version.

### Low

- **P3**: [palette.ts:100](../../packages/core/src/util/color/palette.ts#L100) (core, called from every palette theme)
  appends `, label` per index to the row of `i % colors`. Reached with large counts only when a user picks a 'set' list
  for `particle-index`/`element-index` (their defaults switch to an interpolated scale above the list length) or by
  `model-index` (default set, count = models). `verify/pal10m.ts`: 1M 41 ms, 7.9M chars, 53 MB retained (rope
  strings); 10M 0.6 s, 549 MB. Fix: cap labels per row (first ~10 + `… (+N)`) or build the legend lazily.
- **P4**: `themes.ts` re-run, element/bond ms at 186k atoms: illustrative 20.4/48.5 vs molecule-type 2.0/4.2;
  entity-id 7.1/15.2 vs chain-id 3.5/7.9; operator-hkl 16.4/34.3 (per unit) vs operator-name 2.2/4.8. Fix: memoize
  `lighten` in a `Map<Color, Color>`; per-chain serial `Int32Array` for entity-id/chain-id; `Map<Model, index>` for
  entity-source; `Map<unit.id, serial>` for operator-hkl.
- **P5**: `position()` runs per directed edge, so the bonded-H search runs twice per H-bond with a new closure each time.
  Fix: cache donor-H/acceptor points per undirected contact; module-level callback. Not measured; small next to P1.

## Bugs

| Id  | Title                                                                    | Where                                                                                                                      | Severity | Verdict   |
| --- | ------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------- | -------- | --------- |
| B1  | intra-unit `includeParent` filter reads `members[offsets[i]]`            | [interactions-intra-unit-cylinder.ts:135](../../packages/graphics/src/props/computed/representations/interactions-intra-unit-cylinder.ts#L135) | medium   | confirmed |
| B2  | entity-source indexes polymer-only `sequence.sequences` by entity index  | [entity-source.ts:81](../../packages/graphics/src/theme/color/entity-source.ts#L81)                                        | medium   | confirmed |
| B3  | `ColorTheme.Palette` 'set' with offsets writes `undefined` (black)       | [color.ts:128](../../packages/graphics/src/theme/color.ts#L128)                                                            | low      | confirmed |
| B4  | bridge visual rebuilds geometry on every update (size theme by reference) | [interactions-bridge-cylinder.ts:276](../../packages/graphics/src/props/computed/representations/interactions-bridge-cylinder.ts#L276) | low      | confirmed |

- **B1**: the loops `for (i = offsets[f]; i < offsets[f + 1]; ++i)` already walk member positions, yet read
  `elements[members[offsets[i]]]` at :135, :144, :150. For `i >` feature count `offsets[i]` is `undefined`, so the
  element is `undefined` and `SortedArray.has` is false. The inter-unit visual (:153/:167/:179) and the intra
  `eachInteraction` marking loop use `members[i]`. Failure: `includeParent` on, child = part of a unit → stub/full/between
  decisions use unrelated atoms; `parentbug.ts` on 1tup with an even-residue child: 67 of 1198 directed edges decided
  differently (3 of 438 with a contiguous child, at its boundary). Fix: `members[i]`; test with a child structure.
- **B2**: `sequences` holds polymer entities only (`fromAtomicHierarchy` skips non-polymers; `basic/sequence.ts` follows
  `entity_poly_seq`); the entity-index-keyed map is `byEntityKey`. Failure: `data/ent.cif` (entity 1 ZN non-polymer,
  entity 2 polymer with an `entity_src_gen` row): `sequences = ['2']`, factory throws `Cannot read properties of
  undefined (reading 'sequence')`; `getParams` also calls `getMaps`, so the theme UI and the carbon/illustrative
  'entity-source' style fail too. Any file with a non-polymer before a polymer, or a src row for a branched entity,
  throws; standard PDB order (polymers first) works. Also :179: `label_seq_id` 0 or past the array gives `NaN` →
  `palette.color(NaN)` → `undefined` (black). Fix: `byEntityKey[entityIndex]`, skip rows without one; treat a missing src
  as `DefaultColor`; replace `models.indexOf` per element (:173) with a `Map`.
- **B3**: `count = min(10^maxPrecision, 4096)` is 1 for integer offsets, so `normalize(i, 0, 0)` is `NaN`/`Infinity`;
  and for `t` past the last stop `findPredecessorIndex` returns `off.length` and `src[j]` is `undefined` (the interpolate
  branch clamps, :111). `palette.ts` re-run: set `[[blue, 0], [red, 1]]` → `[undefined, undefined]`; set
  `[[blue, 0.2], [red, 0.6]]` → last 4 of 10 `undefined`. Reached through `volume-value` (:77, :100) when a user picks a
  set list with offsets (the default is interpolate). Fix: `count = max(2, …)`, clamp `j` to `[0, src.length - 1]`;
  extend `theme/_test/palette.test.ts`.
- **B4**: `setUpdateState` adds `newTheme.size !== currentTheme.size`; `Theme.create` on every structure-representation
  update ([representation.ts:151](../../packages/plugin/core/src/state/transforms/structure/representation.ts#L151))
  makes a new size theme, so a color or alpha change rebuilds the bridge mesh. `ComplexMeshVisual` already rebuilds on
  `!SizeTheme.areEqual` ([complex-visual.ts:461](../../packages/graphics/src/repr/structure/complex-visual.ts#L461)).
  Fix: drop the clause. Cost is small (bridges are few).

## Checked and dropped

- Nothing refuted. Corrections: P6 measured (340 ms at 10M, `particle-hierarchy` only; `particle-entity` is 40 ms) and
  raised to medium; P3 retains 53 MB per 1M labels as rope strings, not ~1.4 MB per 100k; P5 left plausible (unmeasured,
  small); P2 gain shrinks to ~1.6× after slice 08 P5 and applies only to the fallback after P1.
- The ASA theme on coarse units is a known issue, not repeated.

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/ac6cef1a-5954-43a5-ae83-c22fef35f065/scratchpad/perf-audit/graphics-theme-props/`,
bundled with `node ../build.mjs X.ts X.mjs` (esbuild, `molstar-src` condition), run with `node X.mjs`.

| Script            | What                                                                                   |
| ----------------- | -------------------------------------------------------------------------------------- |
| `themes.ts`       | factory and per-element/bond eval of every color theme on 1tup-2x2x2 (186k atoms)      |
| `inter.ts`        | intra-unit interactions: dashed mesh vs impostors on the same builder props, `RS=8 node inter.mjs` |
| `parentbug.ts`    | B1: ignore decisions with `members[offsets[i]]` vs `members[i]` on 1tup + child        |
| `entsrc.ts`       | B2: entity-source factory on `data/ent.cif` (non-polymer entity first)                 |
| `palette.ts`      | B3: 'set' palettes with integer offsets and a last stop below 1                        |
| `palbench.ts`     | P3: `getPalette` 'set' legend time and size, 1e3–1e6                                   |
| `verify/pal10m.ts`| P3: legend time and retained heap at 1M and 10M (`node --expose-gc`)                   |
| `verify/ph.ts`    | P6: `particle-hierarchy` / `particle-entity` factory on 10M synthetic particles        |
