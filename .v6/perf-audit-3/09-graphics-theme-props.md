# Slice 09: graphics theme and props (v3)

Scope: the color and size themes, palettes and layers in [`theme/`](../../packages/graphics/src/theme); the interactions
representation (intra-unit, inter-unit, bridge visuals), the ASA, interaction-type, SIFTS and cross-link themes and
volume streamlines in [`props/`](../../packages/graphics/src/props). HEAD `16277d822` (2026-10-25), Node 26.12, Apple
M4. The machine was **very loaded** during this run (load average 22–34, not 4–5): wall times are 3–10× v2's and only
in-process or back-to-back A/B ratios, counts and allocation are used below. A/B base: `git archive 0944823d4` of
`packages/{core,io,model,graphics}`, bundled with the same bench (`build-base.mjs`, an `@molstar/*` → base `src`
resolver). v2: [../perf-audit-2/09-graphics-theme-props.md](../perf-audit-2/09-graphics-theme-props.md). No finding
touches Canvas3D or `gl`.

PRs since v2 that touched the slice: #2032 (P16, B2, B7), #2033 (B3), #2034 (P6), #2042 (cached chain-id, entity-id and
operator-hkl serial maps), #2044 (P18), #2045 (B5, B6), #2048 (P4 part: illustrative `lighten` memo); #2025 (P15) in
`repr/.../link.ts`, read for this slice. #2043 (Model Archive PAE) has no file in these paths.

## Hot paths now

| Path | Cost (fixture) | Where the time goes |
| --- | --- | --- |
| `particle-hierarchy` per `Theme.create`, > 16 compartments | 17 × 100 entities: 127–156 ms every update; 20 × 500: 1.4–3.4 s every update (10M or 100k particles alike) | `distinctColors` rejection sampling (`getSamples`, color-space conversions ~95 % of the profile); the 16-entry LRU of #2046 thrashes on one key per compartment (P19) |
| `particle-compartment` per `Theme.create` | 1M: 20–90 ms; 10M: 0.3–0.65 s + a 40 MB `Int32Array` per update | O(N) dense-index pass, not cached as #2034 did for its siblings (P20) |
| inter-unit interactions on assemblies | asm-32: 41,864 edges, 96 % = 1,258 within-copy edges × 32 operators | unchanged since v2 (P17) |
| `repr:interactions` fresh build (impostors, flat-32) | f32 allocated 106 MB for 82.8 MB kept (v2 before #2025: 436 MB) | P15 fixed; remaining: builder, dashes (slice 08) |
| interactions `sizeFactor` update (flat-32) | 8.5–13 ms, alpha update 4–8 ms on the same load (v2: 64–66 vs 2.2–2.4) | P16 fixed: no geometry rebuild |
| theme `color(location)` per 186k elements (flat-32) | entity-id 23/63–94 ms elem/bond, illustrative 25/67–72, chain-id 9–10/22–36, molecule-type 5–7/9–27 | per-element template key + `Map` (P4) |
| `getPalette` 'set' legend | 1M: 251 ms / 53 MB retained; 10M: 5.8 s / 549 MB (loaded) | O(count) string concatenation in core (P3) |
| pipeline `theme-eval:*`, `theme-apply:*` (7qpd) | alloc: theme-eval 0–3 MB, theme-apply cartoon/ball-and-stick/spacefill 0–1 MB, gaussian-surface 4–14 MB (v2 7.8) | times not comparable (load 27); allocation as v2 |

## v2 → now

14 entries (v2's open, partly fixed and new items of this slice): **fixed 9, partly fixed 1, still open 4, regressed 0,
dropped 0**. No new bug or cost introduced by a fix since v2 was found (two small notes under the table).

| v2 id | Title | Status | PR | Note |
| --- | --- | --- | --- | --- |
| P15 | Impostor interactions: cylinder count ×6 twice, first chunk 6×, copied | fixed | #2025 | `link.ts:392` `create(cylinderCount(), max(1024, linkCount/4))`; `interbench` flat-32 f32 allocated 106.2 MB (v2 base 435.8, v2 patched 106.2), kept 82.8 MB |
| P16 | Impostor interactions rebuild geometry on `sizeFactor` | fixed | #2032 | `sizeFactor` moved to `setMeshUpdateState` in all three visuals; flat-32 update 8.5–13 ms vs alpha 4–8 ms, 0 MB f32 |
| P17 | Inter-unit interactions not instanced on assemblies | still open | | code unchanged; `interop` asm-32: 41,864 edges, 40,256 within copy (1,258 distinct × 32), 1,608 between copies (as v2) |
| P6 | `particle-hierarchy`/`particle-entity` rebuild the O(N) index per `Theme.create` | fixed | #2034 | `WeakMap<ParticleList, Map<by, …>>` keyed by `entityInfo.size`; 10M (8 compartments × 50): 2nd/3rd call 0 ms (base 1.1–1.4 s each), entity 0–9 ms (base 0.74–1.3 s). Left: P19 for > 16 compartments; sibling `particle-compartment` not covered (P20) |
| P4 | Per-element `lighten` and string keys | partly fixed | #2048, #2042 | A/B flat-32, base → head: illustrative 50–60/130–152 → 25/67–72 ms (elem/bond); serial maps cached (entity-id create 0.09–0.12 → 0.03–0.07 ms, operator-hkl 0.39–0.46 → 0.04). Open: per-element keys in entity-id ([entity-id.ts:144](../../packages/graphics/src/theme/color/entity-id.ts#L144)), multi-model chain-id ([chain-id.ts:55](../../packages/graphics/src/theme/color/chain-id.ts#L55)); illustrative's default style is entity-id, so it is now entity-id + ~10 % |
| P5 | H-bond anchor search per directed edge, closure per call | still open | | `intra:115/:126`, `inter:103/:115` unchanged |
| P2 | Mesh fallback: dash caps, estimate ~5× low | still open | | `link.ts:38` `dashCap: true`, `link.ts:144` estimate unchanged |
| P3 | `getPalette` 'set' legend by O(count) concatenation | still open | | core `palette.ts:100` unchanged; `pal10m` 1M 251 ms / 53.4 MB retained, 10M 5.8 s / 549 MB (loaded machine; v2 46 ms, 617 ms) |
| P18 | volume-segment `indexOf` per location | fixed | #2044 | `Map` of places; `seg.ts` 50k segments: create 4–28 ms, one color per segment 3–6 ms (v2 476 ms); 25 distinct colors in the first 5k (palette intact) |
| B2 | entity-source indexes `sequences` by entity index | fixed | #2032 | `byEntityKey[getEntityIndex]`, rows without a sequence skipped, out-of-range seq id → default; `entity-source.test.ts` |
| B7 | volume-value `isRelative` against absolute values | fixed | #2032 | `'relative'` interpolator; direct path converts the sigma domain back (`d * sigma + mean`); `volume-value.test.ts` |
| B3 | `Palette` 'set' with offsets writes `undefined` | fixed | #2033 | `count ≥ 2`, `j` clamped; `palette.test.ts`. A related older defect remains: B8 |
| B5 | Mesh fallback recreated on every update without `fragDepth` | fixed | #2045 | `mustRecreate` asks `checkCylinderImpostorSupport` (intra :323/:340, inter :336/:353, bridge :427/:444) |
| B6 | Bridge visual ignores the child with `includeParent` | fixed | #2045 | `getChildLegs` by `parentDisplay` (bridge :129), same feature-in-child rule as inter-unit; child change rebuilds (`complex-visual.ts:217`) |

Notes on the fixes (not bugs):

- #2032 (P16): with no geometry rebuild on `sizeFactor`, the fallback bounding sphere (`intra:198`, `inter:213`,
  `bridge:305`, expanded by `1 * sizeFactor` at build time) keeps the old factor: at most `(new − old) × 1 Å` too small
  (≤ 0.8 Å with `InteractionsParams.sizeFactor` max 1). The sphere is from endpoints anyway
  (`cylinders.ts:150`); negligible.
- #2034 (P6): the cached index is keyed by `entityInfo.size`; with `by: 'function'` the grouping reads the entries'
  `function`, so a same-size edit of an entry would be stale. Nothing edits entries in place today (Simularium only adds).
- #2042: the serial maps are cached in `structure.currentPropertyData` (chain-id and operator-hkl on the root,
  entity-id on the structure itself, its source serials come from its models). `asParent()` proxies share the parent's
  `state` and so its property data, with the same units: correct. Inputs (`Model.Index`, `AsymIdOffset`) are set once
  when the model cell is created, before structures.

## Open performance findings (ranked)

| Id | Status | Title | File | Kind | Estimate | Effort | Frozen |
| --- | --- | --- | --- | --- | --- | --- | --- |
| P19 | new | `particle-hierarchy`: per-compartment distinct colors recomputed every update past 16 compartments (LRU thrash) | [particle-hierarchy.ts:239](../../packages/graphics/src/theme/color/particle-hierarchy.ts#L239) | time | measured: 17 × 100 entities 127–156 ms → 0 per update; 20 × 500 1.4–3.4 s → 0 | S | no |
| P17 | still open | Inter-unit interactions on assemblies: one non-instanced visual | [interactions-inter-unit-cylinder.ts:47](../../packages/graphics/src/props/computed/representations/interactions-inter-unit-cylinder.ts#L47) | memory+time | v2 estimate kept: 1tup asm buffers ÷~17; real capsids ÷1–3 | L | no |
| P20 | new | `particle-compartment` rebuilds its O(N) index per `Theme.create` (missed by #2034) | [particle-compartment.ts:75](../../packages/graphics/src/theme/color/particle-compartment.ts#L75) | time+memory | measured: 1M 20–90 ms, 10M 0.3–0.65 s + 40 MB per update → ~0 | S | no |
| P4 | partly fixed | Per-element template keys + `Map` in entity-id (and illustrative's default style), multi-model chain-id | [entity-id.ts:144](../../packages/graphics/src/theme/color/entity-id.ts#L144) | time | measured: entity-id 23/63–94 ms vs molecule-type 5–7/9–27 (flat-32 elem/bond); `keymb` multi-model key 217 vs typed 3.3 ms per 1M (loaded) | S | no |
| P3 | still open | `getPalette` 'set' legend O(count) concatenation | [palette.ts:100](../../packages/core/src/util/color/palette.ts#L100) | time+memory | measured: 1M 53 MB retained, 10M 549 MB → O(colors) | S | no |
| P5 | still open | H-bond anchor search per directed edge | [interactions-intra-unit-cylinder.ts:115](../../packages/graphics/src/props/computed/representations/interactions-intra-unit-cylinder.ts#L115) | time | v2: ~5 % of the interactions build | S | no |
| P2 | still open | Mesh fallback: dash caps, estimate ~5× low | [link.ts:144](../../packages/graphics/src/repr/structure/visual/util/link.ts#L144) | time+memory | v2: mesh fallback ~2× faster and smaller | S | no |

### P19: distinct group colors past the LRU's 16 entries

- **Mechanism**: `buildHierarchyColors` calls `getDistinctGroupColors(entityCount, baseColors[ci], variability, shift)`
  once per compartment on every `Theme.create` (the particle transforms create the theme on every update,
  `representation.ts:126` in plugin-core). Each call is a `distinctColors` with a hue range around its compartment's base
  color, so one cache key per compartment. #2046 (01/P9) put an `LRUCache` of **16** entries in front of `distinctColors`
  ([distinct.ts:119](../../packages/core/src/util/color/distinct.ts#L119)); cycling through 17 or more keys in order
  evicts each one before its reuse, so every call recomputes. `computeDistinctColors` draws
  `max(1000, count × 100)` samples by rejection inside a 20° hue band: ~7.5 ms per compartment of 100 entities,
  70–150 ms per compartment of 500.
- **Evidence**: `ph2.ts` (synthetic particle list, `NC` compartments × `NE` entities), three factory calls:
  16 × 100: 248, 0, 0 ms; **17 × 100: 314, 156, 127 ms**; 8 × 50 at 10M: 3.4 s (index), 0, 0; 20 × 500 at 100k:
  1.4–3.0 s on every call (base the same). Profile of 20 × 500: `getSamples` 20 %, `xyz_rgb`/`xyz_lab`/`toHcl`/
  `rgb_xyz`/`toLab` 52 %, `checkColor` 5 %, the index 0 %.
- **Fix**: memoize the color part of `buildHierarchyColors` next to the cached index: in `getHierarchyIndex`'s entry,
  keep the last `{ key, colors, baseColors, baseLabels }` with `key = JSON.stringify([compartmentPalette, variability,
  shift])`. Alternatively raise the LRU to ~256 (entries are small `Color[]`), but the theme-level memo also skips the
  `JSON.stringify` and `slice` per compartment. Test: 17 compartments, two creates, `distinctColors` spy called 17 times
  in total.
- **Risk**: none; same colors.

### P20: `particle-compartment` index per update

- **Mechanism**: `buildCompartmentColorIndex` (two passes over `count`, a `Map` lookup per particle, a fresh
  `Int32Array(count)`) runs in every factory call; it is the default theme for compartment-only particle lists
  (`provider.ts:66` in plugin-core) and the transform recreates the theme on every update. #2034 cached the same
  pattern in `particle-hierarchy` and `particle-entity` only.
- **Evidence**: `pc.ts`, 20 compartments: 1M 20–90 ms per call, 10M 295–650 ms per call (load 22), each also
  allocating a 40 MB `Int32Array` (off-heap, not in `total_allocated_bytes`).
- **Fix**: the `WeakMap<ParticleList, …>` cache of `particle-entity.ts:68` (keyed by nothing else: the index depends on
  `compartments` only). Add a case to `particle-colors.test.ts`.
- **Risk**: none for immutable lists (the same assumption as #2034).

### P17, P4, P3, P5, P2

Unchanged from v2 except as noted in the table; the v2 sections hold (mechanism, fix). P4's remaining part: the
entity-id key `${entityId}|${sourceSerial}` per element (and per bond end) — replace by a per-model `Int32Array`
chain-index → serial built with the cached maps (`getCachedSerialMaps` is the natural home now); the same for chain-id
when `models.length > 1`. Illustrative then follows its style theme.

## Bugs

| Id | Severity | Status | Title | File | Failure scenario |
| --- | --- | --- | --- | --- | --- |
| B8 | low | new (pre-existing, not from #2033) | `Palette` 'set': samples at `i / (count − 1)` and `il = max(count, list.length)`; a stop is lost when stops ≥ samples | [color.ts:139](../../packages/graphics/src/theme/color.ts#L139) | volume-value, kind 'set', 11 stops at 0.0, 0.1, … 1.0 → 11 texels `0,1,2,3,4,6,7,8,9,10,10`: stop 5 never drawn, value 0.5 shows stop 6 |

- **B8**: one decimal gives `count = 10`, but 11 stops sit on a 0.1 grid; the loop then runs `max(10, 11) = 11` times
  with `t = i / 9`, so sample 10 is past 1 (clamped to the last stop) and samples 5–9 are one stop ahead. The palette
  texture is read with a nearest filter at `floor(v × width)`, so texel k covers `[k/11, (k+1)/11)` while it was
  computed at `k/9`. Probe `palette2.ts`: stops shown over `v ∈ [0, 1]`: `0,1,2,3,4,6,7,8,9,10`. With fewer stops than
  samples the misalignment is within one sample (e.g. 5 stops on a 0.25 grid, 100 samples: boundaries off by ≤ 1 %).
  Reached through volume-value's color list (`offsets: true`) with kind 'set'. Fix: sample texel centers,
  `t = (i + 0.5) / n` with `n = max(count, 2 × list.length)` (or `10^p + 1` and `t = i / 10^p`), and drop the
  `max(…, list.length)` loop bound; extend `theme/_test/palette.test.ts` with the 11-stop case.

## Already good, re-checked

- The #2042 serial-map caches (correct sharing, see notes), operator-hkl per-instance key (granularity 'instance': one
  call per unit, 416 on flat-32).
- `Theme.create` element-symbol + physical on flat-32: base 0.43–0.87, head 0.24–0.41 ms per call (A/B, 200 calls);
  chain-id + uniform 0.004–0.009 ms.
- Theme tests `packages/graphics/src/theme/_test`: 8 files, 27 tests pass at HEAD.
- `sifts-mapping`: per-structure accession scan, per-key color memo; module-level accession map grows with every
  accession seen (strings only, by design: "same colors for same accessions").
- `particle-attribute`, `particle-index`, `cross-link`: O(1) factories.

## Not measured

- Absolute stage times (`theme-eval`, `theme-apply`, `repr:interactions`): the machine ran at load 22–34 throughout
  (7qpd `repr:gaussian-surface` 1.1–8 s vs v2 375 ms). Allocation matched v2; ratios above are same-process or
  back-to-back. Re-run `node scripts/pipeline-bench.mjs run --fixture 7qpd --runs 5 --only
  'valence-model,interactions,theme-*,repr:*'` on a quiet machine.
- P17's buffer sizes on asm-192 (not re-run: code unchanged, edge counts on asm-32 identical to v2).

## Bench scripts

In `/private/tmp/claude-501/-Users-davinci001-projects-molstar-molstar-v6/31ee827f-4715-4967-bec3-085f5013d11f/scratchpad/pa3/09/`:
`node build.mjs <abs in.ts> <abs out.mjs>` (HEAD, `molstar-src`), `node build-base.mjs …` (v2 base in `base/`),
`build2.mjs` (with patches).

| Script | What |
| --- | --- |
| `ph2.ts` (`node ph2.mjs N NC NE`) | particle-hierarchy/entity factory per call (P6, P19) |
| `pc.ts` | particle-compartment factory per call (P20) |
| `seg.ts` | volume-segment create + eval at 1k/10k/50k (P18) |
| `palette2.ts` | B8 reproduction |
| `interbench.ts` (`COUNT=1`) | interactions build, f32 allocation, sizeFactor vs alpha (P15, P16) |
| `interop.ts` | inter-unit edges within vs between copies (P17) |
| `themes.ts`, `tc.ts`, `keymb.mjs` | per-element theme cost, `Theme.create` A/B, key lookups (P4) |
| `verify/pal10m.ts` | 'set' legend at 1M/10M (P3) |
| `pb2-7qpd.log`, `pb-flat32.log` | pipeline theme/repr stages (loaded machine) |
