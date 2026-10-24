# Verify C: slices 07–12 (graphics repr/geo/theme, Prism viewer/UI, app level)

Adversarial pass over the batch-C findings. HEAD `16277d822`. Load average during the reruns was 4–9, against the
auditors' 20–34, so wall times below are lower than the reports' and are given only to correct the scale. Probes are
the auditors' own scripts in `scratchpad/pa3/<slice>/`, re-run unchanged. No source was changed.

| Id | Verdict | Severity / estimate | Reason |
| --- | --- | --- | --- |
| 07/P17 | CONFIRMED | low; adds ~35 ms to a flat-192 recolor that takes 1.5 s, under load | `updateValues` calls plain `ValueCell.update` on `uLod` and `uClipObjectType`, and the fingerprint watches both. Rerun `fp.mjs` 7qpd: 7/7 members re-checked on each update. The fix is not in frozen code: it goes in `repr/merged.ts` and `geo/geometry/base.ts`, while `getMergeKey`/`canMergeValues` (`gl/renderable/merged.ts`) stay as they are |
| 07/P20 = 09/P17 | CORRECTED | one finding (v2 09/P17), not new; medium; asm-192 265 MB → ~11–20 MB, not ~2–3 MB | The inter-unit visuals are `ComplexVisual`s with one instance. The between-copy edges stay non-instanced: about 4 % of all edges (1,608 of 41,864 on asm-32), roughly 10 MB on asm-192. v2's ÷~17 (~20 MB) holds |
| 08/P15 + P12 | CONFIRMED | as reported (alloc +30–35 MB per cold cartoon build; MC 2 → 7–13 MB) | `--trace-deopt` on small cases. `tube.mjs 1000`: `addTube` "embedded weak objects cleared", and the run after it re-allocates 33.8 MB. `mc.mjs 60`: `interpolate` ×5, `processCell` ×4, `stridedRow` ×4, `get3dOffsetFromEdgeInfo` ×4, plus `add`/`add3`. Not on every major GC (1 drop in 6 gc'd tube runs), but it recurs |
| 08/P18 | CORRECTED | low; 5.5 ms/glyph at quality 3 (load 5), not 7–10; 40–60 glyphs ≈ 0.22–0.33 s | Mechanism confirmed: 2 full EDTs over a ~218×405 cell. `glyphs-head.mjs`: q1 1.3, q2 3.0–3.3, q3 5.5–5.7, q4 9.0–9.2 ms per glyph |
| 08/P19 | PLAUSIBLE | low; gain small | Confirmed 12 `caAdd*` calls per line. Rerun `cylh.mjs`: 300k lines take 7–19 ms, as fast per float as the cylinder fast path (15–20 ms). The claimed "74–350 ms" is load inflation. A P16-style fix is worth perhaps ~1.5×, a few ms per 300k lines |
| 08/P20 | CONFIRMED | low; +48 MB per 50k labels (~1 KB/label), add +15–25 % time | `TextMeasure.layout` does `split('')`, then a `Float64Array` and `previous + grapheme` per pair. Rerun: base 9.3–9.8 MB vs HEAD 57.6–58.0 MB; add 126–148 vs 152–163 ms |
| 08/B10 | CONFIRMED | low | `get` draws first. On overflow it returns `placeholder` without recording `mapped[char]`, so every later `get` redraws (2 EDTs) and warns again |
| 09/P19 | CORRECTED | low; 17×100: 33–45 ms per update (load 5); not caused by #2046 | Thrash confirmed (16 compartments: 0 ms on repeat; 17: 33–45 ms every call). Before #2046 there was no cache at all, so this is a gap in #2046's LRU, not a regression from it |
| 09/P20 | CORRECTED | low; 10M: ~46 ms + 40 MB `Int32Array` per update; 1M: 4–5 ms | The O(N) index is uncached on every `Theme.create` (default theme for compartment-only lists, `provider.ts:66`). The reported 0.3–0.65 s is ~7–14× load inflation |
| 09/B8 | CONFIRMED | low | `palette2.mjs`: 11 stops on a 0.1 grid give texels `0,1,2,3,4,6,7,8,9,10,10`, so stop 5 is lost. Code: `count = 10^1`, `t = i/(count−1)` over `max(count, len)` samples. Reachable only through user colour lists of kind 'set' with offsets |
| 10/P17 | CONFIRMED | low; N=200 2–17 ms per op under load, ~1–3 ms idle (est.) | `refreshSequence` loops over all of `structureRows` and calls `row()` + `sequenceNode()` on each bump. The bump has no scope. N=200 is rare |
| 10/P18 | CONFIRMED | low–medium (long agent sessions); O(E) per append, O(E²) per session | `draftOf` copies 4 Maps + `order`, and `writesOf` scans every key. Rerun `bench-record.mjs`: 285 KB/101 µs at 100 entries → 2.08 MB/831 µs at 4000 (load 5). The fix belongs in `packages/ai/src/threads` |
| 10/P19 | CONFIRMED | low | `groupEye` awaits one `component.visibility`/`measure.visibility` run per row inside `api.group`. Components K is small |
| 10/B7 | CONFIRMED (code) | low | `hasEye` is a getter. The row's watch sources are label/description/isHidden/picked/place. Emptying the only component leaves the group's label, "1 component" and hidden=false (`allHidden([])` is false) unchanged, so the row never re-renders. Not run in a DOM |
| 11/P12 | CONFIRMED | low–medium (only while the assistant runs a long tool call) | `canvasBottom` has the non-passive `wheel` listener (`shell-view.tsx:391`). `StatusLine` renders `props.after` (the turn card) inside it, and that card holds a `ScrollArea` with max-height 420px. The comment "nothing in it scrolls" is wrong |
| 12/P10 | CORRECTED | medium; deferrable ≈ 0.4–0.57 MB, not 877 kB | Prism's default `graphicsBackend` is `'auto'`, so `getAllBackends()` and the WebGPU code are used at start on WebGPU browsers: the 308 kB is needed at start in the common case. Of the extensions, the lazy ones (video, zenodo) are already stubs. Model-archive and assembly-symmetry register props/behaviours before load, so only their bodies could be lazy. That leaves MVS 234 + transpilers 87 + AI 27 + part of the extensions' 222 |
| 12/P12 | CORRECTED | low; `typecheck-full.test.mjs` 2.86 s for 2 tests at load 5 (not 13 s) | The file is untagged and spawns `tsc`, so it belongs in medium or heavy. The `structure-focus-representation` tests sleep 100 ms + 25 ms polls with real timers (`:67`, `:69`, `:131`): slow by waste, as reported |
| 12/P13 | CONFIRMED | low | Each created root model runs 5 loops over the roots: `handleModel`'s `maxModelIndex` and `asymIdOffset`, then `setMaxIndices`' max, model set and structure set. That is O(N) per model, O(N²) per load. The 36 ms at N=200 was measured under load |

## Evidence notes

- **07/P17**: `isMergeRelevantKey` keeps `/^d[A-Z]/`, `*Type`, `uInstanceGranularity`, `uLod` and `lodLevels`.
  `BaseGeometry.updateValues` (`base.ts:250`, `:265`) bumps `uClipObjectType` and `uLod` on every call. It runs from
  `units-visual.ts:397` and `complex-visual.ts:295`. `canMergeValues` never reads `uClipObjectType`, and `getMergeKey`
  compares `uLod` by value. So the fingerprint is stale after every update, and the full check plus the re-record
  runs. Rerun `fp.mjs 7qpd.fw2.cif`: "sorts=7, changed 7/7" for both recolors and for alpha. The fix (a narrower key
  filter, or `uLod` set only when it changed) is outside `gl/`.
- **07/P20 / 09/P17**: v2 has this as 09/P17, verified in v2 at medium with ÷~17 on 1tup and ÷1–3 on real capsids.
  07 lists it as new with ~2–3 MB, which leaves out the between-copy edges that must stay complex. Keep one id (09/P17).
- **08/P15+P12**: `scratchpad/pa3/08/verify/vC-tube.txt` and `vC-mc.txt` hold the trace output. In the tube runs only
  run 4 of 6 re-allocated (33.8 MB, 41.9 ms against 5.6 ms), which matches "after a major GC that clears the maps",
  not every GC.
- **08/P18**: fontSize = 64·(q+1) and radius = fontSize/3, so `lineHeight` is ~405 at q3. The atlas texture is a
  43 MB `Uint8Array` at q3, allocated per atlas (side note, not this finding).
- **08/P19**: MC lines `get` re-adds every line (`builder.ts:110`). Alloc in some runs (15 MB) is the bench caller's
  boxing, not the builder.
- **08/B10**: `font-atlas.ts:301–304`. When a narrow glyph overflows, `currentY` has already moved down a row, so
  later glyphs cannot use the rest of the last row either (minor).
- **09/P19**: `distinct.ts:119` `LRUCache.create(16)`, key `${count}|${JSON.stringify(p)}`, where `p` holds the
  compartment's hue range. One key per compartment, visited in order, so with 17 or more keys every lookup misses.
- **09/P20**: `pc.mjs 1e7 20`: 46 ms ×3. The 40 MB is off-heap and does not show in the alloc counter.
- **09/B8**: `color.ts:129–143`. Even with texel-centre sampling, 10 samples cannot show 11 stops, so the sample
  count has to change as well (the report's fix does both).
- **10/P17**: `sequence-scene.ts:42` sends one `version` for every source. `scene.ts:1608–1630` visits every
  structure.
- **10/P18**: `host.ts:166` `run` → `draftOf` (`state.ts:51`) → `writesOf` (`state.ts:105`), on every append to an
  unbounded record.
- **10/B7**: `scene.ts:331` (getter), `:935` (`emptiedNode`), `:1194` (`allHidden` of no eyed rows is false), `:1325`.
  `scene-tree.tsx:752`/`:1061` sources, `:994` eye drawn from `hasEye` only on a render.
- **11/P12**: `shell-view.tsx:229–233` (`StatusLine … after={props.agent}` inside `canvasBottom`),
  `canvas-chrome.tsx:1141` (`{props.after}` in the status DOM, no portal), `agent-turn.tsx:282` `ScrollArea`.
- **12/P10**: `options.ts:37` default `'auto'`, then `compose-plugin.ts:111`, then `graphics-backends.ts`'s
  `getAllBackends()`. Making WebGPU lazy needs async backend selection in frozen code, and it only saves bytes on
  WebGL-only clients.
- **12/P12**: `npx vitest run scripts/workspace/_test/typecheck-full.test.mjs` took 2.86 s (tests 99 %). The file has
  no tag or `@module-tag`.
- **12/P13**: `structure-info.ts:153–165` (`handleModel`), `:97–102` (`setMaxIndices`), `:222` (called for every
  created root model).
