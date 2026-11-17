# Gaps (§gaps-queue)

## Built but not agreed (the orchestrator's decisions, to confirm)

1. Undeclared typed columns are tagged in JSON and the S-expression printing (`drift/custom-columns.md`).
2. A range with start > end on a cyclic sequence wraps the origin; elsewhere an error (`drift/cyclic-ranges.md`).
3. Queries (`drift/query-decisions.md`): `chain` returns Nodes (a system); `.n` labels quoted; a boundary at a system's
   end belongs to a range ending there; Location equality structural (rects not merged).
4. The protocol and store (`drift/protocol-details.md`): malformed answers not retried; the ranges form `list`; a stable
   `row` identity for ranges sliced across tiles; per-consumer want sets; the oldest stale request cancelled at once
   when view work waits; a narrow `AssetReader` for runtimes.
5. Lane B's `seq.lang.*` codes and letter alphabets (`drift/lang-codes.md`).

6. Configurable defaults added (the author's rule): `ResolveConfig` (`maxPathSteps` 4, `maxPathCandidates` 16, `encoding`
   `cheapest`, `openQueryLength` 2^31−1); `sequence.name` (default: the node's id); `map.table` (NCBI code, default 1;
   1, 2 and 11 built); `followRange`'s `arc` (`shortest`).
7. Paths are cached across resolves keyed by the map entries of every system the search read (safer than the entries
   used; lane E's call).
8. The TP53 exon coordinates in `resolve/_test/p53-geometry.ts` were written from memory (lengths add to 3 × 393 + a
   stop): check them against RefSeq NM_000546 before relying on them.

9. Lane F's `RenderConfig` defaults: density thresholds (device px per position: letters 14, cells 3, dots 1, bars
   1/64), label column 112 and gap 8, inset 8, track gap 4, section gap 8, track heights (letters 20, ranges 36, values
   32), lane 12 and gap 2, letters 9–14 px, mono advance 0.6 em, fonts (`12px sans` labels, `10px sans` items), fade 12
   px, squircle radius (a third of the height, max 8, none under 6 px), outline widths (hover 1.5, selection 1.5, cursor
   2.5), want margin (half the view each side), `followArc` shortest; `renderSvg` (`timeoutMs` 10 s, light theme, no
   assets unless a reader is given, `dpr` 1); a values layer's `domain`, else open (`drift/render.md`).
10. Lane G's `ViewerConfig` defaults: `hover.closeCodons` true, `viewport.maxPxPerPosition` 32,
    `viewport.groupEqualAxes` true, `viewport.onResize` `keep-zoom`, `follow.followerMoves` `ignore`, `follow.cyclicArc`
    `wrap`, `navigate.wheelZoomStep` 1.25, `navigate.keyPanFraction` 0.1, `navigate.coalesceWheel` true; the default
    preset has no click or drag-select bindings until the selection lane (`drift/viewer-models.md`).
11. Lane H's defaults: `canvas.height` `content` (max 2048), `canvas.dpr` `device`, `canvas.settleMs` 150,
    `pointer.rehit` true, `pointer.pinchFromCtrlWheel` true, `pointer.pinchSpeed` 0.01, `pointer.pinchMaxDelta` 50;
    `BrowserAssetPolicy` (https and relative under the page, confined to the base folder, 256 MB, 30 s, no credentials),
    `HeadlessAssetPolicy` reads nothing (`drift/viewer.md`).
12. The theme fix's `theme.source` (`root` when mounted: tokens from the UI root's computed style, following
    `data-msui-theme`; `config` the alternative); a headless viewer uses the config's `theme`.
13. `pnpm check:workspace` is red on `@molstar/ui`'s export `./components/logo/_generated/glass-wordmark.module.css`
    (no source; pre-existing, `packages/ui` unchanged from master): for the UI's owner, not this slice.
14. Lane I's defaults: `StoreSettings.summaries` (`factor` 16, `maxLength` 2^24; `false` off; not yet in
    `ViewerConfig`), `ResolveConfig.maxNestedDepth` 2 (EC-179), `for_each` (`depth` 0, `order` `key`, `-` descending,
    `filter` none, string equality), `delayedSource`/`delayedRuntime` (`delayMs` 0, real clock), node retention (the
    entries that carried a node, its source's coverage, an explicit `hold`) (`drift/regions.md`).

## Agreed but not built

Stages 3 to 7; the follow-ups in `followups.md` (the tree core's pending union, version-aware validate, `checkRefs`,
`composeLanguage`, patches, dedupe; the theme slot; `rect`/`row` extras; byte dedupe across tiles).

## Built but not adopted

- `seqx-core`'s `query` and `data` have no user yet (lane E is their first).

## Waiting for the author

- Items 1 to 5 above.
