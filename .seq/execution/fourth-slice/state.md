# State: the fourth slice

## Ready to start

The third slice is closed (2026-11-20; its record in `.seq/records/third-slice/`): RR merged, the size bounds set
(195 / 780 / 1250 / 880 kB), its fold-back applied, its worktrees cleaned (CE's kept, parked). Start with
`/loop run .seq/execution/fourth-slice/loop.md`.

**Base**: `packages-seq`. Every brief's base is set to its head when the lane starts.

## Decided

- Gate 0 (the author, 2026-11-19): "Accept all" (plan, Gate 0).
- Stage 0's designs agreed and promoted: structure syncing (a library, not a port), track header parts, seq data ports,
  extension openings (DR, SP, XD; plan, Stage 0).
- Lanes at once: two or three, each strong or fast by its work (the author, 2026-11-20: "2-3 lanes, all either heavy
  or light depending on context"); limited by files: three tracks (molecules, formats and data, views; plan, Order and
  parallelism).
- DP kept as a strong lane (nothing else builds the seq data ports); the order of HP1 and ST1 is the orchestrator's
  (the author, 2026-11-20: "order doesn't really matter, decisions as recommended").
- `@molstar/core`'s standard components (ML1): "Accept" (2026-11-20). BW's range reads: the surface in
  `drafts/bw-range-reads.md` "Accept", 64 KiB pages, "fail on change" (2026-11-20). The first fetch of a service's
  fixtures is still asked when SV reaches it.
- A unit's positions count from 0, as every system's ("From 0, as built", 2026-11-20); molecules.md reworded.
- Small findings collect in `fixes.md` for the fix points; the author's findings from checking the third slice go
  there as they come (the author, 2026-11-20).
- A label that does not fit its box: an ellipsis inside before it is dropped, at least 4 characters by default,
  configurable (`labels.minChars`; the author, 2026-11-20; `fixes.md`).
- Module size (§module-size, `.seq/_checks/file-size.ts`): lane SZ splits the four files over 1500 lines after FX1,
  while no other lane is in them (the author, 2026-11-20: "Yes, add all three").
- Size targets are bounds once reset after RR (the author, 2026-11-20: "reduce render first, then move the size
  targets").

## In progress

- 2026-11-20: the loop started. ML1 (strong; `seq/ml1-lane`), FX1 (strong; `seq/fx-1`; fixes.md's open items but the
  orchestrator's looks), AR (fast; `seq/ar-lane`), each in its own worktree from `e2769319e`. ML1 builds the
  `@molstar/core` standard components last, in commits that can be dropped, for the author's acceptance.

## Next

1. On each report: check, look (FX1's and AR's stories), review in its own worktree, fix round, merge, record.
2. Then each track in the plan's order: DP after ML1; SZ after FX1; MC beside ML2 (after ML1's node shapes).

## Findings by lane

One row a merged lane (loop, step 5, record): review findings and look findings by kind.

| Lane | Merged | Bug | Missing design | Look | Notes |
| --- | --- | --- | --- | --- | --- |
| AR | 2026-11-20 | 3 (review: duplicate names past the bounds; a written archive failing the ratio; case-folded names) | 0 | 0 | six low (a Configuration table, host spec bounds, version, exports, the fuzz's asserts, the story's reads); the look clean |
| NG | 2026-11-20 | 3 (review: public data access called missing in the _spec; a partial copy of layerSource; tile-edge ranges by label only) | 0 | 0 | three low (the import scan's dynamic imports, the variant heads' width, the missing values on core's side) |
| ML1 | 2026-11-20 | 1 (a section's range unbounded on units and molecules) | 1 (a unit's positions: from 0, decided) | 0 | the core table duplicated `AminoAcidNames` (moved to `util/sequence-constants`); five low |
| FX1 | 2026-11-20 | 4 (review: the levels fix's wiring untested; a whole level still built when none was drawn; an empty band while building, against EC-169; a second build budget) | 0 | 2 (look: the "2 hidden" chip over the overview; a narrow frame squeezing the alignment, decided for FX2) | 10 items fixed, 4 partly, the rest to FX2 with the decided designs |
| SZ (1) | 2026-11-20 | 0 | 0 | 0 | review: no behavior change; an import cycle grown to 7 modules, the ranges modules unlisted, three nits |
| DP | 2026-11-20 | 6 (review: malformed descriptions crashed; values on the wrong unit; row-order content; tile nodes losing identity; narrowed runtimes compared by object; a removed node's selection unreported) | 0 | 0 | public names trimmed to the design's |
| BW | 2026-11-20 | 6 (review: range reads unused by any host; memory exhaustion on malformed files; If-Range breaking servers; an abort shared across readers; maxPages not a memory bound; offsets past 2^53) | 0 | 0 | core helper reuse left (new public names) |

## How to verify

The loop's Checks on `packages-seq`; the plan's ticks against `git log --merges`.
