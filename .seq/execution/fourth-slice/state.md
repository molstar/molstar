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
  targets"). Moved for the query language: reading a spec 230, headless 820, minimal 920 kB (the author,
  2026-11-20), the default viewer 1,290 kB ("Move to 1290 kB", 2026-11-20).

## In progress

- 2026-11-20: QE merged (d0e354c4c). ST1 (strong; `seq/st1-lane`), HP1 (strong; `seq/hp1-lane`) and AC (strong;
  `seq/ac-lane`, from d0e354c4c) build in their worktrees.

## Next

1. QE's heavy bound (`:index`, `itWithin` 400 ms) runs at stage 1's gate (the lane could not run the heavy lane).
2. ST1 and HP1: check, look (HP1's stories), review, fix round, merge, record.
3. HP2 after ST1 and HP1; MR, CU, MX on the molecules track; FX3 after stage 1's gate.

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
| ML2 | 2026-11-20 | 2 (review: a reference to a unit without a system threw; a ranged unit over a source drew no letters) | 2 (the atom form's meaning, decided: monadic, then `elements`/`groups`; equivalents hooks) | 0 | internals exported; model.ts 1410 lines |
| MC | 2026-11-20 | 5 (review: a crafted operator expression exhausted memory; quadratic lookups and a rethrown error; branch links on the complex; PO4 classed by what was modeled; the fuzz's own crash) | 0 | 0 | the io allowlist narrowed; the re-count made independent; the lazy chunk's imports reported (90 kB) |
| FX2 | 2026-11-20 | 3 (review: the rail not back when a narrow frame hid the overview; a large group's card cut; a possible resize loop) | 2 (no tree toggle; model issues unwired) | 0 | the ruler's allocation; the look clean |
| MK | 2026-11-20 | 2 (review: provenance on the wrong run; overlaps on `to` unchecked) | 2 (the ranking rule and an unused via, decided) | 0 | a weak property test; model.ts 1483 → 1370 |
| SZ (2) | 2026-11-20 | 0 | 0 | 0 | review: no behavior change; a test anchored on shims, viewer-motion's viewports, two nits (to the next fix lane) |
| QE | 2026-11-20 | 3 (review: an empty column printed as everything; a pending property dropped in node filters, EC-169; column lengths unbounded before zipping) | 0 | 0 | `:index` scanned rows × residues (4.2 s → indexed); monadic laws untested; internals exported; PYL; the default viewer 2.7 kB past its bound, moved to 1,290 kB (the author) |

## How to verify

The loop's Checks on `packages-seq`; the plan's ticks against `git log --merges`.
