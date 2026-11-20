# State: the fourth slice

## Not started: waiting on the third slice's close

The loop starts when the third slice is closed (its `loop.md`, Finalize): lane RR (the render reduction) merged and
the size targets reset to its measured sizes with room; the fold-back applied; the record written; the worktrees
cleaned. Then run `/loop run .seq/execution/fourth-slice/loop.md`.

**Base**: `packages-seq`. Every brief's base is set to its head when the lane starts.

## Decided

- Gate 0 (the author, 2026-11-19): "Accept all" (plan, Gate 0).
- Stage 0's designs agreed and promoted: structure syncing (a library, not a port), track header parts, seq data ports,
  extension openings (DR, SP, XD; plan, Stage 0).
- Lanes at once: one strong and one or two fast (the author, 2026-11-20: "1+1-2").
- DP kept as a strong lane (nothing else builds the seq data ports); the order of HP1 and ST1 is the orchestrator's
  (the author, 2026-11-20: "order doesn't really matter, decisions as recommended").
- A public API change in `@molstar/core` (BW's range reads, the standard components) and the first fetch of a
  service's fixtures are asked when a lane reaches them (loop, Stop and ask).
- Small findings collect in `fixes.md` for the fix points; the author's findings from checking the third slice go
  there as they come (the author, 2026-11-20).
- A label that does not fit its box: an ellipsis inside before it is dropped, at least 4 characters by default,
  configurable (`labels.minChars`; the author, 2026-11-20; `fixes.md`).
- Module size (§module-size, `.seq/_checks/file-size.ts`): lane SZ splits the four files over 1500 lines after stage
  1's gate (the author, 2026-11-20: "Yes, add all three").
- Size targets are bounds once reset after RR (the author, 2026-11-20: "reduce render first, then move the size
  targets").

## Next

1. FX1 (fast; `fixes.md`'s open items and the author's third-slice findings) beside ML1 (strong; `briefs/ml1-*.md`).
2. Then the plan's order: DP, ML2, MC, MK, AC, MR, CU, MX, stage 1's gate, SZ, the header and syncing lanes.

## Findings by lane

One row a merged lane (loop, step 5, record): review findings and look findings by kind.

| Lane | Merged | Bug | Missing design | Look | Notes |
| --- | --- | --- | --- | --- | --- |

## How to verify

The loop's Checks on `packages-seq`; the plan's ticks against `git log --merges`.
