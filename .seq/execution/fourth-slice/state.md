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
- A public API change in `@molstar/core` (BW's range reads, the standard components) and the first fetch of a
  service's fixtures are asked when a lane reaches them (loop, Stop and ask).
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

## How to verify

The loop's Checks on `packages-seq`; the plan's ticks against `git log --merges`.
