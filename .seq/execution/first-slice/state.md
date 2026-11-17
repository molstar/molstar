# State: first-slice   (updated 2026-11-17, by the orchestrator)

Branch: `packages-seq`. Lane A running (branch `seq/a-tree-core`, its own worktree, strong model).

## Done

- The plan approved and Gate 0 settled (the author: "approve both"; range reads wait for bigWig, the slice reads its
  FASTA whole). The settled choices are copied into language-description, locations-and-queries and protocol.
- `loop.md` written; lane A's brief written (`briefs/a-tree-core.md`).

## In progress

- Lane A (stage 1, the tree core), started from `54532c3f7` with `briefs/a-tree-core.md`.

## Next

On lane A's report: verify its checks, review (strong, diff against the designs), merge, add the changelog entry, run
the stage 1 gate, report to the author. Then write briefs B, C, D and start stage 2.

## Decisions taken

- Gate 0, as recommended (see `plan.md`).

## Blockers / open questions

None.

## Verify

The scaffold checks and `pnpm format:check` (green at the loop's commit).
