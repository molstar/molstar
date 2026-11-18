# State: second-slice   (updated 2026-11-18, by the orchestrator)

Branch: `packages-seq`. Not started: waits for the first slice's close and Gate 0 (`plan.md`).

## Done

- Plan, loop and state written; the loop rewritten for this slice.
- Before the loop (the author, 2026-11-18): the history designs corrected; the patch and diff design drafted
  (agreed the same day as `design/patches.md`, git's model); Gate 0 aligned with the first slice's record;
  stage 0's lane renamed S0; the first slice's worktrees and branches removed.

## In progress

- None.

## Next

- A focused Fable 5.2 re-review of the interfaces (`design/patches.md`, `drafts/merge.md`, `drafts/immutable.md`), after
  the second review's points 1 to 5 were applied (the author, 2026-11-18: "let's go with 1-5").
- The author reviews `drafts/merge.md` and `drafts/immutable.md`; then the loop starts: Gate 0 to the author (plan.md,
  with the first slice's open decisions and evidence), stage 0 beside it.

## Decided

- History as a Seq\* port, the host owning the undo stack; `@molstar/ui/history` one adapter (the author,
  2026-11-18: "yes, go ahead"). Gate 0's history item is settled; the designs say so (changes.md, tree-core.md, base-library.md,
  viewer-models.md).

## Decisions taken (awaiting the author)

None yet.

## Verify

As loop.md, Checks.
