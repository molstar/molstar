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

- Fable 5.2 reviews `design/patches.md` with the state semantics (agreed 2026-11-18: git's model and the three other
  decisions) and `drafts/immutable.md`
  (`@molstar/core/immutable`, used in `repr`; lane R); then the loop starts: Gate 0 to the author (plan.md, with the first slice's
  open decisions and evidence), stage 0 beside it.

## Decided

- History as a Seq\* port, the host owning the undo stack; `@molstar/ui/history` one adapter (the author,
  2026-11-18: "yes, go ahead"). Gate 0's history item is settled; the designs say so (changes.md, tree-core.md, base-library.md,
  viewer-models.md).

## Decisions taken (awaiting the author)

None yet.

## Verify

As loop.md, Checks.
