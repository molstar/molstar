# State: first-slice   (updated 2026-11-17, by the orchestrator)

Branch: `packages-seq`. No lanes started.

## Done

- The plan approved and Gate 0 settled (the author: "approve both"; range reads wait for bigWig, the slice reads its
  FASTA whole). The settled choices are copied into language-description, locations-and-queries and protocol.
- `loop.md` written; lane A's brief written (`briefs/a-tree-core.md`).

## In progress

Nothing.

## Next

Start lane A (stage 1) per `loop.md`, step 4, with the brief and the meso\* source path in the launch instructions.

## Decisions taken

- Gate 0, as recommended (see `plan.md`).

## Blockers / open questions

None.

## Verify

The scaffold checks and `pnpm format:check` (green at the loop's commit).
