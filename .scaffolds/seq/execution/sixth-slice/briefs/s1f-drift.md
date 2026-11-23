# Lane S1F: stage 1's drift notes resolved

**Model**: Sonnet (each note says what to do). **Base**: `packages-seq` at the orchestrator's head (`git reset --hard
packages-seq` first). Branch `seq/s1f`.

**Why**: the author, 2026-11-23, on stage 1's fifteen drift notes: "All as recommended". Each note in
`.seq/execution/sixth-slice/drift/` has its options and a recommendation; apply the recommendation.

**For each note**:

1. Where the recommendation keeps the code as it is (amend the design): change the design's text where it states the
   old decision (the Decision, Configuration or Defaults section) to what the code does, in one or two sentences, with
   "(amended 2026-11-23: the code's way, the author as recommended)". Keep the design's reasons where they still hold.
2. Its `## In the specs` row: it now points at the spec section that states the code's way (the design agrees with it),
   no longer "design-only: not built as designed".
3. Where the recommendation needs code (blosum62-provenance's picker text, packages-and-examples' folder-graph test,
   layouts' band history step, the TODO's `composeLanguage` entry to drop, unknown-kinds' version-aware validate), the
   row stays design-only ("not built: the close-out's"), and you list it in your report: the orchestrator places it.
4. Append to the note: "**Resolved** (the author, 2026-11-23, as recommended): <what was done, where>".

**Where**: the designs in `.seq/design/`, the drift notes. Nothing else: no spec, no code, no hot file, never
`state.md`, `plan.md` or `loop.md`. Another lane (S1E) edits comments in `packages/` and the core spec at the same time:
do not touch them. Write only in your worktree and your own temp space. No network.

**Checks** (after `git add .seq`): `node .seq/_checks/designs-in-specs.ts`, `node .seq/oz/_checks/principles-ids.ts
--package .seq/oz --instance .seq --allow '^(?!\.seq/)'`, `node .seq/oz/_checks/merge-gate.ts --package .seq/oz
--instance .seq`, `pnpm format:check`.

**Report** (under 200 words): commits; per note what changed; the code items for the orchestrator to place.
