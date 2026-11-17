# State: first-slice   (updated 2026-11-17, by the orchestrator)

Branch: `packages-seq`.

## Done

- Gate 0 settled; `loop.md` and the briefs.
- **Stage 1** (lane A, the tree core in `@molstar/core/repr`): merged `fbc27cde5`; the review's blockers fixed in one
  round; the gate green on the merged tree (73 tests, about 0.3 s; tsc, lint, format, scaffold checks); the changelog
  entry `.v6/changelogs/core-repr.md`; format.md (tagged undeclared columns) and language-description.md (a column's own
  `required`) updated.

## In progress

- Stage 2: lanes B (language, fast model), C (locations and queries, strong), D (protocol and store, strong), from the
  merge commit, briefs in `briefs/`.

## Next

On each report: verify, review (strong), one fix round, merge (exports in `packages/seq/core/package.json` are the
orchestrator's: the lanes ask for them), gate, report. Then stage 3 (lane E).

## Decisions taken

- Gate 0 as recommended.
- Undeclared columns tagged in JSON and the printing (the orchestrator, from lane A's review): to confirm with the author.

## Blockers / open questions

- The author's confirmation of the tagged columns (non-blocking: lanes C and D do not depend on it).

## Verify

`pnpm test:unit packages/core/src/repr packages/seq/`, `npx tsc -b packages/core packages/seq/core packages/seq/viewer`,
`pnpm lint`, `pnpm format:check`, the scaffold checks.
