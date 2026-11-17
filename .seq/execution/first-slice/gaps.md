# Gaps (§gaps-queue)

## Built but not agreed (the orchestrator's decisions, to confirm)

1. Undeclared typed columns are tagged in JSON and the S-expression printing (`drift/custom-columns.md`).
2. A range with start > end on a cyclic sequence wraps the origin; elsewhere an error (`drift/cyclic-ranges.md`).
3. Queries (`drift/query-decisions.md`): `chain` returns Nodes (a system); `.n` labels quoted; a boundary at a system's
   end belongs to a range ending there; Location equality structural (rects not merged).
4. The protocol and store (`drift/protocol-details.md`): malformed answers not retried; the ranges form `list`; a stable
   `row` identity for ranges sliced across tiles; per-consumer want sets; the oldest stale request cancelled at once
   when view work waits; a narrow `AssetReader` for runtimes.
5. Lane B's `seq.lang.*` codes and letter alphabets (`drift/lang-codes.md`).

## Agreed but not built

Stages 3 to 7; the follow-ups in `followups.md` (the tree core's pending union, version-aware validate, `checkRefs`,
`composeLanguage`, patches, dedupe; the theme slot; `rect`/`row` extras; byte dedupe across tiles).

## Built but not adopted

- `seqx-core`'s `query` and `data` have no user yet (lane E is their first).

## Waiting for the author

- Items 1 to 5 above.
- Whether range-set operations (Location's union, intersect, subtract) move to `@molstar/core/data/int`.
- The edge-case outcomes proposed in `.seq/design/edge-cases.md`.
