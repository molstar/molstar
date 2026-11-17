# State: first-slice   (updated 2026-11-17, by the orchestrator)

Branch: `packages-seq`.

**Overnight run** (the author, 2026-11-17: "continue all the way to the end of stage 7 as an overnight run"): the loop
continues through every gate to the end of stage 7 without waiting; decisions a gate would ask for are made as
recommended and listed in gaps.md for the author's morning review; the stop conditions of loop.md (a refused
permission, a check red twice, a public API change beyond what is agreed, a budget a design blocks) still stop it, and
the stage 7 performance targets are measured and reported, not confirmed.

## Done

- Stage 1 (lane A, the tree core): merged `fbc27cde5`.
- Stage 2 (lanes B language, C queries, D protocol and store): merged `77d03dfc1`, `5fa7c43d4`, `6e885495a`; each
  reviewed and fixed in one round; the gate green on the merged tree (233 tests, about 0.3 s; tsc with tests, lint,
  format, scaffold checks). Designs updated: format (tagged columns), language-description (column `required`,
  `seq.query.*`), locations-and-queries (`chain` as Nodes, `:of`, quoted `.n` labels, the boundary rule), protocol
  (`AssetReader`, `form`, `list`, row identity, malformed not retried, `over_budget`), data-store (want sets, `evicted`).

## In progress

- Stage 3: lane E (resolve, strong) and lane C2 (range sets into core, fast), from the commit after this one.
  The author restarted the loop without objecting to the five stage-2 decisions; the loop proceeds on them as
  recommended, still open to change (gaps.md).

## Next

On confirmation: write lane E's brief (the static model implementing `QueryModel` incl. `chain` and annotation
`fields`; maps; paths; projection; the resolved view; binding `SourceEntry` to the language's source node) and start it.

## Decisions taken (awaiting the author)

See gaps.md, "Built but not agreed".

## Verify

`pnpm test:unit packages/core/src/repr packages/seq/`, `npx tsc -b packages/core packages/seq/core packages/seq/viewer`,
`pnpm -C packages/seq/core typecheck:test`, `pnpm lint`, `pnpm format:check`, the scaffold checks.
