# State: first-slice   (updated 2026-11-17, by the orchestrator)

Branch: `packages-seq`.

## Done

- Stage 1 (lane A, the tree core): merged `fbc27cde5`.
- Stage 2 (lanes B language, C queries, D protocol and store): merged `77d03dfc1`, `5fa7c43d4`, `6e885495a`; each
  reviewed and fixed in one round; the gate green on the merged tree (233 tests, about 0.3 s; tsc with tests, lint,
  format, scaffold checks). Designs updated: format (tagged columns), language-description (column `required`,
  `seq.query.*`), locations-and-queries (`chain` as Nodes, `:of`, quoted `.n` labels, the boundary rule), protocol
  (`AssetReader`, `form`, `list`, row identity, malformed not retried, `over_budget`), data-store (want sets, `evicted`).

## In progress

Nothing: the loop waits for the author's confirmation of the orchestrator's decisions (gaps.md) before stage 3.

## Next

On confirmation: write lane E's brief (the static model implementing `QueryModel` incl. `chain` and annotation
`fields`; maps; paths; projection; the resolved view; binding `SourceEntry` to the language's source node) and start it.

## Decisions taken (awaiting the author)

See gaps.md, "Built but not agreed".

## Verify

`pnpm test:unit packages/core/src/repr packages/seq/`, `npx tsc -b packages/core packages/seq/core packages/seq/viewer`,
`pnpm -C packages/seq/core typecheck:test`, `pnpm lint`, `pnpm format:check`, the scaffold checks.
