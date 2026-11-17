# State: first-slice   (updated 2026-11-17, by the orchestrator)

Branch: `packages-seq`.

**Overnight run** (the author, 2026-11-17: "continue all the way to the end of stage 7 as an overnight run"): the loop
continues through every gate to the end of stage 7 without waiting; decisions a gate would ask for are made as
recommended and listed in gaps.md for the author's morning review; the stop conditions of loop.md (a refused
permission, a check red twice, a public API change beyond what is agreed, a budget a design blocks) still stop it, and
the stage 7 performance targets are measured and reported, not confirmed.

**Gate decisions** (the author, 2026-11-17: "also for gate decisions: if there is a way to make something configurable,
then do that as choose your recommendation as the default"): a choice a gate or a lane would make one way is made a
config item (or a param) with the recommendation as its default, and listed in gaps.md.

## Done

- Stage 1 (lane A, the tree core): merged `fbc27cde5`.
- Stage 2 (lanes B language, C queries, D protocol and store): merged `77d03dfc1`, `5fa7c43d4`, `6e885495a`; each
  reviewed and fixed in one round; the gate green on the merged tree (233 tests, about 0.3 s; tsc with tests, lint,
  format, scaffold checks). Designs updated: format (tagged columns), language-description (column `required`,
  `seq.query.*`), locations-and-queries (`chain` as Nodes, `:of`, quoted `.n` labels, the boundary rule), protocol
  (`AssetReader`, `form`, `list`, row identity, malformed not retried, `over_budget`), data-store (want sets, `evicted`).
- Stage 3 (lanes C2 range sets into core, E resolve): merged; gate green (500 tests); `seqx-core` exports `./resolve`
  and `./view`.

## In progress

- Stage 4: lane G (viewer models) merged after one review round (merge, then `4643e6653`: `./models` and `./behaviors`
  entries); the merged tree 547 tests (repr, data/int, packages/seq), about 0.4 s; tsc, lint, format, scaffold checks
  green. Lane F (render, branch `seq/f-render`, checks green when rerun: 334 tests) is under review (strong); its fix
  round, then the merge with the `./render` export and `_test` null.

## Next

- Merge F; reconcile its hit type with G's (`models/hover.ts` `Hit`: `{layer?, section?, system, position?, element?,
  extent?, target?}`); fold `drift/render.md` and `drift/viewer-models.md` into the designs; the stage 4 gate.
- Then lane H's brief (stage 5) and start it.

## Decisions taken (awaiting the author)

See gaps.md, "Built but not agreed".

## Verify

`pnpm test:unit packages/core/src/repr packages/seq/`, `npx tsc -b packages/core packages/seq/core packages/seq/viewer`,
`pnpm -C packages/seq/core typecheck:test`, `pnpm lint`, `pnpm format:check`, the scaffold checks.
