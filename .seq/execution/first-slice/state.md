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
- Stage 4 (lanes F render, G viewer models): merged, each after one review round; `seqx-core` exports `./render`,
  `seqx-viewer` `./models` and `./behaviors`; the viewer's `Hit` extends the render hit; gate green (581 tests).
- Stage 3 (lanes C2 range sets into core, E resolve): merged; gate green (500 tests); `seqx-core` exports `./resolve`
  and `./view`.

## In progress

- Stage 5: lane H merged after one review round (merge, then `d410648d7`: entries `./context`, `./canvas`, `./ports`;
  the fixtures plugin registered). The merged tree: 665 quick tests, typechecks, lint, format, medium lane, `pnpm build`
  green. `check:workspace` fails on a pre-existing `@molstar/ui` export (`./components/logo/_generated/
  glass-wordmark.module.css`, no source; `packages/ui` is unchanged from master): not this slice's.
- The integration gate's visual pass (Storybook on 6007) found two bugs, being fixed before the gate passes:
  1. the Tracks story's canvas draws the Conservation values only on the right third (the SVG draws all): a fixer
     (strong) on branch `seq/fix-values-canvas`;
  2. the viewer and the stories stay light in the dark theme: the browser tokens are not read from the UI root's
     computed style (tokens.ts says they should be); next, after fix 1.

## Next

- Merge both fixes, rerun the checks, screenshots light and dark, pass the gate.
- Then lane I (`briefs/i-regions.md`, written).

## Decisions taken (awaiting the author)

See gaps.md, "Built but not agreed".

## Verify

`pnpm test:unit packages/core/src/repr packages/seq/`, `npx tsc -b packages/core packages/seq/core packages/seq/viewer`,
`pnpm -C packages/seq/core typecheck:test`, `pnpm lint`, `pnpm format:check`, the scaffold checks.
