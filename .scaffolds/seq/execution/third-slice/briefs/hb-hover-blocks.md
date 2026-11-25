# Lane HB: the hover's outline on its own system only

**The bug** (the author, 2026-11-19: "hovering a CDK copy also highlights on cyclin, but only when within range"):
Complex › Blocks › Split (`packages/seq/examples/src/complex/blocks.ts`, `copies: 'each'`): hovering a residue of
CDK2 (C) outlines the same position in CDK2 (A), and also in cyclin A (B) and (D) while that position is within the
cyclin's length. In a section of blocks each block is a strip with its own axis; `hoverCell`
(`packages/seq/core/src/render/marks.ts`) and the letters' hover pass (`render/letters.ts`, `drawHover`) outline the
hovered position in every strip whose axis holds it (`axis.cell(h.position)`), never checking the block's system
against the hit's (`Hit.system`).

**Fix**: a layer's hover pass outlines a position only on a strip whose system is the hit's. Another copy of the same
chain (CDK2 (A) for a hover on (C)) still lights: one Location keyed by system, by design
([layouts.md](../../../design/layouts.md), Copies; `view/project.ts`). Apply the same check to every layer type's hover
pass that can be drawn over blocks (letters, logo, ranges, values, grid, any other using `hoverCell` or a position
outline). Allocation once warm stays 0 (the hover's guard). Minimal code.

**Tests** (fail on the parent): Complex › Blocks in `each`, `row`/`stack` and on a ring: a CDK2 position lights both
CDK2 copies and no cyclin, at a position within the cyclin's length and past it; canvas and SVG; the hover allocation
guard unchanged.

**Not in scope**: anything else.

**Model**: fast (sonnet). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/hb-hover-blocks
<that commit>`.

**Checks**: `pnpm test:unit packages/seq/`; `pnpm test:medium packages/seq/`; `npx tsc -b packages/seq/core
packages/seq/viewer packages/seq/examples`; `pnpm -C packages/seq/core typecheck:test`; `pnpm -C packages/seq/storybook
typecheck`; `pnpm lint`; `pnpm format:check`; `pnpm -C packages/seq/core test:perf`. Green only on exit 0 with a
non-zero test count.

**Headers**: a modified file keeps its header, year range to 2026, `@author David Sehnal <david.sehnal@gmail.com>` added
if missing. No mention of `.scaffolds/seq/` in code.

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`,
`biome.json`, `vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/` (write the entry's text in the report),
every `package.json` and `tsconfig*.json`, anything in `.scaffolds/seq/`.

**Stop when** green, or a check fails twice after a fix.

**Report**: commits; checks with exit codes and test counts; the stories to look at in light and dark.
