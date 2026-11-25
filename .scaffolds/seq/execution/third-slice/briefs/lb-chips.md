# Lane LB: a one-position block is a chip in every layout

**Goal**: [layouts](../../../design/layouts.md), Blocks, "A one-position block is a chip in every layout" (decided, the
author, 2026-11-19: "the ligand rendering really isn't good", "chips is a prolly good way to do these"; "contacts
should be optional and configurable by the app ... we are building a library, not a one-off viewer"). Today a
one-position block (Complex › Blocks' ATP) is a sliver of one position's share at the fit view, its label clipped
("AT⁠P" overlapping its neighbour). Build, as the design says:
- **Tracks**: a one-position block takes a fixed slot its label's width (the text chip's rule: the label and `×N`, 5 px
  either side), never a position's share; zoom does not scale it. The chips follow the polymers' blocks (`row`) or
  share one short row of their own (`stack`); each track holding something at the position draws one pill there (the
  letters' band the label, a feature band its mark), the others leave it blank.
- **Rings**: a fixed-angle arc a chip after the polymers' arcs.
- **Text**: as built (a chip flowing like a word); check it uses the same rule.
- **The row's x is not linear over a chip**, as over a gap (EC-102): the window, the ruler, the navigators and
  `windowExtent` skip chips as they skip gaps. Hit, hover, selection and links read a chip as its one position.
- **Contacts** stay an optional `pairs` layer a builder adds (the Complex's), never drawn by default.
The `_spec`s (core's Layouts › Blocks) say it; the Design page's Blocks line if it changes. Minimal code;
allocation once warm unchanged; the blocks' timing and size guards.

**Tests** (fail on the parent): Complex › Blocks in `row`, `stack`, split copies and on a ring: the ATP chips at their
label's width at the fit view and zoomed, never clipped, canvas and SVG, both themes; a chip's hover, selection and a
link to it; the window and ruler skipping chips; a zoom about a chip.

**Look**: run Storybook from your worktree on port 6009 (`pnpm --filter @molstar/storybook exec storybook dev -p 6009
--ci --no-open`, background; stop it after) and look at Complex › Blocks (all four) in light and dark; say what you saw.

**Not in scope**: anything else. RC (rows), HL (highlights), FS-A2 (records) run beside you.

**Model**: strong (opus). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/lb-chips <that
commit>`. **The language stays 0.4.0.**

**Checks**: `pnpm test:unit packages/seq/`; `pnpm test:medium packages/seq/`; `npx tsc -b packages/seq/core
packages/seq/viewer packages/seq/examples`; `pnpm -C packages/seq/<core|viewer|examples> typecheck:test`; `pnpm -C
packages/seq/storybook typecheck`; `pnpm lint`; `pnpm format:check`; `pnpm -C packages/seq/core test:perf`. Exit codes
explicit (`; echo exit=$?`), non-zero test counts.

**Headers**: a modified file keeps its header, year range to 2026, `@author David Sehnal <david.sehnal@gmail.com>` added
if missing; a new file has the Mol\* header with that author. No mention of `.scaffolds/seq/` in code.

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`,
`biome.json`, `vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/` (write the entry's text in the report),
every `package.json` and `tsconfig*.json`, anything in `.scaffolds/seq/`.

**Stop when** green, or a check fails twice after a fix, or a public API change is needed.

**Report**: commits; checks; what you saw; drift.
