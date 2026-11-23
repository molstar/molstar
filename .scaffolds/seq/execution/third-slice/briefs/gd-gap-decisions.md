# Lane GD: the gaps decided

**Goal**: build what the author decided on `gaps.md` (2026-11-19, "go as recommended"; gap 6 "Yes, unique ids"), each
already written into its design, plus PG's last allocation item:

1. **A lists-free palette module in `@molstar/core`** ([themes](../../../design/themes.md), Sizes): the value-palette
   evaluation moves to a module in `@molstar/core/util/color` that imports no color lists; `value-palette.ts`
   re-exports it with `ValuePaletteLists`, so Mol\*'s and MolViewSpec's code is unchanged. Seq\* imports the new module.
   **Agreed surface**: that module and the re-export; anything past it stops the lane. Measure the sizes before and
   after (the gap expected about 13.6 kB off headless rendering and the viewer).
2. **Every color-blind substitute safe** (themes, "A reader's color-blind substitution"): `classes-safe`, `classes`'
   groups on Okabe–Ito hues, passing the conformance check in both themes; the unsafe protein schemes (`classes`,
   `clustal-x`, `zappo`, `taylor`) substitute to it; a test over the registry that no scheme names an unsafe substitute.
   `classes` stays the default.
3. **Letter cells in dark** (themes, "Letter cells in dark"): a hairline cell edge, 1 device px, in a new token
   `seq-cell-edge` (3:1 against the dark page; none in light), on canvas and SVG.
4. **One motion item** ([stories](../../../design/stories.md), Configuration; ux-foundations section 7):
   `view.motion: platform | cut | animate` (default `platform`), read by jumps and stories; `story.reducedMotion`
   dropped. The viewer `_spec`'s Configuration table follows.
5. **The pyramid's measured tolerance** ([grids-and-rows](../../../design/grids-and-rows.md)): the tests assert the
   recorded bounds (share within 0.05 on average and 0.3 at worst; the letter exact where the winner leads by more than
   0.25); check RW's tests already do, and the `_spec` says so.
6. **Clade names by unique ids** (grids-and-rows, Clades): a clade is named by its node's id where the id is unique in
   the tree and not a number, else `(clade :of [first last])`; a hover, a fold and `clades.keep` on two clades whose
   nodes share a label (`((A,B)95,(C,D)95)`) reach the right one.
7. **`Placement.at` reads the point from an array** (PG's review, A7): the interface and its three implementations, so a
   grid cell's `hitView` allocates 0 B beyond the test's own call; tighten `GridPathBound`.

**Model**: strong (opus). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/gd-gap-decisions
<that commit>`. **The language stays 0.4.0.** Nothing runs beside you in `packages/seq` right now.

**Minimal code**: reuse what exists; no private copies; sizes reported before and after.

**Tests**: each item with a test that fails on its parent (name the command); drawing against the data for 2 and 3 (a
substituted scheme's view; the dark cells' edges on canvas and SVG), both themes.

**Checks**: `pnpm test:unit packages/seq/` and `pnpm test:unit packages/core/src/util/color/` (core changed);
`npx tsc -b packages/core packages/seq/core packages/seq/viewer packages/seq/examples`; `pnpm -C packages/seq/<core|
viewer|examples> typecheck:test`; `pnpm -C packages/seq/storybook typecheck`; `pnpm lint`; `pnpm format:check`;
`pnpm -C packages/seq/core test:perf`; `pnpm -C packages/seq/core size`; `pnpm test:unit packages/mvs/` (it reads value
palettes). Green only on exit 0 with a non-zero test count.

**Headers**: a new file has the Mol\* header `Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file
for more info.` with `@author David Sehnal <david.sehnal@gmail.com>`; a modified file keeps its header, its year range
extended to 2026 and that author line added if missing. No mention of `.scaffolds/seq/` in code; code cites its own `_spec`.

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`,
`biome.json`, `vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/` (write the entry's text in the report),
every `package.json` and `tsconfig*.json` (ask in the report: a new export of `@molstar/core`), anything in `.scaffolds/seq/`.

**Stop when** green, or a check fails twice after a fix, or the work needs a public API past the agreed surface, or a
permission is refused.

**Report**: commits; each check with exit code and test count; sizes before and after; the stories to look at (p53 ›
Themes and Schemes Gallery with `theme.cvdSafe` on, a letters view in dark, a story's motion) and what to look for in
light and dark; drift; followups.
