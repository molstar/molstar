# Lane RF: the rows' numbers and the last row

**Goal**: three findings by the author in Globins › Rows (the 10,000-row views, Scale with navigators; light):

1. **The end numbers in view** (decided, grids-and-rows, The row header): a row's first and last residue numbers in its
   gutters are those of the residues in view, following a pan and a zoom; at the fit view the row's own; a residue
   counted, a gap never; a row with no residue in view shows none. Today they are the row's overall numbers and do not
   change when the view shows columns 50–150.
2. **An end number past the row's range**: GEN_4844/7-147 shows 247 as its last number. Suspect the lazy rows' count
   (`render/rows.ts`, the last residue numbers "counted once a row" from tiles) counting columns or gaps, not residues;
   check held rows too. The numbers must agree with the label's range where the label has one.
3. **The last row cut off**: at the bottom of the rows band the cells stop before the last label's row, so a number
   shows beside no cells (the author: "there is number label, but the rendering is cut off too soon"). A partly shown
   row draws its cells as far as its label and numbers, clipped at the band's foot, on canvas and SVG alike.

Each with a test that fails on its parent (name the command): numbers at the fit view and after a pan and a zoom, held
and lazy; the count against the label's range over the generated 10,000 rows; the last partial row's cells reaching the
band's foot. Draw stays 0 B warm (the allocation tier's rows cases); the lazy frame stays within its target.

**Model**: fast (sonnet). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/rf-rows-fixes
<that commit>`. Lane NV runs beside you (navigators; it may touch `tracks.ts`): keep to the rows layer and its tests.

**Minimal code**: reuse the rows' existing number and tile paths; no private copies.

**Checks**: `pnpm test:unit packages/seq/`; `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`;
`pnpm -C packages/seq/core typecheck:test`; `pnpm lint`; `pnpm format:check`; `pnpm -C packages/seq/core test:perf`.
Green only on exit 0 with a non-zero test count.

**Headers**: a modified file keeps its header, its year range to 2026, `@author David Sehnal <david.sehnal@gmail.com>`
added if missing; code cites its own `_spec` sections; update the core `_spec`'s rows section (the numbers rule).

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`,
`biome.json`, `vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/`, every `package.json` and `tsconfig*.json`,
anything in `.seq/`.

**Stop when** green, or a check fails twice after a fix, or the fix needs a design change.

**Report**: commits; each check with exit code and test count; the stories to look at (Globins › Rows: Scale, Lazy,
Family) and what to look for, light and dark.
