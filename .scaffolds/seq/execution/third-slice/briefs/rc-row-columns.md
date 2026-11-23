# Lane RC: row columns

**Goal**: [grids-and-rows](../../../design/grids-and-rows.md), "Row columns" (the author, 2026-11-19: "an extra column
with some annotation next to each header… and being able to position it say to the right of the sequence data",
"ok, schedule it like that"; the base `packages/sequences`' data columns, Jalview's annotation columns):

1. **The rows layer's `columns`**: a list of `{ field, side: left | right, as: text | number | bar | chip, width,
   title }`, a param within 0.4.0. `left` in the header between the label and the numbers; `right` after the end
   numbers, before a vertical navigator (the rows rail).
2. **Fields** from the rows table: the Stockholm reader keeps `#=GS <seq> <tag> <value>` lines as per-row fields
   (`OS`, `DE`, `AC`, any tag; today read past in `packages/seq/core/src/data/alignment-readers.ts`); per-row measures
   (`length`, `gaps`, `identity` to the reference row, from the held cells and the existing column measures where they
   help); and a records source's fields joined by row key (a records source whose `target` is the rows system; FS-A1's
   runtime is merged). A missing field an empty cell; a field no row has reported once (`seq.render.rows.column`).
3. **Titles** in the band above the rows (beside the fixed tracks' labels); a click on one orders the rows by it
   (`order: { field, direction }`, which exists; a second click flips). Hover names a cell's value; a column's hit is
   its row's header (select, hover as the header does).
4. **Narrow**: columns yield before the label, right ones first, each down to a floor then hidden (EC-107), a hidden
   count chip in the title band. Canvas, SVG and the text figure draw the same columns.
5. **Globins › Rows**: the family view shows species (`OS`) on the left and identity on the right (DG will build on
   it); `Features/Rows/Columns` a story with its headless twin (left and right, each `as`, ordering by a title, a narrow
   width, a joined file).

The `_spec`s (core's Render › rows, Language; the viewer's) say it. Minimal code; allocation once warm unchanged (the
rows' guards).

**Tests** (fail on the parent): `#=GS` fields read; each `as` drawn in canvas and SVG, both themes, left and right;
identity and length against hand-computed values; a joined records source; a title click orders the rows; the narrow
yielding order and the chip; the text figure's columns.

**Not in scope**: anything else. HL (highlights), IZ (navigators) and FS-A2 (GFF3, the matrix read, `dataFilter`,
p53's entry) run beside you.

**Model**: strong (opus). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/rc-row-columns
<that commit>`. **The language stays 0.4.0.**

**Checks**: `pnpm test:unit packages/seq/`; `pnpm test:medium packages/seq/`; `npx tsc -b packages/seq/core
packages/seq/viewer packages/seq/examples`; `pnpm -C packages/seq/<core|viewer|examples> typecheck:test`; `pnpm -C
packages/seq/storybook typecheck`; `pnpm lint`; `pnpm format:check`; `pnpm -C packages/seq/core test:perf`; `pnpm -C
packages/seq/core size`. Exit codes explicit (`; echo exit=$?`), non-zero test counts.

**Look**: run Storybook from your worktree on port 6010 (`pnpm --filter @molstar/storybook exec storybook dev -p 6010
--ci --no-open`, in the background; stop it when done) and look at Globins › Rows › Family and the new story in light
and dark; say what you saw.

**Headers**: a modified file keeps its header, year range to 2026, `@author David Sehnal <david.sehnal@gmail.com>` added
if missing; a new file has the Mol\* header with that author. No mention of `.seq/` in code.

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`,
`biome.json`, `vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/` (write the entry's text in the report),
every `package.json` and `tsconfig*.json`, anything in `.seq/`.

**Stop when** green, or a check fails twice after a fix, or a public API past the param above is needed.

**Report**: commits; checks; sizes before and after; what you saw; drift.
