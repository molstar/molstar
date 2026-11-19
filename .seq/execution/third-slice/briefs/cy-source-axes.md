# Lane CY: cyclic and source-backed axes

**Goal**: the gaps the plan's CY item lists (found by DP, DT, DG and DC, 2026-11-19), each designed already. Read the
item in [plan.md](../plan.md) first; the designs are [axis](../../../design/axis.md), [layouts](../../../design/layouts.md)
(rings, the linear detail), [sources](../../../design/sources.md) and [grids-and-rows](../../../design/grids-and-rows.md).

1. **A length from coverage**: a sequence whose letters come from a source has its length from the source's coverage
   wherever an axis is read: resolve's static model, the ring, a linear detail across the origin, the GFF3 wrap (pairs
   are done), so a cyclic `seq` file works as inline letters do.
2. **Ranges' arrows by item px**: whether an item draws its arrowhead is decided per item's px width, not per layer
   step (a 2.7 kb ring at fit shows its arrowheads, EC-232).
3. **Focus by the shortest arc** for an item across the origin (EC-12): DP's `it.fails` test passes (turn it to `it`).
4. **A ring's overflow as density**: a ranges layer that overflows its lanes on a ring draws its density as a ring
   band (radially), never refused (the 200 kb ring, overlapping sites).
5. **A grid over a source's letters**: its values read after the source opens, not at resolve (the 1EHZ dot plot reads
   its letters from the file); the tRNA arcs' legend of pair kinds; the modified residues' codes as labels where they
   fit.
6. **Blocks over source sequences**: a blocks section gets its axis from coverage (`seq.render.axis.pending` today).
7. **Row keys on a source alignment**: `(row :key …)` resolves on an alignment read from a source
   (`seq.query.row.unknown` today), once the source's rows are held.
8. **The text layout's `chrome`** never throws on an alignment source not yet open (`render/layout/text.ts`, `p.marks`
   undefined).

Then the example views read their data from the files: DP's Plasmid views read the record from its `seq` file (and
draw sites over their recognition sequence), DT's 1EHZ dot plot its letters, DG's Hemoglobin blocks drop their
explicit ranges, DC's Complex its letters (`derived-data.ts`'s inline letters dropped where the files hold them), QS's
row query back on the sourced seed (the inline four rows kept as a second example). Their facts tests stay green.

**Tests** (fail on the parent): one a gap, in core where the gap is; the views' twins over the files.

**Not in scope**: the TW lane's track padding and well shapes (it may run beside you later: `render/tracks.ts`'s
background and emphasis passes are its); new language params (none needed: stop if one is).

**Model**: strong (opus). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/cy-source-axes
<that commit>`. NB (the numbering kind) may be in its fix round beside you (`resolve/model.ts`'s numberings,
`query/numbering.ts`, the ruler): merge carefully, keep out of those.

**Checks**: `pnpm test:unit packages/seq/`; `pnpm test:medium packages/seq/`; `npx tsc -b packages/seq/core
packages/seq/viewer packages/seq/examples`; `pnpm -C packages/seq/<core|viewer|examples> typecheck:test`; `pnpm -C
packages/seq/storybook typecheck`; `pnpm -C packages/seq/core test:perf`; `pnpm lint`; `pnpm format:check`; `pnpm -C
packages/seq/core size` before and after. Exit codes explicit, non-zero test counts. Look at the changed views in
Storybook in light and dark.

**Headers**: a modified file keeps its header, year range to 2026, `@author David Sehnal <david.sehnal@gmail.com>` added
if missing; a new file has the Mol\* header with that author. No mention of `.seq/` in code; code cites its `_spec`.

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`,
`biome.json`, `vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/` (the entry's text in the report), every
`package.json` and `tsconfig*.json`, anything in `.seq/`.

**Stop when** green, or a check fails twice after a fix, or the work needs public API or language past what exists, or
a permission is refused.

**Report**: commits; checks; sizes; what you saw; drift; followups.
