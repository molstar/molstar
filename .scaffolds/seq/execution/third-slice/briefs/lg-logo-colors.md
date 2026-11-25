# Lane LG: the logo's letter colors

**Goal**: [grids-and-rows](../../../design/grids-and-rows.md), The logo, "A letter's color is its class's, made to read
as a glyph" (the author, 2026-11-19: "this letter coloring doesn't really work for the logo plot"). Today
`packages/seq/core/src/render/logo.ts` fills each stacked letter (and each bar below 2 device px) with the letters
channel's color, a scheme's pale cell fill, so the logo reads washed out (Globins › family, `classes`, light theme).
Make each letter's color its class color moved in oklab lightness until it reaches 3:1 against the track's surface,
in each theme: the existing `ink` rule of [themes.md](../../../design/themes.md) (an ink under 3:1 against its fill is
moved in oklab lightness) applied with the surface as the fill. Reuse that rule's code; hue kept. Computed once per
color and theme (the token table's version keys it), never per frame per letter.

**Tests**: every color of every built-in protein and nucleotide scheme, both themes, at 3:1 or more against the
surface as drawn by the logo (fails on the parent); a drawing check of the Globins logo in canvas and SVG, both themes;
the logo's allocation and timing guards unchanged.

**Not in scope**: anything else. Lane FS-A1 runs beside you in `packages/seq/core/src/data` and the source params:
stay out of them.

**Model**: fast (sonnet). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/lg-logo-colors
<that commit>`.

**Checks**: `pnpm test:unit packages/seq/`; `pnpm test:medium packages/seq/`; `npx tsc -b packages/seq/core
packages/seq/viewer packages/seq/examples`; `pnpm -C packages/seq/core typecheck:test`; `pnpm -C packages/seq/storybook
typecheck`; `pnpm lint`; `pnpm format:check`; `pnpm -C packages/seq/core test:perf`. Green only on exit 0 with a
non-zero test count.

**Headers**: a modified file keeps its header, year range to 2026, `@author David Sehnal <david.sehnal@gmail.com>` added
if missing; a new file has the Mol\* header with that author. No mention of `.scaffolds/seq/` in code; code cites its own `_spec`.

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`,
`biome.json`, `vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/` (write the entry's text in the report),
every `package.json` and `tsconfig*.json`, anything in `.scaffolds/seq/`.

**Stop when** green, or a check fails twice after a fix.

**Report**: commits; checks with exit codes and test counts; the stories to look at in light and dark; drift.
