# Lane DX: a double-click restores only its own first click's clear

**The bug** (the author, 2026-11-19: "double clicking the vertical slider selects the last selected column. is that by
design?"; not by design): `packages/seq/viewer/src/input.ts` keeps the selection a click cleared (`_cleared`, set in
`click`) so that a double-click whose first click cleared it gives it back (`_restoreCleared`, called by `doubleClick`
and by a click with `detail > 1`). `_cleared` is reset only on the select path of `click`; the rows rail's and the
ruler rail's click branches and a suppressed click (the click ending a drag) return before it. So: select a column,
click an empty cell (cleared, kept), then double-click the rows rail: the selection comes back. Any later double-click
reached without a select-path click does the same.

**Fix**: what a click cleared is kept for that click's double-click only: reset it on every click that is not a
double-click's second (every early return of `click` included), and on a pointer down that starts a drag. Keep the
intended case working (double-clicking the only selected element: the first click clears it, the double-click gives it
back and zooms; a Shift+double-click). Minimal code.

**Tests** (headless, the viewer's input tests): each path fails on the parent: a clear then a double-click on the
rows rail; on the ruler rail; after a drag; on a tree node; on the filter chip; plus the intended restore still passing.

**Not in scope**: anything else. Lanes FS-A1 and LG run beside you (data, sources, the logo): stay out of them.

**Model**: fast (sonnet). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/dx-double-click
<that commit>`.

**Checks**: `pnpm test:unit packages/seq/`; `pnpm test:medium packages/seq/`; `npx tsc -b packages/seq/viewer`;
`pnpm -C packages/seq/viewer typecheck:test`; `pnpm lint`; `pnpm format:check`. Green only on exit 0 with a non-zero
test count.

**Headers**: a modified file keeps its header, year range to 2026, `@author David Sehnal <david.sehnal@gmail.com>` added
if missing. No mention of `.seq/` in code.

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`,
`biome.json`, `vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/` (write the entry's text in the report),
every `package.json` and `tsconfig*.json`, anything in `.seq/`.

**Stop when** green, or a check fails twice after a fix.

**Report**: commits; checks with exit codes and test counts; the story to try it in (Globins › family) and how.
