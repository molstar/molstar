# Lane CR: the navigators' cursors by part

**Goal**: [navigators](../../../design/navigators.md), "The cursor says the gesture" (the author, 2026-11-19: "when
the action is resize on the sliders, it should get an appropriate icon"). Today a grip shows `grab` like the body
(`ViewerInput.cursor` in `packages/seq/viewer/src/input.ts`, `NavigatorModel.cursor` in `navigator.ts`). Make the
cursor follow the part: a grip (`start`, `end`) `ew-resize` on a horizontal navigator and `ns-resize` on a vertical one
(the rows rail, a vertical slider), held through the grip's drag wherever the pointer goes; the body `grab`, `grabbing`
while dragged; outside `pointer` (a click centres; the view's rails as the mounted navigators do); idle the default.
One rule for the ruler rail, the slim rail, the rows rail, the strip and the miniature. Minimal code; the `_spec`'s
navigator section says the rule.

**Tests** (fail on the parent): per navigator kind and orientation, the cursor over each part and during each drag,
headless through the models; an idle navigator's default.

**Not in scope**: anything else. LS (logo colors) and FS-A1 run beside you.

**Model**: fast (sonnet). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/cr-cursors
<that commit>`.

**Checks**: `pnpm test:unit packages/seq/`; `pnpm test:medium packages/seq/`; `npx tsc -b packages/seq/viewer`;
`pnpm -C packages/seq/viewer typecheck:test`; `pnpm -C packages/seq/storybook typecheck`; `pnpm lint`;
`pnpm format:check`. Exit codes explicit (`; echo exit=$?`), non-zero test counts.

**Headers**: a modified file keeps its header, year range to 2026, `@author David Sehnal <david.sehnal@gmail.com>` added
if missing. No mention of `.seq/` in code.

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`,
`biome.json`, `vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/` (write the entry's text in the report),
every `package.json` and `tsconfig*.json`, anything in `.seq/`.

**Stop when** green, or a check fails twice after a fix.

**Report**: commits; checks; the stories to try (Complex › Blocks › Split's rail, Globins › Rows › Family's rows rail,
a mounted navigator).
