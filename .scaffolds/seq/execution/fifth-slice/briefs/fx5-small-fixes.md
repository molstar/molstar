# Lane FX5: six open items of fixes.md

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/fx-5`. Model:
Sonnet. Follow `fx-fix-lane.md` (one commit per item, a test failing on its parent, a look in light and dark). Each
item's full text, the author's words and its tests are in `fixes.md` (Open) under the heading named here: read it
there and do what it says.

**Beside you**: VC (`runtimes/vcf`) and RG (core's store, `data/`, `render/pairs*`, `grid*`, `layer-data.ts`, rings'
images, the matrix layer). Touch none of those files; if an item needs one, stop that item and say so.

**Items**:

1. **"The navigator's readout covers the control in use"**: the readout never over the window, its grips or the
   pointer during a drag.
2. **"The square navigator's sweep has no feedback while dragging"**: the sweep drawn live, snapped to cells; Escape
   cancels.
3. **"The palette control"**: all four bullets (no Color entry where the layer ignores a palette; one control, the
   menu; the submenu arrow at the right edge; a palette change a patch that undo and the original entry restore
   exactly).
4. **"The vertical overview cut at the right on a wide page"** (story layout).
5. **"The tree's leftmost vertical line is thinner"**.
6. **"Three element classes exported from stories files show as stories"** (move each out, empty `Known`).

Out of scope: the hover card, the empty well, the drag rule (class items, the strong model's), the navigator's window
and zoom strip (a look design first), anything in `@molstar/core` or `@molstar/io`, the language.

**Where**: a worktree on `seq/fx-5`. Headers updated (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot files
never staged (`.scaffolds/seq/`, `.v6/changelogs/`, `packages/seq/TODO.md`, `package.json`s, `tsconfig*.json`): list the edits.
Do not leave a storybook server running; stop any you start before you report.

**Checks** (all, with numbers): `pnpm test:unit packages/seq/` (under 2 s a touched folder), `npx tsc -b
packages/seq/core packages/seq/viewer packages/seq/examples`, `pnpm -s -C packages/seq/{core,viewer,examples}
typecheck:test`, `pnpm -C packages/seq/storybook typecheck`, `pnpm lint`, `pnpm format:check`, both perf tiers,
`pnpm -s -C packages/seq/core size`, the `.scaffolds/seq/_checks` (features-planned with `--plan
.scaffolds/seq/execution/fifth-slice/plan.md`). No `pnpm test:heavy` (another lane runs it); say which items touch a bound.

**Report**: under 250 words: commits per item; checks with numbers; the story ids for the orchestrator's look; items
stopped and why; the `.v6/changelogs/` text.
