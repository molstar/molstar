# Lane FX6: the two decided drifts

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/fx-6`. Model:
Sonnet. Follow `fx-fix-lane.md` (one commit per item, a test failing on its parent, a look in light and dark). Each
item's full text and its tests are in `fixes.md` (Open) under the heading named here; the designs say what is decided
(viewer-composition.md on the export; genome-formats.md, Defaults 2, on a mixed run).

**Beside you**: RG (core's store, `data/`, `render/pairs*`, `grid*`, `layer-data.ts`, rings' images, the matrix
layer) and GB (`runtimes/genbank`, the services story's dropped file). Touch none of their files; if an item needs
one, stop that item and say so.

**Items**:

1. **"The SVG export blank on a spec with an error"**: export what resolves, as the canvas draws it, the issues in
   the picture's `<desc>`. Update `Features/Builder`'s Check story text and twin to the new behaviour.
2. **"A merged run of mixed categories drawn neutral"**: the majority fill, ties by the legend's order; the hover over
   a mixed run lists its categories. Check `seq-next-features-genome-formats-vcf--tp-53-clin-var` reads by category.

Out of scope: the hover card's layout (a class item), anything in `@molstar/core` or `@molstar/io`, the language.

**Where**: a worktree on `seq/fx-6`. Headers updated (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot files
never staged (`.scaffolds/seq/`, `.v6/changelogs/`, `packages/seq/TODO.md`, `package.json`s, `tsconfig*.json`): list the edits.
Stop any Storybook server you start before you report.

**Checks** (all, with numbers): `pnpm test:unit packages/seq/` (under 2 s a touched folder), `npx tsc -b
packages/seq/core packages/seq/viewer packages/seq/examples`, each `typecheck:test`, the storybook typecheck,
`pnpm lint`, `pnpm format:check`, both perf tiers (the ranges layer is on a frame's path: a warm frame allocates
nothing), `pnpm -s -C packages/seq/core size`, the `.scaffolds/seq/_checks` (features-planned with `--plan
.scaffolds/seq/execution/fifth-slice/plan.md`). No `pnpm test:heavy`.

**Report**: under 200 words: commits; checks with numbers; story ids; the `.v6/changelogs/` text.
