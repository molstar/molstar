# Lane FX11: a track's height while loading

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/fx-11`. Model:
Sonnet (a decided design). Follow `fx-fix-lane.md` (a commit per part, a test failing on its parent, a look in light
and dark).

**Why**: the author, 2026-11-22 (Features › Genome formats › GFF3 › TP53 while loading): "another issue is the size when
loading. could be handled better": the Transcripts track shows one lane with "TP53-227" drawn below its band, then
jumps to 240 px. layouts.md's Track sizing already decides it; FX9 left three parts unbuilt:

1. **Labels stay in the band** (layouts.md, "Loading: the labels stay in the band"): a pending track's placeholder lane
   and any band clip labels and leaders to the band's box; nothing is drawn below a track.
2. **Growth is batched to the settle**: a track grows at most once per settled read (a store `settle`, data-store.md),
   so tracks below move once, never tile by tile; today a global version counter relays out per arrival.
3. **The content from the reachable region**: the lanes packed over the data read for the section's reachable region
   (its tiles at the fit's level, as the navigator reads), so the first settle gives the height a pan will keep; never
   shrunk by a pan or a zoom; reset when the data's identity changes.

**Beside you**: FX10 writes the rails (`render/rail.ts`, `viewer/src/input.ts`, `layout/text.ts`) and a tracks
section's scroll (the section viewport's `top`, `section.maxPx`); LA writes the bridge (`layout/bridge.ts`,
`pairs-ribbons.ts`, the viewer's lock). Your files are the ranges layer's sizing (`ranges.ts`, `known.ts`), the store's
settle signal as the viewer and render read it, and `tracks.ts` only where growth is applied: keep `tracks.ts` changes
small and list them.

**Tests**: a headless twin on a fake clock for TP53's GFF3 by region: tiles arriving one by one move the tracks below
exactly once per settle (count layouts that change a track's height); the first settled height equals the height after
a pan across the reachable region; a placeholder with long labels draws nothing outside its band (canvas and SVG);
EC-169 in a title. Each fails on its parent. Both perf tiers (a warm frame while tiles arrive allocates nothing new).

**Where**: a worktree on `seq/fx-11`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot files never
staged (`.seq/`, `.v6/changelogs/`, `packages/seq/TODO.md`, `package.json`s, `tsconfig*.json`): list edits. Stop any
Storybook server you start and close any browser tab you open.

**Checks** (with numbers): the quick lane (under 2 s a touched folder), `npx tsc -b packages/seq/core
packages/seq/viewer packages/seq/examples`, each `typecheck:test`, the storybook typecheck, `pnpm lint`,
`pnpm format:check`, both perf tiers, `pnpm -s -C packages/seq/core size` (headless 923.3 of 930, default 1,464.2 of
1,470, minimal 1,039.3 of 1,040 kB: a bound passed stops the lane with a report), the `.seq/_checks`
(features-planned with `--plan .seq/execution/fifth-slice/plan.md`).

**Report**: under 250 words: commits, the tests that failed on the parent, sizes, perf, the look (TP53 GFF3 and GTF
loading, light and dark), the `.v6/changelogs/` text.
