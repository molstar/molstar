# Lane FX9: track sizing, as decided, and the tree's way back

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/fx-9`. Model:
Sonnet (a decided design). Follow `fx-fix-lane.md` (a commit per part, a test failing on its parent, a look in light
and dark).

**Build** layouts.md's section "Track sizing: content height, a cap, then its own scroll" (amended 2026-11-22, all its
Decided items yes): a ranges track sized to its lanes as known (growing as rows arrive, never a blank band, one loading
lane while pending); `track.maxPx` 240; past it its lanes scroll inside the track (wheel, keys, a slim rail like the
rows', header and labels fixed); a label with no room floats at the data's left edge; an export draws every lane up to
4,000 px (`export.lanes: 'view'` for the crop); a ring compacts past `rings.trackMaxPx`, no inner scroll; the rows band
without a `height` takes the same rule with a 480 px cap; the five-lane content cap retired. Close `fixes.md`'s
"Sizing a track of many lanes: content height, a cap, then its own scroll", and "A tree opened from its toggle has no
way back on screen" (FX7 stopped it on RG's files, now merged: the `LayoutContext` field in `layer-type.ts`, the line
in `tracks.ts`, and the menu entry "Show tree / Fold tree").

**Beside you**: FR in `core/src/runtimes/*` and the matrix layer's legend; keep out of both.

**Tests**: the section's Defaults as twin tests (the GFF TP53 transcripts sized to content; a 40-lane track capped at
240 px with its own scroll, the header fixed while lanes scroll; an export with every lane; a ring compacting; the
rows band capped); EC-new-1 to EC-new-5 from the section, named as the catalog numbers them (ask me if unnumbered);
both perf tiers (a warm frame of a scrolled track allocates nothing); the Configuration rows in the `_spec` tables.

**Where**: a worktree on `seq/fx-9`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot files never
staged (`.scaffolds/seq/`, `.v6/changelogs/`, `packages/seq/TODO.md`, `package.json`s, `tsconfig*.json`): list edits. Stop any
Storybook server you start and close any browser tab you open.

**Checks** (with numbers): the quick lane (under 2 s a touched folder), `npx tsc -b packages/seq/core
packages/seq/viewer packages/seq/examples`, each `typecheck:test`, the storybook typecheck, `pnpm lint`,
`pnpm format:check`, both perf tiers, `pnpm -s -C packages/seq/core size` (headless 917.2 of 920 kB: a bound passed is
a question), the `.scaffolds/seq/_checks`. No `pnpm test:heavy`.

**Report**: under 250 words: commits; checks; story ids for the look (GFF3 TP53, the Globins rows family, a ring, an
export); what was not built and why; the `.v6/changelogs/` text.
