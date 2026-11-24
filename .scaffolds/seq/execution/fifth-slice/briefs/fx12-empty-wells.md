# Lane FX12: no track drawn as an empty well, the ring's links, a timing bound

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/fx-12`. Model:
Sonnet (decided items). Follow `fx-fix-lane.md` (a commit per part, a test failing on its parent, a look in light and
dark).

**Items** (`.scaffolds/seq/execution/fifth-slice/fixes.md`, Open: read each entry whole):

1. **No track drawn as an empty well** (the author, 2026-11-21: "Yes, group them"): the rule as the entry states it: a
   track with no items in what it covers never draws a full-height blank well; by default it collapses to a one-line
   row naming it with a quiet "no items" (`track.empty: 'note' | 'collapse' | 'keep'`, default `note`); pending is not
   empty (EC-169: the placeholder stays); FX11's settle decides "empty". Then its three instances:
   - **A stacked block's tracks** (`arrange: 'stack'`, Complex › Blocks › Stacked, the ATP chip): a block's track with
     no items for that block is left out of that block (`blocks.emptyTracks: 'collapse' | 'keep'`, default
     `collapse`; the letters track kept), and the block gets its heading like the others ("ATP ×2").
   - **The plasmid's strip** (`examples/src/navigators.ts:218`, `plasmidStrip`): the story's strip gets content and
     sits next to the detail it navigates; the library: a strip with no `content` draws its system's ticks, not a
     blank well.
   - **The Promoters track** (`files/plasmid-gff.ts`): drop the empty track and add a Sites track for the file's 14
     `restriction_enzyme_recognition_site` rows (no fixture re-record needed if the rows are there); a story test that
     every declared track of the story has at least one item.
   - A sweep: a test over the stories' headless twins that no settled track draws a blank well of its full height.
2. **The genome ring's links** (GA2's look): the generator's links structured (a few strong inter-chromosome blocks and
   local links, as synteny or SVs), deterministic on core's `PCG`, so `Features › Region reads › Pairs › On a ring`
   reads as structure, not a disc; and find what draws the row of stray dots above the ring at the top left, and fix
   it.
3. **The 3,000-contig drag bound** (FX8's note in the drag-rule entry): `navigator-blocks.quick.ts`'s wall-clock 2 ms
   bound into the heavy lane under `itWithin`; its allocation bound stays in the quick tier.

**Beside you**: RL (rails, track header actions, hover cards on rails: `rail.ts`, `rows-rail.ts`, header hover) and LA
(bridges, the viewer's input for bridges). Keep out of their files; list any touch. Design text for `.scaffolds/seq/` designs
(ux-foundations.md, composition.md, navigators.md) goes in your report for the orchestrator.

**Where**: a worktree on `seq/fx-12`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot files never
staged (`.scaffolds/seq/`, `.v6/changelogs/`, `packages/seq/TODO.md`, `package.json`s, `tsconfig*.json`): list edits. Stop any
Storybook server you start and close any browser tab you open.

**Checks** (with numbers): the quick lane, `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`,
each `typecheck:test`, the storybook typecheck, `pnpm lint`, `pnpm format:check`, both perf tiers, `pnpm -s -C
packages/seq/core size` (headless 930, default 1,480, minimal 1,060 kB), the `.scaffolds/seq/_checks` (features-planned with
`--plan .scaffolds/seq/execution/fifth-slice/plan.md`).

**Report**: under 300 words: commits, the tests that failed on the parent, sizes, the look (Complex › Blocks ›
Stacked, Plasmid › Ring › Detail with a strip, Files › Plasmid GFF3, Pairs › On a ring), the design text, the
`.v6/changelogs/` text.
