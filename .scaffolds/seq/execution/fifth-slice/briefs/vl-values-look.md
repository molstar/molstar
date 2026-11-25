# Lane VL: the values layer's look, variants to pick from

**Base**: `packages-seq` at the orchestrator's head after RL merges (`git reset --hard packages-seq` first). Branch
`seq/vl`. Model: Sonnet (a story of variants; the author picks; a second step builds the pick). Follow
`fx-fix-lane.md` for commits and looks.

**Why**: the author, 2026-11-22 (p53's Conservation, zoomed out): "is there a way to make the 'barcharts' look nicer
too?" At a bar or two of pixels each, the values layer draws hard columns of full ink: a jagged block, not a profile;
the tallest bars touch the well's top; a missing value reads as a grey stub at the floor; the fill's ends meet the
well's rounded corners square.

**Step 1, variants (no default changes)**: a render config key for the values layer's look (internal, for the story;
not in the `_spec` tables until picked), drawn by the layer's one path, beside today's:

- **A, profile**: below about 4 px a value, one stepped area: the fill the series colour at about 50%, a crisp top
  line in the full colour; at 4 px or more, bars.
- **B, profile and bars**: A, and at 4 px or more separate bars with a 1 px gap and slightly rounded tops.
- **C, B quieter**: B with a lighter fill (about 35%) and the top line one step thinner.

All variants: about 10% headroom at the top; a missing value a gap, never a stub at zero; the fill clipped to the
well's rounded corners; a hairline baseline at zero. A story (`Design/Values looks`, beside `Design/Rail looks` or its
successor): each variant and today's on p53's Conservation and a bigWig signal, each zoomed out (a value under a pixel
or two) and zoomed in (bars of 6–10 px), light and dark, a sentence each; a headless twin. Hover works in each (the
hovered value's bar or step lit).

**Step 2**: stop and report the story id; the author picks. (The orchestrator sends the pick back: build it as the
default, delete the others, update the designs' values lines through the report and the `_spec`.)

**Where**: a worktree on `seq/vl`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot files never
staged: list edits. Stop any Storybook server you start and close any browser tab you open.

**Checks**: the quick lane, `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`, each
`typecheck:test`, the storybook typecheck, `pnpm lint`, `pnpm format:check`, both perf tiers (a warm values frame 0 B
in every look), `pnpm -s -C packages/seq/core size` (variants cost bytes: report them).

**Report**: under 200 words: the story id, a sentence per variant, sizes.
