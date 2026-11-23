# Lane RL: the rail's look, variants to pick from

**Base**: `packages-seq` at the orchestrator's head after FX10 merges (`git reset --hard packages-seq` first). Branch
`seq/rl`. Model: Sonnet (a story of variants; the author picks; a second step builds the pick). Follow `fx-fix-lane.md`
for commits and looks.

**Why**: the author, 2026-11-22: "should we iterate on the look of the rail? i think the contrast breaks the view quite
a bit + the ticks breaks the 'rounded look' of the rest of the ui". Today (`core/src/render/rail.ts`, the window pass):
a dark line in `msui-line`, the window's segment, square tick grips in `seq-grip` (`GripRise` 3), the idle hairline;
at the fit a full-width dark line with end ticks; the grips overhang the well's rounded corners. One rail serves the
ruler, rows, lanes and text (layouts.md, A slim rail: one rail for every box that scrolls).

**Step 1, variants (no default changes)**: a render config `navigator.rail.look` (internal, for the story only; not
in the `_spec` tables until picked) with three looks beside today's, drawn by the one rail path:

- **A, capsule**: the window a rounded pill (round caps, radius half its thickness) in a soft fill (the accent at low
  alpha, or `msui-line` at rest); its ends are the grips (no ticks; a resize cursor near them; they grow on hover);
  no groove.
- **B, capsule in a groove**: A over a faint rounded groove the axis's length (a native scroll track), the groove
  barely visible, the window the thumb.
- **C, quiet until hovered**: B at rest in low contrast and thin; on hover or drag the full ink and `grownPx`; at the
  fit only the groove.

All three: ends inset by the well's corner radius (never past its rounded corner); one look for ruler, rows, lanes and
text rails (grips only where the box zooms). A story `Design/Rail looks` (or `Features/Navigators/Rail looks`): each
variant and today's side by side on a ruler rail (one zoomed, one at the fit), an alignment's rows rail and a lanes
band (GFF3 TP53), light and dark, with a sentence each; a headless twin. Hover and drag work in each.

**Step 2**: stop and report the story id; the author picks. (The orchestrator sends the pick back: build it as the
default, delete the others, update layouts.md and navigators.md's rail lines and the `_spec`, the tokens if new.)

**Where**: a worktree on `seq/rl`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot files never
staged (`.seq/`, `.v6/changelogs/`, `packages/seq/TODO.md`, `package.json`s, `tsconfig*.json`): list edits. Stop any
Storybook server you start and close any browser tab you open.

**Checks**: the quick lane, `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`, each
`typecheck:test`, the storybook typecheck, `pnpm lint`, `pnpm format:check`, both perf tiers (a warm rail frame 0 B in
every look), `pnpm -s -C packages/seq/core size` (variants cost bytes: report them; the pick's are what stays).

**Report**: under 200 words: the story id, a sentence per variant, sizes.
