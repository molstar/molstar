# Lane FX35: the novice pass's four design gaps, as decided

**Model**: Sonnet (four small decided changes). **Base**: `packages-seq` at the orchestrator's head, after FX34's merge
(`git reset --hard packages-seq` first). Branch `seq/fx35`. No lane beside you.

**Why**: [novice-findings.md](../novice-findings.md) (C1, C9, A13, B7 and B8, B11) and [gaps.md](../gaps.md), Open;
decided by the author, 2026-11-24: "Agree all four".

**The decisions** (write each into its design, its Storybook page and the spec in the same commit):

- **D1, opening a chromosome on the ring** (`ring-extension.md`, the interaction matrix's "a block" row and The ring
  zoom): a double-click on a chromosome's label or its arc (the hue ring, the cytobands) opens it: the cut, one undo
  step, as the label menu's "Show only"; a double-click on a ribbon or a stem keeps its meaning. Each ring story's text
  names the gestures: Cmd/Ctrl+wheel or a pinch over a chromosome, a double-click on it, the label's "Show only", the
  square, and on The Partner Keys story the `.` and `,` keys.
- **D2, a legend swatch's click** (`themes.md`, Legends): keeps the highlight (CL1's answer A); the swatch's tooltip
  says "Click to highlight · the eye hides the layer" (its accessible description too).
- **D3, every lane from the track's menu** (`track-header-parts.md`, the menu): on a track whose lanes scroll inside
  (its content past its cap), the menu offers "Show every lane"; on a track at `height: 'content'`, "Cap the lanes"
  (back to absent `height`). Each sets the track's `height` by one patch, one undo step. Not offered on a ring's track.
- **D4, a reference line's value** (`values-marks.md`, Decision 4): a line with a label draws "label · value"
  ("genome-wide · 7.3"), the value in the layer's number format; a line without a label draws the value alone. The
  accessible text unchanged ("reference line at 7.3 (genome-wide)"). The values y axis stays a TODO.

**Items**, a commit at each: 1 D1, 2 D2, 3 D3, 4 D4; each with tests failing on its parent (a gesture through the
viewer's model; the tooltip's text; the menu entries and one undo step; the label boxes and the SVG export), the
stories updated (the ring stories' texts; Manhattan shows the values) and their twins.

**How**: nothing per frame; host neutral; module size (new code lazy where the default viewer does without it).

**Where**: a worktree on `seq/fx35`. Headers per the root `AGENTS.md`. Hot files never staged: `package.json`s,
`tsconfig*.json`, the lockfile, `biome.json`, `.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`,
`packages/seq/_measure/sizes.mjs`, the work folder. Never a bare `git stash`. No network. Do not run the heavy lane.

**Checks**, green only on exit 0 with a non-zero test count: `pnpm test:unit packages/seq/ cli/seqx-render/`;
`pnpm test:medium packages/seq/`; `pnpm -s -C packages/seq/core test:perf` (load under 3: `uptime` first);
`pnpm -s -C packages/seq/core size` (every row; name a raise if one is needed); `npx tsc -b`; `pnpm lint`;
`pnpm format:check`; the scaffold checks in `.scaffolds/seq/AGENTS.md` (stage first).

**Stop when** green, or a check red twice after a fix, or a step is refused by the permission system (report it;
don't work around it).

**Report** (under 350 words): commits; each item with its proof; the checks with exit codes and counts; the size rows;
the stories to look at; the changelog text.
