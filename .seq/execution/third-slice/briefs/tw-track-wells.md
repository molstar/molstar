# Lane TW: track padding and squircle wells

**Goal** (the author, 2026-11-19: "Yes, as proposed"), as [ux-foundations](../../../design/ux-foundations.md) §3
(padding) and §5 (the well a surface) and [language-description](../../../design/language-description.md)'s `track`
row say:

1. **`padding` on a track**: px inside its well, vertical only; a number is top and bottom, `(top n)`, `(bottom n)` or
   `(top n bottom n)` one or both sides; default 0. The well's height is the track's height plus its padding; the
   data's band (what every layer lays out, hits and emphasizes in) is inset by it. Rows, pairs, grids, text bands and
   rings honour it (a ring's padding is radial). A render config default (`trackPadding`, top and bottom) for hosts,
   0. Language 0.4.0 (an addition, recorded in its Versions), the builder and the language description updated;
   validation: non-negative numbers, an issue otherwise.
2. **Squircle wells**: `drawBackground`'s well (`render/tracks.ts`, `seq-track`) drawn as the draw lists' squircle
   (radius by the config's squircle rule over the well's height, capped by `squircleMax`; square below `squircleMin`),
   and each track's data clipped to its well with the same path (one clip a track, not a mark; the frame-time and
   allocation tiers hold). A config `trackShape: 'squircle' | 'square'`, default `squircle`. Canvas and SVG draw the
   same curve (the clip a `clipPath` in SVG). Emphasis drawn per track (a selection box a track) follows the well's
   shape; tints under the data stay inside the clip.

3. **Every emphasis a squircle** (the author, 2026-11-19: "Yes, into TW", then "Yes, round all"; ux-foundations §5,
   "Every emphasis is a squircle"): as the selection's tint and outline already are (`Paint[3] = 1` in
   `render/tracks.ts`): the projected hover (`Paint[3] = 0` today), the highlight stripes in every style (`tint`,
   `band`, `outline`: 0 today), the cell hover of letters, logo and values (`render/marks.ts` `hoverCell`, radius 0)
   and a chip's hover outline, each round where its span really ends, square where cut, square below `squircleMin`,
   an outline's radius the mark's plus its offset (concentric). Tests: each of these in canvas and SVG has the radius,
   a cut end stays square, a cell narrower than the minimum stays square.

**Tests** (fail on the parent): padding moves the data's band and not the well's top; a hit at the padding is no
item; a top-only padding; a ring's radial padding; validation issues; the well's path in canvas and SVG at several
heights (square below the minimum); the clip present for a letters track and a grid track; `trackShape: square`
draws today's rect; the SVG export unchanged except for the path and clip (the twins' recorded SVGs updated). A story
under Features (Themes or Text, wherever tracks are shown plainly) with a padded top, and the Complex Ligands view
padded on top so its arcs clear the section's label.

**Not in scope**: horizontal padding (never), group gaps (§3's 8px between groups), anything else.

**Model**: strong (opus). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/tw-track-wells
<that commit>`. NB (the numbering kind) and CY (source axes) may run beside you; DC's Complex views may still be
merging: if `examples/src/complex/ligands.ts` is not on your base, skip its padding and say so.

**Checks**: `pnpm test:unit packages/seq/`; `pnpm test:medium packages/seq/`; `npx tsc -b packages/seq/core
packages/seq/viewer packages/seq/examples`; `pnpm -C packages/seq/<core|viewer|examples> typecheck:test`; `pnpm -C
packages/seq/storybook typecheck`; `pnpm -C packages/seq/core test:perf`; `pnpm lint`; `pnpm format:check`; `pnpm -C
packages/seq/core size` before and after. Exit codes explicit, non-zero test counts. Look at it in Storybook in light
and dark.

**Headers**: a modified file keeps its header, year range to 2026, `@author David Sehnal <david.sehnal@gmail.com>` added
if missing; a new file has the Mol\* header with that author. No mention of `.seq/` in code; code cites its `_spec`.

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`,
`biome.json`, `vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/` (the entry's text in the report), every
`package.json` and `tsconfig*.json`, anything in `.seq/`.

**Stop when** green, or a check fails twice after a fix, or the work needs public API past the params above
(`padding`, `trackPadding`, `trackShape`), or a permission is refused.

**Report**: commits; checks; sizes; what you saw; drift; followups.
