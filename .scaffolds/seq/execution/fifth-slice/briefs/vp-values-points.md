# Lane VP: values as points, colour by value at full strength

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/vp`. Model:
Sonnet. **The language, alone**: the minor bump to `0.8.0` is this lane's (`SeqLanguageVersion`, its history entry, the
catalogs, the builder's types and docs); no other lane writes in `lang/` while you run. LA3 is in resolve, the store
and the grid layer; SX in `render/layout/` and the rail: stay out of them.

**Why** (the author, 2026-11-22): after the values look C merged (VL), "should we use the same look for the bars for
the plddt etc too?", then "do we have a simple line plot somewhere too? both single color and colored points", then
"Yes, one lane" to the orchestrator's proposal:

1. **Colour by value at full strength**: look C tints the series colour (bars at `BarAlpha` 60%, the profile's area
   at `FillAlpha` 35%). That suits one series colour. Where the colour comes from the value (`col`, `perBar`: a scheme
   or a ramp, pLDDT's four bands), the colour is the data: tinted, it no longer matches its legend and pLDDT's "very
   high" and "confident" blur. So a value's colour draws at full strength in bars, the profile's area nearly opaque
   (enough to read as an area under its line), a line and its points full; a uniform series keeps C's tints.
2. **A line shows its points when there is room**: from about 6 CSS px a value, `encoding: 'line'` and `'area'` draw a
   small dot (a circle, about 3 px across, the value's colour or the series colour, a 1 px ring in the well's ground)
   at each value; below that a plain line. No new param. A missing value has no dot.
3. **`encoding: 'points'`**: dots only (a scatter of values along the axis), one colour or by value; below a pixel a
   value, the device column's extent as a thin vertical mark (min to max), as the line follows each column's maximum
   today. Headroom, the baseline hairline and the clip to the well as C. Hit and hover as the line's (the nearest
   value in the column), the hover's emphasis a ring on the hovered dot.

**Language**: `points` in the `encoding` enum and its doc; the bump to `0.8.0` with its history line; the builder and
config-keys and catalog tests. Nothing else in the language.

**Stories**: p53's "Values as a line and an area" (`examples/src/views.ts`, `values-line`) gains a points track and a
pLDDT-coloured line with points; the Values looks story (`Features/Values looks`) shows uniform against by-value colour
at both densities; light and dark, canvas and SVG, headless twins on a fake clock.

**Tests**: by-value colours drawn equal their scheme's colours exactly (the legend's swatches), uniform ones keep C's
alphas; dots appear at 6 px a value and not below; a missing value has no dot and no mark; `points` hits the nearest
value; SVG export draws the same dots (circles) as the canvas. **Performance**: a warm frame of each encoding 0 B (an
allocation case each), points' hover 0 B warm.

**Where**: a worktree on `seq/vp`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); commit at each
green step; hot files never staged (`package.json`s, `tsconfig*.json`, lockfile, `biome.json`, `.scaffolds/seq/`,
`.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`): list them. The core `_spec`'s values sections
change in the same commits; design text for `layouts.md` ("The values look") listed for the orchestrator. Anything in
`@molstar/core` or `@molstar/io` stops the lane. A file past 800 lines is a question (`render/values.ts` is near it:
split the SVG or points code out if it grows past). **Sizes** soft (headless 940 at 937.9, default 1,495 at 1,481.3):
trim first, report growth with its cause, do not edit bounds. Timing and allocation bounds are gates.

**Checks**: the quick lane, `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`, each
`typecheck:test`, the storybook typecheck, `pnpm lint`, `pnpm format:check`, both perf tiers to completion,
`pnpm -s -C packages/seq/core size` (read the rows for FAIL), `pnpm test:medium packages/seq/`, the `.scaffolds/seq/_checks`
(features-planned with `--plan .scaffolds/seq/execution/fifth-slice/plan.md`, `defaults-designed.ts`). A fix's test fails on its
parent.

**Report**: under 350 words: commits; checks with numbers; story ids; sizes; the version and its history line; the
design text; hot-file needs; the `.v6/changelogs/` text.
