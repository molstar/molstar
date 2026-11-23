# Lane G2: stage 2's gate findings

**Goal**, three findings at stage 2's gate:
1. **The text layout's gutter labels overlap** (the look, "Examples/p53/Text › Text", canvas and SVG): the track labels
   of the thin bands ("Secondary structure", "Domains") are drawn on top of each other in the left gutter, because a
   band is a few px tall and its label 12 px. Lay the labels so none overlaps (the design: [layouts](../../../design/layouts.md),
   Text): e.g. the band labels stacked at their own line height in the first line's gutter, or only the letters'
   track labelled with the bands named in a legend-like list; pick what the design's text figure intends and say why.
   A test that no two gutter labels' boxes overlap, canvas and SVG.
2. **The domains count strip reads as a heavy dark rule** under nearly every line of the text figure: overlapping
   domains overflow the one-lane band into a count strip drawn in a dark ink. In a thin text band, draw the count
   strip in the density tokens' light end (as the tracks layout's density strip), so it reads as a band, not a rule;
   check light and dark. A test that the strip's color is the density token's.
3. **`render/_test/frame-time.test.ts`'s max bound flakes under load** ("reversed, 150,000 pieces › pan frame: never
   over 16 ms" failed in a full heavy run beside the build, passing 3/3 alone): the TODO item says best of three or a
   percentile. Use `itWithin` (best of three, the harness's rule) for its bounds, or a p99 instead of the max with the
   reason; remove the TODO item.

**Model**: fast (sonnet). **Base**: `packages-seq` at `52d69162f`: `git checkout -B seq/g2-gate-fixes 52d69162f`; check
`git merge-base --is-ancestor 52d69162f HEAD`.

**Touches**: `packages/seq/core/src/render/**` (text layout, the density strip, the frame-time test), the examples'
text records if they change, the core `_spec`, `packages/seq/TODO.md`.

**Checks**: `pnpm test:unit packages/seq/`; `pnpm test:heavy packages/seq/core/src/render/_test/frame-time.test.ts` (no
other heavy run: `pgrep -fl "vitest.*heavy"`); `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`;
`pnpm -C packages/seq/core typecheck:test`; `pnpm -C packages/seq/storybook typecheck`; `pnpm lint`; `pnpm format:check`.

**Headers**: a modified file keeps its header, year 2026, `@author David Sehnal <david.sehnal@gmail.com>` added if
missing. No mention of `packages/sequences` or `.seq/` in code. **Hot files**: never stage `package.json`,
`tsconfig*`, the lockfile, `.v6/changelogs/` (write an entry's text in the report if user-facing), `.seq/`. Fixes land
with a test shown failing on the parent (name the command).

**Report** (under 300 words): the commits, each check with its exit code and test count, what to look at (light and dark).
