# Lane LS: the logo's glyph colors saturated

**Goal**: [grids-and-rows](../../../design/grids-and-rows.md), The logo, "And saturated as a glyph" (the author,
2026-11-19). LG (merged) moves each logo letter's class color in oklab lightness to 3:1 against `seq-track`
(`render/logo.ts`, `fillOf`, through `inkAgainst` in `theme/compile.ts`). In light the pale class fills come out
greyish and muddy. Now: raise the glyph color's oklab chroma, hue kept, to a floor clipped to the sRGB gamut at its
lightness, then apply the 3:1 lightness rule; a color already above the floor keeps its chroma; a grey (chroma near
0: a gap, `missing`, a neutral) stays grey. The floor is a new token per theme, `seq-logo-chroma` (start near 0.12 in
light and 0.10 in dark; tune by eye on Globins › Rows › Family and by the test). Reuse `@molstar/core/util/color`'s
oklab code; no private copies. Memoized as LG's map is (per color and token version).

**Tests**: every built-in protein and nucleotide scheme, both themes: each glyph color at 3:1 or more against the
surface; chroma at or above the floor or at the gamut's limit; hue within a stated tolerance of the class color's;
greys unchanged in chroma. LG's tests keep passing. The logo's allocation guard unchanged.

**Not in scope**: anything else. HB (hover checks, may touch `logo.ts`'s hover pass) and FS-A1 run beside you.

**Model**: fast (sonnet). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/ls-logo-saturation
<that commit>`.

**Checks**: `pnpm test:unit packages/seq/`; `pnpm test:medium packages/seq/`; `npx tsc -b packages/seq/core
packages/seq/viewer packages/seq/examples`; `pnpm -C packages/seq/core typecheck:test`; `pnpm -C packages/seq/storybook
typecheck`; `pnpm lint`; `pnpm format:check`; `pnpm -C packages/seq/core test:perf`. Exit codes explicit
(`; echo exit=$?`), non-zero test counts.

**Headers**: a modified file keeps its header, year range to 2026, `@author David Sehnal <david.sehnal@gmail.com>` added
if missing. No mention of `.scaffolds/seq/` in code; the core `_spec`'s "The logo" updated.

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`,
`biome.json`, `vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/` (write the entry's text in the report),
every `package.json` and `tsconfig*.json`, anything in `.scaffolds/seq/`.

**Stop when** green, or a check fails twice after a fix.

**Report**: commits; checks; the floor values chosen and why; the story to look at.
