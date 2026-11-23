# Lane FX19: genotype codes get tokens of their own

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/fx19`. Model:
Sonnet.

**Why** (the author, 2026-11-23, stage 3's gate: "Own tokens"; drift note `matrix-colours.md`): the `matrix` layer
colours genotype codes 0, 1, 2, 3 and 255 (no copy, some copies, every copy, partly called, uncalled) with tokens borrowed
from the value ramp (`seq-ramp-*`) and `seq-missing`, so a theme changing the ramp repaints genotypes.

**Task**: tokens `seq-genotype-0` to `seq-genotype-3` beside `seq-missing`, in the light and dark token sets, today's
colours exactly (no visual change); the matrix layer reads them; the core `_spec`'s matrix Defaults and token list name
the mapping; a test that every code's drawn colour equals its token in both themes, and that changing a ramp token
leaves genotypes unchanged (fails on the parent).

**Where**: a worktree on `seq/fx19`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot files never
staged (`package.json`s, `tsconfig*.json`, lockfile, `biome.json`, `.seq/`, `.v6/changelogs/`, `storybook/packages.ts`,
`packages/seq/TODO.md`, `packages/seq/_measure/sizes.mjs`). No language change. Anything in `@molstar/core` or
`@molstar/io` stops the lane.

**Checks**: `pnpm test:unit packages/seq/ cli/seqx-render/`, `npx tsc -b packages/seq/core packages/seq/viewer
packages/seq/examples`, `pnpm lint`, `pnpm format:check`, `pnpm test:medium packages/seq/`, the core perf tier.

**Report**: under 120 words: commits, checks, the token names and values, the changelog line.
