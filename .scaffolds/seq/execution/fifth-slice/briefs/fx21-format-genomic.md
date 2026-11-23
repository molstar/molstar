# Lane FX21: `formatGenomic` parses each spec once (@molstar/core)

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/fx21`. Model:
Sonnet. **Scope**: `@molstar/core`'s `util/string-format` (the author approved this change, 2026-11-23, stage 3's gate)
and seq's ruler label cache.

**Why** (PF, fifth slice): `formatGenomic` goes through `formatValue`, which parses its format spec (a regex) on every
call: 13.6 KB a frame on the chromosome pan before seq's ruler cached its labels.

**Task**: in `packages/core/src/util/string-format.ts` (read it and its tests): parse each spec `formatGenomic` uses
(`,.Nf`, N 0 to 6) once and pass the parsed spec to the formatter (a parsed-spec overload or a small cache keyed by
the spec string, bounded); output byte-equal to today for every input (a property test over random numbers and every
N); an allocation case or a count showing no parse per call. Then decide whether seq's ruler label cache
(`packages/seq/core/src/render/ruler.ts`, PF's) is still worth its code: measure the chromosome pan's bytes a frame with
and without it; keep it only if it still saves materially, and say the numbers.

**Where**: a worktree on `seq/fx21`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026; the core file's
range extended); hot files never staged (`package.json`s, `tsconfig*.json`, lockfile, `biome.json`, `.seq/`,
`.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`, `packages/seq/_measure/sizes.mjs`). No public API
removed in `@molstar/core` (an addition only, and only if needed); other callers of `formatValue` unchanged.

**Checks**: `pnpm test:unit packages/core/src/util packages/seq/ cli/seqx-render/`, `npx tsc -b`, `pnpm lint`,
`pnpm format:check`, `pnpm check:workspace`, the seq storybook perf tier (the chromosome cases).

**Report**: under 150 words: the change, commits, the readings with and without seq's cache, the decision, the changelog
line.
