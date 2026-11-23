# Lane GA3: inline data over an assembly, by name

**Base**: `packages-seq` at the orchestrator's head after FR merges (`git reset --hard packages-seq` first). Branch
`seq/ga3`. Model: strong. **The language, alone**: the minor bump to the next language version is this lane's (the
version, its history entry, the catalogs); no other lane writes in `lang/`, resolve or the store while you run.

**Why**: [genome-coordinates](../../../design/genome-coordinates.md) "Annotations over an assembly" (the author,
2026-11-22: "Ok we can use this naming"): inline rows of two ends or axes name their sequences in columns parallel to
maps' `from_name`/`to_name`: `a_name` and `b_name` for pairs, `x_name` and `y_name` for a grid's cells, `name` for a
matrix's items. GA2 read sources over an assembly and inline `ranges`/`values` by their `target` field; inline pairs,
grid cells and matrix items over an assembly have no column to say their sequence, so they are refused or unread today
(`packages/seq/TODO.md`, "Region reads' follow-ups": "inline data over an assembly (GA3)").

**Task**:

1. **The language**: `a_name`/`b_name` on inline pairs rows, `x_name`/`y_name` on inline grid cells, `name` on inline
   matrix items; each required when its side's `target` (`a`/`b`, `x`/`y`, `target`) names an assembly and refused
   (a check issue, catalogued) when it names a sequence; unknown names counted and said once
   (`seq.resolve.assembly.name`, as GA2 says it). The builder's entries and docs beside `from_name`/`to_name`. A
   matrix's `source` typed as the other kinds' (FR removed the cast; check what is left untyped in `language.ts` and
   the builder). One minor bump with its history line; the config-keys and catalog tests.
2. **Resolve and the store**: inline rows split by those columns through the names rule into per-sequence (or
   per-block-pair) rows, the same access plan GA2 gives a source; a pair between two sequences of one assembly is a
   chord between blocks; a property: inline rows over an assembly draw as the same rows given per sequence, joined
   (generated assemblies with aliases, `include`, `order`).
3. **Held pairs over blocks by the pairs level** (region-reads.md Decided 7, GA2's partial): held pairs over an
   assembly past the gate bin by the pairs level, not by device column, so held and BEDPE pairs over a genome are one
   picture (a test that fails on its parent).
4. **Stories**: `Features/Genome Coordinates` gains inline links over the toy assembly (a chord ring written in the
   spec, `a_name`/`b_name`) and an inline matrix over two chromosomes; light and dark, canvas and SVG, headless twins on
   a fake clock. The builder showcase (`Features/Builder`) shows one inline pairs row with `a_name`/`b_name`.

Out of scope: the whole-assembly level cut into block pairs (region-reads.md Default 5: list it if you see the way);
sided pairs on side `b` (LA's); the matrix colour tokens (an open drift, the author's).

**Where**: a worktree on `seq/ga3`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); commit at each
green step; hot files never staged (`package.json`s, `tsconfig*.json`, lockfile, `biome.json`, `.scaffolds/seq/`,
`.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`): list them. Read `packages/seq/core/_spec/README.md`
(Language, Resolve) and GC1's bump for the pattern. Anything in `@molstar/core` or `@molstar/io` stops the lane. A file
past 800 lines is a question, past 1500 a failure. Sizes (reading 260, headless 920 — at 917.2 —, default 1,460,
minimal 1,040 kB): the language's growth is the risk; keep docs and catalogs out of the reading path; a bound passed
stops the lane with a report (sizes are the author's).

**Checks**: the quick lane (under 2 s a touched folder), `npx tsc -b packages/seq/core packages/seq/viewer
packages/seq/examples`, each `typecheck:test`, the storybook typecheck, `pnpm lint`, `pnpm format:check`, both perf
tiers, `pnpm -s -C packages/seq/core size`, `pnpm test:medium packages/seq/`, `pnpm test:heavy packages/seq/` when load
is under 6, the `.scaffolds/seq/_checks` (features-planned with `--plan .scaffolds/seq/execution/fifth-slice/plan.md`,
`defaults-designed.ts`). A fix's test fails on its parent.

**Report**: under 350 words: commits; checks with numbers; story ids; sizes; the public surface added (names, against
the design); the version and its history line; what is left; hot-file needs; the `.v6/changelogs/` text.
