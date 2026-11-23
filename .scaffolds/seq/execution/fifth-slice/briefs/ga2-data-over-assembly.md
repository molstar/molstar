# Lane GA2: data over an assembly, read by sequence

**Base**: `packages-seq` at the orchestrator's head after RG merges (`git reset --hard packages-seq` first). Branch
`seq/ga2`. Model: strong (resolve and the store: alone in them). Nothing else writes in core while you run.

**Why**: [genome-coordinates](../../../design/genome-coordinates.md) ("Annotations over an assembly") says `target`
(ranges, values, pairs, grid) may name an assembly, each row on the sequence its `target` field names through the names
rule. GC2 checked it (`resolve/genome-model.ts` `checkAssemblyData`: a mismatched source, unknown names said once) but
never read it: `resolve/unread.ts:32` says "A … over an assembly" as unread, so every data node over an assembly draws
nothing. The design's "Not built" line lists it. It blocks the genome views the slice is for: an ideogram's bands over
every chromosome (GC3's follow-up in `packages/seq/TODO.md`), a BED, VCF or bigWig over a genome, the genome ring's
bin chords and Hi-C over the genome (region-reads.md Defaults 3 and 5), and GX's Genomes dataset next. RS found it as
gap 1 (`packages/seq/core/src/render/_test/region-gaps.test.ts`, an `it.fails` test).

**Task**:

- **Resolve**: a `ranges`, `values`, `pairs`, `grid` or `matrix` node whose `target` (or `a`/`b`, `x`/`y`) names an
  assembly resolves to an access plan by sequence: inline rows split by their `target` field through the names rule
  (unknown names counted and said once, as `checkAssemblyData` already says); a source read by region per sequence of
  the assembly (each want on the block's sequence, the source's names mapped by the names rule both ways, `17` ↔
  `chr17`); `include` and `order` respected (a sequence left out reads nothing). Remove the unread info for these kinds.
- **The store and render**: every layer over such data draws on each block of the section (rows, a ring of
  chromosomes, the square's blocks × blocks), reads per block bounded by the per-layer bounds (`maxTilesPerLayer`,
  `maxBytesPerLayer` summed over blocks), small blocks by the joint small-blocks rule (genome-coordinates.md and
  region-reads.md, one text); pairs between two sequences of one assembly draw as chords between blocks (the genome
  ring); a grid over the assembly × itself its block pairs (Hi-C over the genome through the whole-assembly level).
- **Held pairs by the pairs level** (RS's gap 6, region-reads.md Decided 7): held pairs past the gate bin by the same
  pairs level as region pairs (one reader), so held and BEDPE pairs are one picture; its `it.fails` test turned on.
- **The ideogram's bands over every chromosome** (GC3's follow-up): the navigator over an assembly reads its BED bands
  by sequence.
- **Stories**: move RS's three single-system stories onto assemblies (`Features/Region reads`: the genome ring at the
  fit, Hi-C over the genome, synteny on two assemblies) and add `Features/Genome Coordinates` stories: a BED over the
  toy assembly's rows, and the cytobands over every chromosome as the ideogram navigator. Light and dark, canvas and
  SVG, headless twins on a fake clock.
- **Tests**: the `it.fails` tests for gaps 1 and 6 pass; a region read over an assembly equals the per-sequence reads
  joined (a property over generated assemblies with aliases, `include` and `order`); a row naming no sequence counted
  once; EC-187, EC-189 and EC-191 in titles where they apply; the per-layer bounds summed over blocks; a warm frame over
  25 blocks allocating nothing (both perf tiers).

Out of scope: BEDPE through the store (RS's gap 2, FR's); the language (no new param: `target` naming an assembly is
already declared); sided pairs on side `b` (LA's).

**Where**: a worktree on `seq/ga2`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); commit at each
green step; hot files never staged (`package.json`s, `tsconfig*.json`, lockfile, `.scaffolds/seq/`, `.v6/changelogs/`,
`storybook/packages.ts`, `packages/seq/TODO.md`): list them. Anything in `@molstar/core` or `@molstar/io` stops the
lane. A file past 800 lines is a question; past 1500 a failure. Sizes: reading 260, headless 920, default 1,460,
minimal 1,040 kB; a bound passed is a question.

**Checks**: the quick lane (under 2 s a touched folder), `npx tsc -b packages/seq/core packages/seq/viewer
packages/seq/examples`, each `typecheck:test`, the storybook typecheck, `pnpm lint`, `pnpm format:check`, both perf
tiers, `pnpm -s -C packages/seq/core size`, `pnpm test:medium packages/seq/`, `pnpm test:heavy packages/seq/` when
load is under 6, the `.scaffolds/seq/_checks` (features-planned with `--plan .scaffolds/seq/execution/fifth-slice/plan.md`). A fix's test
fails on its parent.

**Report**: under 350 words: commits; checks with numbers; story ids; sizes; the public surface; what is left; hot-file
needs; the `.v6/changelogs/` text.
