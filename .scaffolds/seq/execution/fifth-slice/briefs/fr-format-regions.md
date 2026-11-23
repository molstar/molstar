# Lane FR: the formats by region, past level 1

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/fr`. Model:
strong (untrusted input, the store's protocol). **GA2 runs beside you alone in resolve and the store**
(`core/src/resolve/*`, `core/src/data/store*.ts`, `summaries*.ts`): touch none of those. Your work is in the runtimes
(`core/src/runtimes/*`) and what RG's protocol already offers them; a change the store must make is listed in the
report for after GA2 merges, not made.

**Task**: [genome-formats](../../../design/genome-formats.md) "What every runtime shares" (Levels of a tab-line source
by region, a view's tiles together bounded), "VCF" (the matrix), "GFF3 and GTF" (span rows by region), Bounds,
Defaults 3, Consequences' FR line; [region-reads](../../../design/region-reads.md) for the shapes (items at a tile size,
`density`, `bins`, `matrix`, `MatrixSummary` as form `summary`, span rows, `readCost`):

1. **Tabix levels from the index's bytes** for `bed`, `bedgraph`, `gff3`, `gtf`, `vcf` with an `index`: candidate levels
   `1, 16, 256, 4096` up to `maxLevel`; a level listed only when every tile on every system fits `maxRegionBytes`
   (`regionBytes` from the index, never the data); above 1 a `density` (ranges) or `bins` (matrix) **counted exactly**
   from that tile's records; `readCost(request)` from `regionBytes` so the store's `maxBytesPerLayer` bounds a view
   (ClinVar at the fit: "zoom in" until the view is a few hundred Mb wide). Levels from a dense and a sparse generated
   index tested; the recorded ClinVar index's view cost past `maxBytesPerLayer` at the fit ("zoom in", no data read).
2. **GFF's span rows by region**: a tile inside an intron gives the transcript's span row (EC-229); every tiling of a
   generated GFF equal to the whole read clipped (EC-227).
3. **VCF's matrix**: items the ALT alleles in the tile × states the samples (`ShapeInfo.states`, header order, narrowed
   by `samples`), `codes` the genotype class (0 no copy, 1 some, 2 every copy incl. a haploid `1`, 3 partly called,
   255 uncalled), `value` float32 by `RegionRequest.field` (`GT` copies, `DS`, `GQ`; NaN uncalled); coarse levels as
   bins; genotypes parsed only for a `matrix` request and only for the `samples` held; `maxSamples` bound; EC-220,
   EC-222 (a long deletion one item across the SNPs it covers). `samples` given on a sites read stays said once (VC's
   info) only where no matrix is asked.
4. **BEDPE through the store** (RS's gap 2): `BedpeSource` lists the `pairs` shape and serves pairs reads (read whole,
   indexed by both ends, within `summaries.maxBytes`), so a BEDPE file and held links are one picture. If the store must
   change for it, list the change for after GA2.
5. **Stories**: `Features/Genome formats/VCF` gains a genotype matrix (generated, 40 samples, with a long deletion) and
   ClinVar at the fit reading coarser levels (exact counts) before "zoom in"; BED and GFF by region at a coarse level;
   the BEDPE row in `Features › Region reads › Pairs` drawn from the BEDPE runtime. Light and dark, canvas and SVG,
   headless twins on a fake clock that compose exactly as the stories do.

**Performance** (each with its path): a VCF region of 10,000 records decoded under 100 ms (sites, built) and **with 100
samples' genotypes under 200 ms** (heavy, `itWithin`; 1,000 quick); a level's counts from a tile without allocation per
record in a warm decode (an allocation case if it is on a frame's path).

**Tests**: per item above; the conformance suite (`data/_test/source-suite.ts`) at every listed level for each format;
a region read equals the whole read clipped at every level (a property over generated files); fuzzed and truncated
files without a throw; each bound refused before it allocates; quadratic scans checked by counting (VC's review found
three); EC-171, EC-220, EC-222, EC-227, EC-229 in titles.

**Where**: a worktree on `seq/fr`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); commit at each
green step; hot files never staged (`package.json`s, `tsconfig*.json`, lockfile, `biome.json`, `.scaffolds/seq/`,
`.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`): list them. Anything in `@molstar/core` or
`@molstar/io`, the language, or the network stops the lane. A file past 800 lines is a question, past 1500 a failure.
Sizes (reading 260, headless 920, default 1,460, minimal 1,040 kB; each runtime under 30 kB): a bound passed is a
question.

**Checks**: the quick lane (under 2 s a touched folder), `npx tsc -b packages/seq/core packages/seq/viewer
packages/seq/examples`, each `typecheck:test`, the storybook typecheck, `pnpm lint`, `pnpm format:check`, both perf
tiers, `pnpm -s -C packages/seq/core size`, `pnpm test:medium packages/seq/`, `pnpm test:heavy packages/seq/` when load
is under 6, the `.scaffolds/seq/_checks` (features-planned with `--plan .scaffolds/seq/execution/fifth-slice/plan.md`). A fix's test fails
on its parent. Stop at a green commit if you run long and report what is left.

**Report**: under 350 words: commits; checks with numbers; story ids; sizes; the public surface; store changes needed
after GA2; what is left; hot-file needs; the `.v6/changelogs/` text.
