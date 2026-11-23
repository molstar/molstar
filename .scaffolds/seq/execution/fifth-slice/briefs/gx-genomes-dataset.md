# Lane GX: the Genomes dataset, generated, and its views

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/gx`. Model:
strong. Beside you: **RX** (the ring extension: `packages/seq/ext-rondo`, and in `seqx-examples` only
`src/genomes/ring/` with the GRCh38 lengths, the 12,000 links and the ring views), **FX19** (genotype colour tokens: the
VCF genotype classes' theme), **FX20** (`test-harness/`), **FX21** (`formatGenomic` moved to `@molstar/core`). Stay out
of their files; read the genotype classes' colours through the theme as they are.

**Why**: [examples](../../../design/examples.md), "Genomes (the fifth slice; promoted 2026-11-21)", agreed by the author
2026-11-21 with its three decisions. Read it whole first (the views table, the generator outline, sizes, edge cases,
configuration), then the five genome designs it names as you need them (genome-coordinates, region-reads,
genome-formats, across-sections; ring-extension only for what RX leaves you).

**No fetch**: the four real slices are recorded already (`packages/seq/examples/fixtures/genomes/`: cytoBandIdeo,
ClinVar and its `.tbi`, Ensembl GFF3 and GENCODE GTF around TP53, with `sources.json` rows). Never touch the network; a
slice you think missing is a question in the report.

**Task**, a green commit at each step:

1. **The generator**, `packages/seq/examples/src/genomes/generate/` (modules under 800 lines): `generateGenomes(options)`
   returning bytes and truths, steps 1 to 6 and 8 of the outline (assemblies with toy A prefixed and toy B unprefixed,
   the alignment truth, genes in GFF3 and GTF, the bedGraph signal, the VCF, methylation wide and long, the Hi-C and dot
   plot grids). **Step 7 (the links) is RX's.** Randomness only from core's `PCG` and `hash3`, a stream per
   `(seed, concern, index)`; no clock, no `Math.random`. Where a generator already makes the design's data (the
   40-sample VCF of `generated-genotypes.vcf`, the BED12 genes), fold it in with its bytes unchanged, or say why they
   change. The record script's `genomes` case writes the fixtures and their `sources.json` rows
   (`LicenseRef-seqx-generated`); a test that the output is byte-identical to the committed fixtures (checksums), and
   each file under its estimated size.
2. **The views**, builder modules beside the existing `genomes/*.ts`, in the cut order's reverse so a cut drops the last:
   synteny ribbons (an inversion reversed, a translocation, a gap); the dot plot (two levels, the gaps tinted as
   built); signal by region (open and one-sided domains, log); samples' variants as rows with the deletion solid, and
   the genotype matrix; methylation by region and the gene-by-sample heatmap; transcripts in three levels; open on a
   gene by name; the genome by region with its names aliased; the contact map; the liftover (A over A′, built in LA4).
   **Not** the genome ring nor both genomes on one ring (RX's, and the orchestrator's after both merge).
3. **The matrix cells**: each view's test drives the interaction matrix's cells it names (the design's views table and
   the designs' matrix sections) and the edge cases in its row, each EC in a test title.
4. **Stories and twins**: `Examples/Genomes/*` (the dataset's other examples' shape), light and dark, canvas and SVG,
   each with its headless twin on a fake clock; the storybook's `views.ts` registry gains them.

**No language change**: a view that needs one stops with a report (cut that view; the cut order is contact maps, then
the liftover). Anything in `@molstar/core` or `@molstar/io` stops the lane with a report. Core or viewer code of Seq\*
only for a bug a view finds, its test failing on the parent, said in the report.

**Sizes**: the examples' fixture bounds (`maxFileBytes` 1 MB, `maxDatasetBytes` 2 MB): the recorded ClinVar pages are
already about 1.36 MB of the folder; measure the dataset with your fixtures and report it with each file's size. Never
edit a bound; the orchestrator decides (trim first: bgzip where the reader takes it, fewer sites before fewer views).
Viewer bundles must not grow (this lane is examples and stories).

**Where**: a worktree on `seq/gx`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); commit at each green
step; hot files never staged (`package.json`s, `tsconfig*.json`, lockfile, `biome.json`, `.scaffolds/seq/`, `.v6/changelogs/`,
`storybook/packages.ts`, `packages/seq/TODO.md`, `packages/seq/_measure/sizes.mjs`): list what you need. A file past 800
lines is a question.

**Checks**: `pnpm test:unit packages/seq/ cli/seqx-render/`, `npx tsc -b packages/seq/core packages/seq/viewer
packages/seq/examples`, each `typecheck:test`, the storybook typecheck, `pnpm lint`, `pnpm format:check`, both perf tiers
to completion, `pnpm -s -C packages/seq/core size` (read the rows for FAIL), `pnpm test:medium packages/seq/`, the
`.scaffolds/seq/_checks` (features-planned with `--plan .scaffolds/seq/execution/fifth-slice/plan.md`, `defaults-designed.ts`).

**Report**: under 400 words: commits; checks with numbers; the views built and any cut, with why; story ids; the
dataset's size and each file's; the byte-identity test; ECs covered; hot-file needs; the `.v6/changelogs/` text; designs
to amend (examples.md's Genomes: what was built as written and what moved).
