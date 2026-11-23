# Lane RG: pairs and grids by region

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/rg`. Model:
strong, high effort (a design's core: the store). **Alone in core's store and resolve** (plan, Order and parallelism).
Beside you: VC (`runtimes/vcf`, maybe the ranges layer's drawing) and BS (storybook only). Shared files
(`data/protocol.ts`, `render/ranges*.ts`) only grow on either side; list what you add there.

**Task**: build [region-reads](../../../design/region-reads.md) (Decision, all sections; Bounds; Defaults 1 to the
last; the interaction matrix's cells; Configuration; Decided by the author 1–17), as its **Consequences** list them for
`seqx-core`, `seqx-viewer` and Storybook. GC1 declared the language part (`pairs.source`, `pairs.field`,
`grid.source`, `grid.field`, the `matrix` kind and layer type): build what reads it, never change the language (a
param found missing is a question to the orchestrator). In this order, each a green commit:

1. **Splits first, no behaviour change** (§module-size): `data/store.ts` (1,042 lines) and `render/pairs.ts` (1,334)
   grow here; split each by concern before adding to it, a commit of its own titled "split, no behaviour change".
2. **Data**: the answer shapes and forms in `protocol.ts` (`RegionRequest.y`, `.side`, `.items`, `SystemInfo.side`,
   `RangesList.itemStart`/`itemEnd`, `OpenedSource.readCost`, `ShapeInfo`'s members, `MatrixItems`); `keys.ts` (the 2D
   key, the sided part, the items level part, symmetric folding); the store (`Want.y`, keys over the product,
   `maxTilesPerLayer` 1,024, `maxBytesPerLayer` 16 MiB from `readCost`, the margin ring, the new forms' bytes);
   `summaries.ts` (`binPairs`, `coarsenBins`, grid levels, the whole-assembly level); `tables.ts` (`checkAnswer`,
   `sliceTable` with span rows, `joinTables` per form); `seq-format.ts` round trips; the conformance suite's new
   cases (a region read equal to the whole read clipped at every level, every shape; either end found).
3. **Resolve**: access plans for `pairs`, `grid`, `matrix` by source; the states system from the coverage.
4. **Render**: `layer-data.ts` (grid level, pairs levels, `readTiles` with `y`, bins then items); the pairs layer
   (held and region pairs through one reader, a bin pair one arc stroked by count, gated by `pairs.maxArcs`); the grid
   layer (store tiles as image tiles, `fill`, mirrored tiles, `Hit.bin`); images on rings as runs; the `matrix` layer
   on the rows layer's band machinery; the hover texts (core) and the viewer's hover card for a bin and a matrix cell.
5. **Stories**: a Design page "Region reads" and `Features/Region reads/*`, one per Default, light and dark, canvas
   and SVG, each with its headless twin on a fake clock (tiles arriving coarse to fine). Their data: **generate it in
   test and story fixtures yourself** (a symmetric Hi-C with two levels and a one-level twin, sided pairs in the `seq`
   format, 12,000 links held and as BEDPE, a genotype-like matrix with a long deletion, values with missing ones);
   the Genomes dataset is GX's, later: keep generators small and in `seqx-examples` where GX can reuse them.

The VCF's matrix (genotypes) and tabix levels are FR's, after you; give FR a `matrix` shape and `readCost` to fill.

**Tests** (design, Consequences, Tests): the property tests (region equals whole clipped at every level and shape;
pyramid counts sum per level and bin pair; a pair from either end deduped once; a symmetric tile transposed; held and
region pairs equal bins); drawing against the data (every bin pair in view its arc and hit; every cell hit as its bin;
across a 2D tile edge; a long deletion solid). **EC-59, EC-62, EC-108, EC-125, EC-169 to EC-173, EC-206, EC-208,
EC-220 to EC-223, EC-229 (span rows) named in titles.** Performance, each in its design's form: the pairs pyramid of
12,000 under 10 ms (quick) and 600,000 within 300 ms of slices (heavy, `itWithin`); items over 10 Mb within
`maxTilesPerWant`; the hysteresis at the threshold; a region grid's pan frame under 8 ms with 4 × 4 tiles; the
reference view's 2D bytes under 44 MB (counted); a whole-genome Hi-C at the fit in one want; a ring of 40 matrix
states under 8 ms; a warm frame allocating nothing in both perf tiers (pairs, grid, matrix).

**Two consumers** for the shared pieces: Seq\*'s viewer and a host adapter (`seq-data-ports`) serving the new shapes
through the conformance suite, in one test.

**Configuration**: the design's table copied into core `_spec`'s, every default present (`config-table.test.ts`).
**Docs**: core `_spec` Data (Protocol, Tables, Summaries, Store), Resolve, Render (Data, pairs, grid, the matrix layer,
rings' images, `bin`, Configuration). The other designs the Consequences amend (`.seq/design/`) are hot: list each
edit for the orchestrator.

**Sizes**: shapes, keys and summaries about 6 kB in headless, the matrix layer about 5 kB; bounds reading 240,
headless 870, default 1,400, minimal 1,000 kB. A bound passed is a question, not a raised number.

**Stop and ask**: anything in `@molstar/core` or `@molstar/io`; a language change; a budget that needs a design change;
a design question (write it in the report as a drift note: what, where, options, a recommendation). If the lane runs
long, stop at a green commit after step 3 or 4 and report what is left rather than rushing the rest.

**Where**: a worktree on `seq/rg`. Headers updated (`@author David Sehnal <david.sehnal@gmail.com>`, 2026). Hot files
never staged (`package.json`s, lockfile, `biome.json`, `.seq/`, `.v6/changelogs/`, `storybook/packages.ts`,
`packages/seq/TODO.md`): list exact edits. A file grown past 800 lines is a question in the report, past 1500 a
failure.

**Checks**: as `fx1-checks-and-carried.md` lists them; the quick lane under 2 s per touched folder; both perf tiers;
`pnpm -s -C packages/seq/core size`; `pnpm test:medium packages/seq/`; `pnpm test:heavy packages/seq/` (alone: not
while another lane runs its heavy tests; say when you run it); the `.seq/_checks` and `features-planned.ts --plan
.seq/execution/fifth-slice/plan.md`. A fix's test fails on its parent; a failure is the base's only when run there.

**Report**: under 400 words: commits; checks with numbers; story ids for the orchestrator's look; the public surface
added; sizes; drift notes; what is left if stopped early; hot-file needs; the `.v6/changelogs/` text.
