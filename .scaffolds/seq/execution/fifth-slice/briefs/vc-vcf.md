# Lane VC: VCF sites

**Base**: `packages-seq` at `0d05f5232` or the orchestrator's later head (`git reset --hard packages-seq` first).
Branch `seq/vc`. Model: strong (a file is untrusted input). Runtimes only, never the language: GC1 declared
`source.samples`; this lane reads only what is already declared. The only agent running (one at a time, Decided).

**Task**: build [genome-formats](../../../design/genome-formats.md)'s `runtimes/vcf`, **sites**, whole and by region at
level 1 (Decision: "What every runtime shares", "VCF"; Bounds, the VCF rows; Configuration rows `maxLineBytes`,
`maxHeaderBytes`, `maxInfoFields`, `maxSamples`, `maxAlleles`), on BZ's io (`@molstar/io/common/bgzf`,
`common/tabix`) and `runtimes/internal/` (text, lines, issues, whole-lines). Read `runtimes/bed*.ts` (BE: whole and by
region at level 1, `bed-indexed.ts`) and `runtimes/gff*.ts` (GG) as merged and follow their shape; reuse, never copy:

- `readVcfHeader` within `maxHeaderBytes`: `##fileformat`, `##contig` (names, lengths: EC-188), `##INFO` and
  `##FORMAT` (`Number`, `Type`), the sample names (at most `maxSamples`); a header without `#CHROM` malformed;
- **a row per ALT allele** (EC-196), every allele of a line with **the line's extent by the index's rule**: `POS − 1`
  to `INFO/END` when INFO has `END=` (whatever its alleles, an SNV's too), else `POS − 1 + len(REF)`; the anchor base
  kept (Decided 12). `id` `pos:ref>alt`, `ref`/`alt` over 32 characters as `<length>:<FNV-1a hex>` (core's hash, not
  a new one), `#2`, `#3` … on a repeated key in file order; `ID` the field `name`, the label;
- fields: `ref`, `alt`, `qual` float32, `filter` category, `category` (`SNV`, `MNV`, `insertion`, `deletion`, `indel`,
  `structural` for `<…>`, `breakend` for `]`/`[`); declared INFO keys typed by the header (`Integer` int32, `Float`
  float32, `Flag` and `String` categories; `Number=A` split by allele, `Number=R` the reference's and the allele's,
  other numbers their text), at most `maxInfoFields`, the source's `fields`/`types` first, then header order, the
  rest said once with their count; a key clashing with a fixed field said once; ALT `.`, `*`, `<NON_REF>` no item
  (counted, info); past `maxAlleles` the record left out, said;
- **genotypes are not parsed** in this lane: the matrix, `samples` and genotype classes are FR's (plan, stage 3).
  `samples` given on a VCF source is accepted and said once as not yet read (`info`), never an error;
- whole (plain or bgzipped, within `limits`) and by region at level 1 through tabix with `index` (EC-186: named,
  never `<uri>.tbi`); one line parser either way; a row past a known length clipped, counted (EC-51); issues
  `seq.source.vcf.*` by the shared rule; every bound refused before it allocates;
- `runtimes/vcf` an opt-in entry, in no bundle that does not list it, its own size row under 30 kB in
  `packages/seq/_measure/sizes.mjs` (about 15 kB with io's BGZF and tabix); the viewers' bounds unchanged;
- **the default picture** (design, Defaults 2 and 3): a variant a mark at its REF extent coloured by `category` with
  the ranges layer as built; a long deletion a bar (EC-222's drawing); past level 1 or the read cost, the built
  too-many-tiles state, never an estimate. If drawing needs a ranges-layer change, keep it small, tested, and say it.

**Tests**: the conformance suite (`data/_test/source-suite.ts`) on a generated fixture and on the recorded ClinVar
slice; **a region read equals the whole read clipped** (a property over generated bgzipped files with random regions,
items deduped across tiles); **the extent equal to `tbx_parse1`'s per line** (with `END=` on an SNV line, a symbolic
deletion, a breakend); ids of duplicates and long SVs; INFO types and each `Number`; multi-allelic lines; fuzzed and
truncated files without a throw; each bound refused before it allocates. EC-51, EC-186, EC-188, EC-190, EC-195,
EC-196, EC-222 (as one item across the SNVs it covers, sites only) named in titles. Performance: **a VCF region of
10,000 records decoded under 100 ms**, sites (heavy, `itWithin`; 1,000 quick). No genotypes test (FR).

**Real slice** (design, "Real slices to record", #2), with BE's script (`scripts/record-pages.ts`,
`examples/src/recorded-pages.ts`): **download approved by the author, 2026-11-21 ("Fetch all four")**: NCBI ClinVar
GRCh38 archive `clinvar_20240107.vcf.gz` and its `.tbi`; record the `.tbi` whole (about 520 kB), the header's pages
and TP53's chunks (GRCh38 chr17:7,668,402-7,687,550; ClinVar names it `17`, say how the names rule meets it). Into a
fresh scratch folder, never run anything from it; **delete the full file after recording** (never committed). Check
NCBI's terms when fetched and say them in its `sources.json` row. If the archive URL has moved, stop and ask; do not
pick another release.

**Stories**: `Features/Genome formats/VCF` (TP53's ClinVar variants, the file's header summary and the source node
beside the canvas and the SVG; a generated file with a long deletion over SNVs; past the read cost "zoom in"; a
misread file and its issues), light and dark, with the headless twin on a fake clock (no real-turn waits). Hover and
click an allele per the design's matrix row "a variant allele" (select its REF extent, double-click focuses it).

**Docs**: core's `_spec` Data section for VCF and the Configuration table rows (the config-keys check); the design's
"Not built" line and other `.scaffolds/seq/` amendments as hot-file needs.

**Where**: a worktree on `seq/vc`. Headers updated (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); commit at
each small green step; hot files never staged (`package.json`s, lockfile, `biome.json`, `.scaffolds/seq/`, `.v6/changelogs/`,
`storybook/packages.ts`, `packages/seq/TODO.md`): list them with the exact edit (the `./runtimes/vcf` export line).
Public API changes allowed in `packages/seq`; anything in `@molstar/core` or `@molstar/io` stops the lane and comes
back as a question. A file grown past 800 lines is a question in the report; past 1500 a failure.

**Checks**: as `fx1-checks-and-carried.md` lists them; the quick lane under 2 s for each folder touched;
`pnpm -s -C packages/seq/core size` (reading 240, headless 870, default 1,400, minimal 1,000 kB);
`pnpm test:medium packages/seq/`; `pnpm test:heavy packages/seq/` (alone); both perf tiers if a frame path changes;
the `.scaffolds/seq/_checks` and `features-planned.ts --plan .scaffolds/seq/execution/fifth-slice/plan.md`. A fix's test fails on its
parent. Do not claim a failure is the base's without running it on the base.

**Report**: under 300 words: commits; checks with numbers; the story ids for the orchestrator's look; the public
surface added; the licence found; sizes; hot-file needs; the `.v6/changelogs/` text.
