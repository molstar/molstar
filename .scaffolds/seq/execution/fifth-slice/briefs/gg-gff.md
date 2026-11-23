# Lane GG: GFF3 and GTF

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/gg`. Model:
strong. Runtimes and the ranges layer's gene line; never the language (GC1 declared `source.features` and
`line: 'gene'`). GC3 runs beside you (navigators and the viewer): no overlap expected; shared files only grow, listed.

**Task**: build [genome-formats](../../../design/genome-formats.md)'s `runtimes/gff` read whole (Decision: "What every
runtime shares", "GFF3 and GTF: transcripts", "A name index"; Bounds; Configuration rows marked GG), on BZ's
`runtimes/internal/` and BE's shared parts (read `runtimes/bed` and `runtimes/internal/` as merged: reuse its line
source, issues and recorded-pages tooling; do not duplicate):

- GFF3 and GTF read whole: `assembleTranscripts(rows)` exported and tested alone (transcript types by
  `transcriptTypes`, pieces, genes, a prokaryote's gene → CDS said once, a shared exon a piece of each parent, cycles
  and depth past 16, attributes past 256); each transcript an item with `thick_start`/`thick_end`, `cds_phase`, its
  fields (`gene`, `gene_name`, `type`) and its own extent (`itemStart`, `itemEnd`); GTF's attributes (quoted and
  unquoted, repeated keys joined) with and without `transcript` lines (`seq.source.gtf.extent`); `features:
  transcripts | genes | rows`; `##sequence-region` lengths (EC-188), a whole-sequence `region`/`chromosome` feature
  describing the system, not an item; `##FASTA` ignored but counted; the fixed conventions said once if written;
  issues `seq.source.gff.*`, `seq.source.gtf.*`;
- **the name index** (EC-228): `findNames` on a GFF3 or GTF read whole (GFF3's `Name`, `ID`, `Alias`; GTF's
  `gene_name`, `gene_id`, `transcript_name`, `transcript_id`), a map built once at open, exact then prefix, at most
  `max`; `Coverage.names`. GC3 wires the locus box's fallback to `findNames`; your index must satisfy its protocol;
- **render**: the ranges layer's `line: 'gene'` (`LineStyles`): pieces as boxes, thick between `thick_start` and
  `thick_end`, thin outside, a hairline across the item's extent with strand chevrons; at 0 bytes a warm frame (both
  perf tiers, B6);
- by region is FR's (after RG): read whole only here, and an `index` on a GFF source says so once;
- `runtimes/gff` an opt-in entry, its own size row under 30 kB in `packages/seq/_measure/sizes.mjs` (about 12 kB);
- tests: the conformance suite on a generated fixture and on the recorded real slices; `assembleTranscripts` (three
  levels, a shared exon, CDS-only, GTF with and without transcript lines, a cycle, a prokaryote's gene, a type outside
  the list); fuzzed and truncated files without a throw, every bound refused before it allocates; EC-188, EC-227,
  EC-228, EC-229 (its whole-read part) named in titles; 100,000 GFF3 lines assembled under 500 ms (heavy, `itWithin`;
  10,000 quick);
- **real slices** (design, "Real slices to record", #7 and #8), with BE's recording script
  (`scripts/record-pages.ts`, `examples/src/recorded-pages.ts`): **downloads approved by the author, 2026-11-21
  ("Fetch all four")**: Ensembl GRCh38 `.chr.gff3.gz` and GENCODE `annotation.gtf.gz`, each cut to TP53's lines (GRCh38
  chr17:7,668,402-7,687,550; keep `##sequence-region` for chr17 and the GFF3 header) by the script, so the committed
  fixture is a few kB. Download into a fresh scratch folder, never run anything from it, and **delete the full files
  after cutting** (40 to 50 MB each; never committed). Check each licence when fetched (Ensembl's and GENCODE's terms)
  and say it in its `sources.json` row;
- stories: `Features/Genome formats/GFF3` and `/GTF` (TP53's transcripts as gene models from the recorded slices, the
  `features` switch, a name found by the locus box if GC3 has merged, else by `findNames` shown; a misread file), light
  and dark, canvas and SVG, with the headless twin on a fake clock (no fixed real-turn waits);
- core's `_spec` (Data: GFF and GTF, names; Render: the gene line); `files.md`'s GFF3 tension closed and other `.seq/`
  amendments as hot-file needs.

**Where**: a worktree on `seq/gg`. Headers updated (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot files
never staged (`package.json`s, lockfile, `biome.json`, `.seq/`, `.v6/changelogs/`, `storybook/packages.ts`,
`packages/seq/TODO.md`): list them with the exact edit. Public API changes allowed in `packages/seq`; anything in
`@molstar/core` or `@molstar/io` stops the lane and comes back as a question.

**Checks**: as `fx1-checks-and-carried.md` lists them, both perf tiers, `pnpm -s -C packages/seq/core size` (viewers'
bounds: reading 240, headless 870, default 1,400, minimal 1,000 kB), `pnpm test:medium`, `pnpm test:heavy
packages/seq/`, the four `.seq/_checks` and `features-planned.ts --plan .seq/execution/fifth-slice/plan.md`. A fix's
test fails on its parent. Do not claim a failure is the base's without running it on the base.

**Report**: under 300 words: commits; checks with numbers; the story ids for the orchestrator's look; the public
surface added; the licences found; sizes; hot-file needs; the `.v6/changelogs/` text.
