# Plan: the fifth slice

**Reconciled 2026-11-21** with the promoted designs (the fourth slice's finalize, step 5): Gate 0 decided
(2026-11-20), stage 0's six drafts agreed and promoted (2026-11-21), the whole-designs review's twelve recommendations
and the fold-in's three taken (`state.md`, Decided). The loop is [loop.md](loop.md) (the fourth slice's with its
fold-back, B1–B11, applied). The first draft of this plan (2026-11-20) is in the git log.

Work id `fifth-slice`. Version 1's fourth milestone, **Genomes** ([version-1](../../design/version-1.md), Milestones
and Deferred): assemblies and chromosomes on one axis, two assemblies on one axis or ring, pairs and grids read by
region, layers across sections (synteny, dot plots, contact maps, a liftover), three genome format families with BGZF
and tabix, the Genomes dataset and its views, and the showcase, **a genome ring like Rondo's, bundled and
count-weighted, as an extension package** (the author, 2026-11-19: "Rondo-like, extension"), on the `<id>:<name>` kinds
the fourth slice built (EX, CN).

## Where this is heading

| Slice | Milestone | In one line |
|---|---|---|
| first to third (closed) | the slice, 1, 2 | the tree core to figures |
| fourth (closed) | 3, integration | molecules and units, structure syncing, track headers, seq data ports, services, bigWig, archives, extension kinds |
| **fifth (this plan)** | 4, Genomes | the designs below, built |
| after | version 1's close-out | the planned rows no milestone took (Gate 0, 9), and what this slice deferred: the Manhattan family (EC-214, 216, 217, 219), `split`, EC-224's join, the designed-only formats when a use asks |

## The designs (promoted 2026-11-21)

| Design | Lanes |
|---|---|
| [genome-coordinates](../../design/genome-coordinates.md) (GA) | GC1, GC2, GC3 |
| [region-reads](../../design/region-reads.md) (RR) | RG |
| [genome-formats](../../design/genome-formats.md) (GF) | BZ; BE, GG, VC; FR (its Consequences, Lanes) |
| [across-sections](../../design/across-sections.md) (LX) | LA |
| [ring-extension](../../design/ring-extension.md) (RN) | RX |
| [examples](../../design/examples.md), Genomes (GD) | GX |

Each lane amends the existing designs its design's Consequences name and that promotion did not (addressing,
locations-and-queries, maps, building-blocks, molecules, data-store, seq-data-ports, performance, security,
viewer-models, ux-foundations, themes, features' rows, edge-cases' outcomes), with their Design pages (loop, step 9).
Promotion already edited axis, extension-openings, files, sources, examples, layouts, navigators, grids-and-rows,
protocol and base-library.

## How this slice runs (from the fourth slice's fold-back)

- **Two lanes at once** ([loop.md](loop.md), Lanes at once). Every agent counts (reviewers, measurements), listed with
  start times in `state.md`'s In progress before it starts (B3).
- **Files limit what runs together**: never two lanes on the language, on core's store or resolve, or in one file;
  the heavy lane alone, on a quiet machine (`uptime` under 3, B5).
- **The orchestrator looks itself**, in a visible pane, with real pointer events for every gesture of a part's matrix
  row (B1, B2, B4); a reviewer checks each `EC-<n>` its brief names in a test title (B8).
- **Sizes are bounds, moved by the author only**; every runtime, the across-sections chunk and the extension load as
  entries or lazily, never in the minimal viewer.
- **Downloads asked before they are made**: the four real slices (genome-formats, Real slices to record: cytobands,
  ClinVar, Ensembl GFF3, GENCODE GTF), each with its licence checked at the ask; tests never fetch.
- **The quick lane's folder budget**: a lane over 2 s a folder moves its own slow tests before its review.

## Done when

- **Genome coordinates**: an assembly's chromosomes as blocks with a direction and a scale each (axis.md's tension
  narrowed); two assemblies on one section, axis or ring (EC-210); small blocks as runs; switching chromosomes and
  core's navigator over several systems with the square (EC-184); the locus and block address forms; sections over a
  source's complex reflowing on its changes (HD's drift).
- **Pairs and grids by region**: the `pairs`, `matrix`, `grid` and `rows` shapes; `RegionRequest.y`; 2D tiles; the
  pairs level and grid summaries with the whole-assembly level; bounds in bytes; images on rings as runs (EC-208,
  EC-223).
- **The genome formats**: BGZF and tabix in `@molstar/io`; `runtimes/bed` (BED, bedGraph, BEDPE, chrom.sizes),
  `runtimes/gff` (GFF3 and GTF as transcripts, both names indexed), `runtimes/vcf`, each with its conformance suite on
  generated and recorded files; the `seq` runtime's sided pairs and a map's runs. bigBed, 2bit, FASTA by region, PAF,
  PIF and chain designed only (no lane).
- **Layers across sections**: the `bridge` (synteny ribbons, the lock, align, select across), contact maps and dot
  plots over two sections, categorical grid axes (EC-209); the liftover (cut last).
- **The Genomes dataset and its views** (examples.md, Genomes): generated byte-identically, the four real slices
  recorded as pages; the ring, synteny, both genomes on one ring, the dot plot, the contact map, signal by region,
  variants and the genotype matrix, methylation, transcripts, open by name, names aliased, the liftover.
- **The ring extension** `packages/seq/ext-rondo`: stems, ribbons wide by count, stubs, hues, cytobands, the ring zoom
  by `section.include`, the square picture, its fallback; 23 chromosomes and 12,000 links at frame time.
- **The fourth slice's carried items** placed: FX1's checks and items, M4's measurements in M5, the novice design
  questions answered at stage 1's gate.
- **Measured** (M5).

## Gate 0 (decided)

**Decided (the author, 2026-11-20: "the recommendation make sense")**: every question as recommended.
1. **Cut order** if the slice runs long: MAF (cut 2026-11-21), then contact maps, then the matrix layer, then the
   liftover. The ring, synteny and the formats are the milestone.
2. **BGZF and tabix** in `@molstar/io`, beside its readers (genome-formats: the surface agreed; no `.gzi`).
3. **The ring extension's package**: `packages/seq/ext-rondo`, in the inventory.
4. **Bundling**: core's pairs level, which the extension reads (region-reads).
5. **Ring zoom**: a cut (`section.include` in place; a second section as the builder's option).
6. **Cytobands**: synthetic in tests (`syntheticBands`), UCSC's in the showcase if its licence allows, asked first.
7. **Real data**: generated first, one real slice a built format, each fetch asked (four fetches).
8. **A navigator over several systems**: core's (genome-coordinates).
9. **Version 1's other planned rows**: the close-out slice after this one.

**Public API** (the author, 2026-11-20): a lane changes `packages/seq`, `core/repr` and `packages/ui` where its design
needs it, recorded in its changelog entry; `@molstar/io`'s BGZF and tabix and core's `ungzip` fix stop at the surface
genome-formats agreed.

**Taken since (2026-11-21)**: three format families built, the rest designed only; the twelve recommendations
(`reviews/whole-designs.md`: no `gzi`; BEDPE indexed by both ends after a whole read in a byte budget; no doubled
BEDPE; a per-view byte bound for tabix levels; a host's bigWig the whole-genome remedy; the Manhattan family, `split`
and EC-224 to the close-out; MAF cut, EC-221 to the matrix, EC-226 to themes; fetch 3 cut; RN's gene detail on GFF and
`coarsenBins`); fetches 4, 5, 6, 9 and 10 dropped; GTF names indexed beside GFF3's.

## Testing strategy

- **Every lane**: quick-lane tests beside the code; a story with a headless twin on a fake clock; drawing checked
  against the data at the fit and across a tile edge (and, here, across a block gap, inside a run, on reversed blocks);
  a fix's test fails on its parent; edge cases named by id in test titles.
- **Formats, offline**: a conformance suite per runtime on a generated file and its recorded real slice (pages);
  fuzzed and truncated files without a throw, every bound refused before it allocates; BGZF block boundaries; a
  region read equal to the whole read clipped at every level.
- **Coordinates**: property tests (a block's direction reversed twice the identity; hit-testing the inverse of
  placement over mixed directions and scales; switching blocks keeps the selection's addresses; the names rule never
  maps an alias twice).
- **By region and across sections**: every shape read by region equal to the whole read clipped; the pairs pyramid's
  counts summing to the items'; each pair drawn once across a bridge's two sides; the twist by geometry equal to GA's
  rule.
- **Performance representatives**: each design's Consequences list them with their paths (quick representative,
  heavy under `itWithin`, the allocation tier for a frame's path: a ribbon frame, a bridge frame, a ring frame).

## Stages and lanes

Each lane: its brief in `briefs/` (written by the orchestrator before it starts, loop step 3), its `_spec` and
Configuration table kept, a changelog entry in its report, an independent review and one fix round, the orchestrator's
look in a visible pane in light and dark before the merge.

### Stage 0: designs (done)

- [x] GA, RR, GF, LX, RN drafted, reviewed and revised overnight; GD drafted and reviewed (2026-11-20 to 21).
- [x] GATE (2026-11-21): the drafts agreed ("generally reasonable and we can fold in the reviews"), the twelve
      recommendations and the fold-in's three taken, promoted to `.seq/design/` (uncommitted at writing).

### Fix points (fast; [fixes.md](fixes.md); `briefs/fx-fix-lane.md`, copied from the fourth slice's)

The fix phase (the author, 2026-11-21: "lets have a 'fix phase' for slice 5 too where we accumulate these smaller
things?"): small findings that need no design collect in `fixes.md`; an FX lane takes its open items at each point. A
fix point with nothing open is skipped. A finding that needs a design is a drift note, never an FX item.

- [x] **FX1** (merged 2026-11-21, `2afb32e91`) (first, before stage 1's lanes; beside BZ when their files differ):
  - the fourth slice's fold-back checks: **B6** a design promoted from now on whose Consequences name a layer, a layout
    or a part has a `Defaults` heading (§defaults-designed; the five designs above are the first it reads); **B7** a
    `features.md` row still "planned" for a lane `plan.md` marks merged fails; both in `.seq/_checks/`, added to the
    loop's Checks;
  - `packages/seq/TODO.md`'s carried items: the services' structure section at the fit (its ruler by author numbering
    at least, or the follower taking the leader's mapped window); the flake watch list (the services story test,
    Prism's `commands.test.ts` and `collapse.test.tsx`: run again alone, a second failure fixed);
  - **a story that reads a bigWig** (M4's gap: no story reads one), so M5 can measure a bigWig pan;
  - `fixes.md`'s open items at its start.
- [ ] **FX2** after stage 1's gate; **FX3** after stage 2's; **FX4** after stage 3's; **FX5** after stage 4's;
      **FX6** before M5 (the novice pass's and the use gates' findings that need no design). Taken early on the
      author's findings (2026-11-21): FX2 (`e678d2dbc`), FX3 (`ded361691`), FX4 (`5bfebf366`); the next fix lanes are
      numbered on from FX5. **FX5** (2026-11-22, `briefs/fx5-small-fixes.md`, Sonnet) beside VC and RG: six items clear of RG's files.

### Stage 1: coordinates (strong, in order; the language one lane at a time)

- [x] **GC1** (merged 2026-11-21, `e3b3a98d0`) (`briefs/gc1-genome-language.md`) the language and resolve, the slice's one minor bump: the `assembly`
      kind, `section.of` an assembly or a list, `section.include`, `block.direction`, `block.scale`,
      `navigator.switch`, the locus and block forms, the builder's `genome()`, `assembly()`, `colinearBlocks()`; the
      names rule; derived sequences from the coverage; a map between assemblies as one map. It **declares the other
      designs' params** with their docs and checks, no reader: GF's (`source.features`, `samples`, `map.source`,
      `line: 'gene'`), LX's (`bridge`, `pairs.shape`, a map as pairs data, `follow_y` a second section, a rows system
      on a grid axis, `member_order`, `member_order_y`) and RR's (`pairs.source`, `grid.source`, the `matrix` kind and
      layer type: the author, 2026-11-21, "Yes, in GC1").
- [x] **SC** (merged 2026-11-21, `2fdc44e76`) (`briefs/sc-complex-sections.md`; strong, after GC1: resolve) **sections over a source's complex** (HD's
      drift, decided 2026-11-21, `packages/seq/TODO.md`): a section `of` a dynamic complex, resolved when the source
      opens and reflowed on its `changes`; the static-model rule's exception for opened sources (seq-data-ports' matrix
      row "a section of a complex after an update: reflows"); its `_spec` and the TODO line closed.
- [x] **GC2** (merged 2026-11-21, `db2364d09`) (`briefs/gc2-genome-axis.md`; after GC1, beside SC: render only) the axis and layouts: blocks with their
      own direction and scale (rulers, the direction mark, `×N`), the gap share on rows and rings, runs of small blocks,
      group gaps and labels, `section.center`, the `block` and `system` fields' domain, export attributes.
- [x] **GC3** (merged 2026-11-21, `61e71891b`) (`briefs/gc3-genome-navigators.md`; after GC2 and SC) navigators and the viewer: windows over blocks, the
      square over blocks × blocks, `navigator.switch`; `navigate.acrossBlocks`, the cursor across blocks (EC-154
      amended), `block.next`, `.previous`, `.reverse`, `section.include` commands, the label menu's entries; history
      by `{ system, position }`; the locus box's fallback to `findNames` (after BZ's types).
- [x] **GATE** (passed 2026-11-21: the five questions as recommended, the author: "Yes, all five"; the author's look done 2026-11-21: every story walked, "these are the things i noticed, so i guess
      we can consider this the gate for stage 1"; its findings in `fixes.md`, five as class items with sweeps), with the
      **novice design questions** from the fourth slice's pass (`packages/seq/TODO.md`), each a yes
      or no with a recommendation: F3 (the fake's "Within 5 Å" on a selection already the chain), F9 (a marked state
      for the track the structure is coloured by), F11 (an empty section saying its issue on the view, beside SC), F14
      (the "n filtered" card naming the filter, its click saying it clears), F25 (a navigator window at the fit drawn
      as one).

### Stage 2: formats (runtimes, not the language; beside stage 1)

- [x] **BZ** (merged 2026-11-21, `167d24491`) (`briefs/bz-bgzf-tabix.md`; strong, first) BGZF and tabix in `@molstar/io` (`common/bgzf`,
      `common/tabix`, `_spec/bgzf.md`), core's `ungzip` over multi-member streams, `seqx-core`'s io allowlist (hot:
      `biome.json`), `runtimes/internal/` (text, lines), `findNames`' types and `data/names.ts`. Stops at the surface
      genome-formats agreed.
- [ ] Then beside each other, after BZ (and after GC1 for a new param):
  - [x] **BE** (merged 2026-11-21, `98028d3f6`) (`briefs/be-bed.md`) the BED family: BED, bedGraph, BEDPE (read whole, indexed by both ends, held rows
        in a byte budget), chrom.sizes; whole and by region at level 1.
  - [x] **GG** (merged 2026-11-21, `9b8cf675b`) (`briefs/gg-gff.md`) GFF3 and GTF read whole: `assembleTranscripts`, `source.features`, `findNames`
        over both (GFF3's `Name`, `ID`, `Alias`; GTF's `gene_name`, `gene_id`, `transcript_name`, `transcript_id`);
        the ranges layer's `line: 'gene'`.
  - [x] **VC** (merged 2026-11-22, `364acfb5e`) (`briefs/vc-vcf.md`) VCF sites at level 1, extents by the index's rule, the anchor base kept.
  - [x] **GB** (merged 2026-11-22) (added 2026-11-21, the author: "Add a GenBank lane"; `briefs/gb-genbank.md`; strong, untrusted input;
        after VC, alone in the formats track) **GenBank read whole**: a short design section first (genome-formats.md,
        "GenBank", agreed in the brief: what the runtime answers and what it refuses), then `format: 'genbank'`: a
        file of one or more records, each a sequence (LOCUS's name, length and topology: `circular` a cyclic system;
        ORIGIN's letters, the molecule type its alphabet) and its features as ranges (the key the category, the
        qualifiers fields: `/label`, `/gene`, `/product`, `/note`); the location grammar: `complement`, `join` and
        `order` (pieces of one feature, as GFF3's parts), `<` and `>` partial ends (kept as flags), a `join` over a
        circular record's origin as one item across it, remote references (`J00194.1:100..202`) and `^` sites refused
        per feature with an issue, never the file; `findNames` over `/label`, `/gene`, `/locus_tag`. Fixtures: the
        plasmid's derived features written as GenBank by `record.ts` (as its GFF3 is: a record NCBI-shaped, bla's
        `join` across the origin, lacZα `complement`), and the real L09137.2 flat file (its only feature `source`; a
        download asked in the brief); fuzz and malformed-input tests as BE's and GG's. Its own size row (under 30 kB
        like the other runtimes). The story side: "Try pUC19.gb" loads it (fixes.md), and a plasmid story reads the
        GenBank file beside the GFF3 one, the two views equal.
- [x] **BS** (merged 2026-11-22, `5bf984d5b`) (added 2026-11-22, the author: "keep it in core for now, but add a storybook example to features section
      that showcases it"; `briefs/bs-builder-showcase.md`; fast, storybook only, beside any lane) **the builder
      showcased**: `Features/Builder`, each story the builder code shown beside the view it builds; the builder stays
      in `seqx-core/lang` (no package of its own for now).
- [ ] **GATE**, with `pnpm build` and `pnpm check:workspace` (`@molstar/io` changed); each runtime its own size row.

### Stage 3: by region and across sections (after GC1; RG alone in the store)

- [ ] **RG** (`briefs/rg-region-reads.md`; strong) the shapes, `RegionRequest.y`, 2D keys and tiles, `side`, items at
      a tile size, `readCost` and `maxBytesPerLayer`, the pairs level and `binPairs`, `coarsenBins`, grid summaries and
      the whole-assembly level, the `matrix` layer's drawing (its params declared by GC1), images on rings as runs (EC-208, EC-223).
- [ ] **GA2** (added 2026-11-22: GC2 checked data over an assembly but never read it, RS's gap 1;
      `briefs/ga2-data-over-assembly.md`; strong, after RG, alone in resolve and the store) ranges, values, pairs,
      grid and matrix over an assembly read by sequence; held pairs by the pairs level (RS's gap 6); the ideogram's
      bands over every chromosome; RS's single-system stories moved onto assemblies.
- [ ] **FR** (`briefs/fr-format-regions.md`; strong, after RG and stage 2's lanes) the formats' by-region rest: tabix
      levels from the index's bytes with exact counts and `readCost` (BED, GFF, VCF), GFF's span rows by region, VCF's
      matrix and genotype classes; the `seq` runtime's sided pairs from either end, `map.source` by either side.
- [ ] **LA** (`briefs/la-across-sections.md`; strong, after RG and FR) the `bridge`: placement, ribbons drawn once by
      ownership, the closed-form hit, the lock, align, select across, the cursor across; contact maps and dot plots
      over two sections; categorical axes; the liftover last (EC-209, EC-212).
- [ ] **GATE**.

### Stage 4: the dataset and the ring (beside each other)

- [ ] **GX** (`briefs/gx-genomes-dataset.md`; strong, after stages 2 and 3) the Genomes dataset generated
      (`generateGenomes`, byte-identical, under the dataset bound) and its views (examples.md, Genomes), the four real
      slices recorded as pages (each fetch asked), stories and twins; the views' matrix cells.
- [ ] **RX** (`briefs/rx-ring-extension.md`; strong, after RG and GC3; on EX and CN) `packages/seq/ext-rondo`: the
      descriptor, `rondo:ribbons`, `rondo:cells`, the `rondo:square` picture, `rondo:hues`, commands and the binding,
      `rondoRing`; the entry's exports and the `LayerType` hooks in core; the fallback; its stories (hot: the new
      package's files, inventory row, root reference).
- [x] **NS** (merged 2026-11-21, `190dfc135`; a Sonnet subagent) (`briefs/ns-nightingale-showcase.md`; fast, any time, alone: storybook only) **a Nightingale showcase**
      (the author, 2026-11-21: "yes, we can do a dev dependency version. maybe add a phase for this to slice 5?";
      reverses the fourth slice's "dropped" for stories only): Nightingale's web components a dev dependency of the
      Seq* storybook alone, a few stories feeding its tracks from the resolved view and the store, the NG consumer
      beside them; no package of Seq* depends on it.
- [ ] **USE GATE** (not blocking): the author tries the genome views and the ring; findings fixed or placed.
- [ ] **GATE**.

### Stage 5: use, measure, finalize

- [ ] **NOVICE PASS** (`briefs/novice-pass.md`, the fourth slice's shape) after the last lane that draws (RX or GX,
      whichever merges last; B10): a task script per new view (find a gene by name, compare two assemblies across a
      bridge, follow a contact on the ring, switch chromosomes on the square), with a pointer path onto each menu and
      submenu; findings counted by kind and placed (`fixes.md` for FX6, a drift note, `TODO.md`).
- [ ] **M5** (alone, quiet machine; `measurements.md`): the heavy lane and the opt-in tiers; a whole human genome by
      region (by chromosome); a VCF region of 10,000 records decoded; the ring of 12,000 links at full scale in the
      browser at dpr 2; a bridge frame; the bundles (each runtime its own row, the extension its own, the minimal
      viewer not grown). **M4's carried measurements**: a bigWig pan (FX1's story), real-vsync frames with the pane in
      view, the 200-copy split patch re-measured.
- [ ] **FINALIZE** (loop, Finalize): reconcile, fold back, the record, the next plan (version 1's close-out),
      worktrees cleaned.

## Order and parallelism

Three lanes at once (loop.md, the author, 2026-11-22). Three tracks, and the builder showcase (BS, storybook only)
beside any of them:

- **Coordinates** (the language, resolve, render, the viewer): FX1 → GC1 → SC and GC2 beside each other → GC3.
- **Formats** (`@molstar/io` and runtimes): BZ beside FX1 or GC1 → BE, GG, VC (after GC1 for their params), two at a
  time beside whatever runs in other files; GB after VC.
- **Views** (the store, render, examples, the extension): RG after GC1 (alone in the store; never beside SC, which
  touches resolve) → FR → LA; then GX and RX beside each other.

Never two lanes on the language (GC1 alone owns the bump; a later param goes through it), on core's store or resolve
(GC1, SC, RG one at a time), or in one file. Fix points run beside whatever runs in other files.

## Hot files (the orchestrator's only)

- `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json` (seqx-core's io
  allowlist gaining `common/bgzf`, `common/tabix`; ext-rondo's restricted imports), `vitest.config.ts`;
- `storybook/packages.ts`, `.v6/changelogs/`;
- `packages/seq/*/package.json` and `tsconfig*.json` (seqx-core's `./runtimes/bed`, `./runtimes/gff`,
  `./runtimes/vcf` entries and `"./runtimes/internal/*": null`), `packages/seq/ext-rondo/package.json` and its
  tsconfig, `packages/io/package.json`, `packages/core/package.json`;
- `packages/seq/*/_spec/README.md` sections shared by lanes: each lane edits its own; GC1 merges the Language section;
- `.seq/`.

## Known risks

- **A direction per block changes what a position is** for paths, hover, selection, follow and export: GC1 and GC2
  land it with a twin per consumer.
- **Region reads in two dimensions move the store's keys**: RG runs alone in the store.
- **The language changes once**: GC1 declares every design's params; a param found late goes through the bump's
  owner, not a second bump.
- **Format licences and sizes**: real slices recorded as pages, small, each asked before it is fetched.
- **The default viewer is at its bound**: runtimes, the bridge chunk and the extension are entries or lazy.
- **Length**: about 20 lanes and six fix points; Gate 0's cut order decides what goes first.

## Out of scope

Mol\* adapters and a live 3D view (the new mol viewer); BAM and CRAM; `.hic` and `.mcool` (the Hi-C example generated
in the `seq` format); the designed-only runtimes (bigBed, 2bit, FASTA by region, PAF, PIF, chain); MAF; the Manhattan
family, `split`, the keyed join (the close-out slice); settings layers; agent tools; a React adapter; writers
beyond SVG, PNG and archives; the compute port, unless M5's profile of the whole genome shows the main thread is the
limit.

## Briefs to write (the orchestrator, before each lane starts)

`fx-fix-lane.md` (copied from the fourth slice's, its items replaced), `gc1-genome-language.md`,
`sc-complex-sections.md`, `gc2-genome-axis.md`, `gc3-genome-navigators.md`, `bz-bgzf-tabix.md`, `be-bed.md`,
`gg-gff.md`, `vc-vcf.md`, `gb-genbank.md`, `rg-region-reads.md`, `fr-format-regions.md`, `la-across-sections.md`,
`gx-genomes-dataset.md`, `rx-ring-extension.md`, `novice-pass.md` (the fourth slice's shape). Stage 0's six design
briefs stay as the record.
