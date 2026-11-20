# Plan: the fifth slice

**Initial draft, 2026-11-20, written while the fourth slice runs** (its stage 2 and later remain). Not agreed: Gate 0
(below) asks the author; "What the fourth slice changes in how this one runs" is provisional until the fourth slice's
`foldback.md` is decided, and is rewritten from it then. The fourth slice's finalize reconciles this file (its loop,
Finalize, step 5).

Work id `fifth-slice`. Version 1's fourth milestone, **Genomes** ([version-1](../../design/version-1.md), Milestones
and Deferred): the deferred list (layers across sections, synteny, contact maps, MAF over a genome, the Genomes
dataset), the genome formats beyond bigWig (2bit, BigBed, VCF, PAF, GFF, chain files; BED and region reads with them),
two assemblies on one axis, and its showcase, **a genome ring like Rondo's, bundled and count-weighted, as an extension
package** (the author, 2026-11-19: "Rondo-like, extension"), on the `ext:` kinds the fourth slice builds (EX).

## Where this is heading

| Slice | Milestone | In one line |
|---|---|---|
| first to third (closed) | the slice, 1, 2 | the tree core to figures |
| fourth (running) | 3, integration | molecules and units, structure syncing, track headers, seq data ports, services, bigWig, archives, extension kinds |
| **fifth (this plan)** | 4, Genomes | assemblies and chromosomes on one axis, the genome formats, pairs and grids read by region, layers across sections (synteny, contact maps), the Genomes dataset and its views, the Rondo-like ring as an extension package |
| after | version 1's close-out | the planned rows no milestone took (Gate 0, question 9) |

## What the fourth slice changes in how this one runs (provisional)

From the fourth slice's state and loop, to be checked against its fold-back:
- **Lanes**: the author's current limit, one line in `loop.md` (the fourth slice ran one, then two). What runs together
  is limited by files: never two lanes on the language, on core's store or resolve, or in one file; the heavy lane
  alone, on a quiet machine (the fourth slice's first stage-1 heavy run at load 22 failed ten time bounds that pass at
  load 3).
- **Sizes are bounds, moved by the author only**: the fourth slice moved them four times, each with the feature behind
  it; the default viewer ended a few hundred bytes under its bound. Lazy chunks are allowed for parts used rarely or
  late (viewer-composition.md); this slice's formats and the ring extension load as runtimes or entries, never in the
  minimal viewer.
- **The quick lane's folder budget**: the fourth slice's gate passed with three folders over 2 s and handed them to a
  fix lane; here a lane over the budget moves its own tests before its review (the third slice's fold-back, held).
- **A lane's report is checked, not trusted**: the fourth slice's checks found failures lanes reported green (a race,
  a story's twin passing on the parent, an example that drew as before); the look and the "fails on its parent" check
  stay steps of every lane.
- **Downloads asked before they are made**: real genome data (cytobands, gene models, a VCF slice) is fetched once with
  the author's yes, its licence recorded in `sources.json`; tests never fetch.

## Done when

- **Genome coordinates** ([axis](../../design/axis.md), [navigators](../../design/navigators.md)): an assembly's
  chromosomes as pieces of one axis with a direction per piece (an inversion between genomes drawn, axis.md's tension
  closed); two assemblies on one axis or ring (EC-210); switching chromosomes and a navigator over several systems
  (EC-184; extension-openings.md, Decided 7); scaled segments and a zoomed ring decided (layouts.md, Open).
- **Pairs and grids by region** ([grids-and-rows](../../design/grids-and-rows.md), Data;
  [protocol](../../design/protocol.md)): the `pairs`, `matrix`, `grid` and `rows` answer shapes; `RegionRequest.y`;
  2D tiles and the store's grid summaries; a pairs level (bin-pair counts) read by region (EC-208, EC-223); dense
  matrices on a ring decided (`seq.render.ring.image`).
- **The genome formats**, each a runtime an application adds (sources.md), opt-in entries like `runtimes/mmcif`, each
  with a conformance suite on recorded or generated fixtures: BED (`base: 0`, BED12 blocks as pieces), bigBed, 2bit,
  FASTA with its index by region, GFF3 and GTF (the `Parent` hierarchy into transcripts, files.md's tension closed), VCF
  with tabix, PAF, chain files; BGZF and tabix in the base library if Gate 0 says so (base-library.md).
- **Layers across sections** ([grids-and-rows](../../design/grids-and-rows.md), [layouts](../../design/layouts.md)):
  a layer drawn between two sections (synteny ribbons between two assemblies; a liftover by chain file); contact maps
  as grids across sections (EC-209); ribbons and categorical grid axes.
- **The Genomes dataset and its views** ([examples](../../design/examples.md), Genomes): two toy assemblies generated
  (an inversion, a translocation, a gap, aliases, an HLA contig name), three-level transcripts, a name index,
  multi-sample signal, GWAS points and LD, a variant matrix by sample, methylation, keyed results, FASTA by region; the
  views: synteny ribbons, a dot plot, Manhattan, a liftover, the genotype matrix; MAF if Gate 0 keeps it.
- **The ring extension** (a package of its own, on `SeqExtension`): `rondo:ribbons` (wide by count, trumpets as
  cluster stems) and the bin-pair level (core's pairs level or `rondo:binpairs`, Gate 0), cytobands, a hue per block,
  outward links to other chromosomes as stubs, the square navigator; 23 chromosomes and 12,000+ links at frame time.
- **Measured**: a whole human genome by region (performance.md's fixture; heavy, by chromosome), a region of a VCF and
  of a bigBed decoded within bounds, a ring frame of 12,000 links, the bundles (each runtime its own row, the extension
  package its own, the minimal viewer not grown).

## Gate 0 (the author, before any build lane)

Each with a recommendation; "Accept all" takes them.
1. **Cut order** if the slice runs long: MAF first, then contact maps, then the liftover. Recommended: as listed; the
   ring, synteny and the formats are the milestone.
2. **BGZF and tabix**: in the base library (`@molstar/io`, a public API addition, its surface agreed in stage 0) or
   private to Seq\*'s runtimes. Recommended: `@molstar/io`, beside its readers, since Mol\* reads BGZF too.
3. **The ring extension's package**: a new package (`packages/seq/ext-rondo`, in the inventory) or an entry of the
   examples. Recommended: a package, so an extension is shown as an outside author would ship one.
4. **Bundling** (research-genome-ring.md, Q2): a `pairs` level read by region (core) or a layer param over a derived
   table (the extension). Recommended: the core's pairs level, since pairs by region is owed anyway (EC-208).
5. **Ring zoom** (Q5): an animated story channel or a cut to a second section. Recommended: a cut to a second section
   first, the channel only if the use gate asks for it.
6. **Cytobands** (Q6): UCSC's (licence to check, fetched once) or synthetic. Recommended: synthetic in tests, UCSC's
   in the showcase if its licence allows, asked before the fetch.
7. **Real data**: which small real slices to record (a gene region's GFF3, a VCF slice, a bigBed track), each asked
   before its fetch; the toy assemblies generated (examples.md). Recommended: generated first, one real slice a format
   for its conformance suite.
8. **A navigator over several systems** (extension-openings.md, Decided 7: designed here): in core's navigators or the
   extension's square navigator only. Recommended: core's, since switching chromosomes needs it in every genome view.
9. **Version 1's other planned rows** (the fourth slice's Gate 0, question 11: track groups, search, patches on
   streamed things, keep this region, the remaining query forms): here or a version-1 close-out slice. Recommended: a
   close-out slice after this one.

## Testing strategy

- **Every lane**: as the fourth slice (quick-lane tests beside the code; a story with a headless twin on a fake clock;
  drawing checked against the data at the fit view and across a tile edge; a fix's test fails on its parent; edge cases
  named by id).
- **Formats, offline**: a conformance suite per runtime on a generated fixture and one recorded real slice (counts and
  fields recorded from the file); malformed and truncated files fuzzed against the security bounds; BGZF block
  boundaries; tabix and bigBed region reads equal to a full read clipped. No test reaches the network.
- **Coordinates**: property tests (a direction per piece reversed twice is the identity; a location carried across two
  assemblies by a chain file and back covers what it started from; switching chromosomes keeps the selection's
  addresses).
- **By region**: pairs and grids read by region equal to the whole read clipped, at every level; the bin-pair pyramid's
  counts summing to the pairs'.
- **Performance representatives**: a whole genome by region (heavy, by chromosome); a ring of 23 blocks and 12,000
  links (quick at a small scale, heavy at full); a VCF region of 10,000 records decoded (heavy, `itWithin`); the
  allocation tier for a ribbon frame.

## Stages and lanes

Each lane: a brief in `briefs/`, its `_spec` and Configuration table kept, a changelog entry in its report, an
independent review and one fix round, the look in a visible pane in light and dark before the merge.

### Stage 0: designs (drafts, decided at Gate 0)

- [ ] GA (strong) genome coordinates: assemblies and chromosomes as pieces, a direction per piece, two assemblies on
      one axis, switching chromosomes, a navigator over several systems, scaled segments and the zoomed ring (axis.md,
      navigators.md, layouts.md amended).
- [ ] RR (strong) pairs and grids by region: the answer shapes, `RegionRequest.y`, 2D tiles, grid summaries, the pairs
      level and its pyramid (protocol.md, grids-and-rows.md amended).
- [ ] GF (strong) the genome formats: one runtime each, their entries and bounds, BGZF and tabix's surface if Gate 0
      puts them in `@molstar/io` (files.md, sources.md, base-library.md amended).
- [ ] LX (strong) layers across sections: what a layer between two sections is in the language, its layout and hit
      testing; synteny, liftover, contact maps (grids-and-rows.md, layouts.md amended).
- [ ] RN (strong) the ring extension: its package, kinds, params, the square navigator, Gate 0's answers 3 to 6
      (research-genome-ring.md as input; extension-openings.md amended).
- [ ] GD (fast, read-only) the Genomes dataset's generator and views: what each view needs from the designs above, the
      fixtures to record, the licences.
- [ ] GATE 0: the drafts agreed and promoted; the questions decided.

### Fix points (fast; `fixes.md`)

- [ ] FX1 before stage 1 (what the fourth slice's close leaves open); FX2 to FX4 after stages 1 to 3's gates.

### Stage 1: coordinates (strong, in order)

- [ ] GC1 the language and resolve: assemblies, pieces' direction, two assemblies on one axis (a language bump).
- [ ] GC2 the axis and layouts: reversed pieces drawn, scaled segments, the zoomed ring.
- [ ] GC3 navigators: switching chromosomes, a navigator over several systems (EC-184, EC-210).
- [ ] GATE.

### Stage 2: formats (beside stage 1; a formats track, fast or strong by format)

- [ ] BZ (strong) BGZF and tabix (in `@molstar/io` or Seq\*'s runtimes, by Gate 0).
- [ ] BD BED and bigBed; [ ] SQ 2bit and FASTA with its index; [ ] GG GFF3 and GTF (transcripts); [ ] VC VCF by
      region; [ ] PF PAF and chain files.
- [ ] GATE (with `pnpm build` and `check:workspace` if `@molstar/io` changed).

### Stage 3: by region and across sections (after stage 1)

- [ ] RG (strong) pairs and grids by region: the shapes, 2D tiles, summaries, the pairs level (EC-208, EC-223).
- [ ] LA (strong, after RG) layers across sections: synteny ribbons, the liftover, contact maps (EC-209).
- [ ] GATE.

### Stage 4: the dataset and the ring (after stages 2 and 3)

- [ ] GX (strong) the Genomes dataset generated and its views (examples.md, Genomes), the stories and twins.
- [ ] RX (strong) the ring extension package: `rondo:ribbons`, the bin-pair level, cytobands, outward stubs, the square
      navigator, its stories.
- [ ] USE GATE (not blocking): the author tries the genome views and the ring.
- [ ] GATE.

### Stage 5: use, measure, finalize

- [ ] NOVICE PASS: a task script per new view (find a gene, compare two assemblies, follow a contact on the ring).
- [ ] M5: the heavy lane and the opt-in tiers, a whole genome by region, the ring at full scale in the browser at dpr 2.
- [ ] FINALIZE: reconcile, fold back, the record, the next plan (version 1's close-out), worktrees cleaned.

## Order and parallelism

Three tracks: **coordinates** (GC1 → GC2 → GC3; the language one at a time), **formats** (BZ first, then the others
beside each other and beside the coordinates track: they touch runtimes, not the language), **views** (RG after GC1,
LA after RG, then GX and RX beside each other). The lane limit is `loop.md`'s.

## Hot files (the orchestrator's only)

As the fourth slice's, plus the new package's (`packages/seq/ext-rondo/package.json`, its tsconfig, its inventory row
and root reference) and `@molstar/io`'s `package.json` if BGZF and tabix land there.

## Known risks

- **A direction per piece changes what a position is** for paths, hover, selection and export: GC1 lands it with a
  twin per consumer, as ML2 did for units.
- **Region reads in two dimensions** may move the store's tile keys: RG runs alone in the store.
- **Format licences and sizes**: real slices kept small, each asked before it is fetched; generated fixtures first.
- **The default viewer is at its bound**: every runtime, the extension and the across-sections layout load lazily or
  as entries; a growth past a bound comes to the author with its feature.
- **Length**: about 25 lanes; Gate 0's cut order decides what goes first.

## Out of scope

Mol\* adapters and a live 3D view (the new mol viewer); BAM and CRAM (reads; a later milestone if a host asks);
GenBank; settings layers; agent tools; a React adapter; writers beyond SVG, PNG and archives; the compute port, unless a
profile of the whole genome shows the main thread is the limit (version-1.md, Deferred).
