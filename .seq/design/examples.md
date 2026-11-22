# Examples: few datasets, many views

## Intent

The author, 2026-11-17: "we could have few example datasets, but should have more views. i would like to keep the
richness, the editing experiences etc."

A handful of datasets, each rich enough to carry many views. A view is one spec over its dataset (or a snapshot of
one), showing one way of looking or one editing experience, and exercising named edge cases
([edge-cases.md](edge-cases.md)). New views cost a spec, not new data.

The exact data, decided 2026-11-18 ("your recommendations for gate 0 make sense generally", the author): small public
entries with permissive licences, fixtures under 1 MB each, the large variants generated in code. For each dataset
this design names the entries, their source and licence, the committed size, what is generated in code, and the views
each feeds, with what a user does in each view and its editing view. How fixtures are read and where views live stays
[examples-and-stories.md](examples-and-stories.md)'s.

## Precedents

### Data

- **The current Seq\*'s fixtures**: Pfam PF00042 (32 rows, HBA_HUMAN first, CC0), 1BRS barnase–barstar, 1IEP
  ABL–imatinib, a synthetic AFDB-shaped PAE, P04637's UniProt entry and Proteins API variants (CC BY 4.0) trimmed to
  what the readers take; recorded **by hand by a script that fetches and trims, never in tests**, each file's URL, date
  and licence in a README beside it. Taken: the record script, the trim to what is read, the provenance beside the file,
  a synthetic file said to be synthetic. Not taken: provenance as prose only (here it is data a test checks), fixtures
  shared by path across packages, a 353 KB variants file kept whole and of mixed licences.
- **This codebase's p53 views**: the protein is real, its domains illustrative, the gene and the conservation generated
  from a seed; a generated file that is also committed (`tp53-gene.seq`) is compared with its generator, and
  `SEQ_RECORDS=update` writes it again. Taken as the rule for every generated-and-committed file.
- **Jalview's example files** (`uniref50.fa`, the ferredoxin project with a structure and a tree): one family, its tree
  and a structure linked to its rows, in one bundle. Taken: a family with its tree and a structure that is one of its
  rows (Globins). Not taken: a project file of a tool's own.
- **Pfam and Rfam seeds** (Stockholm, CC0): curated, small, with `#=GC` lines. Pfam seeds have no lowercase inserts and
  usually no `RF` line, so they do not carry the insert-column case (EC-38); Rfam seeds carry `SS_cons`. Taken whole
  when under the bound.
- **SAbDab** (CC BY 4.0, its API's OpenAPI document): an entry's numbered segments, IMGT only since SAbDab2 (found by DA,
  2026-11-19: the Chothia files are gone). Taken for the wild type's numbering.
  Kabat is not among its files.
- **ANARCI** (BSD-3): numbers any sequence in Kabat, Chothia and IMGT, against HMMs built from IMGT germlines (IMGT's
  own terms). Taken only to make the generated design's test oracle, run by hand (the data's decisions, 3).
- **Addgene and SnapGene plasmid files**: the depositor's or SnapGene's terms, not a stated permissive licence; not
  taken, and nothing of their annotations either (coordinates included). The sequence comes from the INSDC record; its
  features are found in code.
- **GtRNAdb**: its licence is not stated; not taken. The tRNA comes from the PDB, its family from Rfam.
- **PDBe SIFTS** (EMBL-EBI's terms of use, `LicenseRef-EMBL-EBI-terms`, free reuse with attribution): taken for the
  services showcases' chain mappings (2OCJ, 1FIN), as published (SV, 2026-11-21).
- **RCSB/PDBe entries** (CC0) and **AlphaFold DB** (CC BY 4.0): taken, converted at recording to what the views read,
  never the coordinates whole.

### Interactions (what a user does with these data in the precedents)

| Precedent: what a user does | Here |
|---|---|
| Jalview: click a row's name to select it, Shift and Cmd to extend and add rows | planned, [grids-and-rows.md](grids-and-rows.md) (a row of the rows layer: select, extend, add) |
| Jalview: click the tree to cut it into groups, collapse a clade, sort by tree | planned, [grids-and-rows.md](grids-and-rows.md) (a click on a tree node folds its clade, session state; a summary row; tree order; `clades.collapseTo` by a cut); cutting the tree into colored groups: planned as a theme rule over a clade's rows ([themes.md](themes.md)): it colors, never hides |
| Jalview: select columns by the ruler, hide columns | select: built (a range on the axis); hiding columns: waits, no command this slice (decided) |
| Jalview: edit gaps (slide residues), realign | refused: data is not edited in the viewer (ux-foundations, the matrix); the Antibody design is a second sequence, not an edit |
| Jalview: color by Clustal, by conservation threshold | planned, [themes.md](themes.md) (residue schemes, rules) |
| Jalview: a structure linked to its rows (hover and select both ways) | built (the fourth slice, on a fake structure viewer: Features › Structure link); here the link is a map from the chains to their rows |
| ESPript: a figure of the alignment, secondary structure above, frames on conserved columns | planned, [layouts.md](layouts.md) (no interaction: a figure) |
| SnapGene, Benchling: click a feature to select it; select across the origin; double-click to zoom to it | planned, [layouts.md](layouts.md) (a feature across the origin by its shortest arc, EC-12) |
| SnapGene: rotate the map, set a new origin | planned, [layouts.md](layouts.md) (the linear view's start on a cyclic system, a section param, patched) |
| SnapGene: choose which enzymes are shown, cut-site lines | which enzymes: the view's builder input and a filter on the sites layer's data (the query sublanguage's filters in `layer.data`); "digest" and cloning simulations refused (an application's tool, not a view) |
| SnapGene: linear and circular maps linked, a window on the ring | planned, [layouts.md](layouts.md) (the linear detail follows the ring's window, EC-108, EC-114) |
| forna, R2DT: drag the 2D diagram's nodes | refused for v1: the 2D diagram is not one of the four layouts; arcs and the dot plot carry the pairs |
| ViennaRNA dot plots, AFDB's PAE plot: drag a box to select the residue pairs, which lights both in 3D | planned, [grids-and-rows.md](grids-and-rows.md) (the rectangle selects two ranges, EC-145) |
| abYsis, ANARCI: numbering tables in several schemes side by side | planned, the Antibody views (two rulers over one axis, Kabat and IMGT) |
| PDBe-KB, PLIP: the interface residues and a ligand's contacts listed and highlighted | planned, the Complex views (contacts as a `pairs` layer, hovered end to end; each ligand's contacts a ref the builder writes, highlighted when run) |
| AlphaFill's site: transplants listed with their scores, one shown at a time | planned, the Complex views (the ligands as a block, the scores as data; hiding by score a filter in `layer.data`) |
| Mol\*'s sequence panel: chains as blocks, select residues | built (selection), blocks planned ([layouts.md](layouts.md)) |
| The current Seq\*: clade collapse, lazy 10,000 rows with a cursor, CDR bookmarks, ring and linear detail, crosshair dot plot | each planned in its design ([grids-and-rows.md](grids-and-rows.md), [navigators.md](navigators.md), [layouts.md](layouts.md)); the selection's bookmark command, a bookmark's label and removal (`refs.remove`) planned |

## Decision

### The author's decisions on the data

The author, 2026-11-18: the data's open questions "ok as recommended"; in general: "a lot of this is "configurable
behavior" and where appropriate, the code should allow the configuration via spec ports etc."; on sizes: "don't overdo
it on the limits (they were sort of self-imposed and estimate), but obviously try to minimize the amount of code".
On the review of the decisions: "ok" to its recommendations.

1. AlphaFill's transplants are committed only if a CC BY 4.0 or looser licence is stated by recording time; else the
   ligands come from CC0 PDB entries of CDK2.
2. The PAE between CDK2 and cyclin A is generated from 1FIN's interface, said to be synthetic, unless a permissive
   prediction is found at recording.
3. ANARCI's labels-only output is committed as `LicenseRef-anarci-numbering`.
4. A source's attribution and licence are two optional params of a source node, carried into exported files (a
   language change).
5. `@molstar/io` is a devDependency of the examples package for the record script; [examples-and-stories.md](examples-and-stories.md) amended.
6. p53 variants are trimmed to UniProt's and ClinVar's.
7. The views follow the other designs: clade folds are session state with `clades.keep` (a click on a tree node folds);
   copies are `copies: merge` or `each`; hiding columns waits; the selection's bookmark command is planned beside `refs.remove`; an
   interaction with no default gesture goes through a mechanism (a `pairs` layer, a ref the builder writes) or is gone.
8. The author, on the remaining questions: "others as recommended": p53's generator moves to core's `PCG` with the rest, so
   `tp53-gene.seq` and the p53 records change once, in the Globins commit (`SEQ_RECORDS=update`, reviewed as a diff). A drag
   on a ruler: "i think drag on a rule should pan. or rather .. the original solution solved this by showing the rail
   as well ... maybe that's the way to go?": the Globins family view's ruler selects columns in its upper part, its
   lower 40% (the ruler rail) pans.

### Datasets and views

- **Six datasets in version 1, a seventh later**, each a fixture in `packages/seq/examples/fixtures` (small files, read through the source protocol like any
  other data; the large-scale variants generated in code, never stored):

  | Dataset | Holds |
  |---|---|
  | **p53** | human TP53: the protein (UniProt P04637) with features and variants; its gene on chromosome 17 (minus strand, 11 exons, introns), and the gene's region as a small FASTA with its index and a bigWig conservation signal, read by region; a structure's chains (the DNA-binding domain in four copies, author numbering, unobserved residues); an ortholog alignment (about 12 species) |
  | **Globins** | a globin family alignment (about 50 rows, with an insert column and a reference row) and its tree; hemoglobin (α2β2, two molecules, four units, heme ligands); a generated 10,000-row variant for scale |
  | **Antibody** | a Fab: heavy and light chains numbered Kabat and IMGT (insertion codes, IMGT out of index order), CDRs, a designed variant with indels against the wild type |
  | **Plasmid** | a cyclic construct (about 2.7 kb): features crossing the origin on both strands, a CDS across the origin, restriction sites |
  | **tRNA** | a tRNA with base pairs (pairs on one system, nested and crossing), modified bases, a covariance alignment |
  | **Complex** (the author, 2026-11-18: "a protein complex, with PPIs; and perhaps a ligand as well (maybe something from AlphaFill?)"; CDK2 with cyclin A confirmed the same day) | CDK2 with cyclin A (UniProt P24941, P20248): a crystal structure of the complex with ATP (two copies of the complex, so copy groups; unobserved loops); a predicted model of the complex with per-residue confidence and its predicted aligned error (PAE) between the chains; the interface (residues per side, residue–residue contacts across the chains); the ligands: ATP from the crystal and the ligands AlphaFill transplants into the AlphaFold DB model of CDK2 (each with its source entry and fit scores), their contacting residues; features (the PSTAIRE helix, the activation loop, the T160 phosphosite) |
  | **Genomes** (the fifth slice, [version-1.md](version-1.md); the author, 2026-11-19; reconciled 2026-11-21, [below](#genomes-the-fifth-slice-promoted-2026-11-21)) | two small toy assemblies with synteny (an inversion, a translocation, a gap), name aliases and an HLA-style contig name; their alignment and a liftover in the `seq` format; genes with transcripts (three levels, GFF3 and GTF), a name index; a bedGraph signal, a VCF of 40 samples with a long deletion, methylation with missing values; a Hi-C and a dot plot grid; 12,000 links over GRCh38 with real cytobands; FASTA and values read by region (the GWAS points, LD table and keyed results deferred to the close-out slice) |

- **Views per dataset** (the first list; each names the edge cases it carries):

  | Dataset | Views |
  |---|---|
  | p53 | text figure with feature underlines and variant marks (EC-94 to 99); tracks with features, values and secondary structure; the protein on its gene through codons, introns collapsed, minus strand (EC-15, 18, 29, 31, 117, 144), the gene reading along the protein (reversed) and the protein and its gene on one axis, the protein's letters a codon wide (EC-15 to 18, 28 to 31); structure chains linked, unobserved gray, copies merged and split (EC-65, 118); the gene's region by region with its conservation signal (EC-169 to 175, EC-201); the minus-strand gene's features in both directions; orthologs aligned, with hover across rows; queries and search with results stepping and bookmarks (EC-161, 162); a story in snapshots (EC-84, 85) |
  | Globins | family alignment with conservation logo and rows by tree (EC-38, 43, 106); collapsed clades with a summary (EC-100, 103); ESPript-style text figure (EC-97); hemoglobin complex with heme chips and copy groups (EC-64, 65); the 10,000-row alignment, lazy rows and keyboard cursor (EC-41, 44, 154) |
  | Antibody | numbering schemes side by side, Kabat vs IMGT labels (EC-2 to 8); design vs wild type with indel carets (EC-29); CDR selection and editing as bookmarks with undo (EC-162, 164) |
  | Plasmid | ring with a linear detail following its window (EC-108, 114), features as arrows by strand (the builder sets `arrow`, EC-232); a feature across the origin selected and focused by shortest arc (EC-9 to 14); the rotated linear view |
  | tRNA | arcs nested and crossing (EC-129); base-pair matrix as a dot plot with crosshair hover and rectangle selection (EC-109, 145); covariance alignment with pairs |
  | Complex | the chains as blocks with the interface residues and the contacts as arcs across the chains, folded by copy group (EC-35, 62, 65); the PAE as a dot plot over the same window as the tracks, crosshair hover lighting both residues, a dragged rectangle selecting two ranges (EC-109, 145); confidence and secondary structure tracks, unobserved residues gray (EC-118); the ligands as blocks after the polymers, ATP as one position and per heavy atom, AlphaFill's transplants listed with their scores, and the residues contacting a chosen ligand highlighted (EC-64, 101); a text figure with ligand chips flowing like words (EC-99); copies merged into one block labeled ×2 and split again (EC-65); the interface bookmarked as a ref, with undo (EC-162) |
  | Genomes (fifth slice) | a genome ring like Rondo's (rondo.app): chromosomes as blocks with their cytobands, crosslinks bundled into ribbons wide by their count, one chromosome as its own ring, built as an extension package (the author, 2026-11-19: "Rondo-like, extension"; research in the third slice's `research-genome-ring.md`); synteny ribbons with an inversion (direction on segments, EC-206); a dot plot of one genome against the other; both genomes on one ring (EC-210); signal by region with an open domain and with a one-sided domain (EC-169 to 175, EC-201, EC-202); signal by region on a log scale (EC-215); samples' variants as rows with the deletion solid (EC-222); the genotype matrix; methylation matrix by region (EC-220); transcripts in three levels (EC-227, EC-229); a contact map; open on a gene by name (EC-228); a liftover over a map in the `seq` format (EC-212; cut last); the genome by region from FASTA with its index, names aliased (EC-186 to EC-189); deferred to the close-out slice: many signals from one file, the Manhattan family, transcripts by joined results (EC-214, EC-216, EC-217, EC-219, EC-224) |

- **Showcases at the application layer**, in the storybook package (they compose what an application would, and are
  not core features):
  - **an application composing services**: UniProt, PDBe (SIFTS) and AlphaFold DB as runtimes, from recorded
    fixtures, a switch to go live;
  - **predictions onto a structure**: a lab's per-position predictions on a sequence indexed from 1, with no numbering
    of its own, mapped onto an mmCIF structure through `label` and `auth` numbering (an offset, insertion codes,
    unobserved loops, an expression tag), colored on the sequence and in a fake 3D viewer (the Mol\* view waits for the new mol viewer, decided 2026-11-19), hover and selection both ways
    (EC-3 to EC-5, EC-25, EC-118);
  - **addressing**: a text with links that focus, highlight and select in the view and in the (fake) structure viewer; "copy link"
    from the view back to text; the same location in several numberings; rects from a dot plot and an alignment; refs
    as bookmarks ([addressing.md](addressing.md)).
- **Editing experiences are views too**: each dataset has at least one view that is edited (a bookmark added, a track
  hidden or reordered, a theme changed, a snapshot recorded), with undo and redo, a preview that ends in one patch or a
  cancel, and the printed spec updating beside it.
- **The previous implementation's cases** each have a home (its 22 cases, [sources/seqstar](../sources/seqstar/README.md)):

  | Its case | Here |
  |---|---|
  | p53 text and tracks; UniProt features; queries and bookmarks; protein on its genome through codons; reverse-strand gene; streamed chromosome | p53 (the streamed chromosome: p53's chromosome) |
  | protein family; PF00042 alignment; family triage | Globins |
  | antibody numbering; indels in wild-type numbering | Antibody |
  | DNA construct; plasmid ring | Plasmid |
  | RNA base pairs | tRNA |
  | imatinib in ABL; barnase–barstar; AFDB complex; designed dimer | Complex |
  | protein over coding DNA | p53 (the protein on its gene) |
  | genotype matrix | Genomes (fifth slice) |
  | story drafts | p53 (a story in snapshots) |

- **Every view runs three ways**: in Storybook (canvas and SVG side by side), in a quick-lane test (validates,
  resolves, lays out with no unexpected issues, renders to SVG; its printed spec and draw-list record kept), and as a
  spec file a reader without Seq\* code could take (the JSON and its data files).

### One rule for every fixture

- **Real entries, read as published, converted only where no reader exists.** A fixture is in its own format when the
  core reads that format (FASTA, Stockholm, Newick: the Stockholm and Newick readers of [grids-and-rows.md](grids-and-rows.md); TSV, CSV, JSON,
  MessagePack and GFF3 through the records runtime of [files.md](files.md), agreed 2026-11-19; a JSON fixture trimmed,
  never reshaped, so its paths read the live response too), else in `seq` (the
  protocol's shape) converted at recording (mmCIF, GenBank, AlphaFold DB's JSON, AlphaFill's metadata, SAbDab's
  numbered files). Why: a reader without Seq\* can take a published file; a format the core will not read is not a
  reason to commit megabytes it would never open.
- **Recorded by hand, never in tests**: `packages/seq/examples/scripts/record.ts` fetches each entry, trims it to what
  the views read, converts it, and writes the fixture and its provenance. Tests never reach the network
  ([security.md](security.md), the test adapter). Built with the first new dataset (Globins), with the `sources.json`
  test and the seeded generator on core's `PCG` (below); the later datasets add their entries to it. Why: the precedent's rule;
  recording is the only step that trusts a remote host.
- **A real slice of a large indexed file is recorded as pages** ([genome-formats](genome-formats.md), Real slices to
  record; the fifth slice): core's 64 KiB pages the runtime read, keyed by file and page, so the container written by
  the real tools is what the test reads; a read the recording did not make fails naming its page, and a change to a
  runtime's reads is a re-recording. Why: a trimmed file would test the trimming, not the BGZF blocks and the index.
- **Provenance is data**: `fixtures/<dataset>/sources.json` lists every file of the folder: the entry and its version
  (PDB id, Pfam and Rfam release, AFDB model version, INSDC accession.version), the URL, the date recorded, the SHA-256
  of the download and of the fixture, the licence as an SPDX id or a `LicenseRef-` id, the attribution line, what was
  trimmed, and which parts are **derived** (computed by the script or in code, not in the source: annotations found in
  a sequence, contacts, a converted numbering). Why: a test can check it; prose cannot.
- **Licences allowed**: `CC0-1.0` (PDB, Pfam, Rfam), `CC-BY-4.0` (UniProt, AlphaFold DB, SAbDab, and AlphaFill if
  a licence is stated: the data's decisions, 1), `LicenseRef-NCBI-public-data` (INSDC records and ClinVar, on which NCBI places no restriction), and
  `LicenseRef-seqx-generated` (what the package generates: its own, under the repository's licence). Anything else is
  refused at recording. A CC BY dataset's attribution is shown on its Storybook index and carried with an exported spec
  file (the data's decisions, 4).
- **Under 1 MB a file** (decided, 2026-11-18), **under 2 MB a dataset**: checked by a test. Why: the repository carries them and the
  Storybook ships them as assets.
- **Generated, never stored**: every large variant and every computed table (a tree when none is published, contacts
  when they are cheap to compute from what is committed, covariation, restriction sites, sequence features found by
  search, a design) is made in code from a seed (`seeded`, a shared module of the package on core's `PCG`,
  `@molstar/core` `data/util/hash-functions`; no private generator). A generated file that a view reads by location (a `seq` source) is committed and compared with its generator,
  as `tp53-gene.seq` is; `SEQ_RECORDS=update` writes it again.
- **Facts are checked against the data** ("What is drawn is checked against the data"): each dataset has anchor tests
  stating known facts of its entries (a length, a conserved residue, a pair, a feature's span), so a recording that
  changed the data fails a test, not a picture. Facts marked *to verify when recording* below are checked by the record
  script's author before the anchors are written, and corrected here if wrong.
- **Every view's builder takes its dataset as an argument** (the fixture by default), and **each view runs over the
  instance that carries its edge case**: where the published entry does not carry a case (EC-38 in a Pfam seed, EC-7 in
  trastuzumab's IMGT numbering, EC-11 at pUC19c's own origin), the view's story and test run over the generated variant
  that does, and say so. Why: the self-check (below), and an edge case claimed by a view that never meets it is not
  kept.

### The datasets

Sizes are expectations to be confirmed by the record script (it prints each); the bound is the test's.

#### p53 (kept; the structure entry added)

| What | Source | Licence | Size | Form |
|---|---|---|---|---|
| the protein, domains, conservation, the toy gene and its region | as built (the protein real, the rest generated) | P04637 CC-BY-4.0; the rest `LicenseRef-seqx-generated` | 18 KB committed | `seq`, generated |
| the structure: **PDB 2OCJ**, the p53 core domain without DNA, four chains (A–D) at 2.05 Å (Wang, Rosengarth, Luecke 2007) | RCSB/PDBe mmCIF | CC0-1.0 | about 10 KB | `seq`: each chain's `label` and `auth` numbering (an offset onto P04637), its observed residues, its secondary structure (`struct_conf`, `struct_sheet_range`), the zinc ligands as one position each |
| variants on the tracks | the Proteins API, `variation/P04637`, **trimmed to the sources with known licences**: UniProt's own (CC-BY-4.0) and ClinVar's (`LicenseRef-NCBI-public-data`); COSMIC, TCGA and the other large-scale sources dropped (the data's decisions, 6) | per source, recorded per row's source in `sources.json` | under 150 KB | `seq` ranges: position, letters, consequence, clinical significance, source; single-residue changes inside 1–393 |
| secondary structure per residue for the protein's track | UniProt P04637's Helix, Beta strand and Turn features | CC-BY-4.0 | about 5 KB | `seq` ranges |

Why 2OCJ: four copies of the DNA-binding domain with no DNA chains, as the views above ask (copies merged and split,
EC-65; unobserved ends gray, EC-118; author numbering, EC-3). Generated: nothing new. Anchors: each chain's `auth`
numbering equals P04637's; R175, R248, R273 inside every chain's observed range; one zinc a chain (C176, H179, C238,
C242).

Views it feeds (the existing views keep theirs): structure chains linked, unobserved gray, copies merged and split;
secondary structure and variants on the tracks (the small layers, [themes.md](themes.md)). Interactions: hover a chain's residue lights it on the protein
and on the other three copies (copy rule, EC-25); click selects, projected to the protein; double-click focuses;
merging the copies (`copies: merge`, the layouts design's default) and splitting them (`each`) is the editing view's
toggle (a patch, one undo step each);
hiding variants by source or consequence is a filter in `layer.data` (a patch).

#### p53's chromosome (added 2026-11-19)

The author, 2026-11-19: "we should also add a "genome-level" example(s)"; "let's do chromosome scale this slice; and
the whole genomes in appropriate slice" (the Genomes dataset: the fifth slice, [version-1.md](version-1.md)).

Chromosome 17 at its GRCh38 length (83,257,441 bp, a fact, not data), every byte generated in code by region and never
stored, with TP53 at its real place (7,668,402–7,687,550, minus strand) carrying the p53 gene the dataset already has.
It is performance.md's "a whole human genome by region" at one chromosome, and the previous implementation's
"streamed chromosome" case.

| What | Source | Licence | Size | Form |
|---|---|---|---|---|
| the chromosome's letters | generated by region: core's `PCG` seeded by the block's index, so a region's letters never depend on what was asked before; TP53's span from `tp53-gene.seq` | `LicenseRef-seqx-generated` | 0 committed | a source by region (the source protocol) |
| genes | generated by region from one seed (about 1,200, both strands, exons), TP53 placed as built | as above | 0 | features by region, each gene a dynamic node with its exons |
| a signal (GC-like, smooth with noise) | generated by region, summarized by bin-width levels | as above | 0 | values by region and level |

Views it feeds: **the chromosome** (a ruler in Mb, the genes, the signal, the letters when they fit; opened fit to the
whole chromosome) with an overview navigator showing the window; zooming from the whole chromosome to TP53 (a
double-click on the gene focuses it) and down to its letters, every level drawn from the store as regions arrive;
TP53's gene opened as the p53 gene view (the dataset's existing gene, linked by its place). Interactions: wheel and drag
zoom and pan, the navigator's window dragged, hover on a gene (its name, strand and span in the hover card), click
selects a gene, 0 fits. Not here (the Genomes dataset's): open by name (EC-228), gene hierarchies (EC-227), a second
assembly.

What it checks: a frame costs what is on screen at 83 Mb as at 300 residues (performance.md's flat frame time,
measured at the fit view, at a density between two levels and at the letters); regions requested only for the view and
its margin, evicted past the store's bound, never the view (EC-169 to EC-175); the ruler's Mb, kb and bp labels
(`@molstar/core/util/number`).

#### Globins

| What | Source | Licence | Size | Form |
|---|---|---|---|---|
| the family: **Pfam PF00042 (Globin) seed** at the current release | InterPro's Pfam seed, Stockholm | CC0-1.0 | 15–30 KB | Stockholm, as published; a seed over 60 rows is cut to 60 by the record script (HBA_HUMAN, HBB_HUMAN and the rows nearest them by identity kept, the cut recorded) |
| its tree | InterPro's seed tree when served (*to verify when recording*); else computed in code (neighbor joining on identity) and committed, compared with its generator | CC0-1.0 / `LicenseRef-seqx-generated` | about 3 KB | Newick |
| hemoglobin: **PDB 4HHB**, human deoxyhemoglobin (Fermi, Perutz, Shaanan, Fourme 1984): α chains A and C, β chains B and D, four hemes | RCSB/PDBe mmCIF | CC0-1.0 | about 15 KB | `seq`: two molecules (HBA, HBB), four units, copy groups; each chain's helices (A–H, from the entry's `struct_conf`), for the ESPript figure's secondary structure; each heme one position and per heavy atom (EC-64); each heme's contacting residues (heavy atoms within 4 Å, derived) |
| the 10,000-row variant | generated | `LicenseRef-seqx-generated` | none | in code |

- **Reference row**: HBA_HUMAN (the α chain's row); the units of 4HHB link to their rows by an alignment map. A chain
  whose row is not in the seed links by a pairwise alignment computed in code with core's `algorithms/sequence`,
  reported as computed. *To verify when recording*: HBA_HUMAN and HBB_HUMAN are rows of the current seed.
- **Generated, 10,000 × 300** (performance.md, Scale): rows sampled from the seed's per-column profile (each column's
  residue frequencies and gap rate), seeded; the seed's columns first, then columns sampled from its profile to reach
  300; lowercase inserts with an `RF` line every 25 columns: **this variant, not the seed, carries the insert-column
  case** (EC-38), so the family view's insert story runs over it; row names `GEN_<n>/<start>-<end>` with start numbers
  (EC-39); a balanced tree over the rows from the same seed (for collapsed clades at scale). Rows are produced on demand
  by index, never held all at once (EC-41; performance.md, lazy rows): a row is a function of the seed and its index.
- **Numbering in the anchors**: Pfam rows number by UniProt, and P69905 (HBA) and P68871 (HBB) count the initiator Met,
  which the mature chains of 4HHB do not: UniProt = 4HHB's `auth` + 1. So the distal and proximal histidines are His58
  and His87 in 4HHB's α chain, His59 and His88 on the HBA_HUMAN row (and β's His63 and His92 in 4HHB, His64 and His93
  on HBB_HUMAN). The offset is itself an anchor of the alignment map (EC-3, EC-39).
- **Anchors**: as above; every heme contacts its chain's proximal histidine; 4HHB's chains have 141 (α) and 146 (β)
  residues; each chain has eight helices, A to H (*to verify when recording*: the entry's `struct_conf` count).

Views:

| View | Instance | What a user does | Editing view |
|---|---|---|---|
| family alignment, conservation logo, rows by tree (EC-43, EC-106) | the seed; EC-38's inserts on the generated variant | hover a cell: the column lights in every row and in the logo; click a row's header: the row selected; Shift and Cmd on headers: extend and add rows; drag on the ruler's upper part: a column range (its lower 40% is the ruler rail: pans); Cmd/Ctrl+click a tree node: its clade's rows selected | change the residue scheme (a theme patch); hide the logo track |
| collapsed clades with a summary (EC-100, EC-103) | the seed | click a tree node: its clade folds, its summary row shown (session state, no undo step); click it again: unfolded; hover the summary: its columns' consensus; a focus on a row inside opens its clade | `clades.keep` writes the folds into the spec (one patch); pin a row above the tree's order (a patch) |
| ESPript-style text figure (EC-97, EC-100) | the seed with 4HHB's helices on the HBA row | select a column range in text: the same range in the tracks view if linked; frames are drawn, not acted on | none of its own: the figure follows the family view's edits |
| hemoglobin with heme chips and copy groups (EC-64, EC-65) | 4HHB | hover α chain A: chain C lights (copy group) and its row in the family; hover a heme's contact (a `pairs` layer of the derived contacts): the heme and its residue lit; run the "heme pocket" ref (written by the builder: the contacting residues) from the bookmarks list: they are highlighted; merge the copies to `×2` and split again | merge/split copies (`copies`, a patch); bookmark the heme pocket from the selection (the selection's bookmark command) |
| the 10,000-row alignment, lazy rows, keyboard cursor (EC-41, EC-44, EC-154) | generated | arrows move the cursor across rows never held; Page keys by a screen of rows; `]` and `[` zoom the rows; drag the vertical overview's window ([navigators.md](navigators.md)); Enter focuses the cursor's residue; Space on a header folds the row's clade (session state) | pin a row found by the cursor (a patch) |

#### Antibody

| What | Source | Licence | Size | Form |
|---|---|---|---|---|
| the Fab: **PDB 1N8Z**, trastuzumab (Herceptin) Fab with the HER2 extracellular domain (Cho et al. 2003): light chain A, heavy chain B; HER2 (chain C) not taken | RCSB/PDBe mmCIF | CC0-1.0 | about 5 KB | FASTA of the two chains as in the entry |
| IMGT numbering of both variable domains | SAbDab's numbered segments of 1N8Z (889, 890) | CC-BY-4.0 | about 6 KB | `seq` numbering tables: labels per residue, insertion codes |
| Kabat numbering | derived in code from IMGT (SAbDab serves no Chothia): residues counted between the conserved anchors (L23, L35, L88, L98; H22, H36, H92, H103), insertions at Kabat's places (L27, L95; H35, H52, H82A–C, H100); a framework whose IMGT occupancy differs is refused, not guessed; committed as derived, compared with its generator and with trastuzumab's published Kabat numbering | CC-BY-4.0 (from SAbDab's) | about 3 KB | `seq` numbering table |
| CDRs by both definitions | in code (Kabat's and IMGT's CDR ranges are definitions, not data) | `LicenseRef-seqx-generated` | none | in code |
| the design | generated: substitutions in CDR-H2 and CDR-H3, a two-residue insertion in CDR-H3, a one-residue deletion in CDR-L1 | `LicenseRef-seqx-generated` | none | in code |
| the design's numbering oracle (not yet: ANARCI was not installed when DA recorded, 2026-11-19; the design's two rules, Kabat dropping 27 and IMGT dropping 36, are marked unchecked until it is) | ANARCI run by hand on the design, its output (labels only, both schemes) committed for the tests to compare with the generator's numbering; tests cannot run ANARCI | `LicenseRef-anarci-numbering` (the data's decisions, 3) | about 6 KB | `seq` numbering tables |

- **Why trastuzumab**: one of the best known Fabs, its heavy chain carries Kabat insertion codes (expected 52A, 82A–C
  and 100A–C, *to verify when recording*), its sequence is in a CC0 entry and its numbering in SAbDab.
- **What the wild type does not carry**: trastuzumab's IMGT CDR-H3 is 13 residues (105–117), so it has no IMGT CDR3
  insertions; **IMGT's out-of-order labels (111.1, 112.1: EC-7) are on the design only**, whose two-residue H3
  insertion makes them. The numbering views run over the wild type for EC-2 to EC-6 and EC-8 and over the design for
  EC-7.
- **Why the design is generated**: no public design with indels against this wild type has a permissive licence (the
  licence rule decides it). The generator numbers its insertions by the schemes' rules (Kabat: letters after 100; IMGT:
  symmetric about 111 and 112), and the committed ANARCI output is the oracle it is tested against.
- **Anchors**: Kabat L23 and L88, H22 and H92 are cysteines; IMGT 23 and 104 are cysteines, 41 a tryptophan, 118 a
  phenylalanine or tryptophan, in both chains. *To verify when recording*: 1N8Z's chain lengths (214 and 220 as
  deposited).

Views:

| View | Instance | What a user does | Editing view |
|---|---|---|---|
| numbering schemes side by side, Kabat and IMGT rulers (EC-2 to EC-8) | the wild type; the design for EC-7 | hover a residue: both labels shown; a query by label in either scheme (`(at :label "100A")`) selects; a range written backwards in IMGT is an issue with its hint (EC-7) | switch which scheme the ruler shows (a patch of the ruler layer's `numbering` param, the themes design's; insertion labels 100A, 111.1 as minor ticks) |
| design against wild type, indel carets (EC-29) | both | hover an inserted residue: a caret at the next pair on the wild type; select a CDR on the design: the matching range on the wild type | none of its own |
| CDR bookmarks with undo (EC-162, EC-164) | the wild type | select a CDR, bookmark it (the selection's bookmark command), rename, recolor, remove (`refs.remove`); undo and redo each | the view is the editing view |

#### Plasmid

| What | Source | Licence | Size | Form |
|---|---|---|---|---|
| **pUC19c, INSDC L09137.2** (2,686 bp, circular; the record's own title); its only feature is `source` | ENA/NCBI | `LicenseRef-NCBI-public-data` | about 4 KB | `seq`: the cyclic sequence and its name |
| its features, **derived in code**, none taken from another annotation: bla as the open reading frame found in the sequence (on either strand, the one encoding TEM-1, MSIQHFRVAL…); lacZα likewise; the multiple cloning site as the span from its first to its last site (EcoRI to HindIII); the pMB1 origin located by a short signature written in the module, cited to the literature (*to verify when recording*) | in code | `LicenseRef-seqx-generated` | none | in code; bla's codon map with it |
| restriction sites (the MCS's enzymes and a few unique cutters, a list in the module) | computed in code on the cyclic sequence (a site straddling the origin found) | `LicenseRef-seqx-generated` | none | in code |
| the rotated variant: the origin moved inside bla, so the CDS crosses it, **and a generated plus-strand feature across the same origin** (a primer-binding site, named as generated), so one record carries the crossing on both strands (EC-11) | generated | `LicenseRef-seqx-generated` | none | in code |
| the reverse complement (pUC19's usual orientation) | generated | `LicenseRef-seqx-generated` | none | in code, for the property test |
| the ring's scale variant: a 200 kb cyclic sequence with 2,000 features on both strands, a tenth across the origin | generated | `LicenseRef-seqx-generated` | none | in code |

- **Why pUC19c**: the commonest cloning vector, about 2.7 kb as the dataset table asks, an INSDC record with no licence
  question. The record has no feature annotations (bla appears only in a comment, with other coordinates), and the
  commercial maps' annotations are not ours to take; finding the features in the sequence is what plasmid tools do, and
  it is checked. It has no feature across its published origin, so the crossing is made by rotation, the operation a
  user performs in plasmid tools.
- **Anchors**: bla is 861 nt, 286 residues and a stop, starting MSIQHFRVAL; EcoRI and HindIII each cut once, in the
  MCS; every derived feature keeps its length under every rotation (a property test over all 2,686 offsets, on the
  record and its reverse complement).

Views:

| View | Instance | What a user does | Editing view |
|---|---|---|---|
| ring with a linear detail following its window (EC-108, EC-114) | the rotated variant | drag the ring's window (or its handles) to move the linear detail; wheel on the linear detail moves the window (the follower moves its leader); click a feature on either: selected on both | hide a strand's track; choose the enzymes shown (a filter in the sites layer's `layer.data`, a patch) |
| features across the origin selected and focused by the shortest arc (EC-9 to EC-14) | the rotated variant (both strands) | click bla: one item, two runs; double-click: focus by the shortest arc, not the whole ring (EC-12); Shift+drag across the origin extends without wrapping the long way | bookmark bla from its start to its stop |
| the rotated linear view | the published record | set the linear view's start (a section param: a patch, undone) | the rotation is the edit |

#### tRNA

| What | Source | Licence | Size | Form |
|---|---|---|---|---|
| the molecule: **PDB 1EHZ**, yeast tRNA-Phe at 1.93 Å (Shi, Moore 2000), chain A, 76 nt, its modified residues (2MG, H2U, M2G, OMC, OMG, YYG, PSU, 5MC, 7MG, 5MU, 1MA) | RCSB/PDBe mmCIF | CC0-1.0 | about 6 KB | `seq`: the sequence with each modified residue's parent letter and its code; pairs with kind (Watson–Crick, wobble, non-canonical) and whether tertiary |
| its base pairs | the entry's base-pair category (`ndb_struct_na_base_pair`) when present (*to verify when recording*); else a curated table of the published pairs, cited | CC0-1.0 / facts | in the above | `seq` pairs |
| the family: **Rfam RF00005 (tRNA) seed** at the current release | Rfam, Stockholm with `SS_cons` | CC0-1.0 | 100–300 KB (about a thousand rows, *to verify when recording*) | Stockholm, as published |
| covariation on 1EHZ's positions | computed in code directly on 1EHZ's positions: for each pair of 1EHZ positions, mutual information over the seed's columns they align to (through 1EHZ's computed alignment map), with its correction for small counts; a grid layer draws only on its own systems (grids-and-rows.md), so the grid is made on 1EHZ's system, not projected | `LicenseRef-seqx-generated` | none | in code |
| 1EHZ's place in the family | computed in code: a pairwise alignment of 1EHZ's parent letters to the seed row of highest identity (core's `algorithms/sequence`) | `LicenseRef-seqx-generated` | none | in code |
| the dot plot's scale variant: a 4,000 × 4,000 grid (performance.md, Scale) | generated | `LicenseRef-seqx-generated` | none | in code |

- **Why 1EHZ**: the reference tRNA structure, high resolution, modified bases, and tertiary pairs that cross the
  cloverleaf's (EC-129: nested and crossing on one system).
- **The second instance is real and already committed**: every row of the Rfam seed, with its pairs projected from
  `SS_cons`, is a tRNA with a secondary structure; the arcs and dot plot views run over a seed row besides 1EHZ (nested
  pairs only, no modified bases). No second structure is taken: a second tRNA's tertiary pairs would add nothing a view
  draws differently, and a download and a licence each.
- **The dot plot** (the base-pair matrix): one grid over 1EHZ's positions on both axes: below the diagonal the
  structure's pairs (a value per kind), above it the family's covariation computed on 1EHZ's positions (above).
  Why: one figure shows the structure and the evidence for it, both from committed data.
- **Anchors**: acceptor stem 1–7 : 66–72, D stem 10–13 : 22–25, anticodon stem 27–30 : 40–43 (the entry's base-pair table, `ndb_struct_na_base_pair`, has no A31·Ψ39 and has
  Cm32·A38, non-canonical; kept as published, found by DT 2026-11-19), T stem 49–53 : 61–65;
  tertiary 8–14, 15–48 (Levitt), 18–55, 19–56, 26–44, 54–58; anticodon Gm34-A35-A36; the pairs from `SS_cons`
  projected onto 1EHZ land on its stems.

Views:

| View | Instance | What a user does | Editing view |
|---|---|---|---|
| arcs nested and crossing (EC-129) | 1EHZ; a seed row (nested only) | hover an arc: both ends lit, and the pair's cell in the dot plot; click: the pair selected (an item); pan: arcs move continuously, an off-window partner reached past the edge | hide tertiary pairs (a filter by kind in `layer.data`, a patch) |
| base-pair matrix as a dot plot, crosshair hover, rectangle selection (EC-109, EC-145) | 1EHZ; the generated grid for its frame target | hover a cell: a crosshair, both positions lit on the arcs and the sequence; drag a rectangle: two ranges selected; wheel zooms both axes about the pointer | change the covariation's palette (a theme patch) |
| covariance alignment with pairs | the seed | hover a column: its partner column lit by the `SS_cons` arc; hover a 1EHZ residue in the linked sequence section: its cell lights on the seed row of highest identity (the computed alignment map binds 1EHZ to that row) | bookmark the anticodon stem from the selection (the selection's bookmark command) |

#### Complex (CDK2–cyclin A)

| What | Source | Licence | Size | Form |
|---|---|---|---|---|
| the crystal: **PDB 1FIN**, cyclin A–CDK2–ATP at 2.3 Å (Jeffrey et al. 1995): CDK2 chains A and C, cyclin A chains B and D (a construct of P20248), an ATP on each CDK2 | RCSB/PDBe mmCIF | CC0-1.0 | about 40 KB | `seq`: two molecules, four units, two copies of the complex (copy groups); observed residues; ATP one position and per heavy atom; contacts within 4 Å across the chains per copy (the interface, derived), and ATP's contacting residues (derived) |
| features: the PSTAIRE helix (45–56), the activation loop, T160; cyclin A's cyclin box | UniProt P24941, P20248, trimmed | CC-BY-4.0 | about 5 KB | `seq` ranges |
| per-residue confidence and the PAE within each chain: **AF-P24941-F1** (CDK2) and **AF-P20248-F1** (cyclin A2) | AlphaFold DB, the model version current at recording | CC-BY-4.0 | about 400 KB (the two PAE blocks quantized to 0.5 Å, one byte a cell) | `seq`: values and grids |
| the PAE between the chains | generated (the data's decisions, 2) from 1FIN's interface (low across it, high elsewhere, said to be synthetic in its name and in `sources.json`) unless a CC BY prediction of the complex is found | `LicenseRef-seqx-generated` | none | in code |
| AlphaFill's transplants into **AF-P24941-F1** (the AlphaFill entry of the same id): ATP and up to nine other ligands, the best by AlphaFill's transplant clash score, each with its source PDB entry, global and local RMSD and score, and its contacting residues | alphafill.eu, the entry's metadata and model | CC BY 4.0 or looser if stated (the data's decisions, 1) | about 20 KB | `seq`: ligands as positions, scores as data, contacts as pairs |
| the copies' scale variant: 200 copies of CDK2 (performance.md, Scale) | generated | `LicenseRef-seqx-generated` | none | in code |

- **Why 1FIN**: the author's complex, with ATP and two copies of the complex in the asymmetric unit (copy groups,
  contacts folded by copy group, EC-62, EC-65). Unphosphorylated T160: the phosphosite is a feature, not a modified
  residue.
- **The crystal and the prediction differ in extent** (cyclin A's construct against the full chain): an offset map
  between them, its pieces from the entry's `auth` numbering (EC-3).
- **Anchors**: K33, E51 (PSTAIRE) and D145 (DFG) in CDK2; ATP contacts K33 and the hinge (E81, L83) in each copy;
  T160 inside the activation loop; the interface involves the PSTAIRE helix in both copies.

Views:

| View | Instance | What a user does | Editing view |
|---|---|---|---|
| the chains as blocks, interface residues, contacts as arcs across the chains, folded by copy group (EC-35, EC-62, EC-65) | 1FIN | hover a contact (a `pairs` layer): both residues lit, in the copy hovered and its twin; click a contact: the pair selected (both ends) | bookmark the interface as a ref, undo (EC-162) |
| the PAE as a dot plot over the same window as the tracks (EC-109, EC-145): one track section a chain; a grid section whose x follows CDK2's section and y cyclin A's (grids-and-rows.md) | the predictions | crosshair hover lights both residues on the tracks; drag a rectangle: two ranges selected, one on each chain | none of its own |
| confidence and secondary structure tracks, unobserved residues gray (EC-118) | both | hover: the value; select a low-confidence range by dragging | change the confidence palette's bins (a theme patch) |
| ligands as blocks after the polymers (EC-64, EC-101) | 1FIN, AlphaFill | hover a ligand's contact arc: the ligand and its residue lit; run a ligand's "contacts" ref (written by the builder) from the bookmarks list: its contacting residues highlighted; ATP as one position or per heavy atom (this slice: two sequences in the fixture with no map between them, a hover not crossing; the fourth slice joins them, [molecules](molecules.md) gap 14; a patch switches); AlphaFill's scores on hover | hide transplants under a score (a filter in `layer.data`, a patch) |
| text figure with ligand chips (EC-99) | 1FIN | select across a chip; a chip's hover as a ligand's | none of its own |
| copies merged to `×2` and split (EC-65) | 1FIN; the 200 copies for the scale | hover the merged block: every member lit | merge, split: one patch each |

#### Genomes (the fifth slice; promoted 2026-11-21)

Lane GD's design (the fifth slice, stage 0; agreed by the author 2026-11-21, promoted from its draft in
`execution/fifth-slice/drafts/genomes-dataset.md`), on the five genome designs: [genome-coordinates](genome-coordinates.md)
(GA), [region-reads](region-reads.md) (RR), [genome-formats](genome-formats.md) (GF),
[across-sections](across-sections.md) (LX), [ring-extension](ring-extension.md) (RN), and the decisions of 2026-11-21
(three formats built: BED, GFF, VCF; the alignment views read the `seq` format; the Manhattan family, `split` and the
keyed-results join deferred to the close-out slice; MAF cut). A design owner is named per need.

One dataset, `genomes` in `packages/seq/examples`, that every view of the milestone reads: committed fixtures a few
hundred kB each (the examples' `fixtures.maxFileBytes` 1 MB and `maxDatasetBytes` 2 MB, kept), generated by a seeded
module and written by the record script; one real slice for each built format's conformance suite (Gate 0, 7); and a
performance fixture that is never stored (a whole human genome by region, heavy). The views are builder modules as in
the other datasets; each view's test drives the cells of the interaction matrix it names.

##### The views and what each needs

Cut order (Gate 0, 1): contact maps, then the liftover (MAF is cut).

| View | Data (generated unless marked) | Read through | Designs' pieces | Edge cases |
|---|---|---|---|---|
| **Genome ring** (Rondo-like: blocks with cytobands, bundled ribbons wide by count, one chromosome or a pair as its own ring, the square navigator, outward stubs) | GRCh38's 23 chromosomes typed as `assembly.sequences`; 12,000 seeded links (60 % within a chromosome with distance decay, 40 % between, a dozen hub pairs), held and as BEDPE (one picture); **cytobands: real** (UCSC `cytoBandIdeo`, fetch 1), `syntheticBands` (RN) in tests | `runtimes/bed`: BEDPE read whole and indexed by both ends, BED with `fields: { stain: '4' }` | GA `assembly`, `system` field, `section.include`, `navigator.switch`, the square; RR pairs level, `binPairs`, `coarsenBins`, `PairsReader`; RN `rondo:ribbons`, stems, stubs | 108, 169, 171, 206, 210, 218, 225, 235 |
| **Synteny ribbons** (toy A over toy B: an inversion drawn reversed, a translocation, a gap) | the two toy assemblies; their alignment as sided pairs of equal-length runs, written in the `seq` format from one alignment truth in code | the `seq` runtime (sided, either end, GF) | GA direction per block, one map over sequence pairs; LX `bridge`, ribbons coloured by the `a` side's block; RR sides | 206, 207, 208, 209, 218 |
| **Both genomes on one ring** | the toy assemblies, one section `of` both | as above | GA two assemblies on one section, group gap; RR pairs on a ring; RN stubs | 210 |
| **Dot plot of one genome against the other** | the alignment as pairs, and as grid tiles (binned from it, two levels) | the `seq` runtime: pairs, and a grid (as the Hi-C example, region-reads, Decided 10) | LX grid over two sections of blocks (`follow_y` a second section); RR 2D tiles, grid summaries; GA's square | 109, 145, 209, 223 |
| **Contact map** (Hi-C over the toy genome, symmetric) | the links binned into two levels, and a one-level twin for the summaries | the `seq` runtime: grid | RR 2D tiles, symmetric, whole-assembly level; LX between two chromosome sections | 223 |
| **Signal by region** (open domain, one-sided domain), log scale | one bedGraph over both toy genomes, 100 bp bins, a TPM-like column 0 to 10^4, a gap of missing values | `runtimes/bed` (bedGraph), by region with tabix | values layers (built); themes' log (built, EC-215) | 169 to 175, 201, 202, 215 |
| **Samples' variants as rows, long deletion solid**, and the **genotype matrix** | a VCF, 40 samples, 2,000 sites, one 600 bp deletion over 12 SNPs, a haploid call, a no-call, phased and unphased calls | `runtimes/vcf` with tabix, BGZF | GF VCF (genotype classes), RR `matrix` shape and layer; `samples` narrows | 220, 222 |
| **Methylation matrix by region**, and a **gene-by-sample heatmap** | 20 samples by 250 CpG, `value` with NaN missing, wide and long forms; the same layer over 150 genes | wide: the `seq` format's matrix; long: records, pivoted by RR's records transform at open | RR `matrix` (`value` float32, NaN missing) | 220, 221 |
| **Transcripts, three levels** | 150 genes on toy A, 1 to 4 transcripts each, 3 to 12 exons, both strands; GFF3 and GTF of one model | `runtimes/gff` (`features: transcripts`), whole and by region | GF `assembleTranscripts`, span rows, `line: 'gene'` | 227, 229 |
| **Open on a gene by name** | the GFF3's `Name` and `ID` attributes | `findNames` over the GFF3 (or GTF) read whole | GF name index; GA's locus fallback | 228 |
| **Genome by region, names aliased** | letters of both toys from an in-memory adapter on the seed; lengths as `chrom-sizes`; one small FASTA (toy A's `chrM`, 16,569, and `HLA-A*01:01:01:01`, 3,503) for the built FASTA reader; an alias table (`1`, `chr1`, an `NC_`-style accession; `MT` against `chrM` refused) | the adapter (`withAdapter`), `runtimes/bed` (`chrom-sizes`), the built `fasta` | GA names rule, `match: chr`, `(block :name …)` quoted, `seq.resolve.assembly.lengths` | 186 to 189, 200 |
| **A liftover** (cut last) | toy A's map to a revised toy A′ (an insertion, a deletion, a duplicated run), runs in the `seq` format, and positions with truths computed in code | the `seq` runtime (`map.source`, read by either side) | LX liftover (a bridge over the map, A's genes drawn on A′); GA projection by region | 211, 212 |

**Deferred to the close-out slice** (the decisions of 2026-11-21): the **Manhattan family** (points at a height, a
reference line, whole genome by chromosome, LD colouring, a coarse bin keeping its one significant SNP: EC-214, 216,
217, 219), **many signals from one file** (`split` on values, EC-214), and **transcripts coloured by joined results**
(the keyed join, EC-224; RR's `rows` shape is its data). Their generator steps go with them.

##### The generator: outline

`packages/seq/examples/src/genomes/` (modules under 800 lines, §module-size), one `generateGenomes(options)` returning
bytes and truths, one record script case `genomes` (`scripts/record.ts`, run by hand) writing the fixtures and their
`sources.json` rows with licence `LicenseRef-seqx-generated`, and the real slices as their own rows. Randomness only
from core's `PCG` and `hash3` (as `p53/chromosome-data.ts`), a stream per concern seeded by `(seed, concern, index)`, so
adding a concern never moves another's bytes; no clock, no `Math.random`. Output byte-identical on every run, checked by
a test comparing checksums to the committed fixtures (as the other datasets' fixture tests).

Steps (each a module, each exporting its truth for the tests):
1. **Assemblies**: toy A: `chr1` 300,000, `chr2` 200,000, `chr3` 120,000, `chrM` 16,569 (cyclic), `HLA-A*01:01:01:01`
   3,503, a 20,000 gap of N inside chr1; toy B: the same names unprefixed (`1`, `2`, `3`, `MT`) so the alias rule is
   exercised, an 80,000 inversion on chr1, a 60,000 translocation from chr3 to chr2; toy A′ for the liftover. Letters
   by region from `hash3(seed, sequence, tile)` (GC share by region, a few repeats; aligned runs share a base sequence
   with 1 % substitutions and small indels), served by an in-memory adapter; lengths written as `chrom-sizes`.
2. **Alignment truth**: runs (a, b, strand, gaps), written as sided pairs in the `seq` format (rows sharing `row` the
   blocks of one alignment), and the dot plot's grid; the liftover's runs and truth positions.
3. **Genes**: 150 on A, 3 levels (gene, transcript, exon), both strands, names (`TOY0001`...); GFF3 and GTF forms of one
   model, bgzipped and tabix-indexed beside the plain files.
4. **Signal**: a bedGraph with a TPM-like column, smooth waves plus noise, a gap of missing values; bgzipped and indexed.
5. **Variants**: 40 samples by 2,000 sites, the deletion, the haploid, the no-call; bgzipped and tabix-indexed.
6. **Methylation**: 20 samples by 250 CpG, 8 % missing, wide (`seq` format) and long (records).
7. **Links** (the ring): 12,000 BEDPE from the typed GRCh38 lengths, generated at load (about 700 kB as a file: over
   the dataset's share), a 1,000-link BEDPE committed for the reader's tests.
8. **Hi-C and dot plot grids**: the links and the alignment binned into two levels in the `seq` format.

Sizes (committed, estimated; the generator's test asserts each under bound): small FASTA 20 kB; `chrom-sizes` and
aliases 1 kB; alignment and liftover in the `seq` format 25 kB; grids 60 kB; GFF3 plus GTF 90 kB; bedGraph (bgzipped,
tabix) 120 kB; VCF bgzipped with `.tbi` 60 kB; methylation 70 kB; BEDPE of 1,000 links 60 kB; cytobands 4 kB (real).
Total about 0.5 MB, under the 2 MB dataset bound. The real slices are recorded as pages in the formats' own `_test`
fixtures, each its own row; their sizes asked with the fetch.

Heavy performance fixture, never stored: a whole human genome by region (the plan's M5, performance.md): GRCh38's 24
lengths summing to 3.1 Gb as an adapter in the style of `p53/chromosome-data.ts` generalised from one chromosome:
letters a 64 KiB tile at a time from `hash3(seed, chrom, tile)`, about 20,000 genes from the same stream, a GC signal
with exact bin means and a 12,000-link set by region. Tests per chromosome (heavy, `itWithin`), a whole-genome fit under
`maxAssemblyBlocks` (GA), and a quick representative at 3 Mb.

##### Real slices to fetch (each asked before the fetch, Gate 0, 7; none made)

Around TP53 (GRCh38 chr17:7,668,402-7,687,550), recorded as pages (GF). URLs, sizes and licences are checked at the ask.

| # | Format | URL (to confirm at the ask) | Expected download and kept | Licence to check |
|---|---|---|---|---|
| 1 | cytobands (BED) | `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBandIdeo.txt.gz` | a few kB, whole | UCSC conditions of use for database tables; NCBI's public-domain ideogram if refused |
| 2 | VCF sites | `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2024/clinvar_20240107.vcf.gz` and `.tbi` | ranged reads of a file of about 170 MB; kept: the `.tbi` whole (about 520 kB: the per-view byte bound tested at the fit from it), the header, TP53's chunks | NCBI, public domain |
| 7 | GFF3 | Ensembl GRCh38 `Homo_sapiens.GRCh38.<release>.chr.gff3.gz` | about 50 MB whole; kept TP53's lines, `##sequence-region` | Ensembl |
| 8 | GTF | GENCODE `gencode.v<N>.annotation.gtf.gz` | about 50 MB; kept TP53's lines | GENCODE (open) |

Four fetches (the numbers are the first list's): the cytobands serve the showcase (RN), the other three the built
runtimes' conformance suites (GF). Cut: 3 (the 1000 Genomes panel), 4 (bigBed), 5 (2bit), 6 (bgzipped FASTA), 9
(chain), 10 (PAF): their runtimes are designed only, or the generated data covers them.

##### Edge cases this dataset takes

| Edge cases | Where |
|---|---|
| EC-108, 169, 171, 206 to 212, 215, 218, 220 to 223, 225, 227 to 229, 235 | a view above each, named by a test of that view |
| EC-186 to 189, 200 | the names view (aliases, the `chr` rule, the HLA contig, a sequence without a length) |
| EC-214, 216, 217, 219, 224 | deferred to the close-out slice with their views |
| EC-213 | cut with MAF |
| EC-226 | themes' Ranges (deferred there); no view here |

##### Configuration

| Item | Whose | Default | Why | Alternatives |
|---|---|---|---|---|
| `genomes.seed` | option of the generator | one fixed value | byte-identical fixtures | another |
| `genomes.toy.lengths` | option of the generator | 300,000, 200,000, 120,000 | fixtures under 1 MB | larger (heavy) |
| `genomes.links` | option of the generator | 12,000 | RN's frame-time claim | 100,000 (heavy) |
| `genomes.samples`, `.sites` | option of the generator | 40, 2,000 | a screen of rows; a tile of sites | 2,504 samples (heavy) |
| `genomes.whole.length`, `.genes` | option of the adapter | 3.1 Gb in 24 chromosomes, 20,000 | performance.md | none |

##### Decided by the author

2026-11-21, as recommended (the author: "generally reasonable and we can fold in the reviews"; reviewed before promotion (the fifth slice's GD review), the whole-designs review's twelfth recommendation).

1. **Toy B unprefixed, toy A prefixed**, so the alias rule is exercised: **decided yes**.
2. **The 12,000 links generated at load**, a 1,000-link BEDPE committed for the reader: **decided yes**.
3. **The toys' letters from an in-memory adapter**, one small FASTA committed for the built reader and the HLA name:
   **decided yes**.

### What is not here

- **The Genomes dataset's deferred views** (the author, 2026-11-21): the Manhattan family (EC-214, 216, 217, 219),
  many signals split from one file (EC-214) and transcripts coloured by joined results (EC-224) go to version 1's
  close-out slice with their generator steps; MAF is cut (EC-213); the alignment views read the `seq` format, not PAF or
  chain files ([genome-formats](genome-formats.md), Appendix: designed only).
- **Loading by URI in the Storybook** (public files read live): settled with the source protocol
  ([sources.md](sources.md)); every view here reads committed fixtures only.

## The interaction matrix's cells

The cells the datasets' views exercise, with the view that shows each (each view's test drives these).
Cells owned by another design say so; this design adds no gesture of its own.

| | select | extend | add / remove | focus | back | zoom in / out | pan | hide | edit |
|---|---|---|---|---|---|---|---|---|---|
| **an element** | a residue, a feature, a pair (arc), a contact, a ligand, a dot plot's cell (two positions): click; `selection.set` (every dataset) | Shift+click | Cmd/Ctrl+click | double-click; on a cyclic system by the shortest arc (Plasmid, EC-12) | double-click again | — | — | no: hide its track, or filter its layer's data | no: data is not edited in the viewer |
| **a range** | drag on data; across the origin (Plasmid); a rectangle on a dot plot selects two ranges (tRNA, Complex; [grids-and-rows.md](grids-and-rows.md)) | Shift+drag | Cmd/Ctrl+drag, Alt+drag | Enter | Alt+Left | — | — | — | — |
| **the selection** | — | Shift+arrows (Globins 10,000) | Cmd/Ctrl+Shift+arrows | Enter | Alt+Left | — | — | — | a bookmark from it (Antibody CDRs, Complex interface, tRNA stem): the selection's bookmark command, its name and key planned (decided) |
| **a leading section** | — | — | — | — | — | the dot plot zooms both axes ([grids-and-rows.md](grids-and-rows.md)) | the ring by its window ([layouts.md](layouts.md)) | no | the linear view's start on a cyclic system (Plasmid; [layouts.md](layouts.md)) |
| **a follower section** | as leading | as leading | as leading | as leading | as leading | moves its leader (the linear detail moves the ring's window) | moves its leader | as leading | as leading |
| **the view** | — | — | — | 0 fits | Alt+Left / Alt+Right | — | — | — | undo and redo of every editing view |
| **a track or layer** | a row of an alignment by its header (Globins; [grids-and-rows.md](grids-and-rows.md)) | Shift on headers | Cmd/Ctrl on headers | a clade's summary focuses its rows | Alt+Left | — | — | `track.hide`; a clade folded by its tree node (session state); a layer's items hidden by a filter in `layer.data` (tRNA pairs by kind, AlphaFill by score, Plasmid by enzyme, p53 variants by source): the query sublanguage's filters, never a theme | residue scheme and palettes ([themes.md](themes.md)), pins and `clades.keep` ([grids-and-rows.md](grids-and-rows.md)), `copies` merged or split (Globins, Complex, p53) |
| **a bookmark or link** | — | — | `refs.add`, the selection's bookmark command; `refs.remove` (planned) | runs its action (the heme pocket, a ligand's contacts: highlighted) | Alt+Left | — | — | — | its label (Antibody CDRs; planned) |
| **a navigator** | — | — | — | a click centers the view (Globins' vertical overview; [navigators.md](navigators.md)) | Alt+Left | the window's handles | drag the window | the host hides the part | — |

## Configuration

Whose, by the rule the axis design states. None is a viewer setting or a spec param: each is an option of a call (a
builder's, a generator's or the record script's argument, documented in the module that owns it; this table is the one
list, `seqx-examples` having no `_spec`) or a fixed rule of the fixtures, checked by a test.

| Item | Whose | Default | Why | Alternatives |
|---|---|---|---|---|
| `fixtures.maxFileBytes`, `fixtures.maxDatasetBytes` | fixed (the fixtures test) | 1 MB, 2 MB | decided; the repository and the Storybook's assets | larger for a dataset the author names |
| `fixtures.licences` | fixed (the record script and the fixtures test) | `CC0-1.0`, `CC-BY-4.0`, `LicenseRef-NCBI-public-data`, `LicenseRef-seqx-generated`, `LicenseRef-anarci-numbering` | permissive, attributable (decided) | one added by the author's decision |
| `p53.variants.sources` | option of the record script | UniProt, ClinVar | licences known per source (decided) | gnomAD once its licence is checked |
| `chromosome.length`, `.genes`, `.seed` | option of the generator | 83,257,441 (GRCh38 chr17), about 1,200, one fixed seed | p53's chromosome at its real size, TP53 in place | a whole genome's 3.1 Gb by chromosome (heavy only) |
| `globins.seedRows` | option of the record script | at most 60 | the dataset table's "about 50"; a figure that fits a screen | the whole seed |
| `globins.scale.rows`, `.columns`, `.seed` | option of the generator | 10,000, 300, 42 | performance.md's scale; one fixed seed | 1,000 × 2,000, 10M lazy rows (heavy only) |
| `globins.reference` | option of the view's builder | HBA_HUMAN | 4HHB's α chain is its unit | HBB_HUMAN |
| `antibody.design` | option of the generator | the edit list in its module | indels and IMGT insertions certain | another list |
| `plasmid.rotation` | option of the generator | the origin at bla's middle | a CDS across the origin | any offset (the property test sweeps all) |
| `plasmid.enzymes` | option of the view's builder (written into the sites layer's `layer.data` filter) | the MCS's enzymes and unique cutters | what a cloning user looks at | any list |
| `plasmid.scale.length`, `.features` | option of the generator | 200 kb, 2,000 | a BAC-sized ring for the ring's frame target | 3M (the sequence scale) |
| `trna.covariation` | option of the generator | mutual information with a small-count correction | computed from committed data, simple | RNAalifold-style covariation |
| `trna.scale.size`, `complex.scale.copies` | option of the generator | 4,000; 200 | performance.md's scales | — |
| `complex.alphafill.count` | option of the record script | ATP and nine more, by score | a list a user reads; under the size bound | all transplants (over the bound) |
| `complex.pae.between` | option of the generator | generated, from 1FIN's interface (decided) | no permissive prediction of the complex | a ModelArchive prediction; one the author supplies |
| `contacts.cutoff` | option of the record script | 4 Å between heavy atoms | the common definition | 5 Å |
| the seeded generator | fixed: core's `PCG` | — | one generator, no private copy | — |

## Tensions

- **Published files against converted ones.** A Stockholm seed committed as published is what a reader without Seq\*
  takes; a structure converted to `seq` is not, and the conversion is code in a script tests never run. Recommendation:
  as decided (native where the core reads it, converted elsewhere, the download's checksum kept so the conversion can
  be repeated); the fourth slice's structure runtime may later read a trimmed mmCIF instead.
- **The record script needs structure readers.** Contacts, helices and observed residues come from mmCIF, which
  `@molstar/io`'s CIF parser reads; the examples package declares `@molstar/core` and `@molstar/seqx-core` only, and
  examples-and-stories.md says it depends on `seqx-core`. No `noRestrictedImports` rule in `biome.json` covers
  `packages/seq/examples` (only the repository-wide zod and graphics-backend rules apply), so the boundary is the
  package's declared dependencies and that design. Decided (the data's decisions, 5): `@molstar/io` as a devDependency for the record script.
- **Real data moves.** A Pfam, Rfam or AlphaFold DB release changes rows, scores or models; a recording again breaks
  anchors and records. Accepted: the version is in `sources.json`, the anchors say what changed, and a recording is
  a reviewed change like any other.
- **Synthetic and derived data shown as real.** The PAE between the chains, the antibody's design and pUC19c's features
  are generated or found in code; a user may read them as measured or curated. Their names, `sources.json` (derived)
  and the dataset's Storybook index say so; a standing tension while no permissive source for them exists.
- **Edge cases on variants, not entries.** EC-7, EC-11 and EC-38 are carried by generated variants, not by the
  published entries; a view that claims them runs over the variant. Accepted: a real entry carrying each was not found
  under a permissive licence.
- **Attribution travels with a CC BY file.** An exported spec carries data that needs its attribution; the language has
  no attribution field. Decided (the data's decisions, 4): two optional params of a source.

## Self-check

The claim: **this design covers what a user does with the datasets.**

1. **Assume it holds.** Every dataset above has its entries, licences and sizes; every view above
   names what a user does and an editing view; the cells it exercises are named; the large variants are generated;
   the fixtures are read under the policy the tests already use (the Node adapter confined to `fixtures/`, the
   Storybook's plugin serving the same folder, the render CLI with that folder as its root and no network, EC-126).
2. **What it rests on, negated.** It rests on *the data being ours*: fixed, small, committed, known to the tests. The
   negation: **the user brings the data** (their family, their plasmid, their complex), and the inverse of reading a
   fixture, **writing one out** (the exported spec file, the third way a view runs).
3. **The negated world.** A user's alignment has no reference row, no tree, 3,000 rows; a user's plasmid has its
   features across the origin already; a user's structure has no published pairs. A view written against the fixture's
   constants (row 0 is HBA, bla's offsets) fails on any other file. An exported spec carries CC BY data with no
   attribution. And the plural: views tested on one instance each never meet a second one of the same shape, and an
   edge case the one instance lacks is claimed but never met (the review found three: EC-7, EC-11, EC-38).
4. **Both held.**
   - Integrated: every view's builder takes its dataset as an argument, and its test runs it over a second instance of
     the same shape and over the instance carrying its edge case (Globins' generated alignment, the antibody design,
     Plasmid's rotated record with both strands across the origin, the Rfam seed's rows as second tRNAs, the generated
     grid and copies). Anchors state facts only in the anchor tests, never in a view. The record script's checksums
     make a re-recording a visible change.
   - Integrated, planned: the licence, attribution and derived parts of every fixture in `sources.json`, shown on each
     dataset's index; attribution in an exported spec, two optional params of a source (decided).
   - Integrated, from the same negation applied to the views' claims (an interaction claimed that no default gesture
     performs is a claim never met): each view's "what a user does" names a gesture of the matrix or a mechanism (a
     `pairs` layer of contacts, a ref the builder writes and the user runs, a block's `of`); the rest were removed.
   - Not integrated, standing: **a user's own file in the Storybook** (drop a file onto a view): needs readers the
     slice does not build (GenBank, mmCIF) and a host-side upload path; carried to the fourth slice with services.

## Consequences

- **`seqx-examples`**: per dataset a folder `src/<dataset>/` (data, generators, one module per view, `build(data?)`)
  and `fixtures/<dataset>/` (the files and `sources.json`); `scripts/record.ts`, the `sources.json` test and the
  seeded generator on core's `PCG` built with Globins, the other datasets adding to them; `views.ts` lists the new views. Tests
  (quick): each fixture under its bound and listed in `sources.json` with an allowed licence id and its checksum;
  generated-and-committed files equal their generators; the antibody design's numbering equals the committed ANARCI
  oracle; anchors per dataset; every view validates, resolves, lays out and renders with no unexpected issues, over each
  instance it names; drawing checked against the data in each view at the fit view and across a tile edge (every item
  in view has marks, every mark is hit); property tests (Plasmid's rotations over all offsets on the record and its
  reverse complement; the generated alignment's rows equal per index across calls; covariation symmetric). Quick
  representatives of the targets the views carry (the 10,000-row open at 1,000 rows; a ring frame on the scale ring at
  a tenth; a dot plot frame on a 400 grid), the full sizes under `itWithin` in the heavy lane. Medium: every view in
  both encodings (JSON, MessagePack).
- **`seqx-core`**: the Stockholm and Newick readers ([grids-and-rows.md](grids-and-rows.md), with their fuzz tests); nothing else from this design.
- **`seqx-viewer`**: nothing of its own; the views use what the other designs build.
- **Storybook**: `Seq* (next)/Examples/<Dataset>/<View>` for every view, canvas and SVG side by side, a headless twin
  for each story composing more than one consumer (the ring and its linear detail, arcs and the dot plot, the PAE and
  the tracks, hemoglobin and its family); each dataset's index lists its views, sources, licences, attributions and
  what is generated or derived; the Design page "Examples" gains the datasets.
- **`_spec`**: `seqx-core`'s `_spec` gains the readers' sections; the examples package's fixture rules (provenance,
  licences, bounds, record-by-hand, generated-and-compared, the instance per edge case) and its configuration stay in
  this design and the fixtures test, its options in the modules that own them (no `_spec` of its own: built so in the
  third slice).
