# Examples: few datasets, many views

## Intent

The author, 2026-11-17: "we could have few example datasets, but should have more views. i would like to keep the
richness, the editing experiences etc."

A handful of datasets, each rich enough to carry many views. A view is one spec over its dataset (or a snapshot of
one), showing one way of looking or one editing experience, and exercising named edge cases
([edge-cases.md](edge-cases.md)). New views cost a spec, not new data.

## Decision

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
  | **Genomes** (a later milestone, [version-1.md](version-1.md)) | two small toy assemblies with synteny (an inversion, a translocation, a gap), name aliases and an HLA-style contig name; genes with transcripts (three levels), a name index; multi-sample signal in one file, GWAS points, a per-pair LD table, a variant matrix by sample with a long deletion, methylation with missing values, transcript-level results keyed by id; FASTA with an index and values read by region |

- **Views per dataset** (the first list; each names the edge cases it carries):

  | Dataset | Views |
  |---|---|
  | p53 | text figure with feature underlines and variant marks (EC-94 to 99); tracks with features, values and secondary structure; the protein on its gene through codons, introns collapsed, minus strand (EC-15, 18, 29, 31, 117, 144); structure chains linked, unobserved gray, copies merged and split (EC-65, 118); the gene's region by region with its conservation signal (EC-169 to 175, EC-201); the minus-strand gene's features in both directions; orthologs aligned, with hover across rows; queries and search with results stepping and bookmarks (EC-161, 162); a story in snapshots (EC-84, 85) |
  | Globins | family alignment with conservation logo and rows by tree (EC-38, 43, 106); collapsed clades with a summary (EC-100, 103); ESPript-style text figure (EC-97); hemoglobin complex with heme chips and copy groups (EC-64, 65); the 10,000-row alignment, lazy rows and keyboard cursor (EC-41, 44, 154) |
  | Antibody | numbering schemes side by side, Kabat vs IMGT labels (EC-2 to 8); design vs wild type with indel carets (EC-29); CDR selection and editing as bookmarks with undo (EC-162, 164) |
  | Plasmid | ring with a linear detail following its window (EC-108, 114); a feature across the origin selected and focused by shortest arc (EC-9 to 14); the rotated linear view |
  | tRNA | arcs nested and crossing (EC-129); base-pair matrix as a dot plot with crosshair hover and rectangle selection (EC-109, 145); covariance alignment with pairs |
  | Complex | the chains as blocks with the interface residues and the contacts as arcs across the chains, folded by copy group (EC-35, 62, 65); the PAE as a dot plot over the same window as the tracks, crosshair hover lighting both residues, a dragged rectangle selecting two ranges (EC-109, 145); confidence and secondary structure tracks, unobserved residues gray (EC-118); the ligands as blocks after the polymers, ATP as one position and per heavy atom, AlphaFill's transplants listed with their scores, and the residues contacting a chosen ligand highlighted (EC-64, 101); a text figure with ligand chips flowing like words (EC-99); copies merged into one block labeled ×2 and split again (EC-65); the interface bookmarked as a ref, with undo (EC-162) |
  | Genomes (later) | synteny ribbons with an inversion (direction on segments, EC-206); a dot plot of one genome against the other; both genomes on one ring (EC-210); signal by region with an open domain and with a one-sided domain (EC-169 to 175, EC-201, EC-202); many signals split from one file, log scale (EC-214, EC-215); a Manhattan plot of points with a threshold line, colored by LD to a chosen lead (EC-216, EC-217, EC-219); samples' variants as rows with the deletion solid (EC-222); methylation matrix by region (EC-220); transcripts colored by joined results gated on padj (EC-224); open on a gene by name (EC-228); a liftover by chain file (EC-212); the genome by region from FASTA with its index, names aliased (EC-186 to EC-189) |

- **Showcases at the application layer**, in the storybook package (they compose what an application would, and are
  not core features):
  - **an application composing services**: UniProt, PDBe (SIFTS) and AlphaFold DB as runtimes, from recorded
    fixtures, a switch to go live;
  - **predictions onto a structure**: a lab's per-position predictions on a sequence indexed from 1, with no numbering
    of its own, mapped onto an mmCIF structure through `label` and `auth` numbering (an offset, insertion codes,
    unobserved loops, an expression tag), colored on the sequence and in Mol\*, hover and selection both ways
    (EC-3 to EC-5, EC-25, EC-118);
  - **addressing**: a text with links that focus, highlight and select in the view and in the structure; "copy link"
    from the view back to text; the same location in several numberings; rects from a dot plot and an alignment; refs
    as bookmarks ([addressing.md](addressing.md)).
- **Editing experiences are views too**: each dataset has at least one view that is edited (a bookmark added, a track
  hidden or reordered, a theme changed, a snapshot recorded), with undo and redo, a preview that ends in one patch or a
  cancel, and the printed spec updating beside it.
- **The previous implementation's cases** each have a home (its 22 cases, [sources/seqstar](../sources/seqstar/README.md)):

  | Its case | Here |
  |---|---|
  | p53 text and tracks; UniProt features; queries and bookmarks; protein on its genome through codons; reverse-strand gene; streamed chromosome | p53 |
  | protein family; PF00042 alignment; family triage | Globins |
  | antibody numbering; indels in wild-type numbering | Antibody |
  | DNA construct; plasmid ring | Plasmid |
  | RNA base pairs | tRNA |
  | imatinib in ABL; barnase–barstar; AFDB complex; designed dimer | Complex |
  | protein over coding DNA | p53 (the protein on its gene) |
  | genotype matrix | Genomes (later) |
  | story drafts | p53 (a story in snapshots) |

- **Every view runs three ways**: in Storybook (canvas and SVG side by side), in a quick-lane test (validates,
  resolves, lays out with no unexpected issues, renders to SVG; its printed spec and draw-list record kept), and as a
  spec file a reader without Seq\* code could take (the JSON and its data files).

## Provisional

The exact data (which structure for p53, which globins, the toy genomes' sizes; for the Complex, the crystal entry, the
prediction and AlphaFill's entry, checked against their sources the way the TP53 exons were) is chosen when each dataset
is first used, small enough to commit; public files the Storybook may also load by URI are settled with the source protocol
([sources.md](sources.md)).
