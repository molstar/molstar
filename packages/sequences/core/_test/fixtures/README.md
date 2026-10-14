# The shared Seq* test fixtures

Read by the core's, the viewer's and the Storybook's tests and stories. Tests never fetch: these are recorded.

- `1brs.json`, `1iep.json`, `p53.json`, `pf00042.json`: the Seq* test fixtures, made (by the viewer spike's
  `scripts/prepare-fixtures.ts`, removed with the spike) from the Seq* prototype's frozen inputs (`seqstar-dev`, branch
  `prototype-260819`, `prototype/fixtures`; that pack's `metadata.json` has URLs and checksums). `1brs.json`: barnase
  (P00648) and barstar (P11540), sequences, UniProt features, residues observed in 1BRS chains A and D, interface
  contacts (≤ 4 Å heavy atoms), a synthetic confidence per residue (UniProtKB CC BY 4.0, wwPDB 1BRS CC0 1.0).
  `1iep.json`: ABL kinase chain A (225–498) and imatinib (STI, 37 heavy atoms), contacts per residue for the crystal
  pose and AutoDock Vina's 4 docked poses, classified approximately (from `autodock-vina-1iep`, Apache-2.0). `p53.json`:
  P04637 (UniProtKB CC BY 4.0). `pf00042.json`: Pfam PF00042.29 (globin), 32 rows, HBA_HUMAN first (Pfam/InterPro CC0
  1.0); case 22's family and the keyed-sets story's.

- `pf00042-seed-8.aln`: `pf00042.json`'s first eight rows (Pfam PF00042.29, CC0), the rows named `ID/start-end` as the
  File formats story names them, written 2026-10-14 in CLUSTAL W 2.1's layout (clustal_num: blocks of 60 columns, each
  line's cumulative residue count, a conservation line under each block from Clustal's strong and weak groups, `.` gaps
  as `-`) by a script that is not kept (no Clustal program run). The Clustal reader's tests (`model/formats.test.ts`,
  `format/load.test.ts`) and the File formats "Alignments and trees" story read it.

- `af-barnase-barstar-pae.json`: a synthetic predicted aligned error for the fixture model AF-BARNASE-BARSTAR
  (full-length barnase then barstar, 247 residues), in AlphaFold DB's PAE JSON shape
  (`[{"predicted_aligned_error": [[…]], "max_predicted_aligned_error": 31.75}]`, Å to one decimal). Not a prediction:
  made from `1brs.json` by a deterministic formula, low within each folded chain, high to and from barnase's first 47
  residues (its signal peptide and propeptide) and across the chains, but for a lower patch between the residues near
  1BRS's contacts. Case 8 and `@molstar/seq-molstar-poc`'s PAE test read it.

- `variation-P04637-ends.json`: eight of P04637's Proteins API variants
  (https://www.ebi.ac.uk/proteins/api/variation/P04637, fetched 2026-10-14; CC BY 4.0) with the response's `accession`
  and `sequence` (393 residues), each feature trimmed to type, begin, end, wild type, alternative sequence and
  consequence: three inside the sequence, two running past its end (377–400, 385–400) and three starting past it (394, a
  stop lost; 394–398 and 399, frameshifts). The `uniprot` service's test of variants past the sequence's end reads it
  (`format/load.test.ts`).

- `gc.prt`: NCBI's genetic code table, Version 4.6 (27 tables), from `https://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt`
  (fetched 2026-10-09; NCBI data, public domain). The model's `GC_PRT` (`src/model/genetic-code.ts`) copies its
  `ncbieaa` and `sncbieaa` strings; the model's codon test checks every table's 64 codons against this file.
