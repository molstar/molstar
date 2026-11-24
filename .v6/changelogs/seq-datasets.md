## Seq* (next): the datasets

- Plasmid: pUC19c (INSDC L09137.2) recorded; its features (bla, lacZα, the cloning site, the pMB1 origin) found in the
  sequence and its restriction sites computed in code; a ring with its linear detail and features as arrows by
  strand, the rotated variant with features across the origin, the rotated linear view, a 200 kb scale ring, and
  editing views.
- tRNA: PDB 1EHZ's chain A (its 14 modified residues, its 27 base pairs from the entry's own table) and Rfam
  RF00005's seed recorded; arcs nested and crossing, chords on a ring, the base-pair dot plot with the family's
  covariation computed on 1EHZ's positions, a seed row with its `SS_cons` pairs, the seed under its logo, a 4,000 ×
  4,000 dot plot, and editing views.
- Pairs layers place their arcs over a sequence whose letters come from a source; the rows layer's letters legend
  follows the alphabet (RNA lists U); a pairs table's extra category columns can be filtered on
  (`(item :tertiary yes)`).
- Globins: the Pfam PF00042 seed (CC0, cut to the 60 rows nearest 4HHB's chains), UniProt's organisms of its rows (CC
  BY 4.0, the class derived from the lineage) and PDB 4HHB (CC0) recorded; the family by its tree (neighbor joining on
  identity) with the logo, species and identity columns, folded clades, a text figure, 10,000 generated rows,
  hemoglobin as blocks with heme chips, helices and contacts, editing views and the dataset's index. A folded clade's
  summary row keeps its `+n` however long the first label.
- Complex: PDB 1FIN (CDK2–cyclin A–ATP, two copies), six CC0 PDB entries of CDK2 with ligands, UniProt P24941 and
  P20248 and AlphaFold DB's confidence (CC BY 4.0, changes stated) recorded; the chains as blocks (merged ×2, split, on
  a ring), contacts as arcs, ligands as chips with their contacts coloured by ligand, ATP per atom, confidence and
  secondary structure, text, 200 copies, the PAE between the chains and over the whole complex (synthetic), editing
  views and the dataset's index.
- A grid section that follows a section of blocks in a row: one strip a pair of chains, chain boundaries drawn, a
  dragged rect split at them, the crosshair edged in the hover's ink; `gridSummary(view, rect)` (mean, min, max,
  count, each cell once, a triangle layer's half). A track's "n filtered" chip goes under its label when the track is
  tall enough.
- Antibody: trastuzumab's Fab (PDB 1N8Z, CC0) with SAbDab's IMGT numbering (CC BY 4.0, its licence field recorded)
  and Kabat derived from it (a framework of another occupancy refused), read as declared numberings from files; a
  generated design with indels (IMGT's 111.1 and 112.1 before 112; its numbering this generator's rule, ANARCI's check
  later); Kabat and IMGT rulers side by side with hover in both, a ruler-scheme edit, the design against the wild type,
  CDR bookmarks with undo and their list. A headless export opens the source of a numbering a ruler shows.
