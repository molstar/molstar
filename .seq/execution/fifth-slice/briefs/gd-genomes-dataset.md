# Lane GD: the Genomes dataset and its views (stage 0, read-only, fast)

**Task**: write `drafts/genomes-dataset.md`: from examples.md (Genomes, 112, 124, 151, 245, 413-414, 446) and the five
drafts beside you (GA, RR, GF, LX, RN; read them as they land), a table of each view (the ring, synteny ribbons, a dot
plot, Manhattan, the liftover, the genotype matrix, MAF if kept: Gate 0, 1 cuts it first) and what it needs: the
generated data (two toy assemblies with an inversion, a translocation and a gap, aliases, an HLA contig name,
three-level transcripts, a name index, multi-sample signal, GWAS points and LD, a variant matrix, methylation, keyed
results, FASTA by region), the real slices to record (one per format, Gate 0, 7: each asked before its fetch, with its
licence), the formats, the designs' pieces, the edge cases (EC-206 to 224, 227, 228), the performance fixture (a whole
human genome by region, heavy). Then the generator's outline (seeded, deterministic, sizes).

**Model**: fast (sonnet), read only. **Where**: the main checkout, read only except your one file; no code, no commits.

**Report**: under 300 words: the table's gaps (a view no design serves yet), the fetches to ask for, the generator's sizes.
