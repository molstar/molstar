## Seq* (next): formats by region

- BED, bedGraph, GFF3, GTF and VCF files with a tabix index read by region past the finest level: coarser levels from the index's bytes with exact counts, so a view's reads stay bounded and a whole chromosome shows a density before "zoom in".
- A GFF tile inside an intron shows its transcript's span row; a VCF's samples draw as a genotype matrix (no copy, some, every copy, partly called, uncalled), a long deletion one item across the variants it covers.
- A BEDPE file serves pairs by region, so its links and the same links held in the spec draw as one picture.
