## Seq*: GFF3 and GTF

- `runtimes/gff` reads GFF3 and GTF whole as transcripts, genes or rows (`features`), with a name index
  (`findNames`: GFF3's `Name`, `ID`, `Alias`; GTF's gene and transcript names and ids); the ranges layer draws gene
  models (`line: 'gene'`) with their names beside them.
