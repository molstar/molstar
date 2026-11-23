## Seq* (next): version 1's close-out

- Seq\* viewer: a legend's eye hides or shows its layer, one undo step; a hidden layer's row stays in place, dimmed,
  with the eye that shows it again, also for a layer hidden from the start and a bridge's layer (`mountLegend` option
  `layers: false` drops the eye). The track header's menu and the eye name a layer the same way.
- Seq\* core: a bridge's ribbons have an `orientation` field (same or inverted; a data column of that name wins),
  keyed in a legend with the width's meaning; a continuous legend over an open domain says "in view".
- Seq\* genome formats: a transcript tagged `Ensembl_canonical` or `MANE_Select` in GFF3 or GTF is drawn first in its
  gene's stack with a bold label, and its card says so; a matrix's samples (a VCF's genotypes, a methylation matrix)
  are named beside its rows; a genotype card gives the call as the file wrote it (`0|1`, `1/2`, `1`) and what it is; a
  missing value has a colour of its own.
