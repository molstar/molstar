## Seq* (next): one annotation file for several sequences

- A Features › Files story, "Two sequences from one table": a TSV whose `target` column names each row's sequence
  (GFF3's `seqid` works the same), read by one source into a section per sequence, each showing only its own rows; a
  row with an empty target is on the source's own system.
