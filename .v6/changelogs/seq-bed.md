## Seq*: the BED family

- `runtimes/bed`, an opt-in runtime for BED (BED3 to BED12, BED12's blocks as pieces of one item, columns named by
  number through `fields`), bedGraph, BEDPE (indexed by both ends, loaded when used) and chrom.sizes, read whole or by
  region through tabix; BED's name as the label and its strand drawn.
- The ranges layer keeps a list of at most a row per CSS pixel as items at the fit (labels, hits and selection).
- Tests read real files recorded as pages (`seqx-examples`' recording script and `withRecordedPages`).
