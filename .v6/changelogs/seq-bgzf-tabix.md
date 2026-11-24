## IO: BGZF and tabix

- `@molstar/io/common/bgzf` and `common/tabix`: BGZF read whole or by virtual offset over core's range reads; TBI and
  CSI indexes, a region's lines as bytes with extents by htslib's rule; every bound refused before it allocates.
- Fix: core's `ungzip` reads every member of a multi-member stream (BGZF, concatenated `.gz`), skipping trailing zero
  padding; a member cut short throws.
- `@molstar/test-harness/genome/bgzf-writer`: a test-only BGZF and tabix writer.
