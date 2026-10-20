## Performance: large unshared structures

Measured on 1tup copied into one model (`tests/src/pipeline`, PR #1997): 1.1M atoms and 2,496 units (`flat-192`),
5M atoms and 11,232 units (`flat-864`).

- Intra-unit bonds are about 1.4× faster: distance bonds query 3.1 Å (above every bond threshold) instead of 4 Å,
  and dictionary bonds come from the residue's atoms. 5M atoms 5.2 → 3.7 s, 1M atoms 1.1 → 0.76 s; the bonds are the
  same, their edge order may change (#2007)
- Models derived from `entity_poly_seq` without `pdbx_unobs_or_zero_occ_residues` key missing residues per residue,
  not per atom: model creation at 5M atoms allocates 2.7 → 2.1 GB, 2.7 → 2.4 s (#2007)
- Geometry builders trim buffers their size estimate overshot (`ChunkedArray.compactTrimmed`), and meshes reserve 2
  triangles per vertex instead of 3. At 1M atoms the gaussian surface retains 325 → 203 MB, cartoon 381 → 351 MB (#2007)
- `CylindersBuilder.create` and `LinesBuilder.create` read `initialCount` and `chunkSize` as cylinders and lines (they
  counted vertices, 6 and 4 per item, while callers passed items): ball-and-stick bonds allocate their buffers once
  instead of twice (#2007)
