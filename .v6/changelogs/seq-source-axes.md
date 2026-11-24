## Seq* (next): source-backed axes

- A sequence whose letters come from a source takes its length from the source's coverage once open: `openModel(spec,
  store, options)` opens the sources, resolves again with their lengths (`ResolveOptions.lengths`, keyed by system and
  source, so an edited source never keeps the old length) and opens the rest; `renderSvg` uses it. A cyclic `seq`
  file, rows across its origin, rings and details over it, GFF3 rows past its length, grids and blocks work as with
  inline letters.
- A focus on a ring takes the shortest arc (EC-12); arrowheads are decided per item's width; a ring draws a layer's
  density as a band; `(row :key …)` works on an alignment read from a source, waiting (`seq.query.row.pending`) until
  its rows are held; legends include pairs layers; range labels fit narrow marks.
- The example views read their data from the committed files. `@molstar/core/util/map` gains `mapsEqual`.
