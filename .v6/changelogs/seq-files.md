## Seq* (next): files as sources

- A records runtime reads TSV, CSV, JSON and MessagePack, as an array of records or an object of columns alike (a typed
  column kept without a copy when it is already in the protocol's conventions); the `source` node describes the file's
  conventions (`fields`, `records`, `letters`, `system`, `base`, `ends`, `types`, `attribution`, `licence`); row-level
  issues (`seq.source.records.*`) reported with their counts and places; bounds checked before parsing.
- The `seq` format reads MessagePack too, its codec loaded only when a file is one.
- Records sources read GFF3 (`format: gff3`: nine columns, from 1, attributes as fields; a feature ending past a
  ring's length drawn as one item across the origin, clipped with a warning on a linear system) and a matrix
  (`shape: matrix`, an array of arrays checked against its axes' total).
- A layer's data written as `<annotation>/(item …)` (a lone item form, or `and` or `or` of items) is now the track's
  own rows, not a filter: its rejected rows are not counted in the "n filtered" chip, and the viewer's hide and legend
  actions keep it. Only a trailing `(not …)` is the user's filter.
- `Features/Files` stories: UniProt's P04637 entry read as published, the same table as TSV, CSV, JSON and MessagePack,
  a plasmid's GFF3.
