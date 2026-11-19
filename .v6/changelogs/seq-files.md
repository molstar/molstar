## Seq* (next): files as sources

- A records runtime reads TSV, CSV, JSON and MessagePack, as an array of records or an object of columns alike (a typed
  column kept without a copy when it is already in the protocol's conventions); the `source` node describes the file's
  conventions (`fields`, `records`, `letters`, `system`, `base`, `ends`, `types`, `attribution`, `licence`); row-level
  issues (`seq.source.records.*`) reported with their counts and places; bounds checked before parsing.
- The `seq` format reads MessagePack too, its codec loaded only when a file is one.
