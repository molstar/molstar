## IO: faster CIF and BinaryCIF parsing

- Text CIF is tokenized over its UTF-8 bytes: `CIF.parseText` (and `parseCifText`) also takes a `Uint8Array`; a string
  is encoded first. 124 MB, 1.12M atoms: parse 650 → 302 ms (median of 5), retained 320 → 265 MB; model 424 → 425 ms
  with 245 → 174 MB allocated. Field values are unchanged; a UTF-8 byte order mark is skipped; input of 4 GiB or more is
  refused instead of wrapping around.
- `CIF.parse(Uint8Array)` reads text CIF too, told apart from BinaryCIF by the first byte (BinaryCIF was the only
  option; `CIF.parseBinary` stays strict); the plugin's Parse CIF accepts a Binary data cell holding text CIF.
  `cif2bcif` reads text without a string.
- BinaryCIF: a field decodes on its first access, sorted and windowed `atom_site` columns are made on first read, the
  common encoding chains decode in one pass, and coordinates decode straight to Float32. 1.12M atoms: model retained 342
  → 244 MB, allocated 377 → 318 MB (time 391 → 419 ms, within noise); all `atom_site` columns decode in 67.5 ms instead
  of 84.5.
- Behaviour: the coordinates (atomic and IHM coarse) of a BinaryCIF model with quantized coordinates (FixedPoint, the
  usual encoding) are Float32, as a text CIF model's (they were Float64); raw Float64 coordinates stay Float64.
- core: `parseIntBytes`, `parseFloatBytes` (number-parser), `utf8Encode`; `utf8Read` is faster for short ASCII.
- BinaryCIF string columns stay dictionaries (an index per row into their distinct strings; views and windows gather the
  indices): 1.12M atoms, model 419 → ~334 ms, retained 244 → 124 MB, allocated 318 → 124 MB.
- [Breaking] `Column.__array` and `CifField.__array` are the method `storage()`: `{ kind: 'array', array }` or `{ kind:
  'dictionary', indices, values }`. `Column.ofDictionary` is new.
- Fix: `toArray({ start })` of a `Column.window` or `Column.view` of a computed column read from the window's first row.

