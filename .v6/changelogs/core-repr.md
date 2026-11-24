## Core: `@molstar/core/repr`

- `@molstar/core/repr`: a codec-independent tree core for languages described as data. Nodes with `custom` data and
  namespaced kinds; typed columns (eight typed-array types and a dictionary-coded category) identified by element type
  and bytes, every NaN one NaN; the Merkle hash; addresses with query segments, canonical printing and a form-resolver
  table; coded issues with severity. `repr/schema` adds `column` and `table` params, validation reporting every issue,
  normalization on read, defaults and docs. `repr/builder` is a typed builder engine whose methods return addresses;
  `repr/codec` a codec interface with a conformance suite; `repr/codecs/json` and `repr/codecs/sexpr` the JSON and
  S-expression codecs (undeclared columns tagged, long columns elided in printing). Written from meso\*'s tree core.
