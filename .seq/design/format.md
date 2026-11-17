# The format: one tree, JSON and MessagePack

## Intent

A spec is a tree people, agents and other programs read and write. It needs encodings that are lossless and
interchangeable, one of them compact enough for real data (typed columns of numbers), a way to inspect it compactly, a typed builder most specs are written with, and a description of the language that everything else derives from.

## Precedents

- The previous SeqViewSpec: JSON, MolViewSpec-like nodes (`kind`, `ref`, `params`, `children`), parameter tables in
  TypeScript, a generated reference, no published schema ([sources/jbrowse-reply](../sources/jbrowse-reply/README.md)
  point 5).
- meso\*: S-expressions canonical (`.sexpr` files), JSON a lossless second encoding, languages described as
  S-expression data with a TypeScript twin ([sources/meso](../sources/meso/README.md)).

## Decision

- **The tree is canonical; JSON and MessagePack are two lossless encodings of it** (the author, 2026-11-17: "yes, let's
  keep using JSON with s expression printing for easier inspection"; then "maybe msgpack should be the default? it's
  interchangeable with JSON. and we already spent a lot of time getting the builder … understand typedarrays etc?").
  - **MessagePack is the default** for files and transfer: compact, and a typed column of numbers is written as an
    extension type (`0x11`–`0x18`: int8, uint8, int16, uint16, int32, uint32, float32, float64, little-endian), read back
    as that typed array without a copy where alignment allows.
  - **JSON is the readable encoding**: for docs, diffs, hand edits and tools that only speak JSON. A typed column is a
    plain array; a float32 column is written as its shortest decimals; NaN (missing) as `null`.
  - **One data model**: a number means its value (a float32 its shortest decimal); a typed column means the array of
    its numbers, NaN as missing. Encoding a tree in either and decoding it gives an equal tree; a typed array is a
    storage choice, never a different meaning. A reader tells the two apart by the first byte.
  - **The builder keeps typed arrays as given** (a column passed as a `Float32Array` stays one in the tree), and offers
    `float32` and other column kinds for numbers given as lists.
  - There is no `.sexpr` file type; the S-expression spelling is never stored.
- **S-expressions are a printing**: `print(spec)` gives a compact, diffable view for people, docs, tests and agents;
  a reader for the same spelling exists for tests and agent input, and lowers to the same tree. It is never stored.
- **One node shape**: `{ kind, id?, params?, custom?, children? }`. `id` is local (unique among siblings); an address is
  the path of ids from the root with query forms as segments (`/data/p53/(range :from 120 :to 130)`), meso\*'s grammar;
  refs are addresses, and authors name their own in a `refs` container ([addressing.md](addressing.md)).
- **The language is data**: one TypeScript description (`defineLanguage`, typed `as const`) declares kinds, where each
  may sit, the id policy, params with types, defaults, bounds and docs, and refs. From it are derived: validation
  (every issue, stable codes), defaults, the builder's types, the docs pages
  and the agent skill. Nothing is written twice.
- **Versions**: the language has one semver. A spec says its version; a reader migrates older specs with explicit
  tree-to-tree functions, and refuses newer major versions with an issue.
- **Unknown kinds** in a newer minor version, or in an extension namespace the reader lacks, are kept, reported and
  skipped (with the kind's declared fallback when it has one), never a failure of the whole spec.
- **Patches** (`set`, `unset`, `insert`, `remove`, `move`, each by address) are the only way a spec changes; a batch
  is all or nothing and returns its inverse; `diff(a, b)` returns the patch from `a` to `b`, tested as
  `apply(a, diff(a, b)) = b`.

- **The builder is the main way to write a spec** (the author, 2026-11-17: "the main use case will will using the
  "builder" to create the views"): typed from the description, it makes a misplaced kind or a wrong param a type error
  where it is written, and produces the same tree a reader would.
- **No published JSON Schema, initially** (the author, 2026-11-17: "i am not sure we need to publish the JSON schema, at
  least initially. it could end up being bloated"). The description could derive one later; that waits for a consumer
  that needs it (an outside reader or an agent tool that validates without the core).

## Rejected alternatives

- S-expressions canonical (meso\*): rejected by the author for this language; the spelling most readers and tools
  share is JSON.
- JSON Schema or a validation library as the description: one issue at a time, no slots, no builder types; the
  description is the source, and anything else is derived from it.

## Consequences

- `decode(encode(x)) = x` for both encodings, `json(msgpack(x)) = json(x)`, `parse(print(x)) = x` and
  `apply(a, diff(a, b)) = b` are tested properties, typed columns of every kind and NaN included (§round-trip-equality is meso\*'s; here
  they are conformance tests of the core).
- The language reference (generated docs pages) is what documents the JSON for readers without the core.
