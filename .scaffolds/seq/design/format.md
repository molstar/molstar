# The format: one tree, one language, codecs at the edge

## Intent

A spec is a tree people, agents and other programs read and write. It needs one language with one version, a value
model that holds real data (typed columns of numbers) without making equality and hashing slow, encodings that are
lossless and interchangeable, a compact printing for inspection, a typed builder most specs are written with, and a
description of the language that everything else derives from.

## Precedents

- The previous SeqViewSpec: JSON, MolViewSpec-like nodes, parameter tables in TypeScript, a generated reference, no
  published schema ([sources/jbrowse-reply](../sources/jbrowse-reply/README.md) point 5); typed columns written to
  MessagePack as extension types `0x11`–`0x18`.
- meso\*: S-expressions canonical, JSON a lossless second encoding, languages described as data
  ([sources/meso](../sources/meso/README.md)).

## Decision

- **One language, one version** (the author, 2026-11-17: "one"). The spec's levels (`data`, `views`, `refs`,
  `snapshots`) are container kinds of that one language ([composition.md](composition.md)), not languages of their own.
  A spec says its version; a reader migrates older specs with explicit tree-to-tree functions and refuses a newer major
  version with an issue.
- **The tree and its value model are codec-independent** (the author: make the tree core "codec indepedent and about
  abstract patterns"). The tree core ([tree-core.md](tree-core.md)) defines nodes, values and their identity; codecs
  read and write them at the edge and live where their libraries do.
- **Typed columns are declared, not guessed**: a param that holds a column says its element type in the language
  description (`int32`, `float32`, `float64`, `uint8`, a dictionary-coded category). The builder and every codec give
  that param that typed array, so a column decoded from JSON is the same typed array the builder made. Why: equality,
  hashing and printing need no value-by-value comparison across array kinds.
  - **Identity by bytes**: a column hashes by its bytes and element type, cached on its node; two columns are equal when
    their types and bytes are.
  - **Missing is NaN** in float columns, written `null` in JSON; integer columns have no missing value unless their
    type says so.
  - **Undeclared columns are tagged**: a column outside a declared `column` param (in `custom`, an `any` param, under an
    unknown kind) is written with its type, `{"$column": "float32", "values": […]}` in JSON and
    `(:$column float32 :values […])` in the S-expression printing; declared columns are plain arrays when the codec is
    given the language (lane A's review; to confirm with the author).
  - **Printing elides columns** past a short head (`{float32 × 3,000,000: 0.12 0.15 …}`, braces reserved); the S-expression reader
    refuses an elided column.
- **Encodings**, both lossless for the value model:
  - **JSON is the default interchange** (as recommended by the review, accepted by the author: "rest as recommended"):
    readable, diffable, readable by any tool and any other viewer without Seq\* code; columns as plain arrays.
  - **MessagePack is the option for large inline data**: typed columns as extension types `0x11`–`0x18`
    (little-endian), decoded without a copy where aligned. A reader tells the two apart by the first byte. Built (lane
    L, then moved down: the encoder and decoder into `@molstar/core/data/msgpack`, the tree's codec into `@molstar/core/repr/codecs/msgpack`, io re-exporting both) with three more extensions, confirmed by the author (2026-11-18): `0x10` the numbers
    MessagePack drops (NaN, ±Infinity, −0), `0x19` a category column (its dictionary, then int32 codes), `0x1a` a
    bigint.
  - **Large data belongs in sources**, not inline (§bulk-data-not-nodes): a spec carries small inline values; a
    3-million-letter sequence or a dense signal is a file read through the source protocol ([sources.md](sources.md)).
- **S-expressions are a printing**: `print(spec)` for people, docs, tests and agents; a reader lowers the same spelling
  to the same tree (columns aside). Never a file type, never stored.
- **One node shape**: `{ kind, id?, params?, custom?, children? }`. `id` is local; an address is the path of ids with
  query forms as segments ([addressing.md](addressing.md)).
- **The language is data**: one description declares kinds, where each may sit, the id policy, params with types
  (column types included), defaults, bounds and docs, and refs. Validation (every issue, stable codes), defaults, the
  builder's types, the docs and the agent skill derive from it.
- **Unknown kinds** in an extension namespace the reader lacks are kept, reported and skipped, never a failure of the
  whole spec. A bare kind a `0.x` reader does not know is an error (every `0.x` minor may break), and a layer type of a
  newer version is read as its fallback (core's Extensions); before 1.0 an unknown bare kind in a newer minor of the same
  major becomes a warning, kept and skipped with the kind's declared fallback (not built: before 1.0, the package's TODO under `repr`) (amended 2026-11-23: the code's way until 1.0, the author as recommended).
- **Revisions and patches** ([patches.md](patches.md)): a spec's revision is an immutable root named by its hash; a
  change is a commit (a root and its parent revisions), and the patch between two revisions is derived (`diff`).
  Patches (`set`, `unset`, `insert`, `remove`, `move`, `replace`, by addresses of ids) are how a writer asks for a
  change; a batch is all or nothing and returns its inverse; `merge3` merges two lines of edits by node identity, with
  conflicts as data. Columns are replaced whole by a patch, never edited element by element.
- **The builder is the main way to write a spec** (the author: "the main use case will will using the "builder" to
  create the views"), typed from the description; misuse is a type error where it is written.
- **No published JSON Schema, initially** (the author: "it could end up being bloated"); the description could derive
  one when a consumer needs it.

## Rejected alternatives

- MessagePack as the default (proposed by the author as a question, "maybe msgpack should be the default?"): it would
  make the default file unreadable to tools and viewers without a MessagePack reader and Seq\*'s extension types, the
  JBrowse reply's portability point; kept as the option for large inline data.
- A typed array as "only a storage choice", equal by value to a plain array: makes hashing and `diff` of large columns
  slow and subtle.
- A language and a version per level: ten versioned descriptions for one spec.
- S-expressions canonical (meso\*): rejected by the author for this language.

## Consequences

- `decode(encode(x)) = x` for both codecs, `json(msgpack(x)) = json(x)`, `parse(print(x)) = x` (columns aside) and
  `apply(a, diff(a, b)) = b` are tested properties at small sizes.
- The language reference (generated docs pages) documents the JSON for readers without the core.

## In the specs

| Decision | Where |
|---|---|
| One language, one version; the levels are containers | core › Versions |
| A newer major refused with an issue | core › Gate order of checkSeqSpec |
| Older specs migrated by tree-to-tree functions | design-only: not built (pre-release `0.x` moves with no migration, core › Versions; the version-migration hook is in `packages/seq/TODO.md`, under `repr`) |
| The tree and its values codec-independent | repr › repr: the tree core |
| Typed columns declared, not guessed; codecs and the builder give the typed array | repr › Codecs stay schema-free |
| Identity by bytes and element type | repr › Columns by element type and bytes |
| Missing is NaN in float columns, `null` in JSON; integer columns have none | repr › Param types |
| Undeclared columns tagged; declared ones plain given the language | repr › Where a column is plain and where tagged |
| Printing elides columns; the reader refuses an elided one | repr › S-expressions |
| JSON the default interchange | repr › JSON |
| MessagePack for large inline data, its extensions | repr › MessagePack |
| Large data belongs in sources, not inline | core › Bounds |
| S-expressions a printing, never stored | repr › S-expressions |
| One node shape | repr › One node shape |
| The language is data; validation, defaults, builder types and docs from it | repr › A language is plain data |
| Unknown kinds in an extension namespace kept, reported, skipped | repr › An unknown kind in a foreign namespace |
| Unknown kinds of a newer minor kept with the kind's fallback (before 1.0) | design-only: not built: the close-out's (the version-aware validate, `packages/seq/TODO.md` under `repr`; until then a bare unknown kind is an error and `fallback` is not acted on; a layer type is read as its fallback, core › Extensions) |
| Revisions, patches, `diff`, `merge3`; columns replaced whole | repr › Patches |
| The builder the main way to write a spec | core › Language builder |
| No published JSON Schema, initially | design-only: not built (until a consumer needs one) |
| The properties as tests | repr › Tests |
