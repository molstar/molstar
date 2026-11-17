# repr: the tree core

`@molstar/core/repr` (`packages/core/src/repr/`): the codec-independent tree core languages are described and written
on: one node shape, values with identity, languages as data, validation with every issue, addresses with query segments,
a typed builder. Not `@molstar/graphics/repr` (visual representations): `repr` here is a representation of a spec.
Written from meso\*'s tree core (`@molstar/meso-core` in the molstar-ludo repository) with the delta of
`.seq/design/language-description.md`; the designs are `.seq/design/tree-core.md`, `format.md`, `addressing.md`,
`locations-and-queries.md`. This spec is what is built: rules and their reasons.

| Folder         | What                                                                                 |
| -------------- | ------------------------------------------------------------------------------------ |
| `tree/`        | the node, values and typed columns, identity, the atoms' spelling, addresses, issues |
| `schema/`      | languages as data, `checkLanguage`, validation, normalization, defaults, docs        |
| `builder/`     | the typed builder's engine                                                           |
| `codec/`       | the codec interface and its conformance suite; no codec                              |
| `codecs/json`  | the JSON codec, the default interchange                                              |
| `codecs/sexpr` | the S-expression printing and its reader                                             |

**Imports go one way**: `tree` ← `schema` ← `builder`; `codec` imports `tree` and `schema` (the conformance trees are
normalized); `codecs/*` import `tree` and `codec`. Nothing here imports outside `@molstar/core` (`hash-functions` is the
one dependency). Why: a codec needs no language, the tree needs no codec, and MessagePack lives in `@molstar/io` beside
its library.

## The tree

- **One node shape**: `{ kind, id?, params?, custom?, children?, args? }`. `id` is local (unique among siblings), so an
  address is the path of ids. `args` holds a data form's positional values (a query, a theme expression); a node with
  `args` has no id and no children. Why: one recursive form every language is spelled in.
- **`custom`** is extension data on any node: kept, hashed and encoded, never validated against the kind. A codec
  without a language gives its columns back as plain arrays (normalization does not reach it).
- **Kinds** are words, optionally namespaced: `^([a-z][a-z0-9_-]*:)?[a-z][a-z0-9_-]*$` (`ext:track`). Why: extension
  kinds stay apart from the language's and from each other's.
- **Reserved keys**: `kind`, `__proto__` and `$…` in params, custom data and records. `$` is the encodings' (JSON's
  tagged atoms, the printing's `:$custom`).
- **Readers are bounded** at `MaxDepth` (512) levels; trees built in code are hashed, printed and validated by
  recursion.

## Values

- **Atoms** are lossless: float64 numbers with -0, NaN and the infinities, bigints, strings, booleans, `null`.
- **Typed columns** (`tree/column.ts`): `Int8Array` … `Float64Array` (not `Uint8ClampedArray`, whose writes clamp), and
  `CategoryColumn` (`{ dict, codes: Int32Array }`, `-1` missing). Why: real data without a per-value cost.
- **A category is canonical**: made by `categoryColumn` / `categoryFromCodes` only, its dictionary in the order of first
  use with no unused entry. Why: its identity is then the sequence of its values, so an encoding that writes the values
  (JSON's plain array) gives the same bytes back.
- **A column is never mutated** once in a tree: its hash is cached per column object, as a node's is per node.

## Identity

- **`hashNode`**: MurmurHash3 x86 128-bit (core's) over a canonical text: record keys sorted, every atom kind told
  apart, children by their hashes (Merkle), cached per node. Absent `params`, `custom`, `children` hash as empty ones.
  Seed and text are meso\*'s, so a tree without the delta hashes as it did there (a pinned test). Why: equal across key
  order and encodings; an unchanged subtree costs nothing.
- **Columns by element type and bytes** (`hashColumn`, `columnsEqual`): NaN equals NaN, -0 is not 0, an `Int8Array` is
  not a `Uint8Array` of the same bytes, a typed column is not the plain array of its values. Bytes are the platform's
  (little-endian on every engine that runs this).
- **`valuesEqual`** compares atoms by `Object.is`, columns by bytes, lists and records structurally, nodes by hash.

## Atoms' spelling

One spelling for addresses and the S-expression printing (`tree/atoms.ts`), so an address reads the same alone and in a
spec: `true false nil nan +inf -inf`, bigints `123n`, a word that parses as a JavaScript decimal number is a number
(`-1`, `1e5`, `.5`), any other word is a symbol, a string (`94A`, `0x10`, `2nm`). A string prints bare when it reads
back as itself and holds no `/` (else quoted: `"HBA_HUMAN/1-141"`), except a canonical address. Why: residue labels and
quantities need no quotes, and no units table is needed to read them (locations-and-queries.md "Atoms").

## Addresses

- **Grammar** (addressing.md): `address = "/" | ("/" step)+ ["/*"]`, `step = id | form`,
  `form = "(" kind (":" key value)* ")"`, `value = atom | list | form | address`. Read by the tree's own small reader
  (no codec), printed canonically on one line: `parse(print(a)) = a`. Messages never quote the address: an offset says
  where.
- **Selection** (`selectAddress`) is bounded: a required `limit`, `maxDepth`, a `signal` and a `deadline`; more is an
  issue, never a truncation. Kinds and params of forms are checked against the caller's when given.
- **The form-resolver table** (`SelectOptions.resolvers`, by form kind): a resolver receives what the prefix names
  (`Named`: nodes, or an earlier resolver's result) and returns a result or issues. A form without a resolver is
  meso\*'s child selection (`matchesForm`), which applies to nodes only; an id, a resolver-less form or `/*` after a
  result that is not nodes is `repr.tree.address.resolver`. Without resolvers in the pattern, the walk is meso\*'s (one
  path stack, an address copied per match). Why: one grammar, each language deciding what its forms mean.
- **Matching** (`matchesForm`): the same atom (`Object.is`), an atom held by a list param, an equal list, or a node
  fitting a nested form. Defaults are not applied.
- **Relative addresses** (`./part`) are not read: they belong to definitions, which wait.

## Issues

- **One shape**: `{ code, severity, message, hint?, offset?, path?, address? }`. `code` is stable and namespaced
  (`repr.tree.*`, `repr.schema.*`, `repr.schema.language.*`, `repr.codec.*`); `severity` is the code's, from its
  catalog: `error` the input is invalid, `warning` something is lost (kept, skipped, not reported), `info` never changes
  what is drawn. `path` is the step path (`['children', 3, 'params', 'x']`); `address` is beside it where the place has
  one.
- **Catalogs** (`defineIssues`): each code with its severity, description and hint; a test fails on a code used in a
  source and not cataloged. Messages name only the language's vocabulary, never what the input wrote.

## Schema

- **A language is plain data** (`Language`: `id`, `version`, `doc`, `root`, `kinds`); `defineLanguage<const L>` checks
  it (`checkLanguage`) and throws on a broken one: a programmer error. Its type gives the builder's types.
- **Kinds**: `doc`, `parents`, `ids`, `params`, `data` (a data form's positionals), `fallback` (what a reader without
  the kind uses; documented, not acted on yet).
- **Param types**: `string`, `number`, `integer`, `boolean`, `list`, `record`, `map`, `union`, `node` (a slot, its
  language by id, resolved late), `ref` (an address word, or a bare id with `in`), `any`, and the delta:
  - `column`: `{ of, missing?, min-length?, max-length? }`. A float column's missing value is NaN; a category's is
    `null` where `missing` says so; `missing` on an integer column is refused by `checkLanguage` until its validity list
    is built.
  - `table`: `{ columns, key? }`, columns `column` or `list`, each required by its own `required` (as a record's field);
    equal lengths checked; extra columns kept, each a column or a list; `key` values unique.
- **Units are a hook** (`UnitsHook`): with it, `checkLanguage` checks a param's `unit` and validation reads a quantity
  string in a number slot; without it, `unit` is documentation and a string there is a type issue. Why: units stay
  outside the core for now (tree-core.md).
- **Definitions wait**: no `def` type, no uses, no `(arg …)`, no containers of definitions.

## Validation

- **Every issue**, capped by `maxIssues` (default 100, then a `repr.schema.limit_exceeded` warning). Unknown params,
  missing required ones, types, enums, patterns, ranges, ids (missing, forbidden, duplicate among siblings), parents,
  slots (language, kind), refs' spelling, data forms' arity, columns' element types, tables.
- **An unknown kind in a foreign namespace** is `repr.schema.kind.extension`, a warning: kept, reported, skipped
  (format.md "Unknown kinds"). An unknown kind without one is an error.
- **A union** reads a node by the member that takes its kind, anything else by the first member it fits quietly.
- **Refs**: their spelling only (an address that parses, or an id). What they name is a resolve pass's, not built.

## Normalization

- **Codecs stay schema-free; reading types the columns**: `normalize(language, tree)` turns every value of a `column`
  param (a table's declared columns included, inside lists, records, maps, unions and slots) into its declared typed
  array; best effort, a value that does not fit is left for validation to report. `validateAndNormalize` is the read
  pipeline's last step (decode → migrate → validate and normalize). Why: every tree the passes see has typed columns,
  and equality, hashing and printing are defined on normalized trees.
- **Structural sharing**: an unchanged subtree is the same object; a normalized tree normalizes to itself.

## Codecs

- **The interface** (`codec/codec.ts`): `encode`, `decode` (never throws: issues), `typedColumns`, `elides`. The
  guarantees: `decode(encode(x))` normalized is `x` (by hash) unless the codec elides; `encode` is deterministic.
- **The conformance suite** (`codec/conformance.ts`) is data and a function (`checkCodec`), no test framework, so a
  codec elsewhere (MessagePack in `@molstar/io`) runs the same suite from its own test. Its trees cover every atom, odd
  keys, slots, data forms, custom data, extension kinds, columns of every element type with their edge values, a table
  and a deep tree; its columns are short (at most 8 values).
- **JSON** (`codecs/json`): `{ kind, id?, params?, custom?, children?, args? }`; atoms JSON lacks as
  `{ "$number": "NaN" | "Infinity" | "-Infinity" | "-0" }` and `{ "$bigint": "…" }`; a column a plain array, a float
  column's NaN `null` (its infinities and -0 tagged), a category's values strings and `null`. Strict, every issue with
  its path; a byte order mark skipped.
- **S-expressions** (`codecs/sexpr`): a printing for people, never stored.
  `(kind id :param value … :$custom (…) child …)`; one line when it fits in 80 columns, else one argument a line; a word
  opening with `/` is an address, read through its forms. A column prints as a list up to its head (8 by default) and
  elided past it (`[float32 × 3,000,000: 0.5 0.5 … …]`), which the reader refuses (`repr.codec.sexpr.elided`). Data
  forms are read as such when the reader is told their kinds (`dataForms`).

## Builder

- **`createBuilder(language)`**: a kind's method on the bodies of the kinds it sits under, an id as the kind says, its
  params the description's `ParamValues` (a column takes its typed array or a plain array; a number with a unit a number
  or a quantity string). Misuse is a type error where it is written, and thrown at run time where types do not reach (an
  id twice under one parent, a malformed id, a kind under the wrong parent).
- **It builds the normalized tree**: `root` and `forms` return what `normalize` gives, so a built tree and the same tree
  read from JSON are equal. `check` is `validate`. Printing is the codecs' (`printSexpr`).
- **It imports `tree` and `schema` only** (a test): a script that writes specs loads no codec and no language.

## Tests

`_test/` per folder, quick lane: the tree (hash pinned to meso\*'s, columns, atoms), addresses (grammar round trip,
issues, selection, resolvers), the schema (every misuse by code, severity and address; columns and tables;
normalization; units hook; docs), both codecs through the conformance suite and their own cases, the builder (normalized
output, runtime errors, `@ts-expect-error` misuse), and the issue catalog against the sources.
