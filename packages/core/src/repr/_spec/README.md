# repr: the tree core

`@molstar/core/repr` (`packages/core/src/repr/`): the codec-independent tree core languages are described and written
on: one node shape, values with identity, languages as data, validation with every issue, addresses with query segments,
a typed builder. Not `@molstar/graphics/repr` (visual representations): `repr` here is a representation of a spec.
Written from meso\*'s tree core (`@molstar/meso-core` in the molstar-ludo repository) with the delta of
`.seq/design/language-description.md`; the designs are `.seq/design/tree-core.md`, `format.md`, `addressing.md`,
`locations-and-queries.md`. This spec is what is built: rules and their reasons.

| Folder           | What                                                                                 |
| ---------------- | ------------------------------------------------------------------------------------ |
| `tree/`          | the node, values and typed columns, identity, the atoms' spelling, addresses, issues |
| `schema/`        | languages as data, `checkLanguage`, validation, normalization, defaults, docs        |
| `builder/`       | the typed builder's engine                                                           |
| `codec/`         | the codec interface and its conformance suite; no codec                              |
| `codecs/json`    | the JSON codec, the default interchange                                              |
| `codecs/sexpr`   | the S-expression printing and its reader                                             |
| `codecs/msgpack` | the MessagePack codec, on `data/msgpack`; its column extensions as a decoder hook    |
| `patch/`         | revisions, patches, diff and merge; the patch language; the conformance suite        |

**Imports go one way**: `tree` ← `schema` ← `builder`; `tree` ← `schema` ← `patch`; `codec` imports `tree` and `schema`
(the conformance trees are normalized); `codecs/*` import `tree`, `schema` (where a language's declared columns sit) and
`codec`; `codecs/msgpack` also `data/msgpack` (MessagePack's encoder and decoder). Nothing here imports outside
`@molstar/core` (`hash-functions`, `util/array`, `util/utf8` and `data/msgpack` are its dependencies). Each folder is an
entry of the package (`@molstar/core/repr/tree`, …, `@molstar/core/repr/codecs/sexpr`), exported as `./state` is;
`patch/`'s entry is `repr/patch.ts` (as `data/int.ts` is `data/int/`'s), reached through the wildcard export, its suite
`@molstar/core/repr/patch/conformance`. Why: a codec needs no language and the tree needs no codec. MessagePack's
encoder and decoder (`@molstar/core/data/msgpack`) and this codec moved down from `@molstar/io`, which re-exports both
at their old paths: a reader of JSON's binary twin needs no io dependency.

## The tree

- **One node shape**: `{ kind, id?, params?, custom?, children?, args? }`. `id` is local (unique among siblings), so an
  address is the path of ids. `args` holds a data form's positional values (a query, a theme expression); a node with
  `args` has no id and no children. Why: one recursive form every language is spelled in.
- **`custom`** is extension data on any node: kept, hashed and encoded, never validated against the kind. Its columns
  are encoded tagged (Codecs), so they come back typed with no language.
- **Kinds** are words, optionally namespaced: `^([a-z][a-z0-9_-]*:)?[a-z][a-z0-9_-]*$` (`ext:track`). Why: extension
  kinds stay apart from the language's and from each other's.
- **Reserved keys**: `kind`, `__proto__` and `$…` in params, custom data and records. `$` is the encodings' (JSON's
  tagged atoms and `$column`, the printing's `:$custom` and `:$column`).
- **Readers are bounded** at `MaxDepth` (512) levels; trees built in code are hashed, printed and validated by
  recursion.

## Values

- **Atoms** are lossless: float64 numbers with -0, NaN and the infinities, bigints, strings, booleans, `null`.
- **Typed columns** (`tree/column.ts`): `Int8Array` … `Float64Array` (not `Uint8ClampedArray`, whose writes clamp), and
  `CategoryColumn` (`{ dict, codes: Int32Array }`, `-1` missing). Why: real data without a per-value cost.
- **A category is canonical**: made by `categoryColumn` / `categoryFromCodes` only, its dictionary in the order of first
  use with no unused entry. Why: its identity is then the sequence of its values, so an encoding that writes the values
  (JSON's plain array) gives the same bytes back.
- **The tree owns a column it is given** (`toColumn`, `normalize`, the builder): no copy is made, its NaNs are made the
  canonical NaN in place, its hash is cached per column object as a node's is per node, so it is never mutated after.
  Why: a large column costs no copy, and its hash is computed once.
- **One NaN**: every NaN of a float column is the canonical NaN (any payload or sign bit rewritten, `canonicalColumn`,
  once per column, also before hashing and comparing). Why: identity is by bytes, and NaN means missing, one value.

## Identity

- **`hashNode`**: MurmurHash3 x86 128-bit (core's) over a canonical text: record keys sorted, every atom kind told
  apart, children by their hashes (Merkle), cached per node. Absent `params`, `custom`, `children` hash as empty ones.
  Seed and text are meso\*'s, so a tree without the delta hashes as it did there (a pinned test). Why: equal across key
  order and encodings; an unchanged subtree costs nothing.
- **Columns by element type and bytes** (`hashColumn`, `columnsEqual`): NaN equals NaN, -0 is not 0, an `Int8Array` is
  not a `Uint8Array` of the same bytes, a typed column is not the plain array of its values. Bytes are the platform's
  (little-endian on every engine that runs this). Two columns whose cached hashes differ are unequal without a scan.
- **`valuesEqual`** compares atoms by `Object.is`, columns by bytes, lists and records structurally, nodes by hash.

## Atoms' spelling

One spelling for addresses and the S-expression printing (`tree/atoms.ts`), so an address reads the same alone and in a
spec: `true false nil nan +inf -inf`, bigints `123n`, a word that parses as a JavaScript decimal number is a number
(`-1`, `1e5`, `.5`), any other word is a symbol, a string (`94A`, `0x10`, `2nm`). A string prints bare when it reads
back as itself and holds no `/` (else quoted: `"HBA_HUMAN/1-141"`), except a canonical address. Why: residue labels and
quantities need no quotes, and no units table is needed to read them (locations-and-queries.md "Atoms").

## Addresses

- **Grammar** (addressing.md): `address = "/" | ("/" step)+ ["/*"]`, `step = id | form`,
  `form = "(" head value* (":" key value)* ")"`, `head = kind | operator`, `value = atom | list | form | address`. Read
  by the tree's own small reader (no codec), printed canonically on one line, a form's keys sorted:
  `(range :from 94 :to 292)` whatever order it was written in, so equal addresses print alike. `parse(print(a)) = a` up
  to key order. Messages never quote the address: an offset says where.
- **Positional values** (a form's `args`) come before its keyword params, in order: an expression's operands,
  `(filter :test (> (prop b-factor) 60))` (decided by Seq\*'s author, 2026-11-20, "Accept the core change"). Which forms
  take them is the language's that reads the form; the tree reads and prints them. A positional value after a keyword is
  `repr.tree.address.form`. **Operator heads**: `OperatorKinds` (`=`, `!=`, `<`, `<=`, `>`, `>=`) head a form beside a
  kind; never a node's kind in a tree (`KindPattern` is unchanged). Why: a predicate nests deep and reads as MolQL's.
- **A numeric-looking key is a number**: `(chain :id 1)` compares the number 1 and never matches a string id `"1"`;
  quote it to match one, `(chain :id "1")`. The atoms' rule decides, not the param's type (a form is read without a
  language).
- **The child index**: a lookup by id (`resolveAddress`, and inside `repr` `childIndex` and `childById` of
  `tree/address.ts`, not exported by `repr/tree`) scans up to `ChildIndexThreshold` (8) children; past it, a native
  `Map` from id to position, built on the first lookup and cached per children array in a `WeakMap`. The first of two
  siblings with one id wins, as a scan. Why: an edit's path copy makes a new array, indexed again on its first lookup,
  O(fan-out) as the copy was; no persistent collection is needed, since an insert shifts every later position.
- **Selection** (`selectAddress`) is bounded: a required `limit`, `maxDepth`, a `signal` and a `deadline`; more is an
  issue, never a truncation. Kinds and params of forms are checked against the caller's when given.
- **The form-resolver table** (`SelectOptions.resolvers`, by form kind): a resolver receives what the prefix names
  (`Named`: `{ kind: 'nodes', matches }` or `{ kind: 'other', type, value }`) with the caller's `limit`, `maxDepth`,
  `signal` and `deadline`, and returns a result or issues. A form without a resolver is meso\*'s child selection
  (`matchesForm`), which applies to nodes only; an id, a resolver-less form or `/*` after a result that is not nodes is
  `repr.tree.address.resolver`. Why: one grammar, each language deciding what its forms mean.
- **A selected node's address is segments**: ids, or a form where a resolver names a node no id path does
  (`/data/genes/(node :key T1)`, a dynamic node).
- **A resolver's errors fail the selection; its warnings and infos are carried** into the selection's issues. Its
  `complete: false` and `pending` are carried too, reserved for evaluation over data read by region (not acted on).
- **Staged selection** (a pattern with a resolved form) goes segment by segment: the `limit` bounds the final matches,
  not a step's; the signal and the deadline are read between steps and every `CancelEvery` nodes inside them and the
  `/*` walk, which keeps one path stack per root. Without resolvers in the pattern, the walk is meso\*'s (one path
  stack, an address copied per match).
- **Matching** (`matchesForm`): the same atom (`Object.is`), an atom held by a list param, an equal list, or a node
  fitting a nested form. Defaults are not applied.
- **Relative addresses** (`./part`) are not read: they belong to definitions, which wait.

## Issues

- **One shape**: `{ code, severity, message, hint?, offset?, path?, address?, count?, places? }`; `count` and `places`
  (the first ten) are for deduplicated issues, whose dedupe is not built yet. `code` is stable and namespaced
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
  (format.md "Unknown kinds"). An unknown kind without one is an error. Either way its children are checked with the
  language's kinds, with no parent rule for them (what sits under it is the extension's); normalization goes on below it
  too. Why: a spec with one extension node keeps every check of what the language does know.
- **Extension params**: a param key in a namespace other than the language's (`rondo:bundle`) is an extension's: kept,
  carried and not checked, on any kind. A kind's `extensionParam` names the param whose value, when namespaced
  (`type: 'rondo:ribbons'`), makes the node an extension's: its undeclared params are then kept unchecked too (the
  declared ones are checked as ever); so does its `fallbackParam` when given (`fallback`: a newer type's params kept by
  an older reader, which warns elsewhere). Why: an extension's params ride on the language's nodes, checked by the
  extension that owns the namespace; a reader without it keeps them, so a spec prints and patches the same with and
  without it.
- **A union** reads a node by the member that takes its kind, anything else by the first member it fits quietly.
- **Refs**: their spelling only (an address that parses, or an id). What they name is a resolve pass's, not built.

## Normalization

- **Codecs stay schema-free; reading types the columns**: `normalize(language, tree)` turns every value of a `column`
  param (a table's declared columns included, inside lists, records, maps, unions and slots) into its declared typed
  array; best effort, a value that does not fit is left for validation to report. `validateAndNormalize` is the read
  pipeline's last step (decode → migrate → validate and normalize). Why: every tree the passes see has typed columns,
  and equality, hashing and printing are defined on normalized trees.
- **Structural sharing**: an unchanged subtree is the same object; a normalized tree normalizes to itself.
- **`normalizeValue(type, value)`** normalizes one value as `normalize` reads a param of that type: what a patch's `set`
  writes. `validateChanges(language, root, changes)` is validation of what a change introduced (Patches).

## Codecs

- **The interface** (`codec/codec.ts`): `encode`, `decode` (never throws: issues), `typedColumns`, `elides`. The
  guarantees: `decode(encode(x))` normalized is `x` (by hash) unless the codec elides; `encode` is deterministic.
- **The conformance suite** (`codec/conformance.ts`) is data and a function (`checkCodec`), no test framework, so a
  codec elsewhere runs the same suite from its own test. Its trees cover every atom, odd keys, slots, data forms, custom
  data, extension kinds, columns of every element type with their edge values, NaNs with a payload and a sign bit,
  columns in custom data and in an `any` param, a table and a deep tree; its columns are short (at most 8 values).
- **Where a column is plain and where tagged**: given the language (`language`, `languages` options), an encoder writes
  a column in a declared column param (a table's declared columns included; `declaredColumnPlaces`) as a plain list,
  which normalization types on read. Every other column (in `custom`, an `any` param, an undeclared table column, or all
  of them without a language) is tagged, and decodes typed with no language. Why: JSON stays plain where a reader knows
  the language, and no column is ever lost to a plain array nothing types.
- **JSON** (`codecs/json`): `{ kind, id?, params?, custom?, children?, args? }`; atoms JSON lacks as
  `{ "$number": "NaN" | "Infinity" | "-Infinity" | "-0" }` and `{ "$bigint": "…" }`; a plain column a plain array, a
  float column's NaN `null` (its infinities and -0 tagged), a category's values strings and `null`; a tagged column
  `{ "$column": "int32", "values": […] }` or `{ "$column": "category", "dict": […], "codes": […] }` (malformed:
  `repr.tree.value.tag`). Strict: a key given twice in an object is `repr.codec.json.duplicate_key` (JSON.parse would
  keep the last silently), every issue with its path; a byte order mark skipped.
- **S-expressions** (`codecs/sexpr`): a printing for people, never stored.
  `(kind id :param value … :$custom (…) child …)`; one line when it fits in 80 columns, else one argument a line; a word
  opening with `/` is an address, read through its forms. A plain column prints as a list, a tagged one as
  `(:$column int32 :values [1 2])` or `(:$column category :dict [a b] :codes [0 1 -1])`, both up to the head (8 by
  default); past it a column is elided in a brace form, `{float32 × 3,000,000: 0.5 0.5 … …}`, which no value reads as
  (braces are reserved) and the reader refuses (`repr.codec.sexpr.elided`). Data forms are read as such when the reader
  is told their kinds (`dataForms`).
- **MessagePack** (`codecs/msgpack`, `MsgPackCodec`, a `Codec<Uint8Array>`), on `data/msgpack`'s `encodeMsgPack` and
  `decodeMsgPack`/`parseMsgPackValue` and their extension hooks. A node is a map
  `{ kind, id?, params?, custom?, children?, args? }` (keys in that order), a list an array, a record a map without
  `kind`: JSON's shape, so `json(msgpack(x))` is `json(x)`; a reader tells the encodings apart by the first byte (a map,
  never `{`). What MessagePack lacks is an extension (type, payload): `0x10` NaN, ±Infinity, -0 (one byte: 0 NaN, 1
  +Infinity, 2 -Infinity, 3 -0); `0x11`–`0x18` a typed column (the elements' bytes, little-endian: int8, uint8, int16,
  uint16, int32, uint32, float32, float64); `0x19` a category column (u32 dictionary size; each entry a u32 byte length
  and its UTF-8; the int32 codes); `0x1a` a bigint (its decimal digits). Why: `encodeMsgPack` writes a non-finite number
  as nil and -0 as 0, and has no bigint; the specials are one byte, being rare. A column is always tagged
  (`typedColumns: true`).
  - **Zero copy**: the encoder writes a column's bytes where they are; the wire is little-endian, a big-endian host
    byte-swaps a copy on both sides (`swapBytes`). The decoder gives a column as a view into its input where the
    payload's offset is a multiple of the element size, else one copy; never writing into its input, so a float column
    holding a NaN is copied (the tree makes NaN canonical in a column it owns; a foreign encoder may write a payload). A
    tree decoded in place shares its input's memory: the caller does not change the input afterwards.
  - **Reading what is not a tree**: `decode` never throws. Bytes that are not one MessagePack value are
    `repr.codec.msgpack.syntax` with the offset reached; only the parser's stack overflow is `repr.tree.depth_exceeded`;
    a key given twice is `repr.codec.msgpack.duplicate_key`, a key not a string `repr.codec.msgpack.key_type` (a second
    pass finds them, the decoder keeping the last and stringifying), each with offset and path. A decoded value that is
    no tree has the JSON codec's issues (a `bin` is no column: only an extension the decoder made is; a malformed
    extension is `repr.tree.value.tag`). Strings are read as values (a leading U+FEFF kept), the extensions' UTF-8
    strict; encoding refuses a lone surrogate (it throws): UTF-8 cannot hold it.
  - **The column extensions for a plain value** (`msgPackColumnDecoder`, a `MsgPackExtDecoder`): `0x10` and
    `0x11`–`0x19` read as above, a payload that does not read throws, any other type is the decoder's `{ type, data }`;
    maps stay maps, never nodes. Why: a reader of records (a table as MessagePack) keeps typed columns without the
    tree's meaning of `kind`.

## Builder

- **`createBuilder(language)`**: a kind's method on the bodies of the kinds it sits under, an id as the kind says, its
  params the description's `ParamValues` (a column takes its typed array or a plain array; a number with a unit a number
  or a quantity string). Misuse is a type error where it is written, and thrown at run time where types do not reach (an
  id twice under one parent, a malformed id, a kind under the wrong parent).
- **It builds the normalized tree**: `root` and `forms` return what `normalize` gives, so a built tree and the same tree
  read from JSON are equal. `check` is `validate`. Printing is the codecs' (`printSexpr`).
- **Handles are addresses**: a kind's method returns the new node's address (`/data/p53`) when it has one (an id under
  the root or under a node with an address), else `undefined`; the types say which (a slot's body, or a body under a
  node without an id, builds no addresses). A language's builder wraps them.
- **It imports `tree` and `schema` only** (a test): a script that writes specs loads no codec and no language.

## Patches

`repr/patch` (`patch/`; entry `repr/patch.ts`). Git's model: revisions are the truth, patches derived.

- **A revision** is an immutable root named by its `hashNode`. A patch is what `diff(a, b)` says takes `a` to `b`, or
  what a writer asks for; applying one gives a revision, and the revision is what a history records (its patch from a
  parent is `diff(parent, root)` when wanted). The core has no commit type and no undo stack: a holder's history is its
  own. The hash (MurmurHash3) names revisions and is no security boundary.
- **A patch**: `{ ops, base?, result? }`, no metadata (that is its holder's). Checks in order: `result` equal to the
  root's revision is already applied (`changed: false`, the caller's root, no inverse); else a `base` that is not the
  root's revision is `repr.patch.stale` (never rebased automatically); else, with `requireBase`, no `base` is
  `repr.patch.base`.
- **Six ops** by address of ids (`[]` the root; a form or `/*` is never a target): `set` and `unset` one key (`where`
  params by default, or custom; a record param set whole; an equal value or an absent key is a no-op keeping the node; a
  reserved key `repr.patch.key`); `insert` a node with an id unique among its new siblings at a `place` (after a
  sibling, first with `after: null`, or last); `remove` a subtree; `move` to a parent and place, `rename` giving a new
  id (never under itself: `repr.patch.cycle`); `replace` a whole node keeping its id, its kind free (equal by hash: a
  no-op keeping the old object). On the root: `set`, `unset` and `replace` only (`repr.patch.op`). Under a parent with
  any child without an id, `insert`, `remove` and `move` are `repr.patch.unaddressable`: a position after a sibling
  cannot be restored among siblings with no names, so that parent changes by `replace`.
- **All or nothing** (`patch/work.ts`): a draft over the root holds the nodes along the touched paths, each knowing its
  origin (the old root's node and address, or made). Every op is checked whole before it changes anything, and the first
  refused op refuses the patch: the caller's root is never touched. `finalize` shares every untouched subtree. Bounds:
  `maxOps` (10,000; `Infinity` for a diff, so law 1 holds on any diff, `repr.patch.limit`), and depth(parent) +
  height(node), a value's nesting counted, at most `MaxDepth` (`repr.patch.depth`). A materialized parent keeps its kids
  by id in a `Map` (the first of two with one id) and a count of kids without one, so a lookup, a taken id, a place and
  the addressable check are O(1) and many ops on one wide parent are not quadratic (a splice stays O(fan-out)).
- **The inverse** is each op's inverse against the state just before it, in reverse order, with `base` the new revision
  and `result` the old.
- **With a language** (`ApplyOptions.language`, slots through `validate.languages`), in order after the ops:
  - **Normalization** of what the patch made: an inserted or replaced node by `normalize`, a `set` value against its
    kind's param by `normalizeValue`. It costs what the patch made. Without a language nothing is normalized.
  - **Refs follow a move or a rename**, op by op, so a ref an earlier op wrote follows a later move. Walked: values
    typed `ref` (not with `in`) in params, inside lists, records, maps, unions (the member validation would choose),
    tables, data forms' args and slots; not `any`, `custom`, or addresses in strings. A ref matches by its leading id
    run, forms after it kept; addresses nested in its forms follow too. A rewritten ref prints canonically; an untouched
    one keeps its spelling. Without a language nothing is rewritten, and the result says `refs: 'not-rewritten'`.
  - **A stored patch** (a slot of a patch language) is walked by its own rule, not by its params' types: `at` and
    `parent`; the refs of an embedded node, by the target language; a `set` value when the target language types its key
    as a ref in any kind that has it. Typing reads the language alone, never the spec: whether a value is a ref never
    changes under a patch, so a move takes along the same refs forward and back, and no edit walks the steps again. The
    dry run (the spec with the step's earlier ops applied) decides only what a ref names: `at` and `parent` before their
    op, an embedded node's and a `set` value's refs through it. A stored patch inside an op's node or value has its own
    dry run, started from that op's state. So a step that inserts a node and sets a key on it is not dangling, and a
    stored patch follows renames and counts in link integrity as any ref: left dangling or rebound, it refuses. An op of
    a step that no longer applies (its target gone, its sibling renamed, its own move refused) is the step's to report,
    never a refusal of the edit that caused it. `place.after` (a bare sibling id) is not rewritten.
  - **The ref index** (`patch/refs.ts`, internal): `RefIndex` answers `under(prefix)` and `dangling()`; this one is on
    plain `Map`s, cached per root, language and options. Its sites are kept in chunks by node object (a node's own
    sites, a chunk per child that has any): `derive` walks only the paths the ops changed and shares every chunk of a
    node that is the same object at the same position (a subtree without sites is skipped at any position; no site reads
    the root, a stored patch's included; `derive` is by identity and position and reads no `Changes`). `applyPatch`
    leaves the new root's index pending when the old root has one: the nearest index that exists, never another pending
    one, derived when the next apply (or a read) needs it, so a chain of applies builds once and a run of sets that read
    no index walks nothing (always walking a typed-set step cost a one-key set 3 to 30 ms on a 10,000-ref spec with 30
    steps). A pending entry keeps that index alive (its root, chunks and caches): one revision, however long the run.
    `under` groups the sites by first id on the first query with a prefix. The caches hit only for the same `validate`
    object a caller passes (held weakly). Measured on 10,000 refs: building the old root's index was 4.5 ms of an 8.4 ms
    apply, deriving it is about 0.3 ms; the persistent `HashMap` of `data/immutable` is not needed at that size. On a
    p53-size spec with 30 stored patches of 4 ops: a rename about 0.6 ms, 10 reorders 0.08 ms, 10 renames 0.6 ms
    (`scale.test.ts`). A site's `path` holds child indices, so an insert or remove in a wide parent rebuilds the sites
    of its later siblings (nothing reads `path`). An apply has one, derived across its ops (`LiveRefs`): the old root's
    index for what the ops left alone (a site counts while its node is where the moves so far put it, still the old
    root's, and its key untouched), and an overlay, by the draft that holds it, for each key an op set or a rewrite
    changed and each node an op made. The old root's index is built only when an op changes an address (a move to
    another parent, a rename, an insert, a remove, a replace); a patch of sets and reorders builds none.
  - **Validation** (`validateChanges` in `repr/schema`, on the net `Changes`): a made node whole (its subtree, kept
    children under a changed kind included, its parent rule, a root's kind); a node whose params changed, for the
    changed keys and that no unset key is required; a moved node, for its parent rule; a changed child list, for the
    uniqueness of the ids the change brought. A path copy is checked for nothing, so an error present before is never
    reported again and only what the patch introduced refuses it.
  - **Link integrity**, the old root against the new, never op by op. Read through the apply's index: the refs the ops
    wrote, made or rewrote, and the old root's under an address the patch changed (removed, made, moved from or to). A
    stored patch's ref under none of them is not judged, though an edit may change its dry run: that is an op of the
    step no longer applying, the step's to report.
    - Dangling: a ref naming no node that named one in the old root or was written by the patch (made, or a key an op
      set), `repr.patch.dangling`, a warning with `allowDangling`. A ref in made content that was in the old root with
      the same value at the same address and dangled there is not written: a replace that keeps it is not refused.
    - Rebound: a ref the patch did not change whose target changed identity, `repr.patch.rebound` unless `allowRebind`:
      it named no node and names one now, or the node it named was removed or moved away and its address names another
      node. A ref a move rewrote is judged the same way, by what it named before and names now. Identity is lineage: the
      old root's node a draft descends from; a replacement is an edit in place (its node, and below it a node at an id
      the old node had, keep the old identity), an inserted node has none. So a remove then an insert at the address a
      ref names rebinds it, and one below what a ref names passes.
  - **The language's validators** (`validators`), last: `(before, after, changes) => issues`, an error refusing.
- **The result check**: a `result` other than the new revision is `repr.patch.result`.
- **`Changes`** (declared in `repr/schema`, re-exported here) is net: one entry a node whatever ops reached it; `made`
  and `children` in the new root, `removed` in the old (the top of each removed subtree, and a node moved into it from
  elsewhere at its own old address), `moved` from old to new; a made node is not also in `params` or `children`; a node
  removed after other ops changed it is only in `removed`; a reorder is in `children` only; ref rewrites are `params`.
- **The inverse's restoring sets**: every ref the inverse's moves would take along without the forward patch having
  rewritten it gets a `set` of its key's old value at the end of the inverse. Each ref is followed from where the
  inverse has it (the new root, the op restoring its key, or the content a remove or replace restores, by lineage)
  through the inverse's later moves; one that would not end as the old root has it is set back. So a move back never
  drags a ref the move forward did not touch.
- **The allowances turned around** (`inverseNeeds`): a patch that left refs dangling has an inverse that rebinds them
  (`allowRebind`); one that rebound refs, an inverse that needs both. The inverse also needs `allowDangling` when it
  writes a ref that dangled in the old root (a restoring set, or content it restores: an inverse's restored content is
  made, and the old root of the inverse does not hold it, so the "same place" exception cannot see it), each key it
  writes judged as the old root holds it (a restored step may hold an `at` a move rewrote, so its own dry run is not the
  old root's), and `allowRebind` when it makes again content an old root's ref named (from the inverse's side, a node it
  made). Law 3 holds with these only: the patch's own allowances are not the inverse's.
- **Issue paths**: an issue about the spec has the new root's path and address; one about an op the patch tree's
  (`['children', i, 'params', key]`, as `patchToNode` lays it out). Codes: `repr.patch.*` (`PatchIssues`).
- **`diff(a, b)`** (`patch/diff.ts`): `base` and `result` the two revisions. Merkle descent (equal hashes skipped);
  children matched by id within a parent; a kept id with a changed subtree descended into. Fallbacks: a kind change, a
  `replace` of that node; duplicate ids or a child without an id among either side's children, a `replace` of the parent
  (the nearest node with an address); a root of another kind or id, a `replace` of the root. Moves within a parent are
  the fewest that keep the longest run in order (a longest increasing subsequence, O(n log n)). Emission, so every
  `place` is valid in order: removes, then `b`'s children left to right (each not in the kept run placed after its
  predecessor in `b`), then the descents. Changes of parent and renames come out as a remove and an insert.
- **The patch as a tree** (`patch/language.ts`): `patchToNode` lays a patch out as `(patch :base … op …)`, one child per
  op, `at` and `parent` plain address strings, a place as `:after id`, `:first true` or `:last true`, an embedded node a
  slot. `patchLanguage(language)` describes it for a language (memoized; id `<id>.patch`, embedded nodes slots of the
  language); `AnyPatchLanguage` over none (embedded nodes `any`). `patchFromNode` reads it back, a malformed op
  `repr.patch.op` at its path.
- **The conformance suite** (`patch/conformance.ts`, `@molstar/core/repr/patch/conformance`): `PatchCases`, each with
  its options, and `checkPatches(apply, cases?)`: the result and its codes as expected, the root by hash, the `changes`
  and `inverseNeeds` where a case gives them, the caller's root unchanged, and every inverse giving the root back
  applied with its `inverseNeeds` only (the case's allowances stripped). A language runs it with its own cases.

## Merge

`merge3(base, ours, theirs, options)` in `repr/patch` (`patch/merge.ts`): three-way, by node identity.

- **Only `ok: true` has a `root`**, with `changes` the net change from ours to it (`Changes`, so a holder recompiles
  only what the merge changed). `ok: false` has the `conflicts` and a `partial`. Why: a caller cannot commit a merge and
  ignore its conflicts by accident.
- **Keyed by node**: `diff(base, ours)` and `diff(base, theirs)`, each change keyed by the address in the base it
  touches: a set's key on its node, a remove's or replace's node, an insert's new node, a reorder's parent. `diff`
  infers no move across parents and no rename, so every address it writes stands in the base, and in the other side
  unless that side removed or replaced it or above it.
- **Changes that both apply**: to different nodes, to different keys of one node; equal changes (the same value, the
  same node by hash, both removed, the same new order) once. A remove under what the other side removed is in it.
- **Conflicts**, each row both ways round, one a key, node or parent: `set` (a key set to two values, or set and unset);
  `remove` (any change at or under what the other removed, but a remove: a set, an unset, an insert under it, a reorder
  of its children, a replace of it or under it); `insert` (one id under one parent, other content); `replace` (any
  change under it, or another replace of it); `order` (two different new orders of the children on both sides). `move`
  is reserved for a merge of written patches: `diff` gives no move across parents and no rename, so `merge3` never gives
  one. `at` is in the base, but for `invalid`; `ours` and `theirs` the ops of the two diffs, in their order.
- **Orders are compared as orders**, so the result never depends on which child `diff` chose to move: per parent, the
  children on both sides: kept from the base, or inserted on both with equal content. One side's new order of the kept
  children is taken, the same one once; two different ones, or a node both inserted placed two ways, conflict and ours
  is kept. Why: the merge is then symmetric up to the order of siblings one side alone inserted.
- **Placement**: a child only one side inserted (or kept while the other's remove of it conflicts) goes after its
  nearest predecessor on that side that still stands, else first; ours' before theirs' after one sibling. It never fails
  the merge.
- **Independence by node**: a change of theirs depends on theirs' changes that made what it names: the sibling its place
  names when theirs inserted it, and, with a language, the targets of the refs it writes (by their id runs, under an
  insert or a replace of theirs; every run of a ref, those of the addresses nested in its forms too). A stored patch's
  refs make none: they name nodes of its own dry run, and link integrity judges them. A change whose dependency is in a
  conflict is left out with it, listed in that conflict's `theirs`. Other changes to one node do not make a dependency
  (two keys of one node both apply).
- **Theirs as one patch on ours**, under the options (`applyPatch` with `maxOps` unbounded): removes, replaces, sets and
  unsets as they are, then per parent the merged child order by the fewest moves that keep the longest run, as `diff`
  emits them. So normalization, link integrity and the validators judge the merged root against ours: a ref of ours to
  what theirs removed, or two links making a cycle, are refused there. Nothing to apply gives ours, the same object.
- **`partial`** is ours with every change of theirs independent of every conflict, applied the same way; when that is
  refused, ours. **`invalid`**: the merged root refused, `at` the address of the first error in the merged root (or
  `[]`), `ours` empty, `theirs` every change left out, `issues` the refusal's; `partial` ours.
- **Issues**: the apply's warnings (dangling refs kept under `allowDangling`), and without a language
  `repr.patch.unvalidated` (an info): nothing was normalized, rewritten or validated.
- **Rebase**: a stale patch `p` made on `base` is `merge3(base, current, applyPatch(base, p, options).root, options)`,
  the same options on both, so `p`'s moves rewrite refs as they would have. The host's opt-in; `repr` has no helper.

## Tests

`_test/` per folder, quick lane: the tree (hash pinned to meso\*'s, columns, atoms), addresses (grammar round trip,
issues, selection, resolvers), the schema (every misuse by code, severity and address; columns and tables;
normalization; units hook; docs), both codecs through the conformance suite and their own cases, the builder (normalized
output, handles, ownership, runtime errors, `@ts-expect-error` misuse), MessagePack (conformance, every column type,
alignment, malformed input, keys, BOM, surrogates, the plain-value hook; decode throughput in the heavy lane: columns at
least 200 MB/s, a node-heavy document guarded at 20 MB/s), the issue catalog against the sources, and a short seeded
fuzz of the three readers (never a throw, always a tree or issues). Patches: the conformance suite, net `Changes`,
structural sharing, issue paths, the patch tree's round trip, and diff's three laws as seeded property tests (without a
language, and with one on valid trees with dangling and rebinding allowed); the child index. Merge: its laws as seeded
property tests (`merge3(b, x, b)` and `merge3(b, x, x)` are `x`, `merge3(b, b, y)` is `y`, symmetry without conflicts up
to one-side siblings' order, every `partial` valid), every conflict kind both ways round, dependents, placement, orders,
`invalid`, a rebase, and its result consumed by `applyPatch` and `diff`.
