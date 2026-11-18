# @molstar/seqx-core

The SeqViewSpec language and everything that runs on it without a DOM: the language description and builder, the
coordinate model and maps, the source protocol and the data store, queries, the resolve passes, draw lists, and the SVG
and recording backends. It imports `@molstar/core` and an allowlist of other Mol\* entries, each with its reason.

This spec holds the contracts of what is built, each a rule and its reason; code comments cite its sections. A section
is added in the change that builds what it describes.

## Versions

- The language has one version, `LanguageVersion` (`0.3.0`). Why: one language, one description, one number for readers
  to check.
- **Pre-release** (`0.x`): a change to the language moves the version with no migration (`0.1.0`: every child of the
  spec has an id, `layer` included; `0.2.0`: the `refs` container and the `ref` kind; `0.3.0`: the `snapshots` container
  and the `snapshot` kind). A spec in an earlier `0.x` reads as one in this one, and is reported where it no longer
  fits.

## Language

`lang/`: `SeqLanguage`, the description of SeqViewSpec on `@molstar/core/repr` (kinds, parents, params, docs), and the
checks it cannot express.

- **Kinds** (first slice): `spec` (root: `version`, `title`); the containers `data`, `views`, `refs` and `snapshots`
  (reserved ids of the same name); under `data`: `source`, `sequence`, `ranges`, `values`, `map`; under `views`: `view`,
  then `section`, `track` (`hidden`, default false: kept in the spec, its row collapsed; its layers keep their own),
  `layer` (`hidden`, default false: kept in the spec and not drawn, so hiding is a one-key `set`). A `ranges` node's
  `map: 'codon'` makes each item (its rows by `id`) a coding sequence of its own (Overlay); a track's `for_each` makes
  it a template (Resolved view, Templates; its kind param is `node_kind`, since `kind` is reserved, and is `ranges` or
  `values`: what a layer draws). Reference params (`target`, `of`, `from`, `to`, `follow`, `source`, `data`) are `ref`
  params naming the kinds they take; a reference is an address.
- **Every child of the spec has an id**: each kind's ids are `required` in the description (a `layer`'s too), and the
  language's check `seq.lang.id.required` covers what the description cannot, an extension's kind. The bodies of slots
  are no children of the spec (a stored patch's ops have no ids). Why: every node is addressable, so hiding a layer is a
  one-key patch, and `repr/patch` can insert, remove and move under every parent.
- **Positions** are 0-based, ranges half-open. A range with `start > end` wraps the origin where its sequence is
  `cyclic` and is an error elsewhere; a section's `range` follows the same rule. Strands are `+` or `-`. Exons of a
  codon map are in genomic order and do not overlap (`strand` says which way the protein runs); on a cyclic gene the
  last exon may cross the origin (`end < start`, ending at or before the first exon's start, EC-11). Why: one convention
  for every coordinate the language writes.
- **Defaults in params**: a codon map's `table` is the NCBI genetic code, 1 when absent; a sequence's `name` is the
  record its source knows it by, its id when absent. Neither is filled in by normalization, so a spec that leaves them
  out hashes as it did.
- **Gate order of `checkSeqSpec`**: (1) the version: a newer major is `seq.lang.version.unsupported` alone, since this
  version's description says nothing of a newer spec; (2) the description's validation with normalization (`repr.*`,
  columns typed, so a spec read from JSON is checked as one built); (3) the language's checks (`seq.lang.*`), only when
  step 2 found no error, so they read a valid, normalized tree. Why: each gate assumes what the one before settled.
- **What the checks cover**: one of two params (`letters` or `source`, `rows` or `source`, `values` or `source`, `uri`
  or `value`), an `index` only for a file, letters of the alphabet, row intervals, bounds and strands (the issue's path
  carries the row), what a map method needs, exons, a codon map's alphabets, section ranges, that a reference names a
  node of a kind its param takes, that a layer's data is what its type draws (reported once: a kind the param does not
  take is the reference's issue alone; a template's layer without data draws each instance's node), that a template's
  `source` names a source, that `follow` chains end. Extents are checked where the letters are inline.
- **Patches** are `@molstar/core/repr/patch`'s, applied with `SeqLanguage`: its conformance suite runs over the core's
  cases and Seq's own (a layer hidden, a track moved with a section's `follow` following a rename, a dangling ref
  refused, a layer without an id refused; with the checks of a change and the slots, a `seq.lang` refusal, a story step
  following a rename and counting in link integrity, a step with a base refused), and a 10-op patch on the p53 spec,
  applied and diffed back, is held under 1 ms in the heavy lane (`_test/patches.test.ts`). The slots' languages are
  `SeqValidateOptions` (`seqLanguages`: the stored patches' `seq.patch`), which validation, the builder and `applyPatch`
  take.
- **Checks of a change** (`check-changes.ts`, `checkSeqChanges`, `SeqPatchValidator`): the `seq.lang.*` checks as
  `repr/patch`'s validator, on what a patch's `Changes` reach: the nodes it made, set or moved; what names them, or
  named any node of a subtree it removed or moved away (walked in the old root); on through the user refs that name
  those, and through a ref with segments after it (`/refs/p/x` where `/refs/p` names a prefix of the node, at most the
  ref depth deep); a template's layers under a track whose params changed; the sections and refs of a chain a changed
  one closes. What names what is a reverse index of the spec's references (a persistent set of referring nodes by the id
  run they name), built once a spec and derived from the last revision's along its changes. The same checks on the same
  nodes before the change are taken away, by code, node (through moves) and place in the node. Why: a patch is refused
  only for what it introduced (an error already there never blocks an edit beside it), and the cost follows the change,
  not the spec. Held by tests that compare it with the whole check's difference, patch by patch, along a chain of
  revisions and over seeded random patches of ref chains; and in the heavy lane, at 16,000 layers, the first check after
  a load and a fan-out to every node no worse than the whole check, a one-node edit a small part of it.
- **Bounds** (`bounds.ts`, `SpecBounds`, `DefaultSpecBounds`; the core's `bounds` item, `SeqCoreConfig.bounds`): a spec
  comes from a URL, a link or an agent, so what it sizes is refused with a coded issue before the work starts.
  `checkSpecBounds` is the first gate of `checkSeqSpec`, before the version and the description: one iterative walk (a
  20,000-deep spec costs no stack) that allocates no tree and stops at the first bound the whole spec cannot pass.
  `maxBytes` (64 MB: strings by length, numbers 8 bytes, columns by their bytes: `seq.lang.bound.size`), `maxDepth` (256
  levels of nodes and of the values nested in params: `seq.lang.bound.depth`), `maxQueryDepth` (64 nested forms in an
  address string: `seq.lang.bound.query_depth`) and `maxColumnLength` (2^31 − 1: `seq.lang.bound.column`); a spec past
  size or depth gets that issue alone, one past the others the first place of each. The issues kept are bounded as they
  are reported (`IssueCaps`): the first `maxPlaces` (10) of a code, at most `maxNodeIssues` (100) a node, the rest
  counted in one `seq.lang.bound.issues` warning a code (`BoundIssues`; a patch's check keeps what it introduced the
  same way, after taking away what was there). The walk reads a node's children, params, `custom` data and `args` alike
  as values (no query depth in `custom`, the extension's) and charges a fixed overhead (8 bytes) for each node, record,
  array, key and element before it pushes or visits them, so a million empty records cost what they say; a node that is
  no record with a kind, or a `children` that is no list, is `seq.lang.bound.shape` (refused, never thrown on). A
  patch's check refuses the same nesting: `checkSeqChanges` walks the subtrees a change made or moved (to the depth the
  address already has) and the params it set, with `checkSubtreeDepth`, before it reads anything else (a subtree too
  deep to walk is refused alone); `seqPatchValidator(bounds)` is the validator with a composition's bounds, and
  `seqBoundsValidator(bounds)` the bounds alone (a preview, a story step). Its size bound is per patch: the bytes of
  what it makes and of the keys of `params` or `custom` it sets (only those keys are walked) against `maxBytes`; a moved
  subtree is walked for depth, not counted. Why: one rule for the two consumers, and the cost of the incremental check
  stays the change's. The table's other rows: nested sources are `maxNestedDepth` (Resolve), a column's length and a
  block's values the data's (Data, Tables), regular expressions in search and archives wait for the features that read
  them (TODO).
- **Open**: the theme slot of a layer's `color` and a query in its `data` wait for the sublanguages (TODO).

### Refs

`lang/refs.ts`: user-defined names for addresses (bookmarks, regions, sets).

- **A `ref`** sits in the `refs` container: `to` (an address, or a non-empty list of them: a set), `label`, `color` (a
  hex literal or `token:<name>`). Its address is `/refs/<id>`, and it stands wherever an address does: in commands and
  links, in other refs, as a prefix (`/refs/p53/(range :from 1 :to 9)` when the ref names `/data/p53`), as an operand of
  a form, and as a layer's `data` (one address; the resolved view and the checks follow it). Why one kind for a bookmark
  and a set: a bookmark is a ref, so adding one is a patch like any other.
- **Following** (`followRef`): the address's leading `/refs/<id>` is replaced by the ref's target, again until no ref
  leads it, the segments after it kept; a set is followed member by member. Refused, never guessed: a cycle, a chain
  past `DefaultRefDepth` (8) refs, a set or a subtree followed by more segments (they would apply to each member).
- **Checks** (`seq.lang.*`, warnings, so the rest of the spec opens): a ref whose address's leading ids, once the refs
  it starts with are followed, name no node (`ref.unresolved`; forms after them are a query's, and a node data makes
  waits for its data), a chain that comes back or runs past the bound (`ref.cycle`).
- **Adding one is one patch** (`addRefPatch`): an insert at the end of `refs`, or of the container with the ref in it,
  made against the root's revision. Refs follow moves and renames (`repr/patch`); a removal that would leave one
  dangling is refused unless allowed.
- **Not followed yet**: refs in data nodes' params (`target`, `from`, `to`, `source`) and a section's `of` and `follow`:
  a ref there names a node of kind `ref`, which the checks refuse (`seq.lang.ref.target`).

### Story steps

`snapshot` nodes in the `snapshots` container (`/snapshots/<id>`): `label`, and `patch`, a stored patch in `seq.patch`
(`patchLanguage(SeqLanguage)`, `SeqPatchLanguage`) with no `base`.

- **A step applies over whatever the spec is**: entering one is navigation, never an edit (the viewer's `SpecModel`).
  Why no base: a story is written once and shown over every later revision.
- **Its addresses are refs of the spec** (`repr/patch`'s stored-patch rule): a rename or a move rewrites them, and link
  integrity counts them. A sibling named by a `place` is not rewritten: such a step no longer applies after the rename,
  and is shown as not applying.
- The builder writes one with `b.snapshots.snapshot(id, { patch, label })` (a `Patch` or its node; a base throws).

## Language builder

`createSeqSpec(options)`: the way most specs are written; `build()` is the normalized tree, `check()` is `checkSeqSpec`
of it.

- **Handles are addresses, branded by kind**: `Address<'sequence'>` is the string `/data/p53` that says what it names.
  Every call that adds a node takes its id and returns its address, a `layer`'s too (`track.layer(id, params)`).
  Section, track and view builders carry their `address` and a `toString`, and stand where an address of their kind
  does.
- **A reference param takes** a handle (or builder) of a kind its description allows, or an address written out (any
  string); a handle of another kind is a type error. Which params are references is read from the description, in the
  types and when the builder converts a handle to its address. Why: no list of param names to keep in step with the
  language.
- **Misuse** that types do not reach (an id twice under one parent, a malformed id) throws where it is written; the rest
  is found by `check`. Calls are recorded and replayed on the engine of `@molstar/core/repr/builder` at `build`.
- **Type tests** (`_test/types.test.ts`) mark wrong kinds and params with `@ts-expect-error`; they are checked by
  `tsc -p tsconfig.test.json`, which the package build does not run.

## Issues

`seq.lang.*` codes (`LangIssues`, and `BoundIssues`, `seq.lang.bound.*`, the bounds'), each with a severity, a
description and a hint; a test fails on a code a source uses and the catalog lacks, and on one cataloged and never used,
and a misuse for each code is in `_test/lang.test.ts`. The namespaces are `repr.*` (the tree core's), `seq.lang.*`,
`seq.query.*`, and the later `seq.resolve.*`, `seq.source.*`, `seq.store.*`, `seq.view.*`, `seq.render.*`,
`seq.address.*`. An issue's `path` is the step path and its `address` the node's, where it has one.

## Queries

`src/query/` (`@molstar/seqx-core/query`): locations, element keys, the three result types, and the query forms as
`repr` form resolvers over a model of systems, numberings and annotations. Forms and addresses are one language: a form
is a segment of an address, `/data/p53/(range :from 94 :to 292)`.

### Result types

- **Nodes** (exact addresses of spec nodes and of nodes data makes), **Elements** (keys of what data makes: an
  annotation's items) and **Location** (positions on systems). A form gives one of them; they travel as `repr`'s `Named`
  (`{ kind: 'nodes' }`, `{ kind: 'other', type: 'elements' | 'location', value }`). A `Named` whose value lacks its
  type's shape is no query result (`fromNamed` checks it). Why: a hover, a click and a layer each hold the type they act
  on, and no form guesses a conversion.
- **Conversions are forms**: `(extent :of <elements>)` gives a Location, `(items :overlapping <location>)` Elements. A
  form given another type is `seq.query.form.type`.

### Location

- `{ parts: { system, ranges: Int32Array }[], rects: { a: part, b: part }[] }`: one part per system (sorted by address),
  its ranges sorted, merged (adjacent ones too), half-open `[start, end)` index pairs, none empty; a cyclic system's
  wrapped range is two ranges. Rects pair two parts (two axes), sorted, each once.
- **A value**: built normalized (`location`), frozen, never mutated. Equality is structural (`locationsEqual`): equal
  parts and an equal list of rects. Rects are never merged, so two rects covering the cells of one are another value.
  `hashLocation` is equal for equal Locations, cached per value.
- **Positions are int32**: `normalizeRanges` throws on an end that is not an int32 integer (a programmer error, never a
  wrap); the forms check data first and report `seq.query.position.invalid`.
- **Set operations** per system: union, intersection, subtraction. Rects intersect side by side over the same two
  systems; a rect is not subtracted, and its complement is refused (`seq.query.form.unsupported`).
- **Overlap** of an item `[start, end)`: a boundary (zero length) at `p` overlaps `[s, e)` when `s ≤ p < e`; one at the
  system's end (`p = length`) overlaps a range ending there, so `(all)` holds it. O(log n).

### Element keys

- `{ owner, key }`: the owner's canonical address and the data's own key, a string or a number. A number never equals a
  string (`1` is not `"1"`), as in an address. Elements are sorted (owner, numbers before strings) and each once;
  compared and hashed as values.

### Numberings and labels

- **Indices** are 0-based half-open inside; **labels** are inclusive names at the edge, converted once where a form is
  read.
- **Built-ins**: `index` (0-based) and `ordinal` (1-based), on every system; without `:numbering` a form reads
  `ordinal`. A declared numbering is looked up on the system, then on its molecule (`parent`, shifted by the unit's
  start, clipped to the unit): a unit's own beats its molecule's of the same id. Why: the author's numbering of a chain
  wins over the entity's. Not across a codon step: that waits for paths.
- **A label** is an integer with an optional suffix, letters (`94A`) or `.n` (`"100.1"`), or else a name. One total
  order: by integer; no suffix, then letters (shorter first, then alphabetic ignoring case), then `.n` by n; names after
  every number; spelling breaks the last ties. A `.n` label is written quoted: bare, `100.1` is the number 100.1 (the
  atoms' rule) and is refused with a hint to quote it; quoted, `"100.10"` and `"100.1"` are distinct labels. An integer
  outside int32 is refused (`seq.query.label.int32`). A repeated label names its first position.
- **Insertion codes**: `(at :label 94 :ins A)` is `94A`; `""`, `?`, `.` mean none; a code on a range, or on a label with
  a suffix, is `seq.query.label.insertion`.
- **Resolving a range**: each end is its label's position; a missing end snaps inward (a start to the first position
  whose label is at or after it, an end to the last at or before it, `seq.query.label.snapped`, info). Then one rule in
  index order for exact and snapped ends alike, never assuming labels rise with the index (IMGT): `from ≤ to` is the
  range; ends snapped past each other over positions with no label are empty (`seq.query.range.empty`); otherwise a
  cyclic system wraps as two ranges short of a whole lap (`seq.query.range.lap`), and any other is
  `seq.query.range.backwards`. A missing `:from` or `:to` is the system's own end.
- **Bounded label work**: a numbering may answer `snap` and `nearest` itself (`index` and `ordinal` do, in O(1));
  otherwise the positions are scanned, each a visit (signal, deadline). A hint lists the labels nearest an unknown one,
  kept a few at a time (no list of every label), and none past `HintScanLimit` positions.

### Forms

| Form                                           | Prefix            | Result                                                                                                          |
| ---------------------------------------------- | ----------------- | --------------------------------------------------------------------------------------------------------------- |
| `(range :from l :to l :numbering n)`           | systems           | Location                                                                                                        |
| `(at :label l :ins c :numbering n)`            | systems           | Location (one position)                                                                                         |
| `(all)`                                        | systems           | Location                                                                                                        |
| `(rect :a <location> :b <location>)`           | as its sides need | Location (one rect; each side one system; sides never wrap, `seq.query.range.wrap`)                             |
| `(item :id k)`, `(item :<field> v …)`          | annotations       | Elements (every param matches; a list param matches any of its atoms; a field the annotations lack is an issue) |
| `(items :overlapping <location>)`              | annotations       | Elements (every item without `:overlapping`)                                                                    |
| `(node :key k)`, `(chain :id k)`               | nodes             | Nodes, addressed by the form itself (`/data/1tup/(chain :id A)`), so location forms may follow a chain          |
| `(extent :of <elements>)`                      | any               | Location                                                                                                        |
| `(and :of […])`, `(or :of […])`, `(not :of x)` | any               | the operands' one type                                                                                          |

- **Combinators take `:of`**: a form's params are keywords only (the address grammar), so `and`, `or` and `not` read a
  list or one operand under `:of`. Operands of different types are `seq.query.form.type`.
- **`not`** complements within what the prefix names: its systems' positions, its annotations' items, its nodes'
  children with ids. Dynamic nodes are not in that universe yet.
- **Operands**: a nested form is evaluated against the same prefix (its resolver, else `repr`'s child selection); an
  address from the root. The location of `items :overlapping` is evaluated against the annotations' targets; a part on a
  system no target is reports `seq.query.location.unreachable` (projection is not built), a rect is refused. An operand
  naming a ref is refused (`seq.query.ref.unresolved`): addressing follows refs before the forms read an address
  (Addressing), so only a caller of `evaluateQuery` itself meets it.
- **The table**: `queryResolvers(model, extra?)`; an extra kind named as a query form throws (a form's meaning is never
  replaced silently). `evaluateQuery` is `selectAddress` with the table, typed.

### The model the forms read

- An interface (`QueryModel`): `system(address)` (length, cyclic, numberings, parent), `annotation(address)` (target,
  fields, `items(over?)` with what is pending), `node(owner, key)` and `chain(owner, id)` for nodes data makes. The
  static model implements it; tests use a hand-written one.

### Staged and bounded

- **Staged**: over data read by region an evaluation gives `{ result, complete, pending }`, `pending` a Location of the
  regions not loaded; a form unions its operands' pending regions. Across two resolved segments of one address the
  selection keeps the last segment's (the union is `repr`'s, to come).
- **Bounded**: the caller's `limit` bounds the nodes and elements a form gives (more is
  `repr.tree.address.limit_exceeded`, never a truncation); the signal and the deadline are read every `CancelEvery`
  visits, one counter across nested forms and nested addresses (`repr.tree.address.cancelled`, `deadline_exceeded`). A
  Location's size is not bounded yet.

### Issues

- `seq.query.*` (`QueryIssues`): `form.prefix`, `form.param`, `form.type`, `form.unsupported`, `numbering.unknown`,
  `label.syntax`, `label.int32`, `label.unknown` (with the nearest labels), `label.insertion`, `label.snapped` (info),
  `range.backwards`, `range.lap`, `range.wrap`, `range.empty` (warning), `position.invalid`, `location.unreachable`
  (warning), `ref.unresolved`, `node.missing` (warning). A resolver's errors fail the selection; its warnings and infos
  are carried. A test checks the catalog against the codes the sources use.

## Addressing

`query/address.ts`, `query/links.ts` (`@molstar/seqx-core/query`): an address resolved to what it names and back, and
links that carry addresses. Forms and addresses are one language, so this sits with the queries.

### Resolving

- `resolveAddressIn(root, model, address, bounds)`: the leading ref followed (`followRef`), the refs nested in form
  operands followed (`:of`, `:a`, `:b`, `:overlapping`; a set as `(or :of [...])`), then the query forms over `model`; a
  set's members give one type of result, together (else `seq.address.type`), gathered and normalized once (linear in the
  members). Why operands only: a quoted string reads as the same value as an address, and elsewhere a string is a value
  (an id, a label) that must never be rewritten. It answers a status, the result, the refs passed (`via`), the issues,
  and what is pending.
- **Status**: `resolved`; `waiting` when the result is empty and a node data makes is not held (`(node :key k)`), or the
  data is not loaded (`seq.address.waiting`, info: it resolves when they arrive); `missing` when it is empty otherwise
  (`seq.address.missing`, warning); `failed` on an error. An issue's `address` is the address asked, its `places` the
  refs passed and what they led to: an address that no longer resolves is reported with what it named, never re-bound.
- **Bounded**: `limit` (10,000 nodes or elements; more is `repr.tree.address.limit_exceeded`, never a truncation),
  `maxDepth` (64), `refDepth` (8); a signal and a deadline as a selection takes them.
- **As positions** (`locationOfResult`): a Location itself; elements by their items' extents; a system's node whole; an
  annotation's by all its items; a section (and what is under it) by the region it shows. What focus and reveal read.
- **The overlay's query model** (`overlayQueryModel(model, overlay)`, `resolve/query-model.ts`): the static model's
  systems and annotations, then the overlay's (item systems, dynamic annotations), and `(node :key k)` naming a node the
  store holds (a string key: `1` is not `"1"`); one per model and overlay. Why: a ref to a dynamic node resolves when
  its node arrives, through the same forms.

### Back

- `addressFor(root, model, result, { prefer, numbering })`: the shortest stable address, by `prefer` (`['ref', 'ids']`):
  a ref of the spec that names exactly it (equal Location, Elements, or set of node addresses; each ref resolved once
  per root and model, `refResults`), else ids and forms. Positions as `(range :from :to)`, one as `(at :label)`, a whole
  system `(all)`, a cyclic range over the origin one wrapping range, several under `(or :of [...])` (from the root when
  they span systems); a rect `(rect :a … :b …)`; elements `(item :id k)` or `(item :id [k …])`; nodes by their
  addresses.
- **Labels** in the numbering asked (`ordinal` by default; an integer label a number, others strings). The address is
  resolved again and kept only when it reads back to the same result; else `ordinal`, which always does (a repeated
  label, a hole). Why: an address that names something else is worse than a longer one.

### Links

- **Mol\*'s link grammar** (`printLink`, `parseLink`): `!key=value&…`, each key and value percent-encoded as a whole
  (`&`, `=`, `#`, spaces, brackets, parentheses, quotes, `!` and `*` too, so an address with `)` survives in a markdown
  destination; `/` kept readable: `!seq-focus=/refs/dbd`), split before decoding. Keys of the prefix (`seq-` by default)
  are actions (`focus`, `reveal`, `select`, `highlight`, `hover`) and `channel`; other keys are kept in order for the
  viewer that knows them. A link is data: never markup.
- **Refused**: an address that does not parse fails the link (nothing runs from a half-read link), a bad encoding
  (`seq.address.link.syntax`), a link past `maxLength` (8,192) or with more than `maxActions` (16) actions
  (`seq.address.link.limit`); an unknown action of the prefix is left out with a warning (`seq.address.link.action`).
- **Tension**: Mol\*'s own reader decodes the whole link before splitting, so an address holding `&` or `=` (only inside
  a quoted string) reads wrong there; this reader is exact.

### Configuration

| Item                           | Default          | What                                        |
| ------------------------------ | ---------------- | ------------------------------------------- |
| `limit`                        | 10,000           | the most nodes or elements an address names |
| `maxDepth`                     | 64               | how deep a walk goes                        |
| `refDepth`                     | 8                | how many refs a chain passes                |
| `prefer`                       | `['ref', 'ids']` | the forms the way back tries, in order      |
| `numbering`                    | `ordinal`        | the labels the way back prints              |
| link `prefix`                  | `seq-`           | this viewer's keys in a link                |
| link `maxLength`, `maxActions` | 8,192, 16        | a link's bounds                             |

### Issues

`seq.address.*` (`AddressIssues`): `missing`, `waiting` (info), `ref.cycle`, `ref.depth`, `ref.set`, `ref.empty`,
`type`, `link.syntax`, `link.action` (warning), `link.limit`; a test checks the catalog against the sources.

### Performance

A set of 16,000 addresses resolves in about 7.5 times the time of 2,000 (the heavy lane bounds the ratio at 16; a
pairwise union was 62). Resolving 1,000 addresses (ids, forms and refs mixed) on the p53 spec: 6 ms median, 5.6 to 12.9
over 20 runs; the heavy lane bounds one sample at 40 ms (`resolve/_test/addressing.test.ts`).

## Data

`src/data/` (`@molstar/seqx-core/data`): the source protocol, the `seq` format and its runtime, the data store.

### Protocol

- A **format runtime** (`FormatRuntime`: `name`, `formats`, `open(source, { assets, signal })`) opens a `SourceEntry`
  (`address`, `format`, `uri` or inline `value`, an `index` for a format read with one) to an `OpenedSource`: its
  `coverage` and `read(request, signal)`. `assets` is a narrow port (`readText(uri, signal)`). Why: a runtime is tested
  with a fake, never the network. A file that cannot be read is a `SourceOpenError` with `seq.source.unreadable`.
- **Coverage**: systems (unique names, lengths ≥ 0, alphabet, cyclic, aliases), shapes per listed system (with a ranges
  table's fields), integer levels ≥ 1 as bin widths, ascending, at least one (`[1]` without levels), nodes known up
  front (`checkCoverage`).
- **An answer** is tables by shape and the dynamic nodes the region brings (`{ key, node }`, a spec node of its kind),
  or a failure: `failed` (retried with backoff), `malformed` (fails at once, with a cooldown), `not-served` and
  `too-large` (never retried; too large asks a coarser level). Cancellation rejects with an `AbortError`, never an error
  state; an `AbortError` the store did not ask for is the source's failure.

### Tables

- Typed columns, each table with a `form`: `letters` (`start`, a string or uint8 codes with an alphabet); `ranges` as a
  `list` (`start`, `end` int32, `strand` int8, `id` category, `row` int32, `fields` as `{ name, column }`) or `density`
  (`bin` int32, `count` float32, `category`); `values` as `positions` at level 1 (`value` float32 from `start`) or
  `bins` above it (`min`, `max`, `mean`, `count` float32 from `start`, aligned to the level). Why a `form`: one shape
  has two encodings by level. Why fields as a list: any name, `kind` and `$x` too, in memory and in a file.
- **A row's identity** is its `id`, else its `row`: a range over a tile's edge is in each tile it touches, and readers
  dedupe by it.
- `checkTables` checks forms, column types and equal lengths. `checkAnswer(answer, request)` checks an answer against
  its request: a record, a known failure with an issue, nodes with string keys, only the shapes asked, each table within
  the region, ranges with `start ≤ end` and strands −1, 0, 1, codes within their alphabet, at most `maxRows` rows.
- **Bounds**: `checkAnswer` refuses a table whose column is longer than `maxColumnLength` (2^31 − 1) or whose values
  (rows times columns) pass `maxBlockValues` (2^24, a grid or matrix block's bound), before it walks the extent, and the
  store settings carry both; a file's `seq` reader refuses the same (`checkTableBounds`) and a system declared longer
  than a column can be. The dynamic nodes of an answer are at most `maxNodes` (10,000) and weigh at most `maxNodeBytes`
  (16 MB) together by the spec's size estimate, each within the spec's depth and shape (`measureNode`); `answerBytes`
  charges them by that estimate, a flat floor under it. Why: an answer's size is a claim of its source, refused before
  it is cut, joined or kept.
- Cutting (`sliceTable`) keeps letters, positions and bins by offset, and a range or a density bin by its anchor (its
  first position in the region) within the part cut; `joinTables` joins consecutive parts, so each row is held once.

### The seq format

- The protocol written down: a coverage and `{ request, answer }` pairs as a tree value. A file is a `seq` node through
  the JSON codec (typed columns tagged, atoms lossless); inline, the same record. A failure's `kind` is written
  `failure` (a record with `kind` reads as a node). Reading reports the first misfit (`seq.source.malformed`, at its
  path) and never throws. Why a tree value: one lossless encoding, the codec's, and MessagePack with it later.
- **Its runtime** joins the pairs of the request's system (by name or alias) and level that overlap the region, in
  position order, each cut to its part; a gap no pair holds is `not-served`, naming its positions (the format is
  explicit: an empty region is a pair with empty tables). A failure among the pairs is the answer; an unlisted system or
  level is not served. A list without ids gets `row` from its place in the file. Inline data (`inlineSeq`) is one pair
  over the whole system at level 1, its length in the coverage.

### Asynchronous sources

- `delayedSource(inner, { clock, delayMs })` answers as `inner` does, each answer after a delay (a number, or by the
  request) on a `StoreClock`; an abort rejects at once and clears the timer; 0 is no timer. `asyncSeqSource(data)` is
  the synthetic by-region source (`seq` data answered so); `delayedRuntime(runtime, { formats })` gives it a format of
  its own (`seq-async`). Why: the store's states (pending, cancelled, evicted, failed) are driven by a fake clock in
  tests and seen in a story with the real one, through one source.

### Summaries

- A source that lists level 1 alone and serves values on a system of known length gets coarse levels when attached
  (`withSummaries`, the store's `summaries` setting; `false`: none): `[1, f, f², …]` (`factor` 16) up to its longest
  values system. A system up to `summaries.maxLength` (2^24) gets a pyramid: read whole once at level 1 (checked by
  `checkAnswer` with the store's `maxRows`), every level built from the one below (`binValues`, `coarsenBins`: min, max,
  count-weighted mean, NaN not counted), then the level-1 values let go; its bytes (`summaryBytes`) count against
  `maxBytes` and go with the source. A longer system's coarse tile is binned on the fly from a level-1 read of that tile
  alone, refused (malformed) past `maxRows` positions. A caller's abort is noticed at once, even while the shared read
  goes on. Other shapes at a coarse level are read at level 1 over the same region. Why: a values layer zoomed out reads
  a coarse level within the tile bound, a level costs one pass over the one below, and memory stays within the bound.

### FASTA

- `FastaRuntime` (`fasta`): the file read whole through the assets port (range reads wait for bigWig); with an `index`
  (`.fai`: name, length, offset, bases and bytes a line) the records are found by it once it is checked against the file
  (five integer fields, unique names, each record after its header line, full lines of the stated size, ending within
  the file), else refused (`seq.source.malformed`); without one, the index is built from the file (`faidx`, samtools'
  layout). A file that is not ASCII is refused: offsets count bytes. Coverage: a system a record, letters, level 1.

### Store

- **Keys**: `<source>|<system>|<level>|<tile>` and `<source>#<key>`; a part holding `|`, `#` or `\` escapes it, so no
  two keys of other parts are equal. A tile is `tileBins` bins (1,024) at the level's width; a want is clipped to the
  system's length when known and snaps to whole tiles; an unlisted level falls back to the coarsest listed no wider
  (EC-170). A want over `maxTilesPerWant` tiles is `seq.store.too_many_tiles` and not asked: work is bounded before it
  starts. A tile is asked for every shape its system serves: the key has no shapes.
- **Sources** are checked on `attach`: a malformed coverage is an issue and the source unavailable. `detach` cancels its
  requests, drops its entries, releases its holds and clears its reports.
- **Entries** are `pending`, `loaded`, `empty` (no rows, no nodes), `failed` (kind and issue) or `evicted` (wanted,
  dropped for the bound), values of a `KeyedState`, so only what reads an entry is told. Each change has a new version
  from one counter; the store's `version` moves with any; an unchanged entry keeps its object. Pending is never empty
  (EC-169).
- **Want sets**: each consumer has one (`wanter()`: `want`, `whenSettled`, `dispose`; the live view's is `want`), the
  whole set it needs now. They merge by best priority, so an export requests and waits without touching the view's wants
  (EC-125). Tiles are deduped, ordered by best priority (view, margin, prefetch) then first want (EC-172), at most
  `maxInFlight` (6) in flight.
- **Leaving the set**: a request in flight is cancelled after `graceMs` (a jitter back keeps it), or at once, oldest
  stale first, when view work waits for a slot (EC-173); a tile not yet asked for or evicted is dropped; an answered one
  stays in the cache; one in backoff or cooldown keeps its timer and its attempts, so going away and back neither resets
  them nor asks a failing source at once.
- **Failures** (EC-174): `failed` is retried `attempts` times with backoff, the entry as it was, then `failed` with a
  cooldown; `malformed` (also a result that is no promise, or an answer that throws when read) goes to the cooldown at
  once; after it a wanted tile is asked again, the entry failed meanwhile (no flicker). `not-served` and `too-large` are
  never retried. Nothing stays in flight after an answer of any kind.
- **Reports**: each cause once per source, at most `maxIssues`, then `seq.store.limit_exceeded`.
- **Eviction**: the loaded bytes are bounded (`maxBytes`, 64 MB); least recently used first among what is not wanted;
  still over, the wanted prefetch then margin tiles become `evicted` and are not asked again until they leave the set
  and come back, so a planner that wants again on every change does not loop. Never the view, never a request on its
  way. Still over: `overBudget`, reported once each time (`seq.store.over_budget`); a coarser level is the planner's.
- **Settle**: no wanted entry pending (a failed one waiting to be asked again is settled), kept on every transition;
  `whenSettled(timeout?)` resolves `false` on a timeout.
- **Nodes** are held, not cached: by each data entry that carried them, by their source's coverage while attached, and
  by `hold`; with no holder they go quietly. A re-answer takes in its nodes before releasing the old ones. The same key
  with other content keeps the held node (`seq.store.node_conflict`). `nodes` lists the entries held now and
  `nodeVersion` moves only when one comes or goes: the overlay's version.
- **Shared per holder**: what consumers share on one store (a source, a nested source, a node, an overlay source) is
  held per holder and goes only when the last lets go; each is tested with two consumers interleaved, one ending while
  the other reads. Why: the viewer and an export share a store, and one consumer's end must not take what the other
  still shows.
- Time is a `StoreClock` (real timers by default): tests drive a fake one and never wait.

## Resolve

`src/resolve/` (`@molstar/seqx-core/resolve`): the static model, tier 1: pure, over a checked spec, nothing from the
store. `resolveSpec(spec, { runtimes, config })`, or `createResolver` to keep entries and paths between specs.

### Configuration

`ResolveConfig`, each choice with its recommended default (`DefaultResolveConfig`): `maxPathSteps` (4),
`maxPathCandidates` (16), `encoding` (`cheapest`; `runs` or `expanded` hold every bounded map so, for tests and
comparisons), `openQueryLength` (int32's max: what the query forms read on a system of unknown length, so labels resolve
there; layout never reads it), `maxNestedDepth` (4: the levels of nested sources the overlay requests, EC-179).

### Systems, sources, annotations

- **Systems**: each `sequence`, with its alphabet, cyclic flag, the built-in numberings (`index`, `ordinal`), its `name`
  (the record its source knows it by: the param, else its id), and its length when its letters are inline; else the
  length is unknown, pending until the source is opened, never a guess. Why: nothing in tier 1 waits, and nothing lays
  out on a made-up length.
- **Sources**: each `source` node bound to a `SourceEntry` made from its params (`SourceParams`, the description's: one
  definition) and to the first runtime of the caller's registry whose `formats` holds its `format`. No format is
  `seq.resolve.source.format`, no runtime `seq.resolve.source.runtime`; either way the source is unavailable.
- **Access plans**: inline letters, rows and values are a `seq` source of one tile over the whole system (the `seq`
  runtime, whatever the registry), named by the system's address; a node with a `source` reads it by region, its system
  by the system's `name`.
- **Wrapped rows**: a row that wraps a cyclic origin is two rows of one item, in the items and the table (EC-9). On a
  sequence of unknown length it is left out with `seq.resolve.length.unknown` until the length is known.
- **The query forms' model**: `system`, `annotation` (ranges; inline items at once, a source's pending over what is
  asked); `node` and `chain` answer none (no overlay, no structures in tier 1), and the form reports
  `seq.query.node.missing`.
- **Caching**: an entry is keyed by its address and the hashes of the nodes it reads (a map: itself and its two
  sequences); an equal spec gives the same model, its paths and resolved views with it.

### Maps

- **One interface** (`Linear`): the image of `[start, end)` forward and back as runs (`from`, `to`, `length`, `dir`:
  `from + i ↦ to + dir·i`), the target ranges into a caller's buffer (`targets`, no allocation beyond its growth), the
  whole map as canonical runs when bounded, the inverse. At most one image per source position; backward may be
  one-to-many (a fusion) and is not ambiguous.
- **Every position has a direction**, a lone one too: canonical runs join only where a run continues the last in the
  same direction, and the expanded encoding keeps a direction per position. Why: a boundary lands after or before a
  position by its direction, so with directions kept every encoding answers alike, and `linearsEqual` (directions
  included) agrees with the boundaries.
- **Encodings**: `implicit` (`k ↦ d + dir·k` on `[lo, hi)`, modulo a cyclic target's length), `runs`, `expanded` (a
  position each, for a region asked). Every encoding of a map gives the same answers; `runs()` is computed once.
  Composition composes encodings: implicit with implicit is implicit, else runs; nothing is expanded.
- **Offset**: implicit; positions whose image falls off a non-cyclic target map nowhere (EC-21); onto a cyclic target
  the image wraps, split at the origin (EC-11).
- **Codon**: the coding sequence's runs, fine place `c` (residue `c / 3`) to a nucleotide; plus strand in genomic order,
  minus strand from the highest exon's end downward (EC-15); one piece is implicit. On a cyclic gene the last exon may
  cross the origin and is split there (EC-11); on one of unknown length it is left out with
  `seq.resolve.length.unknown`. A coding length that is not whole codons, or not 3 per residue with or without a stop,
  is `seq.resolve.map.codon_length` (EC-17), judged only when no piece is left out.
- **Genetic codes**: a codon map's `table` (NCBI id, 1 when absent); built: 1, 2 and 11 (`GeneticCodes`); another is
  `seq.resolve.map.table` and the map has no code. `codonStatus` reads letters through the code: T and U alike, an IUPAC
  letter stands for its set of bases and the codon pairs with every amino acid of the codons it stands for, an unknown
  residue (X) pairs with any codon, a first codon pairs with M where the code starts there, TGA pairs with U and TAG
  with O, other stops pair with nothing; a minus-strand codon is read on the complement (EC-16).
- A map that pairs nothing is `seq.resolve.map.empty`; a partial one is no issue.

### Paths

- **A path** is one `Linear` between two systems and a `step`: `1` (residues to nucleotides), `-1`, or `0`. Across a
  codon step the residue side is in fine places (`3k + j`), every map on it scaled by 3. Why: one linear map answers
  every question of the path; a protein through a chromosome is as sparse as its coding sequence (EC-27).
- **Search** (`findPath`): by counted length from 0, at most `maxPathSteps` counted steps, unit steps free, at most one
  codon step (EC-19), each system once. Candidates that pair nothing are dropped; equal compositions count once;
  distinct ones are `seq.resolve.path.ambiguous` (EC-23, EC-24); more than `maxPathCandidates` is
  `seq.resolve.path.limit`; none within the bound is `seq.resolve.path.none`. A path issue's `address` is the system it
  starts from and its `places` are both systems. The copy rule is in the search: after climbing out of a unit to its
  molecule a path never steps down to a unit of it, so a longer legal path wins (EC-25); a path ending at a ranged unit
  is clipped by its step (EC-26).
- **Edges**: `edgesFrom(system)` lists the spec's maps from a system, both ways: what a search reads, so the overlay
  searches over them and its own.
- **Kept**: a search is done once per model (the reverse read as the inverse), and carried to the next resolve with the
  map entries of every system it read; it is kept while those are the same objects (entries are kept by hash), so an
  edit elsewhere keeps the path object. Why: a map added or changed where the search looked can make it shorter or
  ambiguous; nothing else can.
- **Projection** (`projectRanges(path, ranges, out?)`, `StaticModel.project`): ranges along the path, into the caller's
  buffer when given; from nucleotides to residues, the residues they touch (EC-18). A part with no path is left out with
  its issue; rects are not projected yet.
- **A point** (`projectPoint`): the continuous image of a fractional position: where element `floor(x)` starts on the
  target plus the fraction times the scale in its run's direction. `pathDir` is the path's direction (its first run's).
- **Boundaries** (`projectBoundary`, EC-28, EC-29): the boundary before element `k` lands before its first image, else
  after element `k − 1`'s, else nowhere; onto residues, one inside a codon lands after its residue.
- **Status** (`pathStatus`): paired when the element's whole image exists (a residue's three nucleotides), else
  unpaired; `conflict` comes from letters (`codonStatus`).
- **Closure** (`close`, EC-30): one step through each codon map among the systems named, along the map itself (not a
  shortest path, which two isoforms of one gene make ambiguous): nucleotides gain their residues and lift to their whole
  codons.

### Issues

`seq.resolve.*` (`ResolveIssues`): `ref.missing` and `source.format` (a spec not checked), `source.runtime`,
`map.empty`, `map.codon_length`, `map.table`, `length.unknown`, `path.none`, `path.ambiguous`, `path.limit`, and the
overlay's `node.invalid`, `node.ref`, `node.depth`. A test checks the catalog against the sources.

## Overlay

`resolve/overlay.ts` (`@molstar/seqx-core/resolve`): tier 2, the store's dynamic nodes over the static model, never part
of it (`buildOverlay(model, reader)`; `OverlayBuilder` keeps one a node version and model; `EmptyOverlay`).

- **A node** is `(source, key)`, addressed `<source>/(node :key <key>)` (`nodeAddress`); its depth is 0 under a spec
  source and one more under each source node between (`root`: the spec source). A node whose source is no longer held
  waits.
- **Checked by the description**: each node is validated and normalized as a node of its kind under `data` (cached per
  node object); an invalid one is left out (`seq.resolve.node.invalid`). Refs are absolute addresses (the description's
  `ref` spelling): a `target` names a spec system or an item system the overlay holds, a `source` a spec source or a
  held source node; anything else leaves the node out (`seq.resolve.node.ref`). Each issue once per overlay.
- **Annotations** resolve as the spec's (`resolveAnnotationNode`): inline rows are a `seq` source named by the node's
  address; a `source` reads by region on the target's name.
- **Item systems** (EC-177): an item of a ranges node with `map: 'codon'` (its rows by `id`, else by row; genomic order,
  the first row's strand) gets `<node>/(item :id <id>)` (`itemAddress`): protein, `floor(coding / 3)` residues (a stop
  counts as one), named by its id, with a codon map onto the target (`resolveMapNode`, so its issues are the static
  ones). Made when something names it or a path is first asked for, never before.
- **Paths**: between spec systems the static model's own; with a dynamic end, a search over `edgesFrom` and the
  overlay's item maps with the same bounds, kept per overlay. Why per overlay: an arrival may add a shorter path; one
  search per version, never on hover.
- **Nested sources** (EC-179): `nested(needs)` lists the source nodes at a depth some template on their root needs,
  bound to the first runtime that reads their format (`OverlayOptions.runtimes`; none: `seq.resolve.source.runtime`);
  one deeper than `maxNestedDepth` is `seq.resolve.node.depth` and not listed. Why on need: a nested source is a
  request, and a view that shows no domains asks for none.
- **Named systems**: `systemNamed(name, near)` is the system a nested source's coverage names: an item system of a node
  its parent source (`near`) brought, else a spec system of the name.
- **Holds**: `holds(addresses)` is the nodes of the given addresses (a node, or an item system: its node) and,
  transitively, the nodes their `target` and `source` name: what is shown holds what it is drawn through.

## Resolved view

`src/view/` (`@molstar/seqx-core/view`): the public, renderer-independent form (`resolveView(model, address)`,
`resolveViews(model)`), built over the static model alone, kept per model and view. No data, no dynamic node, no pixels,
no target.

- **A view**: its address, its `layout`, its sections in order, its links, its issues.
- **A section**: an axis, its tracks, and its follow. A section whose `of` names no system is `seq.view.section.system`;
  its layers are not placed.
- **An axis**: its system, its length unless unknown, cyclic, and the positions shown as `[start, end, …]` pieces in
  display order: the whole system, its `range`, or for a cyclic range that wraps the origin two pieces, from its start
  to the origin then on (`[8, 10, 0, 2]`). Without a known length, showing all of the system or wrapping, `ranges` is
  absent: the layout waits for the length. Why: one shape for a wrapped range everywhere, the order it is drawn in.
- **A follow**: the leader's address and axis, the path from the leader's system to the follower's, the scale of a
  position (3 across a codon step from residues, 1/3 back, else 1), its direction (−1 on a minus-strand gene), and
  `point(x)`, the continuous image of a leader position (`projectPoint`), so a follower moves smoothly. No path is
  `seq.view.follow.unplaced` and no follow. `followRange(section, pieces, { arc })` carries a leader's pieces (the
  axis's shape; a piece over a cyclic origin, `start > end`, is split first) to the follower's: the image's hull, or on
  a cyclic follower the shortest arc holding it (EC-12, `arc: 'shortest'`, the default; `'hull'` takes the hull), in the
  axis's shape.
- **Links**: every ordered pair of the view's sections a static path joins, with the path, composed when the view is
  resolved so hover and selection never compose; the same system is the identity.
- **Projection** (`projectOntoView(model, view, loc, { close?, own?, linked? })`, `project.ts`): a Location on the
  view's sections' systems, one Location keyed by system (two sections on one system show the same): each part carried
  along the view's link from a section on its system, else, for a system no section shows, along the model's static
  path; rects are not carried. `close` closes it first among the view's systems (a hover's: a nucleotide lifts to its
  codon and gains its residue); `own` keeps the parts already on a section's system (a hover's and a highlight's), else
  only what it adds where it has no part (a selection's); `linked: false` keeps the own parts alone. Why keyed by
  system, not a function per section: the renderer reads a part per axis with an indexed loop, nothing called or made
  per draw.
- **A track**: its address, its label (its id when absent), its height, and its `template` when it has `for_each`. A
  hidden track (`hidden: true`) is left out of the view, so its row collapses; the spec keeps it.
- **A layer**: its address (with an id) and index, its type and fallback chain (`LayerTypes`; none has a fallback yet),
  its other params as written, its theme channels as values (`color`), its data and the static path from the data's
  system to the section's (none when the same). Data is the section's letters for a `letters` layer without data, an
  address (`kind: 'address'`), or a query kept as written (`kind: 'query'`) on the system its first node is on, with the
  shape of that node (none for a query on a sequence). A type this version does not draw is `seq.view.layer.type`; data
  naming nothing is `seq.view.layer.data`; data of a shape the type does not draw is `seq.view.layer.shape`; data no
  static path brings to the section is `seq.view.layer.unplaced`. Each leaves the layer not placed.
- **Templates** (`templates.ts`): a template is `{ source, kind, depth, order, filter }` (`depth` 0, `order` `key` when
  absent); its layers without data are not placed (no issue). `templateNeeds(view)` is the depth each template source
  needs. `templateNodes` selects the overlay's nodes of the kind at the depth under the source whose params equal the
  filter's values (as written), ordered by key or a param (numbers by value, else as strings; `-` descending), then by
  address. `expandTemplate(model, overlay, section, track)` gives a track a node, `<track>/(node :at <node>)`, labeled
  by the template's label and its top ancestor's key, its data-less layers drawing the node's annotation with its access
  plan and length, carried to the section along `overlay.path` (not placed without one); kept per overlay. Why at
  layout: the resolved view stays pure, and an arrival only makes a new overlay.
- A layer's `data` naming a ref (a chain of them) to one address is read as that address (Language, Refs).
- A layer's `data` may carry `access` and `length` (an instance's, from the overlay): what reads it takes them before
  the model's.
- `seq.view.*` (`ViewIssues`): `section.system`, `layer.type`, `layer.data`, `layer.shape`, `layer.unplaced`,
  `follow.unplaced`; a test checks the catalog.

## Render

`src/render/` (`@molstar/seqx-core/render`): tier 3, the resolved view laid out for a size, a viewport, a target and the
store's entries, as draw lists; hit-testing on the same layout; the recording and SVG backends; `renderSvg`.

### Configuration

`RenderConfig`, each choice with its recommended default (`DefaultRenderConfig`, `renderConfig(over)`): the inset (8),
one label column (112, gap 8), track and section gaps (4, 8), a track's height by its layers' types when the spec gives
none (letters 20, ranges 36, values 32), the lane (12, gap 2), the density thresholds, the letter sizes (9 to 14) and
the mono advance (0.6 em), the label and item fonts (`12px sans`, `10px sans`), the fade (12), the squircle's ratio, cap
and minimum (1/3, 8, 6), outline widths (hover 1.5, selection 1.5, cursor 2.5), the emphasis's look (`emphasis`: the
selection's tint `over` the data, or `under` it; its opacity, 0.18; a highlight channel's tint opacity, 0.3; what is
projected onto a linked section against its own look, `projected`, 0.5), the want margin (half the view a side), a
follower's arc (`shortest`), the tracks a template makes at most (`maxTemplateTracks`, 200), whether a layout lists the
nodes it shows to be held (`nodes.holdShown`, true), what ranges past the lanes do (`ranges.overflow`: `extra-lanes`),
how values aggregate below a pixel (`values.aggregate`: `minmax`), how coarse a level a values layer reads
(`values.coarser`, 1), whether values on the section's own system are placed directly (`values.direct`, true), the
strings draw lists keep interned across frames (`internedStrings`, 4,096; past it a full redraw starts the table again).

- **A values layer's level** (`values.coarser`, k): the coarsest of the source's levels whose bin is no wider than k
  device pixels' positions; 1 is a pixel's (what other layers read). A bin can still be as narrow as k/F pixels, F the
  factor between the source's levels (the store's summaries' or the source's own), so a device column may read up to F/k
  bins or positions: k = 2 halves that worst case, it does not remove it, and a frame's cost still follows the gap
  between levels. What makes the finest level cheap is placing directly (below): 100,000 positions in view, read one by
  one at k = 1, lay out in about 2 ms: so 1 is the default, and it keeps the picture exact (a bin no wider than a device
  pixel); the factor stays a setting.
- **Placed directly** (`values.direct`): data on the section's own system (no path) goes from its positions to x with no
  mapping and no arrays a position; the lists are those of the general mapping (a test compares them), `false` forces
  the mapping.

- **Density steps** (`densityStep(cssPx, dpr, config)`), shared by every layer type, each threshold in its unit
  (`densityUnits`): `letters` from 10 CSS px, `cells` from 3 device px, `dots` from 1, `bars` from 1/64, else `density`.
  The letters step never starts below the width the smallest letter size needs (its advance and 2 px), so it never shows
  cells without letters. Letters, cells and dots draw a mark a position; bars and density at most one a device column a
  lane. Why two units: legibility is in CSS pixels; what bounds the lists is a mark a device pixel.

### Tokens

- `SeqTokens`: every color the lists use (`msui-*` chrome, `seq-*` data: tracks, ink, missing, pending, conflict,
  residue and nucleotide classes, eight categories, the density ramp), as `0xRRGGBBAA`, and the font stacks of the two
  family roles. Node has a light and a dark table (`LightTokens`, `DarkTokens`); the viewer reads its own from CSS.
- **Its version is a hash of its content**, part of the draw lists' key: a theme switch redraws, an equal table does
  not. Why a hash: pure by contract, so an export from Node and a viewer agree on equal tables.
- **A spec's color** is a hex literal or `token:<name>` (`seq-<name>`, then `msui-<name>`, then `<name>`); anything else
  is `seq.render.color` and the layer's default. **Ink follows its fill**: the dark or the light ink, whichever
  contrasts more.

### Text

- The port is `measure(font, text)` in CSS px with one `TextMeasure` (`@molstar/core/util/text-measure`) a font behind
  it (`createTextMetrics`). A font is shorthand naming a family role (`12px sans`, `10px mono`); a backend writes the
  role's stack. Node measures with the bundled table (`bundledTextMetrics`): Arial-compatible sans advances, a 0.6 em
  mono. Letters are never measured: centered in their cells, one size a row, the largest of the scale that fits.
- The port has a `version`, bumped by `clear` (a web font arrived): part of the draw lists' key.
- Font sizes are read from the config once a layout, never per label per draw.

### Draw lists

- CSS pixels, origin top left; a list per pass, drawn in order: `background`, `data`, `highlight`, `selection`, `hover`,
  `cursor`, `labels` (`Passes`; the last four but labels are `EmphasisPasses`).
- **Primitives** as structs of arrays in growable typed arrays kept frame to frame (`reset` zeroes counts): `rect` (x,
  y, w, h, fill, stroke and width, a squircle radius, 0 square, and which ends are round), `path` (move, line, cubic,
  close in shared arrays; fill, stroke, dash), `glyphs` (a run of letters in one font, an x a letter at its center, one
  baseline, a color a letter or one), `text` (an interned string, x, y, anchor, font, color, its measured width and a
  max width: faded past it, never ellipsized), `image` (RGBA bytes the layer owns), `clip` and `unclip`. A pass's `ops`
  keep the order its primitives were added in. Colors are `0xRRGGBBAA`; `NoColor` (0) is none.
- **Bounded by the screen**: at most one primitive a device column a lane, or one a cell when cells are a device pixel
  or wider; tracks below the page are not laid out.
- **Nothing allocated once warm on the draw and emphasis paths** (a full draw of a laid-out view, a hover, a selection
  or cursor redraw; held by the allocation tier, `frame.quick.ts`): one `PassWriter` a pass, kept by the lists
  (`writer`); the interned strings kept across frames (`internedStrings`); the density strip's bytes kept while the size
  holds. Why the rest: in compiled JavaScript a number passed to (or returned from) a call that is not inlined is boxed,
  an allocation, and which calls are inlined depends on the order code ran in. So each pass of a layer type and of the
  view is its own function, called through one table site (`passDraws`, `drawPass`: no call target known there, so each
  is compiled on its own with the writer's calls inlined into it); a loop calls the writer at one site a primitive kind;
  the emphasis passes' rects go through the writer's `args` (an array's numbers are never boxed) and snap inline;
  nothing in a draw makes a closure, an iterator, a destructured array or a number for a call; the hit's numbers are
  read where they are used.
- **A layout allocates a bounded amount a frame**, not nothing: what it makes once (the wants, the regions, the key's
  text, the tiles' keys) and a little an item in view (ranges: an item's key, the lane maps). On the section's own
  system nothing a position or a bin (a values layer places through a kept placer, `AxisTransform.spanIn` and inline
  stores; marks, merges and sorts in kept typed scratch); through a path, each position or bin placed carries its range
  along it (`toSection`: an array a placement), bounded by the screen and the gap between levels. The bounds a scenario
  are the allocation tier's.
- **Reused buffers are tested against fresh ones**: a layer's layout of an input into buffers warmed by another, smaller
  input gives the same lists as a layout into fresh buffers. Why: purity is by contract, not by new objects, and a
  buffer that grows without keeping its contents passes every test whose input fits the first size.

### Emphasis

`Emphasis` (`layer-type.ts`), drawn by the view (`tracks.ts`) in the emphasis passes, per track of a laid-out section:

- `highlights`: the channels in order, each a Location and a color, tinted at `emphasis.highlightTint` on every system
  its Location has (the `highlight` pass).
- `selection` (a Location): its part on the section's system tinted with the accent at `emphasis.selectionTint` and
  outlined; the tints first, then the outlines, so no tint covers an outline (the `selection` pass). Where it has no
  part, `projected.selection` is drawn in the same look at `emphasis.projected` times the opacity.
- `hover`: each layer's own (its cell, its item's marks); `projected.hover` is drawn as the hover's outline at
  `emphasis.projected` on every section but the hover's own (`SectionLayout.id`, the section's id a hit names), so a
  hover from elsewhere (a structure, no layer) is drawn on every section (the `hover` pass).
- `cursor`: the outline of its cell on every section of its system (the `cursor` pass).
- **The order of emphasis** is fixed, below to above: data, highlight channels, selection, hover, cursor; a higher one
  never hides a lower one entirely (tints are translucent, outlines are lines).
- **The tints' place** (`emphasis.selectionPlace`): `over` (the default: the highlights' tints over the data in the
  `highlight` pass, the selection's above them in the `selection` pass with its outline; a selection or a highlight
  redraws the emphasis alone) or `under` (the highlights' tints and then the selection's in the `background` pass, under
  the data and hidden by opaque cells, the outline still in the `selection` pass; an emphasis redraw whose selection,
  projected selection or highlights changed by identity draws `background` and `data` again, `Frame.base`, and a backend
  redraws its base: about a full draw of the data, a millisecond or two on the p53 views, against the emphasis alone).
  Either way the highlights are below the selection. Why `over` by default: a tint under the data is hidden by opaque
  cells, and a selection then costs a data redraw.
- What a shared mark function needs besides its objects (a fill, a stroke, a width) goes through a typed array, and a
  color at an opacity is made inline: the emphasis passes allocate nothing once warm (`frame.quick.ts`: selection and
  cursor; projected hover, selection and highlights).

### Shapes

- Continuous data (cells, value bars) is square; discrete marks (ranges, at every step, a merged run at `bars` too) are
  squircles, radius a third of the height, capped, square below the minimum, so a mark has the shape its emphasis
  outlines. A merged run is round only where its first and last items really end. An end cut by the view, a piece of the
  axis or a piece of the item's own image (an exon) is square; a mark cut at both ends is square. The squircle is one
  curve, `squirclePath` (each corner a cubic on |x|⁴ + |y|⁴ = 1), that every backend draws. Outlines are concentric.

### Layer types

`LayerType<Prep, Lay>`: `name`, `fallback?`, `prepare(layer, model)` (tier 1, cached per resolved layer),
`layout(prep, ctx, data, overlay, prev?)` (per viewport; `prev` is its last layout, to reuse buffers and keep lanes),
`draw(lay, pass, out, emphasis?)`, `hit(lay, x, y)`, `legend?(lay)`. `LayoutContext`: the track's box, the axis
transform, the section, the layer's key, the device pixel ratio, the config, the tokens, the text port, whether the view
settled, `want` and `report`. A layer draws the first type of its fallback chain the target has; none is
`seq.render.layer.type`.

- `letters`: a section's own letters (through a path: `seq.render.layer.letters`); a cell a position colored by class,
  letters on them, a tick a position at `dots`, one bar a piece in view at `bars` and `density` (presence needs no
  data).
- `ranges`: items in sticky lanes, carried onto the section along the layer's path (a domain on its gene: a piece an
  exon); labels inside marks 24 px or wider, drawn in the `labels` pass. At `bars` a lane's marks touching within a
  device pixel merge into a run; at `density`, or a density table, one heat strip, a device column a pixel. Items past
  the lanes collapse into a density row in the last lane, the rest drawn as items (`ranges.overflow: 'extra-lanes'`), or
  the whole layer becomes the strip (`'layer'`).
- `values`: a bar a position, or at `bars` and `density` a bar a device column over the column's extent from the
  baseline (`values.aggregate: 'minmax'`, so negative extremes show; `'max'` the maximum alone), from positions or bins.
  Missing (NaN) is a 2 px neutral mark at the baseline. The domain is the layer's `domain` param where given, else open:
  computed from the values in view when the view settles, held while it does not.

### Data

- A layer's data is its node's access plan (an inline node's own address as the store's source, else its source node),
  the region in view carried back along the layer's path, at the level that fits a device pixel's positions there
  (`fitLevel`; a values layer's, `values.coarser` pixels'); the view at `view` priority and a margin each side at
  `margin`.
- A tile not asked for, pending or evicted draws as a pending placeholder (EC-169), a failed one as a conflict tint; an
  unknown source (not attached) is pending. Two regions sharing a tile (a cyclic view wrapping a short system) read its
  entry once.
- **An item's key** is its `id`, else its `row` (a null `id` falls back to `row`), else its tile and index; never its
  extent. A key seen in an earlier tile is that range again and is not drawn twice; within one tile a repeated key is
  another piece of one item (a wrapped row).

### Tracks layout

- Sections stacked; each an `AxisTransform`: its pieces (`[start, end, …]` in display order) end to end as display
  units, the viewport a window of them (`{ start, end }`, all by default) across the data box. A span that only touches
  the view's edge is not in view (a boundary on the edge is).
- **Leaders first**: axes are made in the order of `follow`, so a follower placed above its leader follows it; a leader
  not in the view is `seq.render.follow.leader`, a cycle `seq.render.follow.cycle`, and the section shows its own range.
- **A follower follows continuously**: its window is the images of its leader's window edges (`point`), in genomic
  order; where an edge maps nowhere, the leader's view has two pieces or the window wraps a cyclic origin, the hull (or
  arc, `followArc`) of the leader's positions in view (`followRange`).
- A section whose length is unknown takes the length its source's coverage gives once attached, else waits
  (`seq.render.axis.pending`).
- **Templates** are expanded here over the frame's overlay (`expandTemplate`), at most `maxTemplateTracks` instances
  (the first in order; more is `seq.render.template.limit`, so the work is bounded before it starts); a template's
  source region in view, on the section's system as the source names it, is wanted at the level that fits, so the nodes
  it brings arrive. The overlay's issues are the layout's. A template below its source also wants each nested source's
  systems (`systemNamed`) in view, carried along the overlay's path from the section, so a nested source that sends its
  nodes by region is asked for them.
- **Shown nodes** (`nodes.holdShown`, true): the layout lists the nodes it shows (`ViewLayout.shown`: the data and the
  path systems of each laid-out layer, through `overlay.holds`), for its owner to hold in the store; with `false`, none.
  Why: a node stays while it is shown, not only while the region that brought it is cached.
- Cell edges snap to the device grid (`snap`). One label column for every track of a view.

### Lanes and Hits

- **Sticky lanes** (`StickyLanes`): an item keeps its lane while it is in view; a new item takes the first lane free
  over its extent; repacked only when the view settled and a fresh packing needs fewer lanes. Why: items must not jump
  when a neighboring region arrives.
- **Hits** are `{ layer, section, element?, extent?, position?, system }`: `layer` the layer's key (its address, else
  its track's and `:index`), `section` the section's id, `element` an item's `{ owner, key }` with its `extent` (a
  Location on the data's system), `position` on the section's system. Cells by arithmetic on the axis; items by
  `LaneIndex` (a lane's marks by left edge with a running maximum of right edges). In a track, the topmost layer with an
  element under the point wins, else the topmost with a position.

### The key

- `Renderer.frame` keys a frame by the model's hash, the view, the size and device pixel ratio, the viewport, the token
  table's version, the text metrics' version, the overlay's version, whether the view settled, the config, and the
  versions of the entries (and the coverages) the last layout read. An equal key keeps the layout and the lists and
  redraws the emphasis passes alone; anything else lays out again into the same buffers. Same inputs, same lists.
- The check compares the parts as they are (and the reads as arrays), building nothing: a hover allocates nothing. The
  key's text (`Frame.key`, `DrawLists.key`) is made when the lists are rebuilt; an emphasis redraw returns one kept
  `Frame` record, its fields set again. So a frame is valid until the next `frame` call (the lists are reused too): a
  caller copies what it keeps.

### The overlay's sources

- `OverlaySources(store, assets).sync(overlay, view)`: a dynamic annotation's inline rows are attached as a `seq` source
  under the node's address, and the nested sources the view's templates need (`overlay.nested`) through their runtimes;
  each kept while its node is the same store entry, detached when it goes (or changes), so a node's data lives as long
  as the node. One that cannot be opened is an issue (`issues`) while its node is held, opened again at the next sync,
  cleared when it opens or its node goes; `renderSvg` reports those left at its end. `hold(nodes)` holds in the store
  the nodes a layout shows and lets go of those it no longer shows (`dispose` lets go of all). On a shared store the
  overlay sources (a viewer's, an export's reading through it) share what one of them attached, never replace it: each
  claims what it needs, and the last claim let go detaches it, so an export ending never takes away a source the live
  view reads (and cancels its reads). One attached by another consumer (not an overlay's) is left alone.

### Backends

- **Recording** (`recordDrawLists`): a line a primitive, pass by pass, numbers to two decimals, colors `#rrggbbaa`
  (`none` for no color): the committed records' format.
- **SVG** (`drawListsToSvg`): a string, no DOM; passes as groups, clips as clip paths, squircles as the core's path,
  images as runs of rects, fades as gradient masks, metadata as `data-*` attributes. Every text and attribute from a
  spec or data is escaped; a character XML forbids (a control but tab, newline and return, U+FFFE and U+FFFF, a lone
  surrogate) becomes U+FFFD, one for one, so a glyph run stays aligned with its x list. Ids take a prefix (`idPrefix`,
  name characters only), a metadata key is written only as a plain name. Held by hostile labels, ids, titles, fonts and
  metadata parsed back as XML (`_test/svg-escaping.test.ts`).

### Export

- `renderSvg(spec, { size, viewport?, runtimes, … })`: checks and resolves the spec, opens the sources its view reads
  (`attachSources`: inline data through `seq`, a source node through its runtime, a template's source; assets refused
  unless a reader is given, `NoAssets`), wants through a want set of its own (EC-125), each round reads the overlay and
  syncs its sources (let go at the end), lays out until its wants and the store's node version hold, waits for them to
  settle, and draws once more with the view settled. It waits for its view's wants alone, never their margins.
- **One deadline** (`timeoutMs`, 10 s, on `clock`: real timers by default, a shared store's clock or a fake one) covers
  opening the sources and settling; past it the opens are aborted and the export is `seq.render.unsettled`.
- **A selection as an address** (`selection`): resolved by the export itself once its data settled, over its own store's
  nodes (`overlayQueryModel`), and drawn with `emphasis.selection`; what it named is `selection` in the result (absent
  with its issues when it named nothing). Why: an export is a consumer of an address of its own, not a copy of a
  viewer's Location.
- **The selection drawn is projected** onto the linked sections by the export (`projectOntoView`), unless the caller
  gives `emphasis.projected.selection`: the export shows on the gene what a viewer shows there.
- It answers `{ svg, settled, issues, lists }`; the SVG records whether it settled and the open domains it computed, and
  its ids take a prefix (`idPrefix`, by default `s` and a hash of the spec). On a shared store it reports only the
  issues of the sources it attached.

### Issues

- `seq.render.*` (`RenderIssues`): `view`, `color`, `layer.type`, `layer.letters`, `axis.pending` (info),
  `follow.leader`, `follow.cycle`, `source`, `template.limit`, `unsettled`; a test checks the catalog.

## Context

`@molstar/seqx-core/context`: the headless composition. One place lists what a render needs; the viewer composes the
same pieces. Nothing registers on import and there is no global registry.

- `createSeqContext({ registry?, config?, ports? })`: the registry's runtimes, the configuration and the ports, and
  through them `check`, `resolve`, `view` (the resolve passes) and `renderSvg`. `retire()` disposes what the context
  made or was told to dispose.
- **Registry entries** are data: `{ kind: 'format-runtime', runtime }` in this version. The same runtime object again is
  no change; another under one name is an error listing every conflict (nothing registered). Layer types, themes and the
  other kinds join when `renderSvg` takes them (a `types` option); run-time, atomic and undoable registration is the
  viewer's followup. Why only runtimes: a layer type the SVG export cannot draw would make the canvas and the export
  differ.
- **Configuration** is `[item, value]` pairs (`SeqConfigItem`, `configItem(key)`); the core's items are `SeqCoreConfig`
  (`theme`, `render`, `resolve`). Items of other keys are another package's (the viewer's `group.item` keys).
- **Ports** are factories with ownership, each called once; unset, the default: no reads (`NoAssets`) and the bundled
  text metrics. A result's `disposeWithViewer` is true unless `false`: a host that passes its own keeps it.
  - **assets**: either the application's `AssetManager` behind `policyAssets` (`ManagedAssetsPort`, with a policy), or a
    reader of its own (`ReaderAssetsPort`: a test adapter that refuses the network, the Node fixtures' adapter). Why
    two: a suite per port needs a second adapter, and an `AssetManager` is not what a test or a CLI has.
  - **text metrics**: `measure(font, text)`; the context wraps it in the per-font cache (`createTextMetrics`).
  - **structure**: none in this version (`runtime` unset).
- **The asset policy** (`AssetPolicy`): a decision on the location, made on the URL it resolves to, never on the string
  it came in as (the parser drops spaces, tabs and newlines, and reads `//host`, `\\host` and `/\host` as another host).
  A location with a scheme is absolute: its scheme must be allowed (none by default) and, when hosts are listed, its
  host. One without is relative: it must resolve to the base's own scheme and host and, with `confine`, stay under the
  base's folder (`..`, `%2e%2e` and a rooted path resolve out of it; an encoded separator is refused). `maxBytes` counts
  bytes (a text's UTF-8 length) and `timeoutMs` bounds the read. `BrowserAssetPolicy` is `https:`, relative under the
  page, 256 MB, 30 s; `HeadlessAssetPolicy` reads nothing. A refusal rejects with `AssetRefused`, which a source's
  opening turns into an issue on the source (`seq.source.unreadable`). Why: a spec's `uri` is untrusted input.
  - **What a read passes to the manager** (`AssetManager.ResolveOptions`, core's `util/_spec/assets.md` §8; a manager
    honors each or refuses the read): `credentials` (`'omit'` unless the policy's `credentials` says otherwise: a
    composition that needs a login says so), `redirect: 'error'` (the location checked is the location read; the
    policy's `redirect` says `'follow'` for a server that redirects. A browser does not say why `fetch` failed, so a
    refused redirect there is the generic download failure; Node.js says, and the read is `AssetRefused`), a signal (the
    caller's and the timeout's: `timeoutMs` cancels the download), and `maxBytes` counted while the body streams (a body
    over it is a refused read, `AssetRefused`). The policy also checks the size once the body is read and abandons a
    read the manager did not end in time, for a manager that falls short; core's `checkAssetManager` is the suite an
    app's manager passes. Why: the policy's claims hold for any manager the application composes.
- **The headless policy** (`headlessAssetPolicy({ root, hosts })`): files under `root` (an absolute path, a Windows
  drive path or a `file:` URL; relative locations resolve against it and stay under it) and `https:` locations of the
  listed `hosts`, nothing else: no network when no host is listed, no `file:` location without a root, with a host, with
  an encoded separator, userinfo or outside it (`..` and `%2e%2e` are resolved by the URL parser first). Credentials are
  never sent, redirects are refused, and the size and time bounds are the browser's. `AssetPolicy.root` is the same rule
  for any policy; the root is a decision on the location, a link in the folder that leads out of it is the reader's to
  refuse. Listed hosts are normalized as the URL parser reads them (case, default port); a root that is no location
  leaves the policy with no file access and no throw. Why: a headless render or a server reads whatever `uri` a spec
  names (server-side request forgery).
- **Bounds**: `SeqCoreConfig.bounds` (`Partial<SpecBounds>` over `DefaultSpecBounds`, Language, Bounds) is the context's
  `bounds`, which its `check` and its `renderSvg` hold (`RenderSvgOptions.bounds`).
- A port factory that throws leaves nothing made: the ports already made are disposed, last first; a structure port made
  is disposed with the context.
