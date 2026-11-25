# The language description and the tree core's delta

## Intent

The second review (2026-11-17) found the first slice blocked here first: no design says what the language description
looks like, what the slice's kinds and params are, how the builder reads, or how far `@molstar/core/repr` departs from
meso\*'s tree core. This design settles them for the first slice.

## Settled

Accepted as recommended (the author, 2026-11-17, at the first slice's Gate 0: "approve both"): fork the tree core into
`repr` now and migrate meso\* later; the slice's kinds and param names as below. The rest follows from agreed designs ([format.md](format.md), [tree-core.md](tree-core.md), [composition.md](composition.md)).

## Decision

### The tree core's delta from meso\*

| Item | meso\* | `@molstar/core/repr` | Why |
|---|---|---|---|
| node | `{ kind, id?, params?, children?, args? }` | the same, plus `custom?` | extension data on any node ([extensibility.md](extensibility.md)) |
| values | atoms, lists, records, nodes | the same, plus **typed columns** (`Int8Array` … `Float64Array`, and a dictionary-coded category `{ dict, codes }`) | real data without per-value cost |
| column identity | — | by element type and bytes, hashed once and cached on the node | equality and hashing of large columns |
| NaN | `{ "$number": "NaN" }` | the same outside columns; `null` inside a float column in JSON | columns are the common case |
| issue | `{ code, message, hint?, path, address? }` | the same, plus `severity` | one meaning of severity ([tree-core.md](tree-core.md)) |
| kind names | `^[a-z][a-z0-9_-]*$` | an optional namespace, `^([a-z][a-z0-9_-]*:)?[a-z][a-z0-9_-]*$` | extension kinds (`ext:kind`) |
| params | `ParamType` (string, number, integer, boolean, list, record, map, union, node, any, ref, def) | the same, plus `column` and `table` | declared columns |

- **Codecs stay schema-free; reading types the columns**: a codec decodes a JSON array as a plain array; the read
  pipeline is decode → migrate → validate and normalize, and normalizing turns every value of a `column` param into
  its declared typed array (MessagePack's extension types arrive typed already). Equality, hashing and printing are
  defined on normalized trees. Why: a codec needs no language, and every tree the passes see has typed columns.
- **Fork now, one home later** (recommended): `repr` is written from meso\*'s tree core with the delta above;
  meso\* keeps its own until its repository joins this one, then migrates to `repr` (the delta is additive). The
  alternative, changing meso\* first, couples two repositories' work now.

### The description's shape

```ts
const SeqLanguage = defineLanguage({
  id: 'seq', version: '0.0.0', root: 'spec', doc: '…',
  kinds: {
    spec:     { doc, parents: [], params: { version: { type: 'string', required: true, doc }, title: { type: 'string', doc } } },
    data:     { doc, parents: ['spec'], ids: 'required' },          // the id is `data`
    sequence: { doc, parents: ['data'], ids: 'required', params: {
                  alphabet: { type: 'string', enum: ['protein', 'dna', 'rna'], required: true, doc },
                  letters:  { type: 'string', doc },                 // inline; or `source`
                  source:   { type: 'ref', kinds: ['source'], doc },
                  cyclic:   { type: 'boolean', default: false, doc } } },
    ranges:   { doc, parents: ['data'], ids: 'required', params: {
                  target: { type: 'ref', kinds: ['sequence'], required: true, doc },
                  rows:   { type: 'table', columns: {
                            start: { type: 'column', of: 'int32', required: true, doc },
                            end: { type: 'column', of: 'int32', required: true, doc },
                            id: { type: 'column', of: 'category', doc }, label: { type: 'list', of: { type: 'string' }, doc } },
                          doc } } },
    // values, map (offset, codon), source, views, view, section, track, layer, theme forms: the tables below
  },
});
```

- **Per kind**: `doc`, `parents`, `ids`, `params` (each with `doc`, `required`, `default`, bounds), `fallback` (a
  built-in kind or layer type a reader without this kind uses), and `data` for data forms (queries, themes) whose items
  are positional `args`, as meso\*'s.
- **`column`**: `{ type: 'column', of: 'int8' | 'uint8' | 'int16' | 'uint16' | 'int32' | 'uint32' | 'float32' | 'float64' |
  'category', missing?: boolean }`; integer columns have no missing value; `missing` on one is refused (`repr.schema.language.column`) until a column
  needs a validity list, whose value shape (a `{ values, valid }` pair or a sentinel per type) is then decided.
- **`table`**: `{ type: 'table', columns: Fields, key?: string }`: columns of equal length, checked; a column is required
  by its own `required: true`, as a record's fields are (lane A, drift resolved);
  extra columns allowed and kept (custom fields of a table, [extensibility.md](extensibility.md)).
- **Slots**: a param of `type: 'node'` naming a sublanguage (`query`, `theme`) holds a data form of it, read and printed
  in place ([locations-and-queries.md](locations-and-queries.md)).
- **Extensions ride on the one description**: `SeqLanguage` is one description; an extension's layer type is a namespaced
  word on `layer.type` with its params plain beside the built-in ones, keys `<id>:<name>` ride on any node, an unknown
  namespaced kind is kept with a warning, and each extension checks its own (`SeqExtension.check`) where it is composed,
  as registration is ([viewer-composition.md](viewer-composition.md), [extension-openings.md](extension-openings.md)) (amended 2026-11-23: the code's way, the author as recommended).
  No `composeLanguage`: the later extension design superseded it; it returns only when an extension needs kinds of its
  own with builder types.

### The kinds built (first and second slices)

| Kind | Under | Params (required in bold) |
|---|---|---|
| `spec` | root | **`version`**, `title` |
| `data`, `views` | `spec` | containers with reserved ids |
| `source` | `data` | **`format`**, `uri`, `index`, `value` (inline), `fields` (column map), `assembly` |
| `sequence` | `data` | **`alphabet`**, `letters` or `source`, `cyclic` |
| `ranges` | `data` | **`target`**, **`rows`** (table: **start**, **end**, id, label, category, strand) or `source` |
| `values` | `data` | **`target`**, `values` (float32 column, one per position) or `source` (by region) |
| `map` | `data` | **`from`**, **`to`**, **`method`** (`offset` with `by`; `codon` with `exons` table and `strand`) |
| `view` | `views` | `layout` (`tracks`) |
| `section` | `view` | **`of`** (a system's address), `range`, `follow` (a section's address) |
| `track` | `section` | `label`, `height` (px), `hidden` (boolean, default false: kept in the spec, its row collapsed; the second slice, as `layer.hidden` was decided); `padding` (px inside the well: a number top and bottom, or `(top n bottom n)`; the third slice, an addition to 0.4.0) |
| `layer` | `track` | **`type`** (`letters`, `ranges`, `values`), `data` (an address or a query), `color` (a theme form), `hidden` (boolean, default false; the author, 2026-11-18); its id required (the author, 2026-11-18: every child in a spec is addressable, so hiding a layer is a one-key `set`) |
| `refs`, `snapshots` | `spec` | containers with reserved ids (the second slice) |
| `ref` | `refs` | **`to`** (an address, or a non-empty list of them: a set), `label`, `color` ([addressing](addressing.md)) |
| `snapshot` | `snapshots` | `label`, `patch` (a stored patch in `seq.patch`, no `base`: a story step, [stories](stories.md)) |

**0.5.0 additions from the fourth slice's designs** (one bump, unreleased): `spec.extensions` (`{ <id>: '<version>' }`,
the least version a figure needs) and a layer's `fallback` and `legend` ([extension-openings](extension-openings.md));
an extension's params plain on its own namespaced type and as `<id>:<name>` keys in another node's `params`, never in
`custom`; `source.ref` and `values.field` ([seq-data-ports](seq-data-ports.md)); a map's `priority` and a section's
`via` ([structure-syncing](structure-syncing.md)); the molecules kinds ([molecules](molecules.md)).

The language's version moved with each kind added under `0.x`, with no migration: `0.1.0` every child has an id,
`0.2.0` `refs` and `ref`, `0.3.0` `snapshots` and `snapshot` (the code's `SeqLanguageVersion`).

### The builder, by example

```ts
const b = createSeqSpec({ title: 'p53' });
const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53_LETTERS });
const gene = b.data.sequence('tp53', { alphabet: 'dna', source: b.data.source('chr17', { format: 'fasta', uri, index }) });
b.data.map('p53-on-gene', { from: p53, to: gene, method: 'codon', exons: EXONS, strand: '-' });
const domains = b.data.ranges('domains', { target: p53, rows: { start, end, label } });
b.views.view('main', (v) => {
  const protein = v.section('protein', { of: p53 });
  protein.track('letters').layer('letters', { type: 'letters' });
  protein.track('domains').layer('domains', { type: 'ranges', data: domains, color: theme.field('label') });
  v.section('gene', { of: gene, follow: protein });
});
const spec = b.build(); // a normalized tree; print(spec) for the S-expression view
```

Handles are addresses (`p53` is `/data/p53`), so a builder value can stand wherever an address is.

### Issue codes

- **Namespaces**: `repr.*` (tree, schema, patch, codec), `seq.lang.*` (checks the description cannot express), `seq.query.*`,
  `seq.resolve.*`, `seq.source.*`, `seq.store.*`, `seq.view.*`, `seq.render.*`, `seq.address.*`.
- **Resolve codes of the fourth slice** (decided by the author, 2026-11-20): `seq.resolve.system.unnamed` (a section over
  a source sequence names no sequence after the source's system, so its length is unknown and nothing lays out: an
  error, a spec fault, never data-dependent, its hint saying to name the sequence), and `seq.resolve.ref.system` (a
  reference to a unit or molecule that made no system is left out, never a throw). `seq.resolve.length.unknown` stays
  the transient warning of a source not yet open: the viewer holds it back while the address is pending
  ([viewer-models.md](viewer-models.md)).
- **First-slice choices**: `sequence.letters` are checked against the alphabet (protein `A-Z`, `a-z`, `*`, `-`; DNA the
  IUPAC codes and `-`; RNA the same with `U`); `layer.color` is an `any` param and `layer.data` a `ref` until the theme
  and query sublanguages land, then a theme slot and a union of a ref and a query slot (the forms above are the target).
- **A catalog**, as meso\*'s `IssueCodeInfo`: each code with its severity, doc and hint, generated into the docs; a test
  fails on a code used and not cataloged.
- **Dedupe** (EC-86): issues of one code and cause are reported once with a count; the first ten places are kept.
- **On a dynamic node**: `address` is the node's form (`/data/genes/(node :key T1)`); `path`, where the language's checks
  give one, is a place inside the node (`['params', …]`); the overlay's own resolve issues carry no path; never a spec
  path (amended 2026-11-23: the code's way, the author as recommended).

## Consequences

- The first slice's step 1 builds the delta, the codecs and normalization; step 2 writes `SeqLanguage` with the kinds
  above and the builder, tested by the p53 example built, printed, encoded both ways and validated with every misuse
  giving its code.

## In the specs

| Decision | Where |
|---|---|
| The node gains `custom` | repr › custom |
| Typed columns and categories | repr › Typed columns |
| Column identity by element type and bytes, cached | repr › Columns by element type and bytes |
| NaN: tagged outside columns, `null` in a float column in JSON | repr › JSON |
| The issue gains `severity` | repr › repr: the tree core › Issues |
| Namespaced kind names | repr › The tree › Kinds |
| `column` and `table` params; `missing` refused on integer columns; extra table columns kept | repr › Param types |
| Codecs schema-free; reading types the columns | repr › Codecs stay schema-free |
| Fork now, migrate meso\* later | design-only: a plan (meso\* in another repository) |
| The description's shape per kind (`doc`, `parents`, `ids`, `params`, `fallback`, `data`) | repr › Schema › Kinds |
| Slots: a theme slot | core › The theme slots |
| Slots: a query slot in a layer's `data` | design-only: not built (waits for its sublanguage, core TODO) |
| Extensions ride on one description: namespaced layer types and keys, each extension checking its own (no `composeLanguage`) | core › Language › Extensions |
| The kinds built and their params | core › Kinds |
| The 0.5.0 additions; versions moved with no migration | core › Versions |
| The builder; handles are addresses | core › Language builder |
| Issue namespaces | core › @molstar/seqx-core › Issues |
| `seq.resolve.system.unnamed`, `seq.resolve.ref.system`, `length.unknown` transient | core › Systems, sources, annotations |
| Letters checked against the alphabet | core › What the checks cover |
| `layer.color` and `layer.data` until the sublanguages | core › The theme slots |
| A catalog per namespace, a test on uncataloged codes | repr › Catalogs |
| Dedupe (EC-86): one code reported with a count, the first ten places | core › Language › Bounds |
| On a dynamic node: `address` its form, never a spec path | core › Overlay |
| On a dynamic node: `address` the node's form, `path` a place inside the node (none on the overlay's own issues) | core › Overlay |
