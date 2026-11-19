# @molstar/seqx-core

The SeqViewSpec language and everything that runs on it without a DOM: the language description and builder, the
coordinate model and maps, the source protocol and the data store, queries, the resolve passes, draw lists, and the SVG
and recording backends. It imports `@molstar/core` and an allowlist of other Mol\* entries, each with its reason.

This spec holds the contracts of what is built, each a rule and its reason; code comments cite its sections. A section
is added in the change that builds what it describes.

- **Entries carry what their readers need, no more**: an entry a reader imports (`lang`, `resolve`, `query`) holds no
  module only an editor or a viewer needs; editing's checks are `lang/changes`. Why: the packages declare no
  `"sideEffects": false`, so a re-export keeps its module's top-level code, and reading a spec paid 26 kB for editing
  (the size targets, checked at each gate, guard it).
- **A restriction names its reason**: a rule here that refuses or limits (an issue code, "only", "never") says why.

## Versions

- The language has one version, `LanguageVersion` (`0.4.0`). Why: one language, one description, one number for readers
  to check.
- **Pre-release** (`0.x`): a change to the language moves the version with no migration (`0.1.0`: every child of the
  spec has an id, `layer` included; `0.2.0`: the `refs` container and the `ref` kind; `0.3.0`: the `snapshots` container
  and the `snapshot` kind; `0.4.0`: a section's `direction`, `pieces`, `collapsed`, `pad` and a layer's `complement`,
  `code`, the layer types `ruler`, `secondary`, `variants` and the layer params `numbering`, `encoding`, `style`,
  `line`, `arrow`, and the themes: the layer's channel slots `color`, `ink`, `size`, `opacity`, `line`, `arrow`, `shape`
  taking a value or a `seq.theme` form, its `scale`, the `scale` kind under `views`, a view's `export`, a view's
  `layout` `text` or `rings`, a section's `layout` and its text params `group`, `numbers`, `group_labels`, `align`,
  `max_width`, `height`, its rings params `rotation`, `labels`, its blocks params `arrange`, `copies`, the `block` kind
  under a section, a section's `of` optional when it has blocks, the view's `export.window`, the layer type `frame` and
  its `across`, the kinds `alignment` and `tree`, a section `of` an alignment, the layer types `rows` and `logo` and
  their params `tree`, `order`, `collapsed`, `pinned`, `hidden_rows`, `height`, `summary`, `numbers`, `tree_width`,
  `logo_scale`, `correction`, `gaps`, and a source's transform `records`, `letters`, `system`, `base`, `ends`, `types`,
  `attribution`, `licence`, with `fields` read). A spec in an earlier `0.x` reads as one in this one, and is reported
  where it no longer fits.

## Language

`lang/`: `SeqLanguage`, the description of SeqViewSpec on `@molstar/core/repr` (kinds, parents, params, docs), and the
checks it cannot express.

- **Kinds**: `spec` (root: `version`, `title`); the containers `data`, `views`, `refs` and `snapshots` (reserved ids of
  the same name); under `data`: `source`, `sequence`, `ranges`, `values`, `map`, `alignment` (rows over columns: its
  `rows` inline, a `label` and aligned `letters` a row, every row as wide, or a `source`; its `reference` row by label;
  the columns system its address, the rows system `rowsSystemOf(address)`, `<address>/(rows)`: rows are a system,
  EC-35), `tree` (`of` an alignment, `newick` its text; leaves matched to rows when drawn), `pairs` (two locations that
  belong together: `a` a sequence, `b` another or `a`'s; its `rows` inline, `a_start`, `a_end`, `b_start`, `b_end` each
  pair's two ends, never wrapping, at least a position each (EC-10), self-pairs allowed; an open `type` its kind, since
  `kind` is reserved: canonical, wobble, contact; `id`, `label`, `a_strand` and `b_strand` (EC-206)), `grid` (a value
  per pair of positions: `x` a sequence, `y` another or `x`'s, `values` x-major float32, NaN missing, `|x| × |y|` of
  them: `seq.lang.grid.length`); under `views`: `view`, then `section`, `block` (a piece of a section on a system of its
  own: `of`, `range`, `label`), `track` (`hidden`, default false: kept in the spec, its row collapsed; its layers keep
  their own), `layer` (`hidden`, default false: kept in the spec and not drawn, so hiding is a one-key `set`). A
  `ranges` node's `map: 'codon'` makes each item (its rows by `id`) a coding sequence of its own (Overlay); a track's
  `for_each` makes it a template (Resolved view, Templates; its kind param is `node_kind`, since `kind` is reserved, and
  is `ranges` or `values`: what a layer draws). Reference params (`target`, `of`, `from`, `to`, `follow`, `source`,
  `data`) are `ref` params naming the kinds they take; a reference is an address.
- **The axis's params** (the writer's: what the figure is; never inferred from what a section draws, so adding a layer
  never turns a section over): a section's `direction` (`genomic`, `reversed`, `follow`, or
  `{ along: <sequence or map> }`; absent: `follow` for a section with `follow`, else `genomic`), `pieces` (`whole`,
  `coding`, a ref to a `ranges` node or a codon `map`, or a list of refs, their union; absent: `whole`, followers
  included), `collapsed` (the form `section.collapse` switches to; absent: `coding`) and `pad` (positions kept each side
  of a piece; absent: 0); a layer's `complement` (`auto`, `always`, `never`; absent: `auto`) and `code` (`auto`, `one`,
  `three`; absent: `auto`), read by `letters`; a layer's `numbering` (a numbering id; absent: the section's default,
  `ordinal`), read by `ruler`, `encoding` (`bars`, `line`, `area`; absent: `bars`) by `values`, `style` (`blocks`) by
  `secondary`, `line` (`solid`, `dashed`, `dotted`, `wavy`, `double`, `box`, `highlight`; absent: blocks) and `arrow`
  (`none`, `end`, `start`, `both`; absent: `none`) by `ranges`; a grid section's `of_y` (absent: `of`), `follow_y` (`x`,
  `none` or a section; absent: `x`) and `aspect` (a number or `none`; absent: 1); by `pairs` its `side` (`above`,
  `below`; absent: `above`), `reach` (absent: 0.25), `height` (`span`; absent: `span`), `crossing` (`same`, `split`;
  absent: `same`), `across_lines` (absent: false) and `mirror` (absent: false: a pair drawn where its data puts it,
  never read off whether the two axes are one system; the tRNA's builder sets it); by `grid` its `aggregate` (`mean`,
  `max`, `min`; absent: `mean`); by both in a grid `triangle` (`both`, `upper`, `lower`; absent: `both`). Why: the
  figure's form is the writer's, never the host's or the data's: one spec draws the same in the canvas, the SVG and
  another viewer. The refs inside `direction` and `pieces` are checked by kind as top-level refs are, and a change's
  check reaches them (its reverse index walks a union's record member). Why: the composition is the application's.
- **Navigators** (`views/<view>/navigators/<id>`, the `navigators` container's id is `navigators`): a `navigator` node
  is a part over a section's axis, its rows or both, placed by the host: `of` (a section or a list of them, refs),
  `axes` (`positions`, `rows`, `both`), `content` (`axis`, `tracks`, `bands`, `rows`), `tracks` (refs), `bands` (a
  `ranges` ref with a `stain` field), `system` (a sequence: a chromosome under a section on a region), `place` (`top`,
  `bottom`, `left`, `right`: a hint), `size` (px across, or `[width, height]`), `hidden`. Absent params are the form's:
  the builder's four helpers (`stripNavigator`, `ideogramNavigator`, `overviewNavigator`, `miniatureNavigator`) write
  them. The view's `export.navigators` (`omit`, `include`) says whether an export draws them. `content: bands` with no
  `bands` is `seq.lang.navigator.bands`. Why a node: a spec carries its navigators to another reader and an export can
  draw them; why not a section: it has no axis of its own to scroll, its window is a section's viewport.
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
- **Checks of a change** (`check-changes.ts`, `checkSeqChanges`, `SeqPatchValidator`; their own entry,
  `@molstar/seqx-core/lang/changes`, so reading a spec carries no editing and its persistent maps): the `seq.lang.*`
  checks as `repr/patch`'s validator, on what a patch's `Changes` reach: the nodes it made, set or moved; what names
  them, or named any node of a subtree it removed or moved away (walked in the old root); on through the user refs that
  name those, and through a ref with segments after it (`/refs/p/x` where `/refs/p` names a prefix of the node, at most
  the ref depth deep); a template's layers under a track whose params changed; the sections and refs of a chain a
  changed one closes. What names what is a reverse index of the spec's references (a persistent set of referring nodes
  by the id run they name), built once a spec and derived from the last revision's along its changes. The same checks on
  the same nodes before the change are taken away, by code, node (through moves) and place in the node. Why: a patch is
  refused only for what it introduced (an error already there never blocks an edit beside it), and the cost follows the
  change, not the spec. Held by tests that compare it with the whole check's difference, patch by patch, along a chain
  of revisions and over seeded random patches of ref chains; and in the heavy lane, at 16,000 layers, the first check
  after a load and a fan-out to every node no worse than the whole check, a one-node edit a small part of it.
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
- **The theme slots** are the `seq.theme` sublanguage ("Themes"); a query in a layer's `data` waits for its sublanguage
  (TODO).
- **A source's transform** (its params besides `format`, `uri`, `index`, `value`, `assembly`): what the file's columns
  mean, written once beside its location, every one optional (Data, Records): `fields` (a role or any other field of the
  ranges, by column name or dotted path), `records` (the path to the records), `shape` (`records`, absent, or `matrix`:
  an array of arrays a grid reads whole as its values, found by `records`; a param because the format cannot say which a
  document is, and a matrix is no record form), `letters` (the path to the letters of the source's own system), `system`
  (its name; absent, the source node's id), `base` (`0` or `1`; absent, `1` for a file, `0` inline), `ends`
  (`inclusive`, `exclusive`; absent, `inclusive` for a file, `exclusive` inline), `types` (a field's type: `int32`,
  `float32`, `float64`, `category`, `string`), `attribution` and `licence` (carried, not read). Why on the source, not
  the reading node: the conventions are the file's whoever reads it, the store keys entries by source, and several
  readers of one file share one open; an edit of the transform is a patch of the source, which reopens it in place (one
  undo step, no side-by-side preview).

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

`snapshot` nodes in the `snapshots` container (`/snapshots/<id>`): `label`, and one of `patch`, a stored patch in
`seq.patch` (`patchLanguage(SeqLanguage)`, `SeqPatchLanguage`) with no `base`, and `state`, a full state (a `spec` node
of the language itself, its refs resolved inside it); neither or both is `seq.lang.snapshot.content`. Beside them:
`from` (a ref to the step a patch applies over: a delta), `keep` (`full` or `spec`, the encoder's form), `view` (the
view and each section's viewport as positions) and the timing (`duration`, `easing`, `linger`). What they mean is
Stories.

- **A step applies over whatever the spec is**: entering one is navigation, never an edit (the viewer's `SpecModel`).
  Why no base: a story is written once and shown over every later revision.
- **Its addresses are refs of the spec** (`repr/patch`'s stored-patch rule): a rename or a move rewrites them, and link
  integrity counts them. A sibling named by a `place` is not rewritten: such a step no longer applies after the rename,
  and is shown as not applying.
- The builder writes one with
  `b.snapshots.snapshot(id, { patch | state, label, from, keep, view, duration, easing, linger })` (a `Patch` or its
  node; a base throws).

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

## Themes

`theme/` (the entry `theme/index`) and `lang/theme-language.ts`: what a layer's channels show, as data. A channel slot
holds a value (a color, a number, a word) or a form of the `seq.theme` sublanguage; the layer type compiles it once a
layout per token table, held domains and the reader's substitution (`LayerChannels`, kept while they agree) to a
`Channel`: a constant, or `at(datum)` over a mark's datum (a field, its value, its letter, its extent).

- **Forms**: `uniform` (one value), `field` (a measure through a palette or a scale), `rules` (clauses tried in order,
  then a default), `scheme` (a registered scheme by name, its params overriding), `letters` (by letter). Palettes:
  `categorical` (values to colors; unlisted values `missing`, `others` or `hash`: FNV-1a over UTF-8 into the list),
  `direct` (the value is the color), `binned` (edges), `continuous` (stops interpolated in Oklab through a 256-entry
  table). Clauses: `eq`, `ne`, `lt`, `le`, `gt`, `ge`, `one-of`, `in`, `within`. Measures: `field`, `value`, `letter`,
  and the measures a host registers (a column's share, a reference row) offered only where registered. Why words for
  clauses: a kind is a word in the description.
- **Scales**: `linear`, `log`, `sqrt`; a `domain` with an open end follows the data in view, computed when the view
  settles and held while it does not; `non_positive` on `log`: `missing` (before any bin or end color) or `clamp` (the
  domain's low end) (`seq.render.scale.log`). A center (`mid`) makes a number's default palette the diverging one. A
  `scale` node (`/views/<id>`) is shared: every layer and channel naming it reads one domain, its open ends the union
  over them at settle (one relayout).
- **Palettes are Mol\*'s value palettes**, evaluated by `@molstar/core/util/color/value-palette-eval` (scale, mid,
  unlisted hash; the module without Mol\*'s named lists, which `value-palette` re-exports with them), converted both
  ways with MolViewSpec's palettes, their scale with them (`toMvsPalette`; `fromMvsPalette` and `fromMvsScale`: type,
  domain, center, `non_positive`, the colors past the ends). Why: one palette model across Mol\*. A named list is one
  the schemes use (`NamedLists`: viridis, cividis, magma, blues, rd-bu, okabe-ito, tableau10, category10) or a token
  list (ramp, diverging, categories); no other of Mol\*'s lists. Why: a spec reads the same wherever Seq\* runs, and a
  list added to Mol\* is not one here until named. Each is imported by its constant (`ColorListViridis`, …), never
  through `ColorLists`: a bundle ships those eight lists, none of Mol\*'s others.
- **Schemes** (`SchemeRegistry`, `builtInSchemes()`): residue schemes as data from core (Clustal X, Zappo, Taylor,
  hydrophobicity, helix, strand, turn, buried), nucleotide ones, ramps and categories (viridis, cividis, magma, blues,
  rd-bu, pLDDT's bands, Okabe–Ito, Tableau 10), secondary structure, Giemsa; each says its subject, the measures it
  needs, whether it is safe under the three color-vision deficiencies (checked: 0.06 Oklab apart under Machado 2009) and
  the safe one a reader's substitution draws instead (`theme.cvdSafe`). Every substitute is itself safe (a test over the
  registry): an unsafe protein scheme (`classes`, the default, included) substitutes to `classes-safe`, `classes`'
  groups on Okabe–Ito hues, the same in both tables. Why: the substitution must never draw an unsafe scheme. A scheme's
  literal colors are the light table's; the dark table's are derived (lightness mapped into the dark band, Oklab L
  0.30–0.62, hue kept) unless it is safe as it is. The band is bounded by the ink: above about 0.62 neither ink reaches
  4.5:1. Provenance is `theme/provenance`'s (`schemeProvenance`), read by a picker, never by drawing.
- **Colors as spelled** (`theme/spelling`: `parseColor`, `isColorSpelling`, `resolveColor`) need no color space: the
  spec reader's checks import them alone.
- **Legends** are data (`Legend`: swatches, a ramp with ticks placed by its scale (a log scale's decades, a center at
  the middle), steps, sizes on the channel's own scale (a radius's `sqrt`); the vocabulary listed, or the values painted
  where unlisted values are hashed: `seq.render.legend.unlisted`), merged when they say the same (two layers on one
  scale: one legend), and converted to Mol*'s legends. A letters form's swatches name their letters: a letters layer
  passes the letters it laid out (`painted`), so its key lists only those (an RNA's U, never its T; no gap swatch where
  there is no gap), then a color once; before any is drawn, every letter.
- **Builder**: `theme.uniform|field|rules|scheme|letters|categorical|…` and `measure.*`; a layer's `field` form over an
  inline table's text column gets the column's values listed (the legend is the vocabulary); a `ranges` layer over a
  `category` column, a `variants` one over a `consequence` column and a `pairs` one over a `type` column default to a
  field form on it.
- **Checks** (`seq.lang.theme.*`): a form or palette that does not exist (`form`), a color that does not read (`color`),
  a channel a form cannot fill (`channel`), a field of the wrong type (`type`), stops out of order (`stops`), a scale
  that does not exist (`scale`), forms nested past `themeDepth` (4, `depth`).
- **Filters**: a theme never hides; a filter does. An item layer's `data` written as a query
  (`<annotation>/(not :of (item :<field> v))`) leaves out the rows its form rejects as the shared item reader reads them
  (`item` by `:id` or a field, a list param any of its values; `not`, `and`, `or`), compiled once a prepare
  (`compileRowFilter`); the track's label column shows "n filtered" (`TrackLayout.filter`), a hit on it
  `target: 'filter'`. A form not read row by row filters nothing and is `seq.render.filter`. **A leading item form is
  the track's base** (`<annotation>/(item :category Domain)`, alone or before a trailing `(not …)`:
  `ResolvedLayer.data.base`, `LayerSource.base`): its rejected rows are left out and not counted as filtered, the
  trailing `(not …)` is the user's (`data.filter`) and counted. Why: one read of a file serves a layer a feature type,
  and the chip counts what the user hid, not what the track is. Any other lone form is the user's, as before; two forms
  that are not an item form and a `not` are read as neither.

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
  systems; **rects subtract** rects over the same two systems: `A × B` less `C × D` is `(A − C) × B` and
  `(A ∩ C) × (B − D)` (at most two rects, their sides of several ranges; four with one range a side), so the cells stay
  exact and a rect selection's Alt+drag is a Location; rects never merge and parts never subtract rects. A rect's
  complement is refused (`seq.query.form.unsupported`).
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
| `(row :key k)`, `(row :key [k …])`             | alignments        | Location: the rows by key (label), positions of its rows system (`seq.query.row.unknown`): a rect's side        |
| `(item :id k)`, `(item :<field> v …)`          | annotations       | Elements (every param matches; a list param matches any of its atoms; a field the annotations lack is an issue) |
| `(items :overlapping <location>)`              | annotations       | Elements (every item without `:overlapping`)                                                                    |
| `(node :key k)`, `(chain :id k)`               | nodes             | Nodes, addressed by the form itself (`/data/1tup/(chain :id A)`), so location forms may follow a chain          |
| `(extent :of <elements>)`                      | any               | Location                                                                                                        |
| `(and :of […])`, `(or :of […])`, `(not :of x)` | any               | the operands' one type                                                                                          |

- **Combinators take `:of`**: a form's params are keywords only (the address grammar), so `and`, `or` and `not` read a
  list or one operand under `:of`. Operands of different types are `seq.query.form.type`.
- **`not`** complements within what the prefix names: its systems' positions, its annotations' items, its nodes'
  children with ids. Dynamic nodes are not in that universe yet. A rect has no complement; **inside `and`**, a
  `(not :of <location>)` operand is subtracted from the others' intersection instead (the same positions for parts), so
  `(and :of [<rect> (not :of <rect>)])` is a rect less a rect, exact (Location).
- **`row`** reads the model's `row(rows, key)`: an inline alignment's labels; a source's rows are not held there, so a
  key of one is `seq.query.row.unknown`.
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

The viewer's `address.*`, `copyLink.*` and `links.*` items; `DefaultAddressBounds`, `DefaultLinkOptions`.

| Item                           | Default          | Why the default; the alternatives                                                                       |
| ------------------------------ | ---------------- | ------------------------------------------------------------------------------------------------------- |
| `limit`                        | 10,000           | the most nodes or elements an address names: past it an issue, never a truncation                       |
| `maxDepth`                     | 64               | how deep a walk goes: an address comes from a link or an agent                                          |
| `refDepth`                     | 8                | how many refs a chain passes (`DefaultRefDepth`): a cycle or a runaway chain is refused, never followed |
| `prefer`                       | `['ref', 'ids']` | the forms the way back tries, in order: a ref naming exactly it survives edits; `['ids']`               |
| `numbering`                    | `ordinal`        | the labels the way back prints (1-based, every system has it); a declared numbering that reads back     |
| link `prefix`                  | `seq-`           | this viewer's keys in a link, so one href carries Mol\*'s and Seq\*'s                                   |
| link `maxLength`, `maxActions` | 8,192, 16        | a link's bounds: a link is untrusted input                                                              |

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

- A **format runtime** (`FormatRuntime`: `name`, `formats`, `open(source, { assets, signal, limits })`) opens a
  `SourceEntry` (`address`, `format`, `uri` or inline `value`, an `index` for a format read with one, and the source's
  other params, its transform among them) to an `OpenedSource`: its `coverage`, `read(request, signal)` and the `issues`
  its opening said without refusing it (a file's rows left out), which whoever opens it reports at the source's address
  (the viewer's store model, `attachSources`). `assets` is a narrow port (`readText(uri, signal)`,
  `readBytes(uri, signal)`; both required: every adapter is the package's own, and a runtime does not test for a
  method). Why: a runtime is tested with a fake, never the network. A file that cannot be read is a `SourceOpenError`
  with `seq.source.unreadable`.
- **`limits`** (`maxRows`, `maxBytes`): a parse's bounds, checked before it allocates; whoever opens a source passes the
  store's settings (`maxRows`, 1,000,000; `maxBytes`, 64 MB: the viewer's source sync, the overlay's sources,
  `attachSources`), `DefaultOpenLimits` the defaults. Why: `JSON.parse` makes several times a file's bytes, and the
  assets policy's 256 MB is a download's bound, not a parse's.
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
- **A row's identity** is its `id`, else its `row`: a range over a tile's edge is in each tile it touches, with its own
  start, and readers dedupe by identity and start; rows of one id with other starts are pieces of one item (the exons of
  a transcript, read in different tiles).
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
  path) and never throws. Why a tree value: one lossless encoding, the codec's. **A file is JSON or MessagePack**, told
  by its first byte: a MessagePack map's marker (`0x80`–`0x8f`, `0xde`, `0xdf`) is MessagePack, anything else JSON (so
  an HTML page or truncated JSON says JSON's error), read through the tree's codec of that encoding; the runtime reads
  bytes (`readBytes`), refuses a file over `limits.maxBytes` before parsing it (`seq.source.too_large`), and loads the
  MessagePack codec only when a file needs it (a dynamic import: a reader of inline specs does not pay for it).
- **Its runtime** joins the pairs of the request's system (by name or alias) and level that overlap the region, in
  position order, each cut to its part; a gap no pair holds is `not-served`, naming its positions (the format is
  explicit: an empty region is a pair with empty tables). A failure among the pairs is the answer; an unlisted system or
  level is not served. A list without ids gets `row` from its place in the file. Inline data (`inlineSeq`) is one pair
  over the whole system at level 1, its length in the coverage.

### Asynchronous sources

- **A region's answer holds the items overlapping it**, not only those starting in it: the store reads whole tiles, and
  a view wholly past a tile's edge reads that tile alone, so an item from the tile before would have no mark. An item in
  two tiles is in both answers, one range (`id` or `row`); the layout keeps it once (key and start). `sliceTable` and
  the test sources follow it; a tiled reader answers a tile with the rows overlapping it.

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

### Cells

- **`cells`** (an alignment's rows on 1D tiles of its rows system, columns whole): at level 1 `codes` (`start` the first
  row, `columns`, a byte a cell row-major: the letter's ASCII, gaps `-` and `.`) with the rows table (`label` a
  category, a row's key; `first`, its first residue's number, EC-39); above it `majority` (a row bin per column: its
  majority letter byte, 255 none, and that letter's share, 0–255), aligned to the level. Why columns whole: an alignment
  is hundreds to a few thousand columns wide, so a tile of 1,024 rows stays within `maxBlockValues` up to 16,384 columns
  (wider is `too-large`), and the store's keys and tiles stay one-dimensional. Cut and joined by rows (`sliceTable`,
  `joinTables`); not written in the `seq` format yet (no reader needs it).
- **An alignment source** (`openAlignment`, `cellsSource`): two systems, `rows` (a position a row, data order) and
  `columns`; `cells` served on `rows` at level 1, a tile a view of the bytes made once. Inline rows read through the
  `alignment` runtime (`AlignmentRuntime`, which resolve binds to an inline `alignment` node); files through
  `AlignmentFileRuntime` (`aligned-fasta`, `stockholm`; the readers' entry). `cellsSource(rows, columns, rowsOf)` makes
  a tile's rows on demand (a generated alignment: no row held before it is asked for, EC-41).

### Alignment readers

`data/alignment-readers.ts`, an entry of its own (`@molstar/seqx-core/data/alignment-readers`: a reader of a spec never
pays for them). Each reads text into plain data and reports what it cannot read, by line, never throwing (EC-40); held
by fuzz tests (mutated text read without a throw).

- **Aligned FASTA** (`readAlignedFasta`) and **Stockholm** (`readStockholm`: a row's letters joined across blocks,
  `#=GF`, `#=GS`, `#=GR` read past, the first alignment of a file): rows of the first row's width, another width or a
  character that is no letter or gap said and left out; gaps `-`, `.` (`~` read as `-`). A row's first residue number
  from its label (`NAME/27-137`, EC-39), else 1. **Insert columns** (EC-38): `#=GC RF`'s `.` and `-`, else a column
  whose letters are all lowercase (A3M's lowercase is on no column). **Pairs**: `#=GC SS_cons` as column pairs, each
  bracket kind (`<>`, `()`, `[]`, `{}`) a nesting of its own and a pseudoknot an uppercase letter with its lowercase.
- **Newick** (`readNewick`): nodes in pre-order (parent, branch length, label, children), leaves in tree order; labels
  bare or quoted (`''` a quote), `[comments]` read past, underscores kept (a label is matched against row names);
  iterative (a deep tree costs no stack), at most `MaxTreeNodes` (1,000,000) nodes; what does not parse leaves the tree
  empty with its issue, and the alignment is kept (EC-43).
- **Leaves to rows** (`matchLeaves`): a leaf names its row by label, else by the label before its last `/`; a row is the
  first leaf's that names it; -1 for none.

### Records

`data/records.ts`, an entry of its own (`@molstar/seqx-core/data/records`: a reader of a spec never pays for it):
`RecordsRuntime` (formats `tsv`, `csv`, `json`, `msgpack`, `gff3`), registered by the viewer's default spec;
`readRecords` the reading without the runtime (tests, a preset).

- **Read whole at open**, converted once into `seq` data (a pair a system over the whole system at level 1: its letters'
  length when the source holds them, else every position a column can, the target's own length standing) and answered by
  region by the `seq` runtime's join (`openSeq`). Why: these files are small, formats are converted at the edge, and one
  answering path is tested already. Files read by region (tabix, BGZF) wait for the Genomes milestone.
- **The syntax**: TSV and CSV are RFC 4180 over the UTF-8 text (`@molstar/core/data/text/delimited`: a BOM skipped, LF,
  CRLF and CR, quoted cells, empty lines skipped; TSV quoted as CSV is); the first row the header, a name twice
  malformed; a row with another count of cells read with the missing ones empty and the extra dropped
  (`seq.source.records.cells`). Cells are text, never trimmed (EC-92 changed: the RFC's rule); a number or position
  parses past the spaces around it. JSON is `JSON.parse` over the text, a plain value (never the tree's codec: a member
  `kind` would read as a node); MessagePack is `@molstar/core/data/msgpack`'s decoder with the tree codec's column
  extensions (`msgPackColumnDecoder`) and nothing else of it. Members are read own only (`__proto__`, `constructor`
  ordinary names, EC-50; a transform's maps are read own only too). A document that does not parse, or a MessagePack one
  deeper than the stack (its decoder recurses; `JSON.parse` does not), is malformed, never a crash (EC-76). Text holding
  U+FFFD is said once (`seq.source.records.encoding`, info): UTF-8 only.
- **GFF3**, a preset of the same reader (the format fixes its columns, base 1 and inclusive ends: a transform's `base`
  and `ends` are not read): nine tab-separated columns (`seqid` the `target`, `source`, `type`, `start`, `end`, `score`,
  `strand`, `phase`, `attributes`), headerless; blank lines, `#` lines and `##` directives are read past and everything
  after `##FASTA` is ignored; `.` is missing (a `null` cell); `key=value;…` attributes become fields by their key, their
  values percent-decoded (a malformed escape left as written), a multi-valued one as its text, `ID` the `id`; `Parent`
  stays a field, no hierarchy is assembled. A row with another count of cells than nine is read with the missing ones
  empty (`.cells`). `fields` maps as for any source (`{ category: 'type', label: 'Name' }`); `seqid` and `ID` are the
  roles, not fields under their names. Why a preset: NCBI serves an INSDC record's features as GFF3, read as published.
- **A row past the system's end is the row reader's** (`collectItems`, `wrap`; render, items): a source's row ending
  past its target's length becomes two rows of one item (`[s, L)` and `[0, e − L)`, one key, one occurrence each) on a
  cyclic system, and is clipped on a linear one (`seq.render.item.range`, warning, once a layout with the count); a row
  starting at or past the end is left out the same way. Why there, not in the runtime: without letters the runtime knows
  neither the length nor that the system is cyclic. A layer of a cyclic system whose source is not inline also wants the
  system's last position (`LayerSource.wraps`, `cyclic`): any row that runs past the end overlaps it, so its tile brings
  the origin's piece though the view is elsewhere; inline rows were wrapped by the model (EC-9).
- **A matrix** (`shape: matrix` on the source, `readMatrix`, `RecordsSource.matrix(expected?)`): an array of arrays
  (JSON, MessagePack, an inline value) found by `records`, served as a grid's values, row-major as `Float32Array`
  (`null` and any cell that is no number NaN, said once); the source serves no region (an empty coverage). Its size is
  checked against the axes' total the grid passes (`expected`, `rows × cols`): a mismatch, or rows of unequal lengths,
  is `seq.source.records.size` (error, the matrix not served; the issue names both sizes). A typed flat array (a
  MessagePack `Float32Array`, an inline one) is kept as it is when it is a `Float32Array`, its dimensions the expected
  size's (a square when none is given: a flat array that is not one is `.size`); any other flat array is copied. More
  than 2^24 values is refused (`seq.source.too_large`). A table format is no matrix (`seq.source.no_data`). Why a param
  and not the format: a document is a document; whether it is records or a matrix is the spec's to say, and the grid
  reading it in the fourth slice's complex PAE is the one reader.
- **Records and columns**: `records` is a dotted path from the root (a member, or an index when all digits; at most 32
  segments); what it reaches is an array of records (AoS: each an object) or an object of columns (SoA: each member a
  list or a typed column, all of one length, checked before any is read); a table file is AoS. Read the same.
- **Fields**: without `fields`, each role from the member or column of its name, and every other member whose values are
  scalars a field under its name (objects and lists are not fields unless a path names one of their scalars). With
  `fields`, a role or a field from its column or path; a column a mapping reads is not a field again. The roles are the
  protocol's: `target` (the row's system, empty: the source's own, named `system` or the node's id; the own system is
  listed only when a row names no target, the source holds letters, or no row names any), `start`, `end`, `position` (a
  row `[p, p + 1)`), `strand` (`+`, `-`, `.`, `?`, `1`, `-1`, `0`; anything else 0, said once), `id` (a category),
  `value` (float32). Another field's type is `types`', else its typed column's, else a text cell's or a boolean's
  `category`, a JSON number's `float64` (`string` held as a category). A path that reaches an object or a list is `null`
  in that row, said once (`.cell`).
- **Positions**: `base: 1` subtracts one, an `inclusive` end adds one (the file defaults, as people, UniProt, GFF3 and
  VCF count; EC-195); an inline `value` defaults to `0` and `exclusive`, as `ranges.rows`, but an inline string `value`
  is a file's text and takes the file defaults. Converted once, at the edge; an end past int32 leaves the row out
  (`.value`). Rows sorted by start once (stable), each keeping its place in the file as `row`. **Values**: with
  `position` sparse (NaN between, from the first position named; one named twice keeps the last, said once, `.cell`),
  without one a value a row from the first position, each at its row's own index (a row left out is NaN there, never a
  shift); more positions than `maxRows` refused. With `start` and `end`, `value` is a field of the ranges only, no
  values table.
- **Typed columns kept without a copy**: a column already the protocol's type for its role (`Int32Array` a start or an
  end, `Float32Array` a value, a category column an id or a category field, any typed field of its declared type), with
  `base: 0`, `ends: exclusive`, one system, every row kept and already in order, and no end clipped; otherwise one copy.
  MessagePack's decoder gives a view where the payload is aligned, else its own copy; a kept view pins its file's whole
  buffer for as long as the source is attached.
- **Issues** (`RecordsIssues`, `seq.source.records.*`; one a cause with its count and its first ten places, a line or a
  record, EC-86): what cannot be read as its format is `.malformed` (no header, a document that does not parse, a
  `records` path to neither form, SoA columns of unequal lengths, a header naming a column twice, too deep) and a matrix
  of another size than its axes' total or of unequal rows `.size`, and a field found in no record `.field` (naming the
  field and listing the columns there are, EC-6): the source unavailable, its `SourceOpenError`'s first issue. One
  severity a code (EC-87): a cell its field cannot hold, read as missing (an object or a list, a number that is none, a
  strand read as 0), and a position named twice are `.cell` (info, the row kept). A row-level misfit leaves the row out:
  a position or value not a number or an end past int32 (`.value`), a position missing (`.missing`, UniProt's
  `UNKNOWN`), an end before its start (`.order`), a row starting before the first position or past the letters' end
  (`.range`; one running past the end is clipped and counted in the same issue, EC-51), hinting the other `base` when a
  first position 0 meets `base: 1` or every row past the end would fit one lower. The opened source carries them
  (`OpenedSource.issues`), and the viewer's store model and `attachSources` report them at the source's address. Why
  left out, not refused: one bad row in ten thousand should not blank a track.
- **Bounds**, the context's `limits`: a file over `maxBytes` is refused before it is parsed, rows over `maxRows` before
  any column is built (`seq.source.too_large`, hinting a format read by region).
- **Configuration**: the runtime is a registry entry (the default viewer spec registers it; a host leaves it out); the
  transform is the spec's (Language); quoting, line ends and encoding are fixed (RFC 4180, LF/CRLF/CR, UTF-8); a path's
  32 segments, ten places an issue and `limits` (the store's settings) bound the work.
- **Performance**: 100,000 TSV rows read under 300 ms (`records-time.test.ts`, heavy lane, best of three); 10,000 in the
  quick lane.

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

`ResolveConfig`, each choice with its recommended default (`DefaultResolveConfig`):

| Item                | Default     | Why the default; the alternatives                                                                                                                                |
| ------------------- | ----------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `maxPathSteps`      | 4           | counted steps a path may take (EC-22): bounded search                                                                                                            |
| `maxPathCandidates` | 16          | candidate paths kept (EC-22): bounded search                                                                                                                     |
| `encoding`          | `cheapest`  | each map in its cheapest form; `runs` or `expanded` hold every bounded map so, for tests and comparisons                                                         |
| `openQueryLength`   | int32's max | what the query forms read on a system of unknown length, so labels resolve before its source opens; layout never reads it (an axis of unknown length is pending) |
| `maxNestedDepth`    | 4           | the levels of nested sources the overlay requests (EC-179): a spec's sources are untrusted                                                                       |

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
- **Pairs and grids** are held whole (`resolve/held.ts`): a `pairs` node's ends and a `grid` node's values read once
  into typed arrays on the entry (`pairs`, `grid`), on two systems (`target` and `ySystem`), with no access plan: they
  are inline, read whole, and their reductions are their layers'. A grid whose values are not `|x| × |y|`, or more than
  `DefaultMaxBlockValues` (the bound the store checks a block), is dropped with `seq.resolve.grid.size` saying the count
  it expected. Pairs and grids read by region wait for the Genomes milestone.
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
- **A point** (`projectPoint`, exported from `resolve`: the viewer's follows read it too): the continuous image of a
  fractional position, place by place: across a codon step residue `x` is fine place `3x`, so each nucleotide of a codon
  is an edge of its own and a codon split by an intron follows its nucleotides, never the intron; where place `floor(x)`
  starts on the target plus the fraction in its run's direction. A fine place within 1e-9 of a whole one is that one (a
  third is not exact in floating point). `projectEnd` is a window end's image: at a whole position the far edge of the
  place before it.
- **A direction** (`pathDirection`, `pathDir`): the direction of the runs covering most of the path's length, `mixed`
  when it has runs both ways (an inversion). Why the majority: a direction that changed while panning would turn an axis
  over under the hand.
- **Coding ranges** (`codingRanges`): the nucleotide side of a path with a codon step, sorted and merged, either way
  along it: what a section's `pieces: coding` shows.
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
- **A section's layout**: its own `layout`, else its view's (`tracks`, `text`, `rings`, `grid`), with its layout's
  params as written (`params`: `group`, `numbers`, `group_labels`, `align`, `max_width`, `height`, `rotation`, `labels`,
  `arrange`, `copies`, `aspect`); the layout reads them with its defaults. Why on the section: a view mixes layouts (a
  ring and its linear detail; a text figure under a tracks overview). What the figure means (group size, numbers) is the
  spec's; spacing is the render config's.
- **A section**: an axis, its tracks, and its follow; its `blocks` when it has `block` children; a grid section its `y`.
- **Grids**: a grid section's `y` (`GridAxis`) is its second axis, down the plot: `of_y`'s system whole (`of`'s when
  absent), top down, and what it shows: `follow_y: x` (the default) x's window (`followsX`; through the path from x's
  system when `of_y` is another), a section's region through `follow` (the path from its system), `none` a window of its
  own (`viewport.y`). A y that cannot follow (no path) is `seq.view.follow.unplaced` and keeps its own window. A section
  whose `of` names no system (with blocks: its first block's) is `seq.view.section.system`; its layers are not placed.
- **Blocks**: each its axis (its system, its `range`, the section's direction), its first display unit in the section's
  (`offset`: 0 in a stack), its label, the blocks it stands for (`members`), and each of the section's tracks' layers
  resolved onto its system (a path a block: placed on each block its data reaches; a layer that reaches no block says
  why once, by its first block's issue). Copies (the same system and range) are one block labeled `×N` with
  `copies: 'merge'` (the default: a copy group's sequence said once, its cost one block's, EC-65, EC-67), a block each
  with `each`. The section's axis is the blocks' display units in order: on one system in a row (or on a ring), the
  blocks' pieces end to end (their breaks between them, a follow along one axis); else units `[0, total)` on no system
  (its `system` the section's address, `pieces: 'blocks'`), the total the sum in a row and the longest block in a stack.
  Why: a position on them is no block's, so a selection, a focus, a bookmark or a cursor on a system never reads them as
  its own; a block's hits carry its system and its index (`block`), and a consumer maps them through the block's axis
  and `offset`. A block that names no system is `seq.view.section.system` and left out.
- **An axis**: its system, its length unless unknown, cyclic, and the positions shown as `[start, end, …]` pieces in
  display order: the whole system, its `range`, or for a cyclic range that wraps the origin two pieces, from its start
  to the origin then on (`[8, 10, 0, 2]`). Without a known length, showing all of the system or wrapping, `ranges` is
  absent: the layout waits for the length. Why: one shape for a wrapped range everywhere, the order it is drawn in.
- **Its direction** (`dir`): the section's `direction` resolved, never read off its layers (a layer added never turns a
  section over): `genomic` 1, `reversed` −1, `follow` the leader's direction times the follow path's (`pathDirection`,
  the majority of a path both ways, `seq.view.direction.mixed`), `{ along }` the direction of the path from that system
  (a map: its own linear part) to the section's. One that cannot be resolved (`follow` written with no placed leader,
  `along` a system no path reaches) is `seq.view.direction.unresolved` and `genomic`. One direction an axis (a direction
  a piece is a tension until Genomes).
- **Its pieces** (`pieces`, the form as written for an export, and `ranges`): `whole`; `coding`, the coding ranges of
  the follow path, else of the path it reads along; a ranges node's inline items on its system, a codon map's coding
  ranges onto it, or their union; widened by `pad` each side (merged), within `range`. Start below end, in display
  order: on a reversed axis the last piece first, a piece `[s, e)` showing `e − 1 − k` at its offset `k` (EC-32). What
  reaches nothing is `seq.view.section.pieces`, and the section shows `whole`.
- **Offsets and positions** (`axis.ts`, exported for the viewer): a fractional offset `u` in a piece `[s, e)` at offset
  `o` is the continuous position `s + (u − o)` forward, `e − (u − o)` reversed (`positionAt`, `offsetOf`); a cell's left
  edge (`cellOffset`) and the element at an offset (`elementAt`). An offset on a break is the next piece's first edge as
  a window's `start` and the previous one's last as its `end`, so a window ending on a break shows nothing hidden.
  Binary searches over the pieces' first offsets, kept per pieces array (no allocation once made).
- **Breaks** (`axisBreaks`: `[offset, hiddenStart, hiddenEnd, …]` in display order): where two pieces leave positions
  hidden between them (genomic order; `start > end` across a cyclic origin); two contiguous pieces, or contiguous across
  a cyclic origin, have none. Zero units wide, so the cells either side touch and a codon follow on a coding axis is
  exact. `breakHolding` finds the break hiding a position.
- **Scroll groups' key** (`axisKey`, `axesEqual`): system, length, pieces and direction.
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
- **A layer**: its address (with an id) and index, its type and fallback chain (`LayerTypes`: `secondary` and `variants`
  fall back to `ranges`, the rest (a `frame` among them) have none), its other params as written, its theme channels as
  written (`color`, `ink`, `size`, `opacity`, `line`, `arrow`, `shape`), its data and the static path from the data's
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
- **Navigators** (`ResolvedView.navigators`, `resolveNavigatorNode` for a host's definition of the same shape): each
  with its bound sections the view has (`of`, in order; `missing` those it lacks), its params with their defaults
  (`axes` `positions`, `content` `axis`, `place` `top`, `right` for rows), the system it draws (`system`, else its first
  section's; `sectionAxis` when it draws that section's own axis: its range, pieces and direction) and the static path
  from each bound section's system to it (identity, a region into its chromosome, a protein into its gene).
- `seq.view.*` (`ViewIssues`): `section.system`, `layer.type`, `layer.data`, `layer.shape`, `layer.unplaced`,
  `follow.unplaced`, `section.pieces`, `direction.unresolved`, `direction.mixed` (info); a test checks the catalog.

## Render

`src/render/` (`@molstar/seqx-core/render`): tier 3, the resolved view laid out for a size, a viewport, a target and the
store's entries, as draw lists; hit-testing on the same layout; the recording and SVG backends; `renderSvg`.

### Configuration

`RenderConfig`, each choice with its recommended default (`DefaultRenderConfig`, `renderConfig(over)`). The sizes are
the UX foundations' numbers, CSS px, settled by the foundation stories; where no reason is given, none is recorded.

| Item                                                           | Default                                                                                             | Why the default; the alternatives                                                                                                                                               |
| -------------------------------------------------------------- | --------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `inset`                                                        | 8                                                                                                   | the page's inset around everything                                                                                                                                              |
| `labelWidth`, `labelGap`                                       | 112, 8                                                                                              | one label column a view, so the data columns align                                                                                                                              |
| `trackGap`, `sectionGap`                                       | 4, 8                                                                                                | gaps between a section's tracks and between sections                                                                                                                            |
| `heights`                                                      | letters 20, ranges 36, values 32, ruler 22, secondary 16, variants 32, frame 16, pairs 64, grid 240 | a track's height when the spec gives none: the largest of its layers' types'                                                                                                    |
| `laneHeight`, `laneGap`                                        | 12, 2                                                                                               | a lane of ranges and the gap within it                                                                                                                                          |
| `density`                                                      | letters 10, cells 3, dots 1, bars 1/64                                                              | each density step's lowest pixels a position (below)                                                                                                                            |
| `densityUnits`                                                 | letters `css`, the rest `device`                                                                    | legibility is in CSS px; what bounds the lists is a mark a device pixel                                                                                                         |
| `letterSizes`                                                  | 9 to 14                                                                                             | a short scale; the largest that fits the cell is used                                                                                                                           |
| `monoAdvance`                                                  | 0.6 em                                                                                              | a letter's size known from its cell with no measuring                                                                                                                           |
| `labelFont`, `itemFont`                                        | `12px sans`, `10px sans`                                                                            | a family role a backend maps to its stack                                                                                                                                       |
| `fade`                                                         | 12                                                                                                  | how far past its max width a text fades out                                                                                                                                     |
| `squircleRatio`, `squircleMax`, `squircleMin`                  | 1/3, 8, 6                                                                                           | a radius a third of the height, capped; none below the minimum                                                                                                                  |
| `outline`                                                      | hover 1.5, selection 1, cursor 1                                                                    | outline widths                                                                                                                                                                  |
| `emphasis.selectionPlace`                                      | `over`                                                                                              | a tint under the data is hidden by opaque cells and a selection then redraws the data; `under` (to confirm with the author)                                                     |
| `emphasis.selectionShape`                                      | `column`                                                                                            | a selection is of positions: one band a section across its tracks and gaps reads as one (the author's look); `tracks`, a box a track                                            |
| `emphasis.selectionTint`, `emphasis.selectionOutline`          | 0.18, 0.6                                                                                           | a soft tint and one thin outline (the author's look)                                                                                                                            |
| `emphasis.highlightTint`                                       | 0.3                                                                                                 | a highlight channel's tint                                                                                                                                                      |
| `emphasis.projected`                                           | 0.5                                                                                                 | what is projected onto a linked section, a factor of its own look's opacity: the same look, lighter                                                                             |
| `emphasis.hoverOffset`                                         | 1                                                                                                   | the hover's outline outside its mark, concentric, never touching the selection's                                                                                                |
| `emphasis.cursorInk`                                           | `msui-text`                                                                                         | the cursor neutral, never the selection's accent (the author's look)                                                                                                            |
| `axis.directionMark`                                           | `label`                                                                                             | a reversed axis says so in its first label (an arrow and "reversed", or "along" what it reads along), in the picture and the export; `none` for a host whose own chrome says it |
| `margin`                                                       | 0.5                                                                                                 | what a want adds on each side, a fraction of the view (the store's `margin` priority)                                                                                           |
| `followArc`                                                    | `shortest`                                                                                          | a cyclic follower takes the shortest arc holding the image (EC-12); `hull`                                                                                                      |
| `maxTemplateTracks`                                            | 200                                                                                                 | the tracks one template makes; past it `seq.render.template.limit`: a source's nodes are untrusted                                                                              |
| `nodes.holdShown`                                              | true                                                                                                | a layout lists the nodes it shows, for its owner to hold in the store; false: nothing held by showing                                                                           |
| `ranges.overflow`                                              | `extra-lanes`                                                                                       | only the items past the lanes collapse into a density row; `layer`, the whole layer                                                                                             |
| `values.aggregate`                                             | `minmax`                                                                                            | below a device pixel a position, each column's extent, so negative extremes show; `max`                                                                                         |
| `values.coarser`                                               | 1                                                                                                   | the coarsest level whose bin is no wider than k device pixels' positions (below)                                                                                                |
| `values.direct`                                                | true                                                                                                | data on the section's own system placed with no mapping (below); false: the general mapping, the same lists                                                                     |
| `ruler.labelGap`                                               | 48 px                                                                                               | the least between two ticks (the foundations' number); labels clear each other by an eighth of it                                                                               |
| `text.scale`, `text.groupGap`                                  | 12 px a position, half a cell                                                                       | the letter size a text opens at (its zoom); a space between words                                                                                                               |
| `text.lineGap`, `text.bandGap`                                 | 8, 2                                                                                                | a line block reads as one when the line gap is the larger                                                                                                                       |
| `text.collapseEmpty`                                           | true                                                                                                | an empty band separates what belongs together (EC-96); false: uniform lines                                                                                                     |
| `text.maxLines`, `text.maxHeight`                              | 100,000; 600 px                                                                                     | a figure, not a genome; a text in a host page scrolls inside its box past it                                                                                                    |
| `text.pieceBreak`, `text.numberGap`                            | `none`; 8                                                                                           | lines stay aligned, the boundary marks the edge; the gap between a line and its numbers                                                                                         |
| `text.heights`                                                 | letters 18, ranges 4, values 3, ruler 14, secondary 6, variants 10, frame 14, pairs 24, grid 120    | a text's tracks thin: a ranges track reads as underlines, a values track as a strip                                                                                             |
| `rings.hole`, `rings.maxDiameter`                              | 0.35 of the radius; 720 px                                                                          | room for the title, chords readable; a map, not a wall (`RenderConfig` preset; 0 to 0.8; px)                                                                                    |
| `rings.blockGap`                                               | 1°                                                                                                  | a gap between blocks reads as a boundary (Circos)                                                                                                                               |
| `rings.letters`                                                | `upright`                                                                                           | readable at any angle; `tangent`                                                                                                                                                |
| `rings.center`                                                 | `title`                                                                                             | the name and length in the hole; `none`                                                                                                                                         |
| `rings.labelGap`, `rings.leaderBend`                           | 8, 12 px                                                                                            | leaders fan out over the ruler, never through one point                                                                                                                         |
| `rings.window`                                                 | `followed`                                                                                          | a window means only what a follower shows; `always` (idle, quiet, when nothing follows)                                                                                         |
| `pairs.maxArcs`, `pairs.hitPx`                                 | 4,096; 4                                                                                            | bounded by the screen: past it arcs binned by their ends' device columns (`seq.render.pairs.limit`); a curve is thin, so a hit within 4 px                                      |
| `grid.rectMinPx`                                               | 12                                                                                                  | rects with outlines only where cells are few and large enough to read; image tiles below                                                                                        |
| `grid.buildBudget`                                             | 65,536 samples a frame                                                                              | a frame while moving under 8 ms; a tile past it draws a coarser level meanwhile, never empty                                                                                    |
| `navigator.window.minPx`, `.hitPx`, `.touchHitPx`, `.handlePx` | 16, 12, 24, 5 px                                                                                    | a window under a few px cannot be seen or grabbed; a press in the middle still pans; a finger is wider                                                                          |
| `navigator.window.floor`                                       | `centred`                                                                                           | the drawn window says where; `thumb` (a scrollbar's mapping)                                                                                                                    |
| `navigator.trackScale`, `navigator.trackMinPx`                 | 0.25, 2 px                                                                                          | a miniature reads as the view; a thinner track is not seen at dpr 1                                                                                                             |
| `navigator.marks`                                              | selection, highlights, cursor                                                                       | where is what I picked; any subset                                                                                                                                              |
| `navigator.rail.grownPx`, `.grabShare`                         | 6 px, 0.4                                                                                           | a target in use; an alignment ruler's labels select columns, its bottom pans                                                                                                    |
| `internedStrings`                                              | 4,096                                                                                               | strings draw lists keep interned across frames, so a frame drawn again interns nothing; past it a full redraw starts the table again                                            |

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

- `SeqTokens`: every color the lists use (`msui-*` chrome, `seq-*` data: tracks, ink, missing, pending, conflict, the
  break ink (`seq-break`), residue and nucleotide classes, eight categories (chosen apart under the three color-vision
  deficiencies in both tables), the density ramp, a five-stop ramp and a diverging triple, secondary structure
  (`seq-ss-helix`, `-strand`, `-turn`), `seq-variant`, `seq-ruler`, the navigator's window, grip and wash, the
  centromere, a letter cell's edge), as `0xRRGGBBAA`, the logo's chroma floor (`seq-logo-chroma`, a number) and the font
  stacks of the two family roles. Node has a light and a dark table (`LightTokens`, `DarkTokens`); the viewer reads its
  own from CSS.
- **A letter cell's edge** (`seq-cell-edge`): in the dark table a hairline, 1 device px, around each cell drawn with its
  letter (the letters layer's `letters` step, a rows layer's letter cells), at 3:1 against the page; none (transparent)
  in the light table. Why: the dark scheme colors sit at 1.4–1.7:1 against the dark page and cannot go lighter without
  their ink failing 4.5:1; an edge costs no figure's colors, a lighter dark page every figure's.
- **Its version is a hash of its content**, part of the draw lists' key: a theme switch redraws, an equal table does
  not. Why a hash: pure by contract, so an export from Node and a viewer agree on equal tables.
- **A spec's color** is a hex literal, `token:<name>` (`seq-<name>`, then `msui-<name>`, then `<name>`) or
  `light-dark(a, b)`; anything else is `seq.render.color` and the layer's default. **Ink follows its fill**: the dark or
  the light ink, whichever contrasts more.

### Text

- The port is `measure(font, text)` in CSS px with one `TextMeasure` (`@molstar/core/util/text-measure`) a font behind
  it (`createTextMetrics`). A font is shorthand naming a family role (`12px sans`, `10px mono`); a backend writes the
  role's stack. Node measures with the bundled table (`bundledTextMetrics`): Arial-compatible sans advances, a 0.6 em
  mono. Letters are never measured: centered in their cells, one size a row, the largest of the scale that fits.
- The port has a `version`, bumped by `clear` (a web font arrived): part of the draw lists' key.
- Font sizes are read from the config once a layout, never per label per draw.

### Draw lists

- CSS pixels, origin top left; a list per pass, drawn in order: `background`, `data`, `window`, `highlight`,
  `selection`, `hover`, `cursor`, `labels` (`Passes`; `window` to `cursor` are `EmphasisPasses`). `window` is the first
  of the emphasis, under the highlights: a window is no data emphasis, and a highlight channel stays above it. It is the
  layouts' own (their chrome): no layer draws there. The recorder writes an empty pass as its header.
- **Primitives** as structs of arrays in growable typed arrays kept frame to frame (`reset` zeroes counts): `rect` (x,
  y, w, h, fill, stroke and width, a squircle radius, 0 square, and which ends are round), `path` (move, line, cubic,
  close in shared arrays; fill, stroke, dash), `glyphs` (a run of letters in one font, an x a letter at its center, one
  baseline, a color a letter or one; a **placed** run gives each letter its own y, angle and vertical scale: letters
  along a ring, a logo's stacked letters, a row's letters; one form for the three, so one code a backend), `text` (an
  interned string, x, y, anchor, font, color, its measured width and a max width: faded past it, never ellipsized),
  `image` (RGBA bytes the layer owns), `clip` and `unclip`. A pass's `ops` keep the order its primitives were added in.
  Colors are `0xRRGGBBAA`; `NoColor` (0) is none.
- **Bounded by the screen**: at most one primitive a device column a lane, or one a cell when cells are a device pixel
  or wider; tracks below the page are not laid out.
- **Nothing allocated once warm on the draw and emphasis paths** (a full draw of a laid-out view, a hover, a selection
  or cursor redraw; held by the allocation tier, `frame.quick.ts`): one `PassWriter` a pass, kept by the lists
  (`writer`); the interned strings kept across frames (`internedStrings`); the density strip's bytes kept while the size
  holds. Why the rest: in compiled JavaScript a number passed to (or returned from) a call that is not inlined is boxed,
  an allocation, and which calls are inlined depends on the order code ran in. So each pass of a layer type and of the
  view is its own function, called through one table site (`passDraws`, `drawPass`: no call target known there, so each
  is compiled on its own with the writer's calls inlined into it); a loop calls the writer at one site a primitive kind;
  the emphasis passes' rects go through the writer's `args` (an array's numbers are never boxed; `curveArgs` a move and
  one cubic from a typed array, styled from `args`: an arc) and snap inline, and so do texts and letters with their
  colors (`textArgs`, `letterArgs`: a color past a small integer is a number boxed too); a mark drawn only on some views
  (breaks, hidden items' ticks, gap marks, carets, three-letter codes) is a function of its own, called only when it has
  something to draw, so the pass it belongs to stays small enough to have its writer calls inlined; nothing in a draw
  makes a closure, an iterator, a destructured array or a number for a call, or reads a property of an imported constant
  object (`Anchor.start`: read once into a module constant; a site whose feedback was cleared deoptimized a warm draw
  without end), nor an imported binding in a loop (`NoColor`: a property read under a module runner, it boxed a number a
  frame in the window's grips); the hit's numbers are read where they are used. A layout sets a buffer it reuses in
  place, never to a length of 0 (that lets the array's storage go).
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

`Emphasis` (`layer-type.ts`), drawn by the view (`tracks.ts`) in the emphasis passes, a strip at a time in its frame
through its warp (Section layouts), over the strip's **band column** (its bands and the gaps between them;
`SectionLayout.columns`, in tracks `SectionLayout.band`; the shape `column`, the default) or over each band's box
(`tracks`), `emphasis.selectionShape`, on the strip's transform and system. Why `column` by default: a selection is of
positions, not of a track, and one band reads as one selection; why a strip: the selection is one look in every layout
(a band on a line in text, a sector over the rings on a ring). A layer with a second axis (`axes: 2`) paints the parts
on its y systems itself, in the emphasis passes' `draw` (it alone knows its y transform); the x painting stays the
strip's.

- **The order of emphasis** is fixed, below to above: data, highlight channels, selection, hover, cursor; a higher one
  never hides a lower one entirely (tints are translucent, outlines are thin lines that never touch).
- `highlights`: the channels in order, each a Location and a color, tinted at `emphasis.highlightTint` on every system
  its Location has (the `highlight` pass), square.
- `selection` (a Location): a soft accent tint at `emphasis.selectionTint` and one thin outline inside it
  (`outline.selection`, 1 px, the accent at `emphasis.selectionOutline`, 0.6), the tints first, then the outlines;
  squircle ends where a range really ends, square where the view or a piece cuts it (the frame rule) (the `selection`
  pass). Where it has no part, `projected.selection` is drawn in the same look at `emphasis.projected` times the
  opacity. `elements`: a selected item of zero length, which no Location holds, is a 2 px accent hairline over its lane
  (the ranges layer, EC-127).
- `hover`: each layer's own (its cell, at least the device pixel it is in so a far zoom never drops it; its item's
  marks), an outline `emphasis.hoverOffset` (1 px) outside the mark, concentric, so it never touches the selection's
  outline under it, in the hover ink (`seq-hover`: dark on light, light on dark); `projected.hover` is drawn the same
  way at `emphasis.projected` on every section but `projected.hoverSection` (the hover's own section's address), so a
  hover from elsewhere (a structure, no layer) is drawn on every section (the `hover` pass).
- `cursor`: a 1 px outline of its cell (`outline.cursor`) over the band, in a neutral ink (`emphasis.cursorInk`,
  `msui-text`), never the accent (the `cursor` pass). Whether it is shown is the caller's (a viewer shows it once the
  keyboard moved it).
- **The tints' place** (`emphasis.selectionPlace`): `over` (the default: the highlights' tints over the data in the
  `highlight` pass, the selection's above them in the `selection` pass with its outline; a selection or a highlight
  redraws the emphasis alone) or `under` (the highlights' tints and then the selection's in the `background` pass, under
  the data and hidden by opaque cells, the outline still in the `selection` pass; an emphasis redraw whose selection,
  projected selection or highlights changed by identity draws `background` and `data` again, `Frame.base`, and a backend
  redraws its base: about a full draw of the data, against the emphasis alone). Either way the highlights are below the
  selection. Why `over` by default: a tint under the data is hidden by opaque cells, and a selection then costs a data
  redraw.
- What a shared mark function needs besides its objects (a fill, a stroke, a width, the corners, an offset) goes through
  a typed array, and a color at an opacity is made inline: the emphasis passes allocate nothing once warm
  (`frame.quick.ts`: selection and cursor; projected hover, selection and highlights).

### Shapes

- Continuous data (cells, value bars) is square; discrete marks (ranges, at every step, a merged run at `bars` too) are
  squircles, radius a third of the height, capped, square below the minimum, so a mark has the shape its emphasis
  outlines. A merged run is round only where its first and last items really end. An end cut by the view, a piece of the
  axis or a piece of the item's own image (an exon) is square; a mark cut at both ends is square. The squircle is one
  curve, `squirclePath` (each corner a cubic on |x|⁴ + |y|⁴ = 1), that every backend draws. Outlines are concentric.

### Layer types

`LayerType<Prep, Lay>`: `name`, `fallback?`, `prepare(layer, model)` (tier 1, cached per resolved layer),
`layout(prep, ctx, data, overlay, prev?)` (per viewport; `prev` is its last layout, to reuse buffers and keep lanes),
`draw(lay, pass, out, emphasis?)`, `hit(lay, at)` (the point in `at[0]`, `at[1]`: an array, so a hit boxes no number),
`legend?(lay, { title })` (its channels' legends), `domain?(lay)` (the y domain it drew). `LayoutContext`: the band's
box and the strip's transform (in the strip's frame), the section, the layer's key, the device pixel ratio, the config,
the tokens, the text port, whether the view settled, `want`, `report` and `theme` (the schemes, the measures, the scales
and their held domains, the reader's substitution), and per strip `strip`, `span` (the union of the strip's bands: what
a frame with `across` outlines), `place` (the section's placement, for a `section` layer), `labels` (`inside`, or a sink
a layout collects labels into: rings), `bandTop` and `rowScale` (a band that scrolls its own rows). Generic additions,
no layout named: `used?(lay)` (the height it used on the strip, 0 for nothing; the band's height when absent: text
collapses a band no layer used, and only the layer knows its lanes), `scope` (`strip`, the default, laid out a strip at
a time; `section`, once a section, drawn in page coordinates: what crosses strips), `axes` (2: a second axis of its own,
the y path in Emphasis).

- **Channels a type reads** ("Themes"): `letters` `color` (the alphabet's scheme by default) and `ink`; `ranges`
  `color`, `size` (a share of the lane, the mark shrunk about its middle, its lane and hit box kept: EC-137) and
  `opacity`, a dictionary-coded field colored once a code; `values` `color` (uniform, or a bar by its value; a line's
  color follows) and a `scale`; `secondary` `color` (its classes in the `seq-ss-*` tokens, keyed by the `secondary`
  scheme; a form reads `kind` as the class word); `variants` `color` (`seq-variant`) and `size` (the head's radius, on
  `sqrt` unless the form says; its hit box kept); `ruler` none (`seq-ruler`). A channel's issues are said in every
  layout that reads it. A layer draws the first type of its fallback chain the target has; none is
  `seq.render.layer.type`.

- `letters`: a section's own letters, or another system's through a path that pairs positions one to one (identity,
  offset, segments, pairwise, an alignment row: a letter its image's cell; a section position no letter lands on, a gap
  mark; a letter with no position, a caret at its boundary, `projectBoundary`, EC-29) or a residue to a codon (its fine
  places carried as runs, each run a part: a codon split by an intron is two parts, the letter drawn in each part with
  room, EC-117; `code`: `auto` three-letter codes where the narrowest part in view has room, `one`, `three`); any other
  path is `seq.render.layer.letters` with its reason (a nucleotide on a protein's axis is a third of a residue and an
  intron or UTR nucleotide has no place there; a many-to-one path gives a cell several letters). The density step
  follows the layer's own cell (a residue three positions wide). A residue its codon does not translate to is marked by
  a bar along the bottom of its cell in the conflict ink (its cell keeps the residue's fill: a lone conflict reads as a
  mark, never a bar of conflict at the fit view) once the section's letters are held (`codonPairs` with the map's
  genetic code, EC-16); nothing is said before. **Complement** (`complement`): `auto` draws a nucleotide as its
  complement where it is read against its strand on screen (the axis's direction times its run's is −1: a reversed DNA
  axis shows the other strand read 5′ to 3′), `always`, `never` (EC-233's rule, changed); a protein never; letters shown
  never change data (model, status and search read the stored ones; colors by letter read the shown one). A cell a
  position (a part) colored by class, letters on them, a tick a position at `dots`, one bar a piece in view at `bars`
  and `density` (presence needs no data). A hit through a path carries the data's element as `extent` (the residue on
  its protein) and its parts on the section's system as `hull` (its codon); its hover outlines each part.
- `ranges`: items in sticky lanes, carried onto the section along the layer's path (a domain on its gene: a piece an
  exon); labels inside marks 24 px or wider, drawn in the `labels` pass. At `bars` a lane's marks touching within a
  device pixel merge into a run; at `density`, or a density table, one heat strip, a device column a pixel. Items past
  the lanes collapse into a density row in the last lane, the rest drawn as items (`ranges.overflow: 'extra-lanes'`), or
  the whole layer becomes the strip (`'layer'`). In a band under a lane's height (a text line's underline) the strip
  uses only the light `ThinStripShare` (0.3) of the density ramp, so it reads as a band, not as a dark rule.
- `values`: a bar a position, or at `bars` and `density` a bar a device column over the column's extent from the
  baseline (`values.aggregate: 'minmax'`, so negative extremes show; `'max'` the maximum alone), from positions or bins.
  Missing (NaN) is a 2 px neutral mark at the baseline. The domain is the layer's `domain` param where given, else open:
  computed from the values in view when the view settles, held while it does not.

- `frame` (EC-100): an annotation's regions outlined across tracks of its section (`across`: the first and the last
  track's ids, the bands between them too; its own band when absent: `ctx.span`). A strip layer: on each strip a
  region's outline covers the span over the region's x, an end the strip cuts (a line break, the view's edge, a piece's)
  left open, its title (the item's label) in its own band where the region really starts; its outline and title hit as
  the region's item, over the bands they cover (the hit tries a layer with `across` before the band's own). A strip
  without its regions leaves its band unused (text collapses it). No fallback: an outline drawn as bars would read as
  data.
- `ruler`: the section's axis labelled by a numbering (`numbering`: a built-in or one the system declares; unknown is
  `seq.render.numbering` and `ordinal` is used). A tick at each label that is a multiple of a step of 1, 2 or 5 × 10ⁿ,
  the least whose ticks are `ruler.labelGap` px apart; a label that does not fit (between the data box's edges and the
  labels already drawn, 6 px apart at least, majors first) is skipped, never shrunk. A label with an insertion code
  (`100A`, `111.1`) is a minor tick from 4 px a position, labeled only where it fits. It labels the positions where they
  are: a reversed axis counts down, a collapsed one labels its pieces. Where a text section prints a group's end number
  (`group_labels` not `none`) the ruler draws the tick without its label: a position carries one label, the group's
  (only where the ruler's numbering is the one the group labels print). Over a piece of more than 4,096 positions, or on
  a built-in numbering, the labels are looked up by number (a tick a step), not a position at a time. A hit is the
  position under the pointer. Fallback: none (omitted).
- `secondary`: the items of a `ranges` layer (the class from the `kind` field, else `category`: DSSP codes or words:
  helix, strand, UniProt's `Helix`, `Beta strand` and `Turn`, anything else a turn), in one band: a helix a block, a
  strand a block with an arrowhead on its C-terminal end, a turn a thin bar. **The arrowhead** is on the piece holding
  the C-terminal end of the item, found among all its rows (an inline annotation's; a tiled source's rows may lie
  outside what was read, so its items have no head: a head on the wrong piece is worse than none) (the end of its last
  row in reading order: across the origin of a cyclic system, the row that starts at it), at the screen side the
  residues run to (the axis's direction against the layer's path's), and not where the view clips that end (EC-130:
  reversed, across the origin, through a minus-strand codon map). An end cut by the view, a piece of the axis or of the
  item's image is square. The class is the whole code or word (DSSP `H G I`, `E B`, `T S`; mmCIF `HELX_*`, `STRN`,
  `TURN_*`; helix, strand, turn). **Items in a hidden stretch** are the break's tick, as `ranges`'. **Copies** (rows of
  one extent and class) are drawn once, the first hit (EC-67). From the coarse level the source's density bins are one
  strip. Fallback: `ranges`.
- `variants`: the items of a `ranges` layer as sites: a stem from the track's floor and a round head (radius 4 px); the
  heads of one site stand on one stem, spaced 9 px or squeezed to fit; a head and its stem are hit as the item. At the
  `bars` and `density` steps and from density bins a strip counts the variants a device column in the density tokens.
  Fallback: `ranges`.
- `values` with `encoding` `line` or `area` (EC-57): the bars' marks as a path through the cells' centers, held flat
  over the half cell at its ends and at a gap (a hole in the data, a break of the axis), broken by a missing value
  (which keeps its mark at the baseline); a color that changes does so halfway between two positions, the line going on
  to the boundary. `area` fills under each run at 0.35 opacity. At the bars and density steps the line follows each
  column's maximum and the column's extent is drawn under it.
- `ranges` with `line`: each item a line in its lane's middle instead of a block (a box and a highlight are blocks), the
  dash and the wave run from the content's start (the item's first position's edge, on its right on a reversed axis, in
  each piece), so a pan moves them with the data; below a position a pixel a run is a plain line. A value that is none
  of the styles is drawn solid and said once a layout (`seq.render.line-style`). With `arrow`, a block is an arrow, a
  line gets a triangle, at the end (the start) of the item on the piece that holds it, by the rule of the strand's
  arrowhead, never on a cut end (nor where the item's rows are not all known); a line's label has the page's color under
  it, so the line stops at it; the hit boxes are the marks' lanes all the same.

- **`bands`** (`render/bands.ts`, fallback `ranges`): a chromosome's stained bands in one band, an ideogram: each item a
  block colored by its `stain` through the `giemsa` scheme (grays from the text ink; `acen` in `seq-centromere`), the
  centromere's two `acen` bands as triangles meeting at it (the p arm's tip at its end, the q arm's at its start, turned
  on a reversed axis), a hairline where no band is; its hit names the band (`state`: its label).

### Data

- A layer's data is its node's access plan (an inline node's own address as the store's source, else its source node),
  the region in view carried back along the layer's path, at the level that fits a device pixel's positions there
  (`fitLevel`; a values layer's, `values.coarser` pixels'); the view at `view` priority and a margin each side at
  `margin`.
- A tile not asked for, pending or evicted draws as a pending placeholder (EC-169), a failed one as a conflict tint; an
  unknown source (not attached) is pending. Two regions sharing a tile (a cyclic view wrapping a short system) read its
  entry once. The regions a layer reads (the view's, or their image through its path) are merged across gaps of up to 16
  positions (`mergeNear`): one read over a run of short pieces, not one a piece.
- **An item's key** is its `id`, else its `row` (a null `id` falls back to `row`), else its tile and index; never its
  extent. A key seen in an earlier tile is that range again and is not drawn twice; within one tile a repeated key is
  another piece of one item (a wrapped row).

### Section layouts

`render/layout/`: a layout (`SectionLayoutType`: `place`, `sameKey`, `nav`, `wraps?`, `chrome`, `hitChrome?`) cuts a
section's display units into **strips** and places each on the page with a **warp**. Why: a layer type then draws one
way, over a straight axis and a box, in every layout; a layout knows strips, bands and warps, never layer kinds.

- **A strip** is one `AxisTransform` (its pieces, its window, x across a width; groups when the layout has them; its
  system the strip's), a stack of **bands** (a track a band, its top and height in the strip's frame, 0 collapsed) and a
  warp to the page. `Strips` is a struct of arrays kept frame to frame: per strip its transform, block, warp, the
  translation (`dx`, `dy`), its top and height on the page, its first display unit (its lays are kept by it); per band
  of a strip (`strip × bands + band`) its top, height and what its layers used.
- **Warps**: a translation (tracks, text), applied in place to what the strip drew (`TranslateWarp`: a text line's lays
  are kept across a scroll, only its translation changes); polar (rings, `PolarWarp`): x to an angle, y to a radius;
  what the strip drew goes through a scratch pass onto the page: a rect an annular sector (its squircle corners warped
  with it), a path's lines cut every few px so they bend, a glyph run placed letter by letter (upright or along the
  ring), an image refused (its layout reports it), a clip applied in the strip's frame before warping (a clip by path
  would change the lists' form for one case). Nothing allocated once warm.
- **The placement** (`Placement`): the section's box, its whole content's height (text: every line, for scrolling and
  export), its strips; `at(x, y)` the inverse of the warp then of the cut, into a reused `PlacePoint` (strip, band, the
  point in the strip's frame, the kind, the display unit); `offsetAt(x, y)` the nearest unit (drags, zoom anchors);
  `place(offset)` where a unit is drawn (on the page: the middle of its strip's first band with height; `at` of it gives
  the unit back). **Kinds of point**: `data`, `label` (a ring's outside label), `gutter` (the label column, text's row
  labels and numbers), `window`, `gap` (between groups or lines: no position for a hover, the nearest for a drag),
  `empty`. Held by property tests: every position in view placed once, `at` its cell's center gives it back, `offsetAt`
  monotone in reading order.
- **The viewport record** (`Viewport`): `{ start, end, top?, scale?, bands?, y? }`, read by the layout: tracks `start`,
  `end` the window, `top` a vertical scroll, `bands` a band's own `bandTop` and `rowScale` (a grid's `y` its second
  axis). One record for every layout (the field's meaning is the layout's): the viewer's groups, follows and history
  stay layout-neutral.
- **Navigation** (`SectionNav`): `pan` (a scroll's px deltas: positive shows what is right and below), `zoom` (about a
  unit), `fit`, `reveal` (the least move that shows a range: the cursor, a linked selection), `focus` (optional: what a
  focus on a range gives; absent, the range itself), `step` (`left`, `right`, `up`, `down`, `home`, `end`; −1 for none),
  `clamp`, pure, over a `NavContext` (the total, the data box, cyclic, the scale's bounds, the placement), each reading
  the viewport it is given (two moves before a layout add up): what the viewer calls instead of doing arithmetic of its
  own. Tracks: the window by px, clamped at `width / maxScale` within the axis, a reveal panning to the nearer edge,
  `up` and `down` none (one line: the cursor crosses sections); with `cyclic` a window wraps, never clamped at the
  total. The viewer sets `cyclic` for a cyclic axis whose layout `wraps` (draws a window over the origin: a ring); a
  tracks leader's transform ends at its total, so the tracks layout does not (TODO: a leader's wrap as two runs).
- **Chrome**: a layout's own marks a pass (`chrome`: once a section, then once a band before its layers) and their hits
  (`hitChrome`, before the bands' layers): tracks its label column, text its gutters, rings their guides, labels and a
  window.
- **A follower reads its leader's viewport** as the leader's layout reads it (`followWindow`: a ring's window, past the
  total when it wraps the origin; text its lines in view), else the leader's transform's window (tracks); a window over
  a cyclic origin is carried as two runs. A layout with `followWindow` places over the whole axis.
- **Layers a strip or a section**: a strip layer is laid out a strip at a time
  (`PlaceContext.layStrip(placement, strip, labels?)`, the layout's label sink reaching `ctx.labels`), a `section` layer
  once a section after the strips, over the placement's box on the page with the placement as its reader, drawn after
  the strips outside their warps and hit before them. A band's own scroll and row zoom (`viewport.bands`, by track
  address) reach its layers as `ctx.bandTop` and `ctx.rowScale`, and the key compares them by value.
- **The lay cache is bounded by the strips in view**: a layer's lay is kept per strip by its first display unit (tracks
  without blocks: by the layer's key alone), and a strip that leaves the view lets its lays go.

### Text layout

`layout/text.ts` (`TextLayout`, `TextNav`, `textGeometry`):

- **Lines of whole groups** (EC-94): positions in groups of `group` (10; 0 none) with `text.groupGap` (half a cell)
  between them; a line holds as many whole groups as fit the width (less the label column and the number gutters, at
  most `max_width`); the last group goes on the last line, short. Breaking is arithmetic on cell widths and gaps (an
  inline group label widens the gap it sits in, so it takes part in breaking). Groups are the axis's (`groups`: a cell
  skips the gaps, a span bridges them), counted from the section's first display unit, so lines stay aligned across
  piece edges; a piece edge does not start a line (`text.pieceBreak: 'none'`), or starts one (`line`: each piece cut
  into lines of its own, the lines' first units kept in the geometry, `lineStart`, `lineEnd`, `lineOf`). Groups are
  bridged, fixed: an item or a selection is cut by a line break only (no second use for cutting at every group). A block
  (with blocks) is to start a line, a one-position block to flow as a chip: the line loop is where they go.
- **A strip a line in view**, translated down the page; its bands the section's tracks (in track order, their heights
  their own, else `text.heights`: thin, a ranges track 4 px reads as underlines, a values track 3 px as a strip) under a
  row of group labels (`group_labels: 'above'`), `text.bandGap` between bands, `text.lineGap` after a line.
- **Collapse** (`text.collapseEmpty`, EC-96): a band no layer used on a line (`used`) is 0 high there and the bands
  below move up (laid out again there). The lines are laid out lazily, only those in view; each line's height is kept
  per section, full until it is laid out (the text's `end` and a follower's window are estimates: a tension), and the
  anchor holds the reader's place as heights land.
- **Its viewport**: `start` the anchor (the display unit at the box's top, fractional within its line), `scale` px a
  position (the letter size, `text.scale` 12 by default, never below the letters' threshold, EC-116); `top` derived;
  `end` derived (what a follower reads: the lines in view at full heights, `followWindow`). The layout places over the
  whole axis, its viewport its own. Its key (`sameKey`): the anchor, the scale, top; a reflow (the width) is the view's
  part of the key, and starts the kept heights again.
- **Numbers** (`numbers`: `start`, `end`, `both`, `none`): a line's first (left gutter) or last (right gutter) position,
  in the numbering of the first of its layers that names one (`numbering`, as a ruler's), else 1-based on the section's
  system; a ruler band's line numbers sit on the ruler's label baseline (one row of numbers across the line), else at
  the band's middle; group labels (`above`: the group's last position over its last cell; `inline`: in the gap after it;
  `none`). The first line's bands carry their tracks' labels in the label column (faded past it), each at its band's
  middle but never closer than a font size to the label above (thin bands would stack them: the labels flow down, the
  next line's label column is empty). A line's gutter (a band's label, its numbers) hits as the band's row
  (`target: 'row'`, the track's address as its layer).
- **The box**: its `height`, else its content's up to `text.maxHeight` (600 px; a host sizing its canvas to the content
  reads it), past which it scrolls inside, its strips clipped to it; an export lays out every line (`figure`) and its
  box is its content.
- **A bound**: at most `text.maxLines` lines (100,000) are laid out; more is `seq.render.text.lines` (give a `range`).
  Why: a chromosome as text is no figure, and the line heights are kept for the whole section.
- **Navigation** (`TextNav`): a pan scrolls lines (vertical: px to the anchor through the line heights); a zoom steps
  the letter size about a unit, its line kept as many lines below the anchor (the breaking's own arithmetic at the new
  size); `fit` the preset's size from the first line; a reveal scrolls a range's lines in by the least move (above: its
  first line to the top; below: its last line's bottom to the box's); a focus brings its first line to the top, the size
  kept; `up` and `down` step a line, the same column (a shorter line's last unit), none past the first or the last; the
  anchor is clamped where the last line meets the box's bottom, the scale to the letters' threshold and `maxScale`.
- Not built yet: chips (one-position blocks: with blocks).

### Window

`render/window.ts`: the window mark (a ring's window over a linear detail, a navigator's, a ruler rail's), drawn in the
`window` pass: a ring's by its layout's chrome in a strip's frame through its warp, a navigator's and a rail's on the
page. One geometry for all, on `navigator.window` (the navigators design: one geometry, extended): its floor
(`floorExtent`: at least `minPx`, at most half the span, `centred` on the truth or a `thumb`, clamped inside the span; a
ring's frame has no span: it only grows about its middle), its parts (`navigatorPart`; a ring's `windowPart` the same
with no span, no `outside`): a grip as a tab just outside each end, hit over half `hitPx` there (half `touchHitPx` by
touch, at least `handlePx`) and at most a quarter of the window inside it, so a press in the middle pans; the body hit
over at least `hitPx`, centred and inside the span; `outside` the rest of the span. The `window` target with its `part`.
Its drawing (`drawNavigatorWindow`, a ring's `drawWindowArgs`): the body in `seq-window` (tinted, outlined), the grips
in `seq-grip` (`handlePx`; lit per part: an end's grip, the body both), a hairline at the true extent where the floor
drew it wider (at least a device pixel), the wash (`seq-wash`, the background at 0.2, at most 0.3: a test) over the span
outside it (a ring: none; why light: an overview's window of a few floored rows of 10,000 washes nearly the whole
picture, whose shape must still read; the window stands out by its edge and grips); idle (nothing can move) a quiet
outline, not a target; an export's an outline with its wash, no grips. Its look (fill opacities, a grip's share of the
thickness) is fixed. Its numbers reach the pass through arrays (`NavWindow`, `args`, `NavWindowFlags` read once as
numbers): a window move allocates nothing once warm.

### Navigators

`render/navigator.ts` (`NavigatorRenderer`, `navigatorBox`, `railSections`, `rowsRailTracks`, `windowText`) and
`render/rail.ts`:

- **The picture**: a navigator lays its system whole out (the fit window, in the bound section's display order;
  `direction: genomic` up whatever the section's) as a view of its own over the spec's data, resolved like any
  (`resolveViewNode`), through the tracks layout with its own renderer and lists (no inset, label column or direction
  mark; the view's image cache, `images`): `axis` no track, `tracks` the section's tracks (`tracks`, else all but
  rulers) at `navigator.trackScale` of their height and lanes (at least `trackMinPx`) with no labels, `bands` one
  `bands` track filling it, `rows` the section's first rows track with no header filling what the other tracks leave (an
  overview's whole column), its rows fitted to the band (`rowScale` the band's height over the rows shown, the view's
  folds); a picture that fills is made again for another height. The same system with no `system` param: the bound
  section's own axis. Colors are the view's: the picture reads the view's renderer's held domains (`domains`, its
  `layoutState`, a picture key mapped to the view's), never a domain of its own over the whole axis, so a value has one
  color in both. Why the layout interface: a miniature drawn by other code drifts from the view it summarizes.
- **What no summary serves is left out, said once**: a lazy alignment's rows (`seq.render.navigator.rows`: left out), a
  ranges layer on a source with no coarse level whose whole axis passes the store's tile bound (the store's settings,
  `storeSettings`; `seq.render.navigator.layer`); a ranges layer reads its coarse level when the source serves one. Its
  wants are `view` priority only (an export waits for them).
- **The windows**: each bound section's viewport as the view's last layout placed it, carried as positions along the
  path to the navigator's system and placed on the picture (a follower's is its own shown range, so a leader's and its
  followers' windows are drawn together; a scroll group's sections share one), in the picture's display order (a section
  whose order the path turns against the picture's: its start is the picture's right end); over a cyclic origin two
  pieces, a grip at each true end (where a piece meets the axis's end it is body); what the section cannot reach (a
  region within its chromosome) washed, fixed; the floor within the picture, not that span (a region's window on its
  chromosome is seen, its hairline the truth). Rows: the band's rows in view on the picture's rows. A window's ends back
  in its section's own display units (`sectionUnitAt`) through the inverse path; past the picture's ends at the rate the
  path carries the edge (a protein position three nucleotides). `moves` says per section which axes can move (0 idle:
  quiet); absent, an export's outlines.
- **Draw lists**: the base (`background`, `data`) is the picture's, keyed as any frame; the `window` pass draws the
  washes, the windows (the one the pointer last acted on on top) and a sweep's outline; then the marks: the selection's
  and the highlights' ticks at the picture's foot (`navigator.marks`), the keyboard cursor's hairline, the pointer's
  hairline. A navigator never draws in a section's lists.
- **Hits** (`hit`): the window part under a point, its section and axis, the position on the navigator's system, the
  display row, a band's name (the `bands` layer's hit `state`). A navigator's hover is its own: a readout.
- **The ruler rail** (`Rail`, `FrameInput.rails` and `railsActive`, by section address, `key`; the viewer's): drawn on
  the first shown ruler's bottom border, else in a slim row of its own under the section's tracks (`SlimRailPx`, 12,
  fixed: the form's size, the least a rail and its grips need; a hidden ruler too; under, not above: the tracks keep
  their place); the whole axis across the data box, the window a segment with a grip each end (over a cyclic origin two
  pieces, a grip at each true end, the axis's ends body), collapsed breaks a notch; its grab zone the ruler's row (on an
  alignment's columns its lower `grabShare`: above it a drag selects columns), the slim row whole; grown (`grownPx`)
  while hovered or dragged (`Emphasis.drag`); idle (not active) a quiet hairline with no grips and no target (at the fit
  it reads as the line); in a figure (an export) a plain segment, no grips, no idle state. Its hit is the `window`
  target with `navigator: 'rail'` (touch the wider hit: `hitView`'s `touch`), its unit `railUnitAt`. Rails are part of a
  section's key by identity (a slim row takes room); which can move is not: the renderer sets it on the rails each
  frame, an emphasis redraw too (an idle flip redraws the window pass alone). `railSections(view, hosted)`: every tracks
  section no shown navigator binds along its positions: a declared navigator and a host's replace its rail, a hidden one
  brings it back.
- **The rows rail** (`placeRowsRail`, `SectionLayout.rowsRails`; `FrameInput.rails` by the rows track's address): the
  same `Rail` upright (`vertical`: along the page's y, its line and grab zone across on x; one geometry, drawn and hit
  by the same code, the coordinates turned at the rect and at the hit), `SlimRailPx` wide at the right edge (right of
  the last numbers' gutter: the gutter widened by it while any rows rail is on, so the data keeps its width idle or
  live), along the band's scrolled area below its pinned rows, in display rows below the pinned band; its window the
  rows in view. Idle while every row is shown (`whole`, the layout's: the renderer sets `active` only when the viewer
  says it can move and not `whole`). Its hit is the `window` target with `navigator: 'rows-rail'` (`RowsRail`), its
  `layer` the track, `position` the display row. `rowsRailTracks(view, hosted)`: each tracks section's tracks with a
  rows layer shown, unless a shown navigator binds the section along its rows (`rows` or `both`): an overview or a
  miniature replaces it, a hidden one brings it back. Why: the design's "simple slider navigator for vertical
  scrolling", the rail a linear section has, so a tall alignment can be scrolled and zoomed with no declared part.

### Rings layout

`layout/rings.ts` (`RingsLayout`, `RingPlacement`, `RingsNav`):

- **Rings outside in** (EC-108): the first track the outermost ring, each its track's height (scaled down together to
  keep a hole of at least `rings.hole` of the radius), round the section's whole axis; the diameter fits the box's width
  less the labels' room (`height` when the section has one), capped at `rings.maxDiameter`; the origin at `rotation`
  (degrees clockwise from the top, a spec param, 0), display units clockwise (a reversed section's positions run
  counterclockwise).
- **A strip a block, one warp**: the unrolled frame is the outer circumference wide (so a layer reduces per angular
  column at the outer ring; inner rings are finer than their pixels, never coarser), x the arc length from the origin,
  so every strip shares one `PolarWarp`. Blocks are arcs with `rings.blockGap` between them; a whole cyclic section is
  the full circle with no gap, its seam no break. A collapsed piece edge is zero width, its hairline marking it.
- **Labels outside with leaders**: the strips' layers give their items' labels to the layout's sink (`ctx.labels`: text,
  extent on the strip's system, element, the item's hit, order); placed when the ring is laid out (never on a window
  move or a hover): each at its item's middle angle (an item over the origin, its rows meeting there, one arc), per side
  (right: angles in [0, π)) sorted by height and spread a line apart within the box, hugging a circle
  `labelGap + leaderBend` past the rings, its leader out along its angle to the bend then to the label; past what a side
  holds, dropped by order (a block's label first, then the bands outside in) then the shorter item. A placed label hits
  as its item (`hitChrome`, the hit the sink was given: never a second hit at an anchor). `labels: 'none'` gives no
  sink.
- **The hole holds the title** (`rings.center: 'title'`): the system's name and its length (`bp`, `aa`; a section of
  blocks its id and its units).
- **The window** (`window` on the layout type, set on a kept placement too): the section's viewport (a ring that
  follows, its leader's range on its own axis), drawn in the `window` pass through the warp: a faint sector over the
  rings, and on a 12 px rim outside the outer ring (the design's least hit) its body and grips (`render/window.ts`), hit
  there as the `window` target with its `part`, so the data under the window stays hit; a grip is hit across the origin
  from either side. A narrow window grows about its middle, its sector and its rim alike. Shown with
  `rings.window: 'followed'` only while a section follows the ring (or the ring follows), `always` idle otherwise; a
  window over the whole ring (nowhere to go: the fit view) is idle too, quiet with no wash, grips or hover; an export
  draws it only with `windows` (`export.window`, `renderSvg`'s `window`).
- **Out of the key**: `sameKey` is always true (rotation and size are the section's and the view's), so a pan of a
  follower lays out the follower alone and draws the ring from its kept lays; its labels are not placed again.
- **Its viewport** is the window: `RingsNav` pans it by px of the box's width, zooms about a unit, reveals a range round
  the origin (unchanged when the range is in the window across the origin), steps wrapping, clamps by `clampWindow` (a
  cyclic window wraps, never clamped at the total; one as wide as the axis is the fit, from 0: a whole turn from
  elsewhere would put the origin mid-window).
- **An image is refused** on a ring (the warp counts it; `seq.render.ring.image`, said once a section by the frame).
- `at` is `atan2` and the radius through the warp's inverse (a label before the rings; the rim `window` while active;
  the hole and outside `empty`; between blocks `gap`); `offsetAt` the nearest block's edge in a gap; both monotone in
  angle from the origin. Held by property tests over rotations.

### Blocks

- **A block** (`block` children of a section: `of`, `range`, `label`; the resolved view's `blocks`): a piece of a
  section on a system of its own; a strip a block in every layout, its transform on the block's system, its layers the
  block's (a path a block), a hit naming its `block`.
- **Tracks** (`layout/blocks.ts`, `BlockPlacement`): `arrange: 'row'` (the default) lays the blocks end to end along one
  viewport with a gap between two (at most 12 px, the gaps inside the window sharing at most a quarter of the width,
  however many: EC-102); `stack` a block a row of the section's tracks, aligned at their first position, one zoom
  (EC-102, EC-115). Each block's label in a row above its tracks, faded at its strip's width. Lays are kept by the block
  and its first unit, so a pan keeps them (by the block too: a stack's blocks all start at 0).
- **Rings**: a block an arc (above).
- **Follow across systems** is refused: a section of blocks on more than one system neither follows nor is followed
  (`seq.render.follow.blocks`, it shows its own range). Why: which block's map a follow would carry is undecided. Blocks
  on one system follow as one axis.

### Tracks layout

- **A cyclic axis's window past its total** (a linear view panned across its origin, `wraps`): the transform takes the
  whole system twice, so the window's units go on over the origin (two runs, no break between); a follower of it wraps
  as of a ring's window.
- **On the interface** (`layout/tracks.ts`): one strip a section (a block in view with blocks: Blocks), translated by
  nothing (its frame the page), its bands the tracks stacked with `trackGap`; its chrome the label column (each track's
  label, the filter chip, the direction mark). Its draw lists are the same before and after the restatement but for the
  `window` pass's header: the committed records are its guard.
- Sections stacked; each an `AxisTransform`: its pieces (`[start, end, …]` in display order) end to end as display units
  in its direction (a reversed piece `[s, e)` shows `e − 1 − k` at its offset `k`, EC-32), the viewport a window of them
  (`{ start, end }`, all by default; `viewports` by section address before it) across the data box. Both ways: a
  position's cell (`cell`, its x extents left to right), a span's x extents per piece (`spanIn`, `spans`) with the cut
  flags on the screen's sides (on a reversed piece a range's start is on the right), the position under x
  (`positionAt`), the positions in view (`visible`), the breaks in view (`breaks`). Lookups are binary searches over the
  pieces, and a layer's loop over the pieces a range touches is `piecesIn`'s (the pieces in view it may touch): a pan
  over a reversed axis of 150,000 pieces (EC-33) never scans them. A span that only touches the view's edge is not in
  view (a boundary on the edge is). Marks mirror with the axis: an item's real ends are round wherever they land, a
  piece's cut ends square (the frame rule).
- **Leaders first**: axes are made in the order of `follow`, so a follower placed above its leader follows it; a leader
  not in the view is `seq.render.follow.leader`, a cycle `seq.render.follow.cycle`, and the section shows its own range.
- **A follower follows continuously**: its window is the images of its leader's window edges as offsets along its own
  pieces in its direction (`followWindow`, the viewer's too: one implementation); a window over a cyclic follower's
  origin is drawn over the arc's pieces (`followRange`, `followArc`).
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

### The grid layout

`layout/grid.ts` (`GridLayout`, `GridNav`, `GridPlacement`): a section of two axes (`layout: grid`).

- **One strip, the plot**, whose transform is x and which carries a second, `strip.y` (`Strips.ys`, an `AxisTransform`
  down the strip's frame, absent in every other layout), read by a layer as `ctx.axisY`; a layout with `axesY` lays its
  `section` layers per strip (they read `ctx.axisY`). Why a strip property: the layer stays a straight-axis layer; only
  the grid has two axes.
- **Every band covers the plot**: the tracks overlay, later over earlier (`fillBands`; `at` gives the last band holding
  the point, the topmost: the placement's `overlay`). A missing value is transparent, so two grids on one pair of
  systems merge later over earlier (EC-59).
- **y** (`ResolvedSection.y`): following x on the same system takes x's pieces, window and direction; following a
  section through a path maps that section's window (as `follow`); else the viewport's own `y` (`{ start, end }`), read
  by `GridNav` (`pan`'s dy, `zoom` and `fit` on both axes) and by `sameKey`. With `follow_y: x` the field is absent.
- **Size** (EC-109): x is the data box (the view's x, aligned with the sections it shares x with); the height makes a
  cell `aspect` (1) times as wide as tall over y's window, capped by the room, the larger of `heights.grid` (a floor,
  decided) and nine tenths of the view's height below the section's top: cells keep their aspect only where the room
  allows. At least 40 px. `aspect: 'none'` takes the tallest track's height. Why the floor: a view sized to its content
  would otherwise feed its height back into the room and shrink.
- **Chrome**: the label column (the first track's label) and y's first and last positions at the plot's top and bottom,
  their strings made at place (a draw allocates nothing).

### Lanes and Hits

- **Sticky lanes** (`StickyLanes`): an item keeps its lane while it is in view; a new item takes the first lane free
  over its extent; repacked only when the view settled and a fresh packing needs fewer lanes. Why: items must not jump
  when a neighboring region arrives.
- **Hits go through the placement**: the section under the point, its layout's own marks (`hitChrome`), then the point
  placed (`at`, reading the page point from an array as a layer's `hit` does, and the strip-frame x through the axis's
  `inGapOf` and `unitInto`: a number passed to or returned from a call that is not inlined is boxed, so placing a point
  allocates nothing on any axis, grouped too, and a grid cell's `hitView` nothing past its caller's call; a text
  section's `hitView` still allocates in its other hit paths, the package's TODO) and the band's layers, topmost first,
  each `hit` with the point in the strip's frame; a layer's hit box is in its strip's frame, so every drawn mark keeps
  its hit box through the warp. `system` is the strip's; `block` its index in the section when the section has blocks;
  `target` `window` (with its `part`) or `row` (a gutter's row) for a layout's marks.
- **Hits** are `{ layer, section, element?, extent?, hull?, position?, system, ySystem?, yPosition?, target? }` (flat: a
  hit of two axes, a rows cell or a grid cell, is the target `cell` with `ySystem` and `yPosition` set, so it allocates
  nothing; `row` a rows header, `tree-node` a tree's node and `summary-row` a folded clade's row, their element the
  clade): `layer` the layer's key (its address, else its track's and `:index`), `section` the section's id, `element` an
  item's `{ owner, key }` with its `extent` (a Location on the data's system), `position` on the section's system. Cells
  by arithmetic on the axis; items by `LaneIndex` (a lane's marks by left edge with a running maximum of right edges).
  In a track, the topmost layer with an element under the point wins, else the topmost with a position.
- **What a break hides is not dropped** (fixed): a `ranges` layer reads the stretches its breaks in view hide too
  (`dataRegions(src, axis, true)`), and its items wholly in one draw one tick on the break a lane in the break ink,
  whatever their number (`ticks`: x, lane, count, the hidden stretch), hit as the break (its `count`); a Location there
  (a projected hover, a selection, a highlight) lights the break. An end of a range on a piece's edge at a break is cut
  there (square). Why: drawing is checked against the data, an item in view with no mark is a lost item, and one tick a
  lane keeps the marks bounded by the screen.
- **A break** is hit within 3 CSS px of its notch above the band (and on a tick of the items it hides):
  `{ target: 'break', layer: <section address>, section, system, position, extent }`, `position` the section's position
  at its edge, `extent` the hidden stretch (`breakHit`); on data the cells and items beside it win (they touch; the
  hairline has no width of its own).
- **Every drawn mark has a hit box.** A `ranges` layer indexes each item's own marks (every lane, every piece), not the
  drawn ones: at `bars` a run merged from several items is hit by the item under the point, its gaps (under a device
  pixel) by the nearest within a device pixel, and the hover outlines that item's own marks inside the run. An item's
  hit carries its `hull` on the section's system (every piece carried along the layer's path): what a focus on it shows,
  where its extent is on the data's system (a domain on its gene, a transcript's domains from a nested source). Why: a
  run drawn as one mark but named by no item was a box with no hit and no hover.

### The pairs layer

`render/pairs.ts` (`PairsType`, `scope: 'section'`, fallback `ranges`): a `pairs` entry's items, each two ends.

- **Arcs**: each end placed with `ctx.place` (page x), the apex in its band (`side`: `above`, `below`), never flat; the
  height by span against the widest in view (`height: span`), so nested arcs nest and the widest fills the band; on a
  single translate strip an end past the window is taken `reach` (0.25) of the window past its edge, so an arc moves
  continuously on pan (EC-129). On a polar strip the apex goes toward the centre (`side: below` reaches the hole): a
  chord, its depth by span. An end not placed in the view, or the ends on different lines without `across_lines`, draws
  a caret at the placed end, never a dangling arc. A curve is one cubic (`curveArgs`), clipped to the data box.
- **Crossing pairs** (`crossing`): `same` draws them crossing; `split` the largest nested subset on `side`, the rest on
  the other (`splitCrossing`: exact up to `SplitExact`, 2,000 pairs, by an interval DP over the sorted ends; greedy past
  it, fixed: the exact subset is quadratic and no figure asks for another bound). On one system only.
- **Bounded by the screen**: past `pairs.maxArcs`, arcs binned by their ends' device columns, one a bin (widest first),
  the rest a density strip of ends, said once (`seq.render.pairs.limit`, info). Ends on two systems that are not the
  section's own report `seq.render.pairs.systems`.
- **Feet**: on the band's edge toward the sequence (`above`: its bottom, the arcs rising; `below`: its top), inset by
  the stroke so the whole arc is in its band; with `split`, on the band's middle, each side half the band. A view puts
  the pairs track on the letters' side (the tRNA's: over its letters).
- **Dots** (with `ctx.axisY`, in a grid): a rect a pair of single positions at every zoom, at least 4 CSS px, outlined
  in the text's ink (it reads over any cell color in both themes); a diagonal segment a stretch, with its partner when
  both ends are on one strand (`a_strand`, `b_strand`), against it otherwise; `mirror` (false: a pair is drawn where its
  data puts it; the tRNA builder sets it) adds (j, i), a self-pair once; `triangle` keeps a side. Pairs at one cell
  share it side by side (EC-63).
- **Themes and filter**: its `color` channel (`channelForms`, held extents, its legend), by default `field type` (a
  categorical palette fixed by value; the builder lists an inline table's kinds, as for ranges' `category`); a pair the
  channel gives no color takes the ranges' ink. The layer's data filter (`<pairs>/(not :of (item :type …))`) is read
  pair by pair at prepare (`id`, `type`, `label`, the strands and ends as fields): a pair left out is not drawn, split
  or hit, counted on its track's "n filtered" (a `section` layer's too); a form not read pair by pair is
  `seq.render.filter`.
- **Split order**: at one place the others' closing ends, each self-pair's two ends, then opening ends; equal intervals
  nest: a self-pair, a duplicate and a shared end cross nothing (held against a brute force).
- **Hit**: candidates by hull over the point's x, then distance to the curve (`pairs.hitPx`); the pair's element, its
  extent both ends. **Hover**: the hovered pair; as arcs, every pair with an end at a hovered position (a base, or
  either position of a cell) or in the hover projected onto the section; as dots, the pair at the hovered cell.

### The grid layer

`render/grid.ts` (`GridType`, `axes: 2`): a `grid` entry's values over the strip's x and `ctx.axisY`.

- **Only on its own systems**: x and y of the section must be the grid's (`seq.render.grid.systems`); no strip y is
  `seq.render.layer.axes`. Why: a grid is not projected; a generator computes it on the systems it is drawn on.
- **Drawn**: at `grid.rectMinPx` a cell and above, rects with outlines; below, `image` tiles (`ImageTile`, 256 × 256
  samples, never smoothed) from the shared 2D reduction (`render/pyramid.ts`: square levels, factor 2, `aggregate`'s
  reducer `mean`, `max` or `min`; a triangle's over its own cells alone, the others missing, so a coarse sample on the
  diagonal keeps its half and mixes in none of the other), one a resolved layer, at the level whose sample is at least a
  device pixel; tiles in the shared image LRU (`render/images.ts`) by layer, level, tile, the pyramid's and the colors'
  versions and the triangle; a tile's box on the device grid, so neighbours share a whole pixel edge (no seam). A frame
  not settled builds at most `grid.buildBudget` samples, drawing a held coarser tile meanwhile (else the first coarser
  within the budget, else the coarsest); a settled frame builds the rest. An axis reversed flips its side of the plot
  (EC-145).
- **Color**: the value's field theme (`value`), its domain the grid's range; a missing value (NaN) transparent;
  `triangle` keeps `upper` (y before x, above the diagonal) or `lower` (y after x).
- **Hit**: any point on the plot is the cell under it (`target: 'cell'`, `position` on x, `ySystem`, `yPosition`), drawn
  or not (outside a triangle too). Why: the crosshair reads the two positions everywhere on the plot.
- **Emphasis along y** (the crosshair, EC-145): for the hover, the selection, its projected part and every highlight
  channel, a column at each lit x (the strip's x painting) and a row at each lit y (the layer's own, through its y
  transform, across x's extent); the hovered cell outlined, its row and column tinted (`seq-hover`); a selected rect
  drawn as its rect, and on a section showing both its systems never its sides projected back (`rectOver`: the strip's x
  painting and the layer's rows skip them). Linked sections light both positions (the hover's location has a y part).
  Nothing allocated once warm: positions and colors go through typed scratch; a hit reads its point from `at`
  (`AxisTransform.positionOf`).

### Tensions (pairs and grids)

- **Images in SVG**: an `image` is runs of rects, so a zoomed-out large grid exports large; a PNG inside SVG needs an
  encoder in core. Rects in this slice; measured with the PNG export, which decides.
- **No navigator for a grid's own y window**: a grid whose y has its own window has no navigator over y; open until a
  Hi-C view needs one.
- **Arcs between sections** (ribbons, EC-209) stay deferred: `follow_y` covers a y following a section of the same view.
- **Emphasis along y** is the two-axis layers' own, through the y path (`axes: 2`).

### The rows layer

`render/rows.ts` (the layer) and `render/rows-model.ts` (what it computes without drawing: clades, display rows,
measures, consensus, the profile). An alignment's rows over its columns, its own y inside its band (`axes: 2`).

- **Data**: the alignment's `cells` on its rows system. Held when rows × columns is within `maxBlockValues` (2^24: the
  10,000 × 300 alignment is held): every tile wanted at level 1 and joined once (kept by the entries' versions); larger
  ones are lazy: the rows in view only, in data order (another order is `seq.render.rows.order`, info: ordering needs
  every row's field), a row not held labeled `row 5,124` and drawn pending; below a device pixel a row, a lazy source
  with no row bins draws a neutral fill (`seq.render.rows.reduced`, info; never a pending stripe that never resolves).
  Pins and hidden rows are by key, and a lazy row's key comes with its tile: on a lazy alignment neither applies
  (`seq.render.rows.lazy`, info).
- **Order** (`order`): `data`, `tree` (the tree's leaf order, rows no leaf names after them in data order), or a field
  of the rows table (`label`, `first`; missing last, ties in data order). Clades fold only in tree order.
- **Display rows** (`DisplayRows`): what each display row shows (a data row, or `-1 − clade` for a folded clade's
  summary row) and where each data row is shown; pinned rows (`pinned`, by label) first in a band of their own and kept
  in their place too (`rows.pinInPlace`, marked `·`; `placeRow` its place in the order: the tree's leaf, the hover and
  the selection are there too); hidden rows (`hidden_rows`) out, counted in a chip ("3 hidden") at the header's foot; a
  clade's "+n" counts its shown rows, one with none goes, one with one folds to that row. Data order with nothing
  pinned, hidden or folded is arithmetic (no arrays): a million lazy rows cost what ten do. Held by property tests over
  random trees and folds (every shown row placed once, found where it is drawn).
- **Clades** (`buildClades`): a tree's nodes over the rows, each a run of its leaf order; named by the node's label
  where it is unique in the tree (leaves counted) and not support (a number, or numbers joined by `/` or `|`, IQ-TREE's
  `87.3/95`: a reinferred tree changes them), else `first..last` (its first and last leaves' labels), so a name survives
  a reorder and two clades never share one (support values as inner labels, `((A,B)95,(C,D)95)`: a hover, a fold and
  `clades.keep` reach each); a `collapsed` name no clade has is said once (`seq.render.rows.clade`, info), never dropped
  silently; `allClades`, `cladesAtDepth`, `cladesWithin(height)` (every largest clade whose leaves lie within `height`
  of its node by branch length) and `cladesHolding(row)` are what the clade commands fold and open.
- **Folds** are the band's (`Viewport.bands[track].folds`, by name: session state, compared by identity and part of the
  section's key), else the layer's `collapsed` (the opening state).
- **The transform**: `rowScale` (px a row; `RowHeight`, 20, fixed, when the band gives none) and `bandTop` (the band's
  top in display rows below the pinned band, fractional, as a step's view declares it; the layout turns it into px,
  clamped so the last rows fill the band) from the band; the pinned band fixed at the top, the rest scrolled under it;
  display rows in view by arithmetic. The band's height: `height` a number of rows, else `fill` (the section's room left
  on the page, at least `heights.rows`, 240; a figure takes the default).
- **A kept window** (`keyAt`, `topOf`): the row at a band's top by its key (a held row's label, a lazy row's `#n`) and
  the fraction scrolled past it; back to a top under any folds and zoom, clamped: a folded row at its summary row, a
  hidden one at the nearest row shown in the order (`DisplayRows.nearestRow`). What the viewer's history keeps, and what
  keeps the window's top row across a fold. `rowsAxis()`: the rows system as an axis in display order (Shift+click
  extends along it).
- **Drawing**: from 3 device px a cell on both axes, a rect a cell (a run of one fill on a row merged at the cells
  step), letters from `density.letters` (10 px) on both axes, ink by the fill; a cell's fill is the `color` channel at
  its letter and column (the alphabet's scheme by default; the column measures registered for its themes), a gap none. A
  summary row's cells are its clade's consensus (`summary: consensus`, the majority letter opaque by its share; `bar`,
  the non-gap share), computed exactly from the cells when first drawn and kept by the display's version. Below 3 device
  px the shared reduction over the scrolled display rows' theme classes (a class a distinct fill; factor 4, row and
  column levels apart) at the level pair whose sample is at least a device pixel, drawn as image tiles (256 × 256
  samples, opaque by the majority's share; the base level holds none: whole) in the shared image cache, offset from the
  pinned band. A level's majority is the majority of the majorities below it (shares averaged), approximate within the
  measured tolerance, tested: the share within 0.05 on average and 0.3 at worst, the letter exact where the exact winner
  leads by more than 0.25. Why: a majority is not mergeable across levels, and below 3 device px a cell the shade reads,
  not the letter. The pinned rows are drawn on their own over a band that does not scroll (rects past `pinnedFrom`), the
  scrolled rows clipped below it. A letter-class reduction (`rows.reduce: letter`) is not built (the package's TODO).
- **The header**, in the label column: the tree (`tree_width`, 120) and the first residue numbers widen the view's label
  column by their width (the widest of its tracks sections, so every section's data starts at one x), the labels keep
  the column's own; the tree by branch lengths, else depth, a folded clade a wedge at its summary row, its geometry in
  display rows kept by the display's version (a pan makes nothing sized by the nodes); labels in the items' font (a
  row's label is data; a Pfam label, `NAME/1-147`, reads whole) from `rows.labelMinPx` (8 px a row), cut with an
  ellipsis past their room, a summary row "first +n"; the first residue numbers right-aligned before the data box and
  the last after it, in a right gutter the view keeps when a rows layer draws numbers (`numbers: gutters`): the numbers
  of the residues **in view** (the columns the cells are drawn for: a pan and a zoom change them; a gap is never
  counted; a row with no residue in view has none), the row's own at the fit view (kept per data, or per lazy row once).
  The right gutter is clipped to the band vertically, as the labels and cells are: a partly shown last row ends at the
  band's foot on every side.
- **Hits** (the header's through the band's gutter: a layer of two axes is asked there), written into one scratch hit (a
  hover allocates no hit object; a clade's element and extent made once a display): a cell is `cell` with its column
  (`position`) and row (`ySystem`, `yPosition`); a header's label is `row` (its `position` the data row on the rows
  system); a summary row `summary-row`, and `tree-node` the inner node whose branch segment is nearest within
  `TreeHitPx` (6, fixed: decided) by distance to the segment, a wedge included (a leaf's branch hits none), each its
  element `{ owner: <alignment>, key: <clade name> }` and its extent the clade's rows.
- **One rows layer a section** (`seq.render.rows.second`, warning: the cursor has one row order); a rows layer draws in
  tracks and text sections only (`seq.render.rows.layout`, warning, elsewhere).
- **In a text section** (the ESPript figure's lines repeat the rows): the band is every shown row (its height the
  display rows times the layout's `heights.rows`), every clade unfolded (`seq.render.rows.text.unfolded`, info, when
  some are folded), no tree, no scroll; a lazy alignment there is drawn pending. Why: a fold would differ per line.
- **Hover** (a clade): a hit on a tree node or a summary row (`tree-node`, `summary-row`: its element the clade) lights
  the clade: its branches and wedges in the hover ink, one band over each run of its rows in the order (a folded clade's
  summary row alone; a pinned row's band at both its places, the scrolled rows' bands never over the pinned band). Why:
  the element is hovered, here and on linked sections.
- **Emphasis along y**: the hover tints its row's band (a pinned row's at both its places) and its column within the
  band and outlines the cell; the selection's parts on the rows system as row bands, its rects of rows × columns as
  their cells, read over the display rows in view only (a binary search a row: a selection of any size costs the rows in
  view); a summary row draws the partial mark (the tint over the share of its rows selected, the outline when all are);
  the cursor with a row (`Emphasis.cursor.ySystem`) outlines its cell (the strip draws a cursor without one).

### The logo

`render/logo.ts`: a conservation logo over an alignment's columns, a fixed track (`heights.logo`, 48) above its rows.
Each column's letters stacked lowest first by `logo_scale` (`bits`: log2(k) − H, k 20 for protein and 4 for nucleotides,
less the small-sample correction `(k − 1) / (2 ln 2 n)` with `correction: small-sample`; or `probability`), the band's
height the full scale; each letter placed and stretched with the glyph runs' per-letter y and vertical scale (its cap
height, 0.72 of its size). Its profile (letter frequencies a column, `gaps: ignore` or `count`) is computed from the
held cells, joined once and shared with the rows layer on one store (`readHeld`), kept by their version; over a lazy
alignment it is a source's to serve, and none is read yet (`seq.render.logo.profile`, info). Below 2 device px a column,
one bar a device column, the tallest stack there (the shared reduction's one-axis max, factor 2), in its top letter's
color. Its color is the alphabet's scheme by default, each letter's (and bar's) color made a glyph (`glyphAgainst`): its
oklab chroma raised, hue kept, to the floor `seq-logo-chroma` (a number token, `SeqTokens.logoChroma`: 0.12 in light,
0.10 in dark, `--seq-logo-chroma` in CSS) clipped to the sRGB gamut at its lightness, then moved in lightness to 3:1
against the track's surface (`seq-track`) by the `ink` rule; a color above the floor keeps its chroma, a grey (chroma
under 0.02: a gap, `missing`, a neutral) stays grey and only moves in lightness. Why a floor and not a second strong
color a class: a scheme names one color a class, and a pale cell fill has little chroma, so lightness alone is muddy.
Once per color and token table (the table's version, which hashes the floor, keys the memo). Hit: the column, the letter
under the point as the hit's `state`. It reads every row (a logo of the rows shown is not built: the package's TODO).

### The key

- `Renderer.frame` keys a frame by the model's hash, the view, the size and device pixel ratio, the viewport, the token
  table's version, the text metrics' version, the overlay's version, whether the view settled, the config, and the
  versions of the entries (and the coverages) the last layout read. An equal key keeps the layout and the lists and
  redraws the emphasis passes alone; anything else lays out again into the same buffers. Same inputs, same lists.
- **Per section** (`LayoutState.sections`): a section whose layout, view parts (model, size, ratio, tokens, text,
  overlay, settled, config, types, themes, held domains), top, transform, viewport (through its layout's `sameKey`:
  tracks start, end, top and bands) and reads (the entries and coverages it read, at their versions) hold is not laid
  out again: its placement and lays are drawn again into the frame's lists (one list a pass for the frame: a backend and
  every reader stay as they are; drawing kept lays is a draw, not a layout).
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
  images as `<image>`s of a PNG data URL (`render/png.ts`: stored deflate blocks, synchronous, no DOM), fades as
  gradient masks, metadata as `data-*` attributes. Every text and attribute from a spec or data is escaped; a character
  XML forbids (a control but tab, newline and return, U+FFFE and U+FFFF, a lone surrogate) becomes U+FFFD, one for one,
  so a glyph run stays aligned with its x list. Ids take a prefix (`idPrefix`, name characters only), a metadata key is
  written only as a plain name. Held by hostile labels, ids, titles, fonts and metadata parsed back as XML
  (`_test/svg-escaping.test.ts`).

### Export

- **Navigators** (`navigators`, else the view's `export.navigators`, else `omit`): included, each shown navigator is
  stacked by `place` around the view (top and bottom across it, left and right along it, a miniature at its size) as a
  group with its windows' `<title>` (`Window: chr17:7,666,403–7,672,402 of chr17`) and `data-navigator` its id, its
  windows outlines and wash, no handles, hover or idle state; the rails drawn as plain windows (no grips, no idle
  state); the export opens the navigators' sources and waits for their wants, colors by the view's held domains.
  Omitted, no rail either (no slim row: a viewer's export is its picture without its chrome). Why omitted by default: a
  window is session, not a figure.
- **What the layouts place** (`figure`): a text section draws every line, its box its content; `size.height: 'auto'` is
  the sections' content's height (measured from a first layout); a ring's window only with `window` (the view's
  `export.window` when absent, false: an exported map is a figure).
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
- **A viewport in positions** (`viewport: { positions: { <section>: [start, end] } }`): each section's window read
  through its axis (`nearestOffsetOf`, its direction and pieces), so a caller (the render CLI) needs no knowledge of
  them; the frame's `viewports` by section address. Offsets (`{ start, end }`) as before.
- It answers `{ svg, settled, issues, lists }`; the SVG records whether it settled and the open domains it computed, an
  empty `<g data-section data-direction data-pieces data-y data-height>` a section (the draw lists are grouped by pass,
  not by section: these say what each axis is), and its ids take a prefix (`idPrefix`, by default `s` and a hash of the
  spec). On a shared store it reports only the issues of the sources it attached.
- **Folds** (`folds`, by a rows layer's track address): the clades the export folds (a viewer passes its own, as it
  passes the selection), else the layer's `collapsed`; a frame's `FrameInput.folds`, part of its key by identity.
- **A ring's window** (`window`, the view's `export.window` when absent, false): an export is a figure (`figure`), so a
  ring's window is drawn only when asked.
- **Legends** (`legend`: `none` by default, `below`, `right`; a view's `export.legend` when the option is absent): a
  band outside the data box with the view's merged legends; each shared scale's domain as drawn is recorded
  (`data-seq-scale-<id>`). `cvdSafe` draws the reader's substitution.

### Issues

- `seq.render.*` (`RenderIssues`): `rows.reduced` (info), `rows.clade` (info), `rows.text.unfolded` (info),
  `logo.profile` (info), `rows.order` (info), `rows.lazy` (info), `pairs.limit` (info), `pairs.systems`, `grid.systems`,
  `layer.axes`, `rows.layout`, `rows.second`, `view`, `color`, `theme.channel`, `theme.scheme`, `theme.measure`,
  `theme.color`, `scale.log` (info), `legend.unlisted` (info), `filter`, `line-style` (info), `numbering`, `layer.type`,
  `layer.letters`, `axis.pending` (info), `text.lines`, `follow.leader`, `follow.blocks`, `ring.image`, `follow.cycle`,
  `source`, `template.limit`, `unsettled`, `navigator.of`, `navigator.path`, `navigator.rows` (info), `navigator.layer`
  (info); a test checks the catalog.

## Stories

`@molstar/seqx-core/story`: a story's states, a transition planned from the diff of two states, sampled at any time.
DOM-free; Seq* runs no clock (the viewer's animation port is the application's).

- **A step's state** (`stateOf(spec, id, { cache, apply, config })`): a full state is its own tree; a patch applies over
  the spec given (the viewer passes its root with the preview), or over the state of the step it starts `from`. Cached
  by the spec's revision and the step (`createStoryCache`, core's `LRUCache`). Why over the spec by default: any step is
  reached directly and editing one step never changes another. A start missing or in a cycle is `seq.story.step.from` on
  the step and on each step after it, which applies over what it can reach. **A patch that no longer applies whole**
  applies op by op: the state is the ops that do, `applies` false, `skipped` the rest (`seq.story.step.stale`). `chain`
  counts the deltas to the step's keyframe; past `story.maxChain` it is `seq.story.step.chain` on the step. The walk
  back goes at most `maxChain` steps from the step asked (a step further back is read over the spec, `applies` false)
  and what such a walk computed is not cached. Why: a chain of any length is bounded work, and a step read alone is
  never given a state a capped walk made. `stepPatchOf(node)` reads a step's stored patch once a node.
- **The encoder** (`storeSteps(spec, steps, { store, maxChain })`, `withSteps`): `patch` a patch over the spec each,
  `full` full states, `smaller` the smallest of a patch over the spec, a full state and a delta from the step before
  while the chain is under `maxChain` (8): so reaching a step costs at most that many applies. `keep: full` stays a full
  state, `keep: spec` a patch over the spec. Sizes compare the stored forms' JSON. States compare without the story's
  own container (`stripStory`).
- **A plan** (`plan(from, to, options)`): plain data; a channel per change, each with its window (`t0`, `t1`,
  fractions), its easing and its slots in a frame. From each section's viewport both states name (`view`, field by
  field: `start` and `end`, `top`, `scale`, each band's `bandTop` and `rowScale`, a grid's `y`); from the diff of the
  roots in its order: a number set over a number (`number`; a section's `rotation` the short way round, its `labels`
  held `none` while it turns: placed once, at its end), a color over a color (`color`, in `story.color.space`: a hex
  literal, or a token through `resolveColor`, the viewer's token table, the end the token itself), an insert (`enter`)
  and a remove (`leave`): a node with a `height` grows from 0 or shrinks to it, one without comes at the window's start
  and goes at its end (what names it, a track naming entering data, never dangles); a move within its parent (`order`,
  at the cut); anything else (`cut`, at `story.cut.at`), the view shown too. Why a token without a resolver or a missing
  value is a cut: interpolation needs both values, never invented. The plan's `reducedMotion` (a viewer passes its
  `view.motion`'s answer) makes every window `[1, 1]`, each channel a cut at the end, and the plan's `duration` 0.
  `pending(section)` marks a view channel that waits on data: it holds its start and the frame counts it (`waiting`). A
  view's window **pans** (width geometric, centre linear: the content under the eye moves evenly) when the windows
  overlap or touch, else goes **out, across and in** (at the middle the window centred between the two and holding
  both); `cyclic` totals move the end's centre by a turn for the shortest arc. Sizes (`scale`, `rowScale`) move
  geometrically.
- **Sampling** (`sample(plan, t, out?)`): every channel at `t` eased, into `out` or the plan's own frame (kept by plan):
  deterministic and allocation-free warm. `sample(plan, 0)` is the start and `sample(plan, 1)` the end exactly (an angle
  ends its short way round, the same angle). Numbers, each color component, a pan's centre and width are monotone. The
  easings are this loop's own (inlined: a double returned from a call is boxed), one copy. Phases by channel
  (`ChannelPhase`: before, moving, done).
- **A frame's picture** (`frameRoot`, `frameView`): the start's root with each channel's value set (a cut once switched,
  an entering node at its grown height), a section's viewport its channel's fields over the start's (a field only the
  start holds kept to the end); the end's root itself once every channel is done, the start's before any began. Why the
  end's own root: `sample(plan, 1)` shows the `to` state exactly, and a frame never makes a revision.

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
