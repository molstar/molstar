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

- The language has one version, `LanguageVersion` (`0.5.0`). Why: one language, one description, one number for readers
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
  `columns`, `logo_scale`, `correction`, `gaps`, and a source's transform `records`, `letters`, `system`, `base`,
  `ends`, `types`, `attribution`, `licence`, with `fields` read, and the `numbering` kind under a sequence; and a
  track's `padding`, an addition; `0.5.0`: the kinds `molecule`, `compounds`, `compound`, `complex` and `unit`, a
  `numbering` under a unit, `pairs` under a complex or a molecule with the link columns `a_unit`, `b_unit`, `a_atom`,
  `b_atom`, `order`, `occupancy` and its `a` optional, a sequence's `modified`, the data references widened to molecules
  and units and a section's `of` to a complex, a block's `form` and its `of` optional (an override of a derived unit), a
  section's `classes`, a layer's `modified`, a `ranges` layer taking `pairs`, a source's `ref` and a values node's
  `field`; a map's `priority`, its methods `segments` and `structure` with `segments` and `status`, a section's `via`,
  and a numbering's methods `offset` and `segments` with `by` and `segments`; `pairs.a`/`b`, `ranges.target` and
  `map.to` naming an alignment, its columns). A spec in an earlier `0.x` reads as one in this one, and is reported where
  it no longer fits: a `0.4.0` spec reads unchanged (every addition is optional, and every kind a reference took it
  still takes).

## Language

`lang/`: `SeqLanguage`, the description of SeqViewSpec on `@molstar/core/repr` (kinds, parents, params, docs), and the
checks it cannot express.

- **Kinds**: `spec` (root: `version`, `title`); the containers `data`, `views`, `refs` and `snapshots` (reserved ids of
  the same name); under `data`: `source`, `sequence`, `ranges`, `values`, `map`, `alignment` (rows over columns: its
  `rows` inline, a `label` and aligned `letters` a row, every row as wide, or a `source`; its `reference` row by label;
  the columns system its address, the rows system `rowsSystemOf(address)`, `<address>/(rows)`: rows are a system,
  EC-35), `tree` (`of` an alignment, `newick` its text; leaves matched to rows when drawn), `pairs` (two locations that
  belong together: `a` a system, `b` another or `a`'s, `a` required under `data` (`seq.lang.link.target`); its `rows`
  inline, `a_start`, `a_end`, `b_start`, `b_end` each pair's two ends, never wrapping, at least a position each (EC-10),
  self-pairs allowed; an open `type` its kind, since `kind` is reserved: canonical, wobble, contact; `id`, `label`,
  `a_strand` and `b_strand` (EC-206); any other category column rides along as a field a data filter reads,
  `(item :tertiary yes)`), `grid` (a value per pair of positions: `x` a sequence, `y` another or `x`'s, `values` x-major
  float32, NaN missing, `|x| × |y|` of them: `seq.lang.grid.length`); under a `sequence`: `numbering` (Numberings,
  below); under `views`: `view`, then `section`, `block` (a piece of a section on a system of its own: `of` (absent: an
  override of a derived unit), `range`, `label`), `track` (`hidden`, default false: kept in the spec, its row collapsed;
  its layers keep their own; `padding`, px inside its well around its data, vertical only: a number top and bottom, or
  `{ top, bottom }` one or both, each non-negative, a `repr` issue otherwise; absent, the render config's
  `trackPadding`: Render, Tracks layout), `layer` (`hidden`, default false: kept in the spec and not drawn, so hiding is
  a one-key `set`). A `ranges` node's `map: 'codon'` makes each item (its rows by `id`) a coding sequence of its own
  (Overlay); a track's `for_each` makes it a template (Resolved view, Templates; its kind param is `node_kind`, since
  `kind` is reserved, and is `ranges` or `values`: what a layer draws). Reference params (`target`, `of`, `from`, `to`,
  `follow`, `source`, `data`) are `ref` params naming the kinds they take; a reference is an address.
- **The axis's params** (the writer's: what the figure is; never inferred from what a section draws, so adding a layer
  never turns a section over): a section's `direction` (`genomic`, `reversed`, `follow`, or
  `{ along: <sequence or map> }`; absent: `follow` for a section with `follow`, else `genomic`), `pieces` (`whole`,
  `coding`, a ref to a `ranges` node or a codon `map`, or a list of refs, their union; absent: `whole`, followers
  included), `collapsed` (the form `section.collapse` switches to; absent: `coding`) and `pad` (positions kept each side
  of a piece; absent: 0); a layer's `complement` (`auto`, `always`, `never`; absent: `auto`) and `code` (`auto`, `one`,
  `three`; absent: `auto`), read by `letters`; a layer's `numbering` (a numbering id; absent: the section's default,
  `ordinal`) and `units` (`none`, `bases`; absent: `none`), `unnumbered` (`mark`, `fallback`; absent: `mark`) and
  `fallbackNumbering` (a numbering id; absent: `ordinal`), read by `ruler`, `encoding` (`bars`, `line`, `area`; absent:
  `bars`) by `values`, `style` (`blocks`) by `secondary`, `line` (`solid`, `dashed`, `dotted`, `wavy`, `double`, `box`,
  `highlight`; absent: blocks) and `arrow` (`none`, `end`, `start`, `both`; absent: `none`) by `ranges`; a grid
  section's `of_y` (absent: `of`), `follow_y` (`x`, `none` or a section; absent: `x`) and `aspect` (a number or `none`;
  absent: 1); by `pairs` its `side` (`above`, `below`; absent: `above`), `reach` (absent: 0.25), `height` (`span`;
  absent: `span`), `crossing` (`same`, `split`; absent: `same`), `across_lines` (absent: false) and `mirror` (absent:
  false: a pair drawn where its data puts it, never read off whether the two axes are one system; the tRNA's builder
  sets it); by `grid` its `aggregate` (`mean`, `max`, `min`; absent: `mean`); by both in a grid `triangle` (`both`,
  `upper`, `lower`; absent: `both`). Why: the figure's form is the writer's, never the host's or the data's: one spec
  draws the same in the canvas, the SVG and another viewer. The refs inside `direction` and `pieces` are checked by kind
  as top-level refs are, and a change's check reaches them (its reverse index walks a union's record member). Why: the
  composition is the application's.
- **Numberings** (`<sequence>/<id>`, a `numbering` child of a sequence): a label a position of its sequence, its id the
  name the forms read (`:numbering kabat`, Queries, Numberings and labels). `method` is `labels`, `offset` (a position
  plus `by`: UniProt's from 1, an author numbering from 94) or `segments` (runs `start`, `end`, `first`: a position in
  no run unlabelled, an `insertion` code on a run of one, `52A`; runs apart, `seq.lang.numbering.segments`); for
  `labels`, `labels` inline, a category column one a position in index order, `null` no label there; or `source` and
  `field` (absent: `label`): a source's rows, a row a position (its `position`, or its `start`), on the sequence by its
  name (Render, Numberings from a source); `scheme` an open word, lowercase (`^[a-z][a-z0-9_-]*$`: `kabat`, `imgt`,
  `chothia`, `uniprot`, `auth`), that says what the labels follow and changes none of them; any such word is accepted (a
  lab's own scheme too). Labels need no order (IMGT, EC-7) and carry insertion codes (`100A`, `111.1`). Checks: one of
  `labels` and `source`, `field` with a source only (`seq.lang.numbering.content`, which also says a param another
  method takes, or `by` or `segments` missing); an id no built-in has (`seq.lang.numbering.id`: the forms read `index`
  and `ordinal` as theirs, so one declared so would never be read); inline labels one a position of inline letters
  (`seq.lang.numbering.length`, error); a label at two positions said once (`seq.lang.numbering.duplicate`, warning: it
  names the first, the queries' rule). A change to a sequence checks its numberings again (Checks of a change). Why a
  child of the sequence: a numbering is the sequence's, and moving or removing the sequence takes its numberings with
  it.
- **Maps** (`map`, `from` and `to` systems): `method` `offset` (`by`), `codon` (`exons`, `strand`, `table`), `segments`
  (a `segments` table: `from`, `to`, `length` a run, `from + i ↦ to + i`, or `to − i` with `strand` `-`; runs apart on
  `from`, each at least a position; many-to-one on `to` allowed: a fusion onto its part, EC-20) and `structure` (a
  sequence onto a unit or the entity's molecule, SIFTS in abstract form: `segments` with a `provenance` word a run,
  `sifts`, `alignment`, `author`; and a `status` table, `start`, `end` on `to` and `conflict` or `unobserved`: the rest
  of a run paired, what no run holds unpaired; its runs apart on `to` too: a residue of the structure is one of the
  sequence). A param another method takes, or one missing, is `seq.lang.map.params`; bad runs or status rows
  `seq.lang.map.segments`. Any map's `priority` (an integer, 0 when absent, not filled in) ranks it among maps joining
  the same systems (Resolve, Paths: priorities sorted weakest first, compared in turn). Why `unpaired` is not written:
  it is what no run holds, so it cannot disagree with the runs.
- **A section's `via`** (a map, or a list): for that section's paths (its follow, its `direction: { along }`, its
  layers' data, what is projected onto it) the named maps' rivals are left out (Resolve, Paths); refs checked as
  `pieces`' are (`seq.lang.ref.target` at the item). Links carry none: where they land is the sections'. Why per use: a
  figure comparing SIFTS with an alignment shows both, a section each.
- **A source named by a host** (`ref`): a host's name for data with no file (an edited or predicted model), read by a
  host's adapter only. Beside a `uri` it is portable (a reader without the host opens the `uri`); alone it is allowed
  and said once (`seq.lang.source.host_only`, info, EC-185); beside an inline `value` it is `seq.lang.source.content`
  (inline data is never a host's). Why: the same spec opens with the host and without it.
- **A values channel** (`values.field`): with a `source`, the column of its values read (`plddt`, `b_factor`); the
  source's default value when absent; without a source `seq.lang.values.content`. Why: one source, one open and one
  parse serve every channel.
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
- **Systems** a data reference takes (`ranges.target`, `values.target`, `pairs.a`/`b`, `grid.x`/`y`, `map.from`/`to`,
  `block.of`): a sequence, a molecule (its residues) or a unit; `section.of` a complex too (Molecules and complexes).
  `pairs.a`/`b`, `ranges.target` and `map.to` also take an alignment, meaning its columns system (EC-48's "an alignment
  its columns"): a family's consensus pairs (Stockholm's `SS_cons`) drawn over its rows, a chain mapped onto the columns
  through the seed row it is closest to (a `segments` map, a run a stretch of the row with no gap: a stopgap until the
  design's `alignment` map kind, the gapped row itself, is built; TODO.md). An inline alignment's width bounds their
  ends and ranges (`seq.lang.pairs.bounds`, `seq.lang.ranges.bounds`; a change to its rows checks them again). Why these
  four: they are what a family's annotation and a structure's link to it name; `map.from`, `values`, `grid` and
  `block.of` wait for a use (a map from the columns is the same map read back).
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
  the ref depth deep); a template's layers under a track whose params changed; a sequence's numberings when it changed;
  the sections and refs of a chain a changed one closes; a molecule's links (counted against its residues); a unit's
  numberings, what names it and its complex's links, also when it is reached through its molecule; and the links of a
  complex a unit left (removed or renamed: links name units by id, not by reference). What names what is a reverse index
  of the spec's references (a persistent set of referring nodes by the id run they name), built once a spec and derived
  from the last revision's along its changes. The same checks on the same nodes before the change are taken away, by
  code, node (through moves) and place in the node. Why: a patch is refused only for what it introduced (an error
  already there never blocks an edit beside it), and the cost follows the change, not the spec. Held by tests that
  compare it with the whole check's difference, patch by patch, along a chain of revisions and over seeded random
  patches of ref chains; and in the heavy lane, at 16,000 layers, the first check after a load and a fan-out to every
  node no worse than the whole check, a one-node edit a small part of it.
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
  block's values the data's (Data, Tables) and an archive's the archive's (Data, Archives); regular expressions in
  search wait for the feature that reads them (TODO).
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

### Molecules and complexes

`lang/molecule-language.ts` (the kinds), `lang/check-molecules.ts` (their checks). Seq\* holds pointers and labels,
never chemistry: a component code, a SMILES, a GlyTouCan id are strings it matches and shows.

- **A sequence's `modified`**: a sparse table (`position` in order, `component`, `length` 1 when absent) of the
  positions whose component is not their letter's standard one; a position not listed has its letter's standard
  component. The letter is the parent's (MSE `M`, PSU `U`), `X` or `N` where there is none. Several rows at one position
  are alternatives (microheterogeneity), the first giving the letter; a component of `length` > 1 stands for several
  parents (CRO for T, Y, G) over their letters. Checks: positions in order (`seq.lang.modified.order`), within the
  inline letters with their length (`seq.lang.modified.bounds`), a length of at least 1 stopping before the next row's
  position (`seq.lang.modified.length`). Why sparse: a protein of 300 residues has a handful of modified ones; why the
  parent's letter: search, alignment, schemes and translation read letters, and the component changes no arithmetic.
- **A `molecule`** (under `data`): what is not a linear polymer of letters (a ligand, an ion, a glycan, a drug). Its
  positions are its `residues` (a component code each, a category column) or its `source`'s
  (`seq.lang.molecule.content`: one of them); its `atoms` optional, beside inline residues (`residue` non-decreasing and
  a position of the residues, `seq.lang.molecule.atoms`; `name`, `element`; heavy atoms only, a hydrogen is
  `seq.lang.molecule.hydrogen`, EC-64); its `class` an open word (`branched`, `ligand`, `ion`, `water`; `ligand` when
  absent); its `name`; its `identifiers` (`namespace` an open lowercase word, `seq.lang.molecule.identifier`; `value`
  matched exactly). Two systems: `<molecule>` (its residues) and `<molecule>/(atoms)`. A molecule with no identifier is
  valid (an unknown ligand). Why a kind and not a sequence: a sequence's arithmetic (alphabets, schemes, codons,
  aligners) does not apply, and folding both into one kind makes every reader of letters check a flag; why no hydrogens:
  they are often absent or model-dependent, so lengths and per-atom data would not line up across entries.
- **Compounds** (`/data/compounds/<code>`, the `compounds` container's id is `compounds`): a `compound` describes a
  component, its id the code (`name`, `parents` a list of codes, `class` the CCD's `type` word, `identifiers`). A code
  needs no node. Why a container: a compound is reference data, not a system, and its codes collide with the writer's
  ids of systems.
- **A `complex`** (under `data`, `name`): units shown together, in order, with the links between them; it adds no
  arithmetic. **A `unit`** (under a complex): an instance of a sequence or a molecule (`of`), its `range` on it (never
  wrapping, within it: `seq.lang.unit.range`, EC-13; all of it when absent), its `observed` runs (`start`, `end`, in
  order within its positions: `seq.lang.unit.observed`), `chain` (`label_asym_id`), `auth_chain`, `operator`
  (MolViewSpec's `instance_id`), `label`, `group` (a BIRD/PRD instance); its `numbering` children (`auth`) label its own
  positions. A unit is a system whose position i is its molecule's `range[0] + i` (Queries, Numberings and labels; the
  paths' unit steps); its built-in `index` and `ordinal` count its own positions, a declared numbering of its molecule
  is shifted to it. Units of one complex with equal `of` and `range` are copies, derived, never declared (EC-65).
- **Links** are a `pairs` child of a complex (between its units: `a_unit`, `b_unit` the units' ids, each end on its
  unit's positions) or of a molecule (between its residues); a child names no `a` or `b`, its parent is both
  (`seq.lang.link.target`); unit columns only under a complex, each end naming a unit of it (`seq.lang.link.unit`) and
  lying within it (`seq.lang.link.bounds`), each end one position (`seq.lang.link.position`: a link joins an atom of
  each; a site or contact of several positions is a pairs node under `data`); `type` (`covalent`, `disulfide`,
  `glycosidic`, `metal`, any word), `a_atom`, `b_atom`, `order`, `occupancy` (0 to 1, NaN 1: `seq.lang.link.occupancy`).
  Both ends on one unit is allowed (a disulfide in a chain); glycoforms and stochastic conjugation are several rows with
  occupancies. Why a `pairs` node: one shape for links and contacts, drawn by the pairs layer; never a map, so a
  glycan's selection never selects its site.
- **Views of them**: `section.of` takes a complex (a block a unit), `section.classes` the classes shown in order
  (`polymer` a sequence's units, else a molecule's class; `[polymer, branched, ligand, ion]` when absent: water left
  out); `block.form` (`whole`, `residues`, `atoms`; `whole` for a molecule or a unit of one, `residues` for a sequence
  or a unit of one, when absent; `atoms` on a system without them is `seq.lang.molecule.form`, a unit read through its
  molecule, a molecule from a source not judged); a `letters` layer's `modified` (`mark`, `code`, `tint`, `none`; `mark`
  when absent); a `ranges` layer takes `pairs` data (each end a mark at its site). What they draw is Render's.
- **Shared with dynamic nodes**: `checkDataNode` runs a molecule's and a sequence's `modified` checks and a pairs node's
  link checks (`DataNodeTargets.parent`), so a host's molecule is held to the spec's rules; `checkUnit` and `checkLinks`
  take what they read as arguments, for the overlay's complexes.

| Item                   | Default                                                                             | Why the default; the alternatives                                                        |
| ---------------------- | ----------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------- |
| `block.form`           | `whole` for a molecule or a unit of one, `residues` for a sequence or a unit of one | a molecule reads at the fit view; one rule for every non-polymer; `residues`, `atoms`    |
| `section.classes`      | `[polymer, branched, ligand, ion]`                                                  | Mol\*'s and MVS's order (EC-101); waters say nothing on an axis; any list, `water` added |
| `copies`               | `merge`                                                                             | unchanged; `each` (poses: the builder's for a pose set)                                  |
| `modified` (`letters`) | `mark`                                                                              | a modified residue drawn as its parent alone hides it; `code`, `tint`, `none`            |
| `code` (`letters`)     | `auto`                                                                              | the component's code where three letters have room; `one`, `three`                       |
| `class` (`molecule`)   | `ligand`                                                                            | the commonest non-polymer; `branched`, `ion`, `water`, any word                          |
| identifier namespaces  | open words                                                                          | a host's registry ids are pointers too; the documented vocabulary                        |
| `ref` (`source`)       | absent                                                                              | portable by `uri`; a host's name for data with no file                                   |
| `field` (`values`)     | the source's default value                                                          | one source, many channels                                                                |
| hydrogens in `atoms`   | refused (fixed)                                                                     | lengths and per-atom data line up across entries (EC-64)                                 |

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
  does. A sequence's numberings are written in its body: `b.data.sequence(id, params, (s) => s.numbering(id, params))`.
- **A reference param takes** a handle (or builder) of a kind its description allows, or an address written out (any
  string); a handle of another kind is a type error. Which params are references is read from the description, in the
  types and when the builder converts a handle to its address. Why: no list of param names to keep in step with the
  language.
- **Molecules and complexes**: `b.data.molecule(id, params, (m) => m.links(id, params))`,
  `b.data.complex(id, params, (c) => { c.unit(id, params, (u) => u.numbering(id, params)); c.links(id, params) })`
  (links' params are a pairs node's but `a` and `b`), `b.data.compound(id, params)` (the `compounds` container made with
  the first); `modified(rows)` and `identifiers(record or pairs)` write the tables from rows.
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
- **A layer with no color takes the next categorical color** (a ranges layer whose spec gives no `color`): the theme's
  `seq-cat-<n>` by its order in the view (`LayoutContext.colorIndex`), counted by what it shows: one data node shown in
  two sections (a ring and its detail) one color, a template's instances one; a spec's color always wins, and a layer's
  items stay one color. Why: one steel blue for every colourless layer reads as a placeholder (the defaults audit).
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
  there is no gap), then a color once; before any is drawn, every letter. A view's legends (`viewLegends`) are its
  tracks' layers' and its sections' own layers' (a pairs layer's arcs and dots), titled by the track holding the layer.
- **Builder**: `theme.uniform|field|rules|scheme|letters|categorical|…` and `measure.*`; a layer's `field` form over an
  inline table's text column gets the column's values listed (the legend is the vocabulary); a `ranges` layer over a
  `category` column, a `variants` one over a `consequence` column and a `pairs` one over a `type` column default to a
  field form on it.
- **Checks** (`seq.lang.theme.*`): a form or palette that does not exist (`form`), a color that does not read (`color`),
  a channel a form cannot fill (`channel`), a field of the wrong type (`type`), stops out of order (`stops`), a scale
  that does not exist (`scale`), forms nested past `themeDepth` (4, `depth`).
- **Filters**: a theme never hides; a filter does. An item layer's `data` written as a query
  (`<annotation>/(except :of (item :<field> v))`) leaves out the rows its form rejects as the shared item reader reads
  them (`item` by `:id` or a field, a list param any of its values; `not`, `and`, `or`), compiled once a prepare
  (`compileRowFilter`); the track's header shows "n filtered", its `seq:filtered` state (Render, Track headers: a chip
  of `TrackLayout.header`, its `address` the filtered layer, hit `header-state`), laid out and drawn as every state. A
  form not read row by row filters nothing and is `seq.render.filter`. **One rule for the base and the filter**
  (`splitDataForms`, read by the layers and the viewer's `dataFilter` alike): a leading run of item forms (`and` and
  `or` of them too) is the track's base (`<annotation>/(item :category Domain)`; `ResolvedLayer.data.base`,
  `LayerSource.base`): its rejected rows are left out and not counted as filtered; a trailing `(not …)` is the user's
  (`data.filter`) and counted. Why: one read of a file serves a layer a feature type, and the chip counts what the user
  hid, not what the track is. A lone form of another kind stays the filter (and is `seq.render.filter` when no row reads
  it); any other mix of forms is `seq.render.filter` too (`data.unread`), every row drawn. Density tables take neither
  the base nor the wrap: they are counts of what a source served.

## Issues

`seq.lang.*` codes (`LangIssues`, and `BoundIssues`, `seq.lang.bound.*`, the bounds'), each with a severity, a
description and a hint; a test fails on a code a source uses and the catalog lacks, and on one cataloged and never used,
and a misuse for each code is in `_test/lang.test.ts`. The namespaces are `repr.*` (the tree core's), `seq.lang.*`,
`seq.query.*`, and the later `seq.resolve.*`, `seq.source.*`, `seq.store.*`, `seq.view.*`, `seq.render.*`,
`seq.address.*`. An issue's `path` is the step path and its `address` the node's, where it has one.

## Queries

`src/query/` (`@molstar/seqx-core/query`): locations, keys, the three result types, and the query forms as `repr` form
resolvers over a model of systems, numberings, annotations and structures; property access apart from the language
(locations-and-queries.md, The language; properties.md). Forms and addresses are one language: a form is a segment of an
address, `/data/p53/(range :from 94 :to 292)`.

### Result types

- **Nodes** (exact addresses of spec nodes and of nodes data makes), **Keys** (keys of what data makes: an annotation's
  items, a structure's groups; the type was _Elements_, renamed by the author, 2026-11-20) and **Location** (positions
  on systems). A form gives one of them; they travel as `repr`'s `Named` (`{ kind: 'nodes' }`,
  `{ kind: 'other', type: 'keys' | 'location', value }`). A `Named` whose value lacks its type's shape is no query
  result (`fromNamed` checks it). Why: a hover, a click and a layer each hold the type they act on, and no form guesses
  a conversion.
- **Conversions are forms**: `(extent :of <keys>)` gives a Location, `(items :overlapping <location>)` Keys. A form
  given another type is `seq.query.form.type`.
- **Pairs answer `items`**: a pairs annotation is queried as its items, a pair one item with a row per end, each row on
  its own system (`QueryItem.system`; a ranges row's is its annotation's `target`): `(items :overlapping <location>)`
  gives the pairs with an end in the location, on either system; `(item :type …)` reads a pair's fields (`type`,
  `label`, its other category columns); `extent` and a result's location cover both ends. Why: a host asks which pairs
  touch a range (the viewer's `linkedPairs`). A location part on neither target is read against the rows' systems before
  it is said unreachable.

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

### Keys

- `{ owner, key }` (`Key`; `query/keys.ts`): the owner's canonical address and the data's own key, a string or a number.
  A number never equals a string (`1` is not `"1"`), as in an address. Keys are sorted (owner, numbers before strings)
  and each once; compared and hashed as values.
- **Their extent** (`key-extent.ts`): an annotation's items by their ranges; a structure group's (`groups`) by its
  owner: a residues or atoms system's position, a complex's unit by its id, or its entity's units.

### Numberings and labels

- **Indices** are 0-based half-open inside; **labels** are inclusive names at the edge, converted once where a form is
  read.
- **Built-ins**: `index` (0-based) and `ordinal` (1-based), on every system; without `:numbering` a form reads
  `ordinal`. A declared numbering (a sequence's `numbering` child, Language, Numberings; `labelsNumbering`: an index of
  its labels built on first use, a label found by its spelling, else ignoring case, as the order of labels ignores it;
  no `snap`, so its positions are scanned) is looked up on the system, then on its molecule (`parent`, shifted by the
  unit's start, clipped to the unit): a unit's own beats its molecule's of the same id. Only a declared numbering is
  shifted: a unit's `index` and `ordinal` are its own (the ordinal 1 of a unit over `[2, 20)` is its molecule's ordinal
  3). Why: the author's numbering of a chain wins over the entity's. Not across a codon step: that waits for paths.
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

| Form                                                                        | Prefix                  | Result                                                                                                          |
| --------------------------------------------------------------------------- | ----------------------- | --------------------------------------------------------------------------------------------------------------- |
| `(range :from l :to l :numbering n)`                                        | systems                 | Location                                                                                                        |
| `(at :label l :ins c :numbering n)`                                         | systems                 | Location (one position)                                                                                         |
| `(all)`                                                                     | systems                 | Location                                                                                                        |
| `(rect :a <location> :b <location>)`                                        | as its sides need       | Location (one rect; each side one system; sides never wrap, `seq.query.range.wrap`)                             |
| `(row :key k)`, `(row :key [k …])`                                          | alignments              | Location: the rows by key (label), positions of its rows system (`seq.query.row.unknown`, `.pending`): a side   |
| `(item :id k)`, `(item :<field> v …)`                                       | annotations             | Keys (every param matches; a list param matches any of its atoms; a field the annotations lack is an issue)     |
| `(items :overlapping <location>)`                                           | annotations             | Keys (every item without `:overlapping`)                                                                        |
| `(node :key k)`                                                             | nodes                   | Nodes, addressed by the form itself (`/data/1tup/(node :key T1)`), so location forms may follow a node          |
| `(chain :id A :operator 1_555)`, `(chain :auth B)`                          | a complex               | Nodes: its units by chain, operator and author chain, each by its id; elsewhere a chain data makes, by the form |
| `(molecule :identifier "ccd:ATP")`, `:class ligand`, `:entity 3`            | `data`, a complex       | Nodes: molecules (a complex's: its units', by their entity)                                                     |
| `(atoms)`                                                                   | a molecule, a unit      | Nodes: its atoms system, `<system>/(atoms)`                                                                     |
| `(elements <items>)`                                                        | structures (Language)   | Location: the schema selection, residues or atoms by its finest field                                           |
| `(groups :residue-test p … :by level)`                                      | structures              | Keys: a key a group                                                                                             |
| `(filter :test p)`                                                          | a Location, Keys, nodes | the prefix's type, a subset                                                                                     |
| `(expand :to residue \| unit \| molecule)`, `(linked :via <pairs> :type w)` | a Location, Keys        | Location                                                                                                        |
| `(extent :of <keys>)`                                                       | any                     | Location                                                                                                        |
| `(merge :of […])`, `(intersect :of […])`, `(except :of [a b …])`            | any                     | the operands' one type                                                                                          |

- **Combinators take `:of`** (`query/combinators.ts`), renamed from `or`, `and`, `not` (the author, 2026-11-20,
  "Rename"; MolQL's names, `and`, `or`, `not` being predicate logic): `merge` and `intersect` a list; `except` its first
  operand less the others, and with one operand (alone or a list of one) the complement within what the prefix names:
  its systems' positions, its annotations' items, its nodes' children with ids. Dynamic nodes are not in that universe
  yet. A rect has no complement; `(except :of [<rect> <rect>])` is a rect less a rect, exact (Location). Operands of
  different types are `seq.query.form.type`. A layer's data filter is a trailing `(except :of …)` (Themes, Filters).
- **`row`** reads the model's `row(rows, key)`: an inline alignment's labels, else a source's rows the query model holds
  (Overlay); while `rowsPending(rows)` the form waits over the whole rows system (`seq.query.row.pending`, the address
  `waiting`); a key no held row has is `seq.query.row.unknown`.
- **Operands**: a nested form is evaluated against the same prefix (its resolver, else `repr`'s child selection;
  `query/evaluate.ts`); an address from the root. The location of `items :overlapping` is evaluated against the
  annotations' targets; a part on a system no target is reports `seq.query.location.unreachable` (projection is not
  built), a rect is refused. An operand naming a ref is refused (`seq.query.ref.unresolved`): addressing follows refs
  before the forms read an address (Addressing; `:via` is an operand too), so only a caller of `evaluateQuery` itself
  meets it.
- **Composition** (MolQL's, monadic; decided by the author, 2026-11-20): a form maps over every element its prefix
  yields (nodes, a Location, Keys read as their extents) and concatenates what each gives; a prefix element with nothing
  to give (a unit with no atoms, a complex with no such chain) yields nothing, said by a warning, never a failure of the
  whole query. A prefix of the wrong type is still `seq.query.form.prefix`.
- **Molecule node forms** (`query/molecule-forms.ts`; the forms' shared kit `query/form-kit.ts`): codes and identifiers
  match exactly (`ccd:atp` is not `ccd:ATP`); a host's `equivalents` adds identifiers to an asked one. A chain form
  gives a unit by its id, so one system has one address whichever form found it (a chain addressed by its form is the
  fallback where the prefix is no complex). `atom`, `entity` and `component` are gone (the author, 2026-11-20: "we
  should not really have (atom), we should have (elements ....)"): `(elements :atom O1G)`, `(elements :component MSE)`,
  `(elements :standard false)`, `(elements :modified true)`, `(elements :component-class saccharide)`,
  `(molecule :entity 3)`.
- **The table**: `queryResolvers(model, extra?, extensions?)`, the built-ins and the registry's forms; an extra kind
  named as a query form throws (a form's meaning is never replaced silently). `evaluateQuery` is `selectAddress` with
  the table, typed; `options.extensions` the registry's (`SeqContext.query`), none by default.

### Language

The forms by role (locations-and-queries.md, The language; MolQL's split): generators `elements` and `groups`, the
filter `filter`, the modifiers `expand` and `linked`, the combinators, and the expression forms a predicate is written
in (`query/expressions.ts`). Structure forms walk a prefix's **frames** (`query/structure.ts`): a unit of a complex (or
a lone molecule or sequence) with its residues system and its atoms when it has them, and what the prefix holds of them;
over nodes, a complex's units, a unit, a molecule, else the complexes under the node (`QueryModel.complexes`: a
source's).

- **`elements`** (`query/elements-form.ts`): its items (`query/items.ts`) are one item as params
  (`(elements :chain A :seq 33)`), columns as lists (row i across them one item, a scalar the same in every row, lists
  zipped, unequal lengths an error), or records `:items [(fields …) …]` (alone; `fields` never evaluated); `itemsOf`
  (exported `selectionItemsOf`, `readItems` as `readSelectionItems`) reads JSON's object or array of objects. Normalized
  to columns by the short names, typed by each property (a literal read as its type: `:entity 3` is the word `"3"`), a
  field absent in every row dropped, a row given twice once; `elementsForm` prints the shorter layout (one item as
  params, else columns or records by the printed length, ties to columns), so `parse(print(a))` equals `a`
  (`itemsEqual`). Range fields: `beg`/`end` over `seq`, `auth-beg`/`auth-end` over the auth numbering in its order (EC-1
  to EC-8: an end it lacks snaps inward); `seq` with `beg` or `end` in one row, `beg` past `end`, one property under two
  names: errors naming the row. Rows count against the limit.
- **Its depth** is its finest field: none below a unit, every matching unit's whole residues system; residue fields,
  those positions; an atom field, positions on the units' atoms systems. A frame without atoms held (a polymer's unit)
  is tested by residue: an atom field reads the heavy atoms of its standard component (`query/residue-atoms.ts`, PYL and
  SEC included, made once a code; a modified residue by its parent, its letter's), so the result is the residues that
  have that atom, as a structure viewer reads them, and a viewer given them picks the atoms (the author, 2026-11-20:
  "Residues, as a structure viewer reads them"); `atom-index` lacks there. Why: the record says how deep it goes, as
  MolViewSpec's does.
- **Sugar, exactly** (`itemsPredicate`): a field `:f v` is `(= (prop f) v)`, `:beg b :end e` `(in-range (prop seq) b e)`
  (an end alone `>=` or `<=`), a row `and` of its fields, items `or` of their rows; a property test holds the walk equal
  to `filter` with that predicate over the same elements. The walk differs in speed only: rows indexed by their molecule
  and unit fields (a mask of the fields a row gives, then the values), each frame looking its rows up once; a row's
  positions found by index from `seq`, its range or its author number (an index of the frame's author numbers, built
  once a frame); only the rest tested per element, residue fields before atom fields. Rows plus matches, never rows ×
  units (a quick representative counts the units' reads).
- **`groups`** tests each level by its own `:<level>-test` (molecule, unit, residue, atom; a finer level's property in a
  test is `seq.query.property.level`) and gives a key a group by `:by` (`residue` by default): a residue or an atom
  (owner its residues or atoms system, key its 0-based position), a unit (owner its complex, key its id), a molecule
  (owner the complex, key the entity id).
- **`filter :test p`**: per position of a Location (a residue, or an atom on an atoms system), per key of Keys (by its
  extent), per unit of nodes (a unit-level test; a node passes when one of its units does). A test reading finer than
  its elements is `seq.query.property.level`.
- **`expand :to`** `residue` (atoms to the residues holding them), `unit` (whole units), `molecule` (the units of one
  molecule in the complex). **`linked :via <pairs> :type w`**: the prefix's residues and every end of a pair (of `type`
  `w` or a list) with an end in them; never a projection through the pairs (links are never maps). Not built:
  `include-surroundings` (no coordinates; a host's contacts are pairs, reached by `linked`).
- **Expressions take positional operands** (the author, 2026-11-20, "Accept, but allow positional operands"), held in a
  node's `args`: `=`, `!=` (on a set of words, membership), `<`, `<=`, `>`, `>=`, `(in-range x lo hi)`, `(in x [ … ])`,
  `(and p …)`, `(or p …)`, `(not p)`, `(prop name)`, `(value <values annotation> :field f)`. Any other form with a
  positional operand is `seq.query.form.param`. `@molstar/core`'s address grammar reads them (`args`, before keyword
  params) and the operator symbols as heads (`OperatorKinds`; the author, 2026-11-20, "Accept the core change"), so
  `/data/1fin/(elements :chain A :atom CA)/(filter :test (> (prop b-factor) 60))` is text, printed back the same. A
  filter over an empty prefix tests nothing and refuses nothing.
- **Types are checked when a test is compiled**: a literal reads as the other operand's type (a number as a word where
  the property is one); an operator over a word and a number, a test that is not boolean, logic over a value:
  `seq.query.form.type`. An element whose test reads no value (a property it lacks) fails the test.

### Properties

properties.md (`query/properties.ts`): a typed symbol read per element of one level, apart from the language.

- **Built in**, the short name printed and the long one read: molecule `entity` (`label_entity_id`), `molecule-class`,
  `identifier` (words); unit `chain` (`label_asym_id`), `auth-chain` (`auth_asym_id`), `operator` (`instance_id`);
  residue `seq` (`label_seq_id`, 1-based on the molecule), `auth-seq` (`auth_seq_id`), `ins` (`pdbx_PDB_ins_code`),
  `index` (`residue_index`), `component` (`label_comp_id`, words: the alternatives), `component-class` (the words of the
  components' classes), `standard`, `modified` (non-standard, or under a `modified`-type annotation item: `type`, else
  `category`, `modified`, `glycosylation` or `lipidation`); atom `atom` (`label_atom_id`), `element` (`type_symbol`),
  `atom-index` (`atom_index`); `b-factor` (`B_iso_or_equiv`) and `occupancy`, values channels read where their channel
  lies (an atoms channel an atom's, a residues channel, as the mmCIF runtime serves them, a residue's: residue-level, so
  a residue's test may read them; the author, 2026-11-20). A bare unknown name is `seq.query.form.param`, its hint
  listing the names (a near spelling first); a namespaced one nobody registered `seq.query.extension.missing`.
- **How each is read** (`FrameReads`, made once a frame and a property): the units (`QueryUnit.entity`, `chain`, …), the
  molecules, the components and the modified index, the numberings (`auth`), the atoms (`QueryAtoms.names`, `element`,
  `residue`, `atomsOf`); a residue property on an atom reads its residue's. `b-factor`, `occupancy` and `(value …)` are
  values channels the `QueryModel.property(system, name, field?)` hook answers (on the atoms system first, then the
  residues) (`resolve/query-values.ts`: a `values` annotation on the system from the store's tiles, `'pending'` while
  one is not loaded, the frame then staged with `seq.query.property.pending`); a polymer's atoms read their residues'
  standard components (Language, `elements`).
- **An element that lacks a property** fails its test and yields nothing; one warning a property,
  `seq.query.property.absent`, naming it and the count. Why: monadic composition.

### Extensions

`query/extensions.ts`; extension-openings.md, the `query-form` and `query-property` entries.

- `{ kind: 'query-form', name, role, params, result, resolve(ctx) }` for a generator, filter, modifier or combinator
  (`ctx`: the model, the form, the prefix as a result, `evaluate`, `tick`, `limit`);
  `{ kind: 'query-form', name, role: 'operator', result: 'value', type, resolve(args) }` over values;
  `{ kind: 'query-property', name, level, type, aliases?, read }`, `read` a function over the model and a system or
  `{ values, field }` (a values annotation's field, by its address or its id on the system).
- **Composed** by `queryExtensions(entries, builtIns)` (the context's `queryExtensionsOf(registry)`): every name and
  alias `<id>:<name>`, none a built-in's, two entries under one name refused; every problem in one error
  (`seq.context.extension.namespace` for a bare one). **Checked by role** when evaluated: a result of its declared type;
  a filter's a subset of its prefix (`seq.query.form.type`). A namespaced form nobody registered is
  `seq.query.extension.missing`, before the selection runs; no fallback. A form's params print by the grammar, sorted
  (core's printer), not in declared order.

### The model the forms read

- An interface (`QueryModel`): `system(address)` (length, cyclic, numberings, parent), `annotation(address)` (target,
  fields, `items(over?)` with what is pending, a row's `system` when not the target: a pair's ends), `node(owner, key)`
  and `chain(owner, id)` for nodes data makes; for the molecule forms, optionally, `components(system)`
  (`QueryComponents`: a sequence's modified rows beside its standard components, a molecule's codes),
  `compoundClass(code)`, `annotationsOn(system)`, `units(complex)` (`QueryUnit`), `molecule(address)` (`QueryMolecule`),
  `equivalents(identifier)`, `atoms(system)` (`QueryAtoms`: names, elements, each atom's residue, a residue's atoms),
  `complexes(owner)` (the complexes at or under an address) and `property(system, name, field?)` (values channels,
  Properties). The static model and the overlay's query model implement it; tests use a hand-written one.

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
  `numbering.pending` (info), `label.syntax`, `label.int32`, `label.unknown` (with the nearest labels),
  `label.insertion`, `label.snapped` (info), `range.backwards`, `range.lap`, `range.wrap`, `range.empty` (warning),
  `position.invalid`, `location.unreachable` (warning), `ref.unresolved`, `node.missing` (warning), `property.level`,
  `property.absent` (warning), `property.pending` (info), `extension.missing`. A resolver's errors fail the selection;
  its warnings and infos are carried. A test checks the catalog against the codes the sources use.

### Examples

The stories `Features/Query` (`packages/seq/storybook/src/_stories/features/query/`) are the examples of the language: a
query box on p53 (what is typed is parsed, checked and resolved by `parseAddress`, `printAddress` and the viewer's
`AddressModel`, its result selected or highlighted, its canonical form, kind, size and issues with hints shown), and one
worked story each for items by category, ranges and slices, a result on the gene (projection is the view's, drawn from
the protein's Location), `merge`/`intersect`/`except`, and an alignment's rows and columns. Each has a headless twin in
`_test/` checking the printed form and the result's positions per system. The structure forms' examples are tests
(`query/_test/elements.test.ts`, `elements-pdb.test.ts`: 4HHB, 1FIN, 3SGJ, 1EHZ) until a structure story lands (ST2); a
form of the design not built yet (`search`, `onto`) has no story.

## Addressing

`query/address.ts`, `query/links.ts` (`@molstar/seqx-core/query`): an address resolved to what it names and back, and
links that carry addresses. Forms and addresses are one language, so this sits with the queries.

### Resolving

- `resolveAddressIn(root, model, address, bounds)`: the leading ref followed (`followRef`), the refs nested in form
  operands followed (`:of`, `:a`, `:b`, `:overlapping`, `:via`; a set as `(merge :of [...])`), then the query forms over
  `model`; a set's members give one type of result, together (else `seq.address.type`), gathered and normalized once
  (linear in the members). Why operands only: a quoted string reads as the same value as an address, and elsewhere a
  string is a value (an id, a label) that must never be rewritten. It answers a status, the result, the refs passed
  (`via`), the issues, and what is pending.
- **Status**: `resolved`; `waiting` when the result is empty and a node data makes is not held (`(node :key k)`), or the
  data is not loaded (`seq.address.waiting`, info: it resolves when they arrive); `missing` when it is empty otherwise
  (`seq.address.missing`, warning); `failed` on an error. An issue's `address` is the address asked, its `places` the
  refs passed and what they led to: an address that no longer resolves is reported with what it named, never re-bound.
- **Bounded**: `limit` (10,000 nodes or elements; more is `repr.tree.address.limit_exceeded`, never a truncation),
  `maxDepth` (64), `refDepth` (8); a signal and a deadline as a selection takes them.
- **As positions** (`locationOfResult`): a Location itself; elements by their items' extents; a system's node whole; an
  annotation's by all its items; a section (and what is under it) by the region it shows. What focus and reveal read.
- **The overlay's query model** (`overlayQueryModel(model, overlay, rows?)`, `resolve/query-model.ts`): the static
  model's systems and annotations, then the overlay's (item systems, dynamic annotations), and `(node :key k)` naming a
  node the store holds (a string key: `1` is not `"1"`); one per model, overlay and rows source. With a `RowsSource`
  (`{ data, want? }`: the store, and how whoever resolves asks for tiles; the viewer's and an export's), `(row :key k)`
  on an alignment read from a source names the first row of the label in its rows' tiles (indexed label → row once a
  table); while any is not loaded the form is pending and asks for the whole rows system through `want`, so no layer
  needs to have read them. Why: a ref to a dynamic node, or to a source's row, resolves when its data arrives, through
  the same forms.

### Back

- `addressFor(root, model, result, { prefer, numbering })`: the shortest stable address, by `prefer` (`['ref', 'ids']`):
  a ref of the spec that names exactly it (equal Location, Keys, or set of node addresses; each ref resolved once per
  root and model, `refResults`), else ids and forms. Positions as `(range :from :to)`, one as `(at :label)`, a whole
  system `(all)`, a cyclic range over the origin one wrapping range, several under `(merge :of [...])` (from the root
  when they span systems); a rect `(rect :a … :b …)`; elements `(item :id k)` or `(item :id [k …])`; nodes by their
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
- **A data adapter is a format runtime**: the core's runtimes are the default adapters of a family, "seq data ports"; a
  host's adapter (a structure viewer serving what it holds in memory) is another runtime for the same formats, put in
  place of the default by `withAdapter` (Context, Registry entries). Why: reading by region, waiting, issues and node
  identity are already the protocol's; a second mechanism would duplicate each. A host's adapter for a source it does
  not hold hands it to the default (its fallback); a source's `ref` (a host's name for data with no file) rides on the
  entry for it, never read by a default adapter (EC-185).
- **Two optional members, every runtime unchanged**: `RegionRequest.field`, the values channel asked (a `values` node's
  `field`: `plddt`, `b_factor`; absent, the source's default value; a source serving one value ignores it), and
  `OpenedSource.changes(fn) → unsubscribe`, the adapter telling that its data changed (a host's edit; whoever attached
  it reopens it in place: Store, Reattach); a change since the source was opened and before the subscription is told
  too, so an edit while it is being opened or attached is never lost. A source read once (a file) has no `changes`. Why:
  one source, one open and one parse serve every channel; a host's edit is a new reading of the same source, not a new
  source.
- **Conformance**: every adapter passes one suite (`data/_test/conformance.ts`): a well-formed coverage whose nodes the
  language takes; each listed shape read whole at level 1 an answer in the protocol's shape, with a channel asked or
  not; an aborted read an `AbortError`; `changes` stopping once unsubscribed. Why: a host is code Seq\* does not own.

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
- **Nodes of the molecule kinds** (a molecule, a `compounds` container, a complex with its units and links) are carried
  as any dynamic node, their tables typed through both encodings; the overlay validates each under `data` and runs the
  language's rules on it. Why: the fixtures' form is the protocol's, one reading for a file and a host.
- **Its runtime** joins the pairs of the request's system (by name or alias) and level that overlap the region, in
  position order, each cut to its part; a gap no pair holds is `not-served`, naming its positions (the format is
  explicit: an empty region is a pair with empty tables). A failure among the pairs is the answer; an unlisted system or
  level is not served. A list without ids gets `row` from its place in the file. Inline data (`inlineSeq`) is one pair
  over the whole system at level 1, its length in the coverage; a ranges list alone over a system longer than a tile is
  served too at the least power of 16 whose one tile covers the system, the same rows with an explicit `row` each (one
  identity at both levels). Why: a layer zoomed out on a long system reads that one tile, never the system a tile at
  level 1 (81,307 at 83 Mb), and a small inline payload needs no level between. Ranges lists only: inline values get
  their coarse levels from the store's summaries, letters are read where they fit, and other shapes stay at level 1.

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
  category, a row's key; `first`, its first residue's number, EC-39; `fields`, the rows' own fields, a column each:
  Stockholm's `#=GS` tags, an inline table's columns past `label`, `letters` and `first`, numbers when every value is
  one; a column not one value a row left out and said, `seq.source.rows.field`, the alignment kept); above it `majority`
  (a row bin per column: its majority letter byte, 255 none, and that letter's share, 0–255), aligned to the level. Why
  columns whole: an alignment is hundreds to a few thousand columns wide, so a tile of 1,024 rows stays within
  `maxBlockValues` up to 16,384 columns (wider is `too-large`), and the store's keys and tiles stay one-dimensional. Cut
  and joined by rows (`sliceTable`, `joinTables`); not written in the `seq` format yet (no reader needs it).
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
  `#=GS <name> <tag> <value>` a row's fields by tag (`OS`, `DE`, `AC`; a tag given twice joined by a space, a line
  without a value said; lines naming a sequence with no row said once with the names, `seq.source.rows.annotation`,
  info, by the runtime), `#=GF`, `#=GR` read past, the first alignment of a file): rows of the first row's width,
  another width or a character that is no letter or gap said and left out; gaps `-`, `.` (`~` read as `-`). A row's
  first residue number from its label (`NAME/27-137`, EC-39), else 1. **Insert columns** (EC-38): `#=GC RF`'s `.` and
  `-`, else a column whose letters are all lowercase (A3M's lowercase is on no column). **Pairs**: `#=GC SS_cons` as
  column pairs, each bracket kind (`<>`, `()`, `[]`, `{}`) a nesting of its own and a pseudoknot an uppercase letter
  with its lowercase.
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
  stays a field, no hierarchy is assembled (rows sharing an `ID` are pieces of one item). A key written twice in a row
  keeps the last, and an attribute named as a column (`type`) leaves the column its own: each said once, `.cell` (info).
  A row with another count of cells than nine is read with the missing ones empty (`.cells`). `fields` maps as for any
  source (`{ category: 'type', label: 'Name' }`); `seqid` and `ID` are the roles, not fields under their names. Why a
  preset: NCBI serves an INSDC record's features as GFF3, read as published.
- **A row past the system's end is the row reader's** (`collectItems`, `wrap`; render, items): a source's row ending
  past its target's length becomes two rows of one item (`[s, L)` and `[0, e − L)`, one key, one occurrence each) on a
  cyclic system; a row longer than one lap is one lap (counted); a start at or past the end moves back by laps
  (`[s − L, e − L)`); on a linear system a row is clipped, and one starting at or past the end left out, both counted
  (`seq.render.item.range`, warning, once a layout). Why there, not in the runtime: without letters the runtime knows
  neither the length nor that the system is cyclic; a source that holds its own letters is clipped by the runtime
  (EC-51, `.range`) before the reader sees a row. A layer of a cyclic system whose source is not inline also wants the
  system's last position (`LayerSource.wraps`, `cyclic`): any row that runs past the end overlaps it, so its tile brings
  the origin's piece though the view is elsewhere. Cost: one more tile a layer read, whose rows enter the layer's rows
  as any tile's do; it is not gated (a coverage does not say whether a source has rows past the length), and a row
  starting past the end is in no tile (a want is clipped to the length), so that rule is for a source that serves one.
  Inline rows were wrapped by the model (EC-9).
- **A matrix** (`shape: matrix` on the source, `readMatrix`, `RecordsSource.matrix(expected?)`): an array of arrays
  (JSON, MessagePack, an inline value) found by `records`, served as a grid's values, row-major as `Float32Array`
  (`null` and any cell that is no number NaN, said once); the source serves no region (an empty coverage). Its size is
  checked against the axes' total the grid passes (`expected`, `rows × cols`): a mismatch, or rows of unequal lengths,
  is `seq.source.records.size` (error, the matrix not served; the issue names both sizes). A typed flat array (a
  MessagePack `Float32Array`, an inline one) is kept as it is when it is a `Float32Array`, its dimensions the expected
  size's (a square when none is given: a flat array that is not one is `.size`); any other flat array is copied. A flat
  array is a typed array or one of numbers and nulls. The document is parsed once at open and the source keeps the array
  the path reaches, not the text or bytes; the size is checked on each call, and every issue carries the source's
  address. More than `DefaultMaxBlockValues` (2^24) values is refused (`seq.source.too_large`). A table format is no
  matrix (`seq.source.no_data`). Why a param and not the format: a document is a document; whether it is records or a
  matrix is the spec's to say, and the grid reading it in the fourth slice's complex PAE is the one reader.
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
  values table. **Keyed records**: an `id` with no `start`, `position`, `value` or `letters`, on an alignment's rows
  system (every record's `target`, else the source's own system, ends `/(rows)`), is a record at its place in the file
  (`[k, k + 1)`, its system as long as its records): a rows system's fields, joined by the row's key (Render, The rows
  layer, Row columns); on any other system such a file is refused as before (`.field`). Why: a table of row annotations
  names its rows, never a position, and a sequence's records without a position are a mistake to say.
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

### bigWig

`runtimes/bigwig.ts` (`@molstar/seqx-core/runtimes/bigwig`, opt in: `BigWigRuntime`, `openBigWig`): a format runtime
(`bigwig`) of `values` by region and zoom level, read by range through the assets port.

- **Reading**: opening reads the header (little- or big-endian, by its magic), the zoom headers and the chromosome B+
  tree; a region reads its level's R-tree (the full data at level 1, a zoom level's records at its reduction), the
  blocks it points to in one call, inflated (`DecompressionStream`, bounded by the header's `uncompressBufSize`) and
  decoded. Why one call: core's range reader batches the pages of all the blocks, adjacent ones in one request.
- **Coverage**: a system a chromosome (its name and size), `values` on each, levels `1` and the zoom levels' reductions
  as they are (protocol, Levels). A level the file does not have, another shape, an unknown chromosome: `not-served`.
- **Answers**: level 1 per position (NaN where the file has no data); a zoom level per bin of its width from a multiple
  of it, `count` the bases with data, `mean` the sum over them by the count, as `binValues` summarizes. A record that
  does not start on the level's multiples (UCSC's start where the data starts) is shared by its overlap with each bin. A
  region past the chromosome's end is empty, not failed.
- **Assets**: `AssetReader.openRange(uri, signal)` when the host has it (`policyAssets`: core's `openRanges` under the
  policy), else the file whole through `readBytes`, bounded by `limits.maxBytes`. A range read that fails (core's
  `AssetRangeError`: a changed file, a mangled `Content-Range`) fails the region (`failed`), never wrong values.
- **Bounds** (security.md: a file is untrusted input): blocks of one region over `limits.maxBytes`, bins over
  `limits.maxRows` and blocks inflating past `maxInflatedBytes` together are `too-large`; a tree node listing more than
  `maxNodeItems`, a tree level wider than `maxLevelNodes` (its offsets de-duplicated), more than `maxIndexBytes` of
  nodes for one open or region, a tree deeper than 16 and a node past its size are `malformed`. Blocks are inflated one
  after another into what they hold, never a buffer of the header's `uncompressBufSize` each. Decoded blocks: the last
  16 kept for a pan back; pages are core's.
- **Configuration**: `bigWigRuntime(over)` (`BigWigRuntime` is `bigWigRuntime()`), `BigWigLimits`,
  `DefaultBigWigLimits`:

  | Item               | Default | Why the default; the alternatives                                                  |
  | ------------------ | ------- | ---------------------------------------------------------------------------------- |
  | `maxNodeItems`     | 4,096   | UCSC writes 256 items a node; 4,096 × 32 bytes is a node of 128 kB at most         |
  | `maxLevelNodes`    | 1,024   | a region of one tile touches a few nodes a level; a crafted level is refused       |
  | `maxIndexBytes`    | 8 MiB   | the nodes one open or region reads: a real file's region reads kilobytes           |
  | `maxInflatedBytes` | 64 MiB  | a region's blocks inflated together: the store's byte bound, so one answer fits it |

- **Tests**: the shared source suite (`data/_test/source-suite.ts`) through core's manager by range and through whole
  bytes; files from a test-only writer on core's `PCG` (`runtimes/_test/bigwig-writer.ts`: every section type, both zoom
  levels, edges, a gap, a chromosome without data, big-endian, uncompressed); its synthetic twin is
  `delayedRuntime(BigWigRuntime)` on a fake clock.
- **Performance**: a region at the zoom level that fits, from a local file, decoded under 50 ms (heavy, `itWithin`).
  Size: its own row, under 30 kB, in no bundle that does not list it.

### mmCIF

`runtimes/mmcif.ts` (`@molstar/seqx-core/runtimes/mmcif`, opt in: `MmcifRuntime`, `mmcifRuntime`, `openMmcif`,
`MmcifIssues`): a format runtime (`mmcif`) of a structure's nodes and its units' values (molecules.md, Reading).

- **One converter**: the file is read whole (gzipped or not, by its magic; `DecompressionStream`, bounded by
  `limits.maxBytes` before and after inflating), parsed by `@molstar/io`'s text CIF parser (an allowlisted entry; the
  lint's restricted imports), its categories read by field name into a `StructureDescription`
  (`runtimes/mmcif-description.ts`), and that into nodes by `structureNodes` (Seq io), never a converter of its own. Why
  no mmCIF schema: its size; the description names the few fields read. Why the parser and the reading by dynamic
  import: a viewer that lists the runtime and opens no mmCIF source pays its own code only. BinaryCIF is not read (its
  decoder's size).
- **What is read**: entities (`entity`, `entity_poly`'s type and `pdbx_seq_one_letter_code_can`); components a position
  (`entity_poly_seq`, `pdbx_entity_branch_list` by `num`, `pdbx_entity_nonpoly`; a water one `HOH` a water of its
  longest unit); identifiers (a non-polymer's `ccd`, a branch's `wurcs` and `linucs` from
  `pdbx_entity_branch_descriptor`, `prd` from `pdbx_molecule`); units (`struct_asym`, or the assembly asked from
  `pdbx_struct_assembly_gen`, each chain under each operator, named as Mol\*'s instance ids: `ASM-1`, `ASM-X0-5`);
  residues a unit from `atom_site`'s first model (a polymer's by `label_seq_id`, a branch's by `auth_seq_id` through
  `pdbx_branch_scheme`'s `pdb_seq_num`, a water's in order), so `observed` is what was modeled and the `auth` numbering
  `auth_seq_id` and `pdbx_PDB_ins_code`; links from `struct_conn` (`covale`: `covalent`, or `glycosidic` with a
  glycosylation role; `disulf`, `metalc`) within one copy, and each branch's `pdbx_entity_branch_link` its molecule's
  `links` (`entityLinks`, once an entity, not once a unit); compounds of the entities read from `chem_comp`, a parent
  from `mon_nstd_parent_comp_id` else `pdbx_struct_mod_residue`; a molecule's atoms from its first unit (natural order
  of `label_asym_id`).
- **Left out, said once** (`seq.io.mmcif.*`): an atom name's second `label_alt_id` (the first read; the others dropped:
  a position is one atom), a `struct_conn` row with an end on a symmetry mate or on no unit read (a water when waters
  are not read). A polymer's backbone `covale` (a modified residue to its neighbor) is the chain, not a link, and a
  branch's own `struct_conn` rows are its list's; neither is said. An assembly the file lacks refuses the source
  (`seq.io.mmcif.assembly`); bytes past the bound are `seq.source.unreadable`, text that does not parse
  `seq.source.malformed`. Hydrogens are `structureNodes`' to refuse (EC-64).
- **Values by `field`**: per residue, over its heavy atoms read: `b_factor` (`B_iso_or_equiv`'s mean), `occupancy` (the
  mean), and `plddt` from ModelCIF's `ma_qa_metric_local` when a metric of type `pLDDT` is there; the default field the
  first of `plddt`, `b_factor`. Level 1 only; a unit, level or field it lacks is `not-served`.
- **Configuration**: `mmcifRuntime(options)` (`MmcifRuntime` is `mmcifRuntime()`); a source's own `assembly` param wins
  over the runtime's:

  | Item       | Default             | Why the default; the alternatives                                                           |
  | ---------- | ------------------- | ------------------------------------------------------------------------------------------- |
  | `assembly` | the deposited model | what the entry is; an assembly id of `pdbx_struct_assembly_gen`                             |
  | `water`    | `false`             | waters say nothing on an axis (molecules.md); `true`                                        |
  | `name`     | `description`       | the entity's own words (`pdbx_description`); `composition` for a branch (`NAG2 BMA`)        |
  | `limits`   | `structureNodes`'   | the description's bounds (Seq io); an assembly's units are made up to `maxUnits`, then said |

- **Tests**: the conformance suite on PDB entries recorded by hand (`examples/scripts/record-structures.ts`, CC0: 1EHZ
  quick, 1FIN, 4HHB and 3SGJ medium), each against the counts recorded from its file and counted again from its
  categories by the parser alone, and valid by the language; a small hand-written file for what no entry has (a polymer
  without canonical letters, an altloc, a hydrogen, two operators, waters, pLDDT, the answers and refusals).
- **Size**: its own row (the runtime and Seq io's converter, 19 kB) under 30 kB; the chunk on first open (the parser and
  the reading of categories) apart, with every import it brings past the default viewer's: 90 kB, half of it
  `@molstar/core`'s linear algebra, which `data/db`'s `Column` (the parser's fields) imports for its tensor schema.

### Archives

`data/archive.ts` (`readArchive`, `writeArchive`, `archiveReader`): a spec with its data files and an index as one zip,
over the assets port, with `@molstar/core/util/zip` (never a private zip). Not a port of its own: functions over
`AssetReader`.

- **Layout**: `index.json` (`{ archive: 'seqx', version: 1, spec, files }`), the spec's entry (`spec.json`, the JSON
  codec's text) and the files the index lists, by the path the spec's sources name them under. An entry the index does
  not list is never inflated. Why: the index is what the writer meant, so a reader never trusts the archive's own
  listing for what is data.
- **Written**: each `uri` and `index` of a `source` is read through the given reader (its policy decides what may be
  read), once a location, and saved under its own path when that is a path under the root, else `data/<n>-<name>`; the
  written spec names the entries. `keep(uri)` leaves a location its host's runtime serves (`generated:`) as written and
  unread. A file that cannot be read, or a bound, stops the write with issues and no bytes: an archive that lacks a file
  its spec names is not the view saved. Why: a partial archive opens to a broken view with no sign of it. Stored, not
  deflated: `.msgpack`, `.gz`, `.bgz`, `.zip`, `.bw`, `.bigwig`, images.
- **Read**: never throws on bytes (`seq.archive.*` issues: `unreadable`, `bounds`, `name`, `index`, `spec`, `data`); a
  spec the reader opens is as written, not yet checked against the language (`checkSeqSpec` is the viewer's). The reader
  serves the files by name (`archiveReader`; `./` dropped; a location not an entry is `ArchiveEntryMissing`, so an
  archive's spec reads nothing of the network or the disk). Why: the archive is the only place its view reads from.
- **Bounds** (EC-91, `ArchiveLimits`; a host raises them for a trusted archive; security.md): `maxEntries` 1,000,
  `maxBytes` 256 MB (the archive's own bytes and the entries' declared sizes together), `maxEntryBytes` 64 MB (the bound
  of a file parsed whole), `maxRatio` 200 for an entry over `ratioFloor` (64 KB), `maxNameLength` 256, the spec within
  `SpecBounds.maxBytes` and its checks. All but the spec are read from the central directory before anything is inflated
  (core's `unzip` limits, and the listing), so a bomb is refused by its declared sizes; a declared size that lies is
  caught by core (the inflated length and CRC-32 must match). Why: the declared sizes are what allocates.
- **Configuration** (`ArchiveLimits`, passed to `readArchive` and `writeArchive`; `DefaultArchiveLimits`; the host's,
  never the archive's):

  | Key             | Default | Whose                                       | Bounds                                         |
  | --------------- | ------- | ------------------------------------------- | ---------------------------------------------- |
  | `maxEntries`    | 1,000   | the host's (security.md, EC-91)             | entries, directories included                  |
  | `maxBytes`      | 256 MB  | the host's                                  | the archive's bytes, and the entries' together |
  | `maxEntryBytes` | 64 MB   | `OpenLimits.maxBytes` (a file parsed whole) | one entry's expanded size                      |
  | `maxRatio`      | 200     | the host's                                  | expanded over compressed, over `ratioFloor`    |
  | `ratioFloor`    | 64 KB   | the host's                                  | an entry this small is not held to the ratio   |
  | `maxNameLength` | 256     | the host's                                  | an entry's name, characters                    |

  The spec's bounds are the host's `SpecBounds` (`specBounds`, `DefaultSpecBounds`). The writer keeps to the same
  bounds, storing an entry that would pass the ratio compressed (a megabyte of `N`) as it is, and refuses the write for
  any other bound with the reader's issue: what is written opens. Why: a saved view that cannot be opened is lost.

- **Names**: two central entries of one name are refused (the end record's count is not the listing's); a path under the
  root only: no empty name, control character, backslash, drive, absolute path, empty, `.` or `..` segment; two entries
  that differ in case or normalization only are refused as one name (a reader on a case-folding disk would write one
  over the other). A trailing `/` is a directory.
- **Known limits**: a runtime that reads files its source does not name (a profile in a `value`) is not saved by the
  writer; two central entries of one name keep the last (core's `unzip` returns a map).

### Store

- **Keys**: `<source>|<system>|<level>|<tile>` and `<source>#<key>`; a part holding `|`, `#` or `\` escapes it, so no
  two keys of other parts are equal. A tile is `tileBins` bins (1,024) at the level's width; a want is clipped to the
  system's length when known and snaps to whole tiles; an unlisted level falls back to the coarsest listed no wider
  (EC-170). A want over `maxTilesPerWant` tiles is `seq.store.too_many_tiles` and not asked, and `keysOf` lists none of
  its tiles (a reader reads it as pending, never as empty): work is bounded before it starts, a key a tile included. A
  tile is asked for every shape its system serves: the key has no shapes.
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
- **Reattach** (`reattach(address, opened)`: a new reading of an attached source, a host's edit): nodes of a key it
  still brings keep their identity (`<source>#<key>`, so their addresses and systems) and take its content, a new
  version; a key it no longer brings goes at once, whoever holds it (the reading says it is gone: as when a patch
  removes a node; its holders' releases then count down to nothing); every tile of it and of its channels is dropped and
  the wanted ones asked again; its reports cleared; the old reading disposed. A key the old reading's tiles brought
  keeps its identity while held: served again by a tile it takes the new content (no conflict), and it goes when the
  tile that brought it is read again without it (`store-sources.ts`). A malformed coverage is an issue and the old
  reading stays. Why: EC-176's "a conflict keeps the held node" is right for two answers of one reading, wrong for a new
  reading.
- **Channels** (`channels.ts`, internal to the package: not in the `data` entry, read by the resolver and the context
  through `data/channels`): a values node's `field` is read under its own address, `channelAddress(source, field)` (the
  source's and the field's joined by a control character no address holds), made when a tile of it is first read, from
  the source as attached (`fieldChannel`: each read asks the field, no nodes, its own summaries), and gone with its
  source; reading its coverage or its keys makes nothing (the source's own); its failures are reported at the source's
  address. Why its own address: tiles, summaries and entries are the field's, and every reader of a layer's source
  (render, describe, the query model) reads it unchanged.
- Time is a `StoreClock` (real timers by default): tests drive a fake one and never wait.

## Seq io

`src/io/` (`@molstar/seqx-core/io/structure`, `@molstar/core` only): a structure's host-neutral description and the one
converter from it to the nodes molecules are read as (Language, Molecules and complexes), shared by every adapter of a
structure (a file's runtime, a host holding its model). Neither imports `@molstar/model` or `@molstar/io`. Its public
names: `StructureDescription` and its table types, `structureNodes` (with `StructureNodes`, its options and
`StructureLimits`, `DefaultStructureLimits`) and `StructureIssues`.

- **`StructureDescription`**: tables of columns named in mmCIF's words, rows referencing other tables' rows by index:
  `entities` (`id`, `type`, the polymer's type, `name`, canonical letters), `components` (a row a position of an entity,
  several at one its alternatives, the first giving the letter), `identifiers`, `chains` (the units: `label_asym_id`,
  `auth_asym_id`, entity, `operator`), `residues` (a chain's positions with `auth_seq_id`, `ins_code`, `observed`),
  `links`, `entityLinks` (a molecule's own: a glycan's bonds), `compounds` (`parents`), a molecule's heavy `atoms`,
  named per-residue `values`. Why columns and indices: a host fills them from its own model without a CIF frame, and a
  description is format-neutral.
- **`structureNodes(description, { source, water?, complex?, limits? })`**: a sequence per polymer entity (letters
  given, else each position's first component's standard letter, else its one parent's, else `X` or `N`, said once:
  `seq.io.structure.letters`; `modified` the positions whose component is not their letter's standard one and every row
  of a position with alternatives), a molecule per other entity (`class` `branched`, `water`, else `ion` when its one
  component is one atom by its compound's `formula` (without one, one atom described), else `ligand`), the `compounds`,
  and a complex of units (`of` the entity's node, `chain`, `auth_chain`, `operator`, `observed` runs when a position was
  not, an `auth` numbering from `auth_seq_id` and `ins_code`) with a `links` child; a molecule's `entityLinks` its own
  `links` child (a polymer's left out: a sequence takes none); water only when asked. It also gives the units' systems
  (their ids and lengths) and `values(field)`, a named value over each kept unit's positions from that unit's own
  residue rows: what a source serves for a `values` node's `field`. Why on the result: the values follow the units the
  reading kept, by the same keys.
- **Checked and bounded** (`structure-read.ts`; a host is code Seq\* does not own): every table and every cell is read
  as its type before anything is sized by it; what does not read is said at the source's address and left out, never
  thrown. A table missing, of another type or with a column of another length is left out (`seq.io.structure.columns`);
  a row with a cell of another type, naming no row of another table, a position past its unit, or a position listed
  twice in a unit (no row wins by order) is left out (`seq.io.structure.rows`); rows, an entity's positions and units
  past their limits are left out (`seq.io.structure.bounds`); hydrogens are left out (`seq.io.structure.hydrogen`,
  EC-64). A seeded fuzz holds it.
- **Keys from identifiers, never indices or row order**: an entity by its `id`, the compounds `compounds`, the complex
  `complex` (`options.complex`); a unit's id its chain and its operator joined by `-` (`A-1_555`; the chain alone
  without one). A chain holding `-` or a character no id takes, an operator holding one, or a key given twice (an
  entity, a unit, a compound) is refused, every row of it (`seq.io.structure.id`): never made into another id, never the
  first row kept, so no two descriptions' rows can meet under one key. Entities, compounds and units are in their ids'
  natural order, links by their ends, identifiers by namespace and value, atoms by residue and name; only a position's
  alternatives keep their rows' order (it is data). Why: every adapter of one structure gives equal nodes, so a
  selection, a bookmark or a ref made with the file holds with the host; and a removed chain, or rows in another order,
  change nothing else (a seeded property over every table).

| Item               | Default    | Why the default; the alternatives                                                   |
| ------------------ | ---------- | ----------------------------------------------------------------------------------- |
| `limits.maxLength` | 1,000,000  | positions of one entity: titin's 35,000 with room; a description is untrusted input |
| `limits.maxUnits`  | 10,000     | units of the complex: a large assembly's chains under its operators                 |
| `limits.maxRows`   | 10,000,000 | rows of one table: a ribosome's residues under its operators with room              |
| `water`            | `false`    | waters say nothing on an axis (molecules.md); `true`                                |
| `complex`          | `complex`  | the complex's id and key; any id                                                    |

## Structures

`src/structure/` (`@molstar/seqx-core/structure`): the translation between Seq\*'s Locations and a structure's units,
for a host wiring a structure viewer beside a sequence view. A library, not a port: nothing registers, nothing is
composed, no structure viewer's types (the vocabulary is the host's adapter's, which alone knows both sides). Plain data
in and out; pure over a model (`StructureModel`: the static model, or it with the overlay a source's complexes are in);
Node-tested. Its public names: `structureFacts`, `toUnits`, `fromUnits`, `unitRows`, `colorRuns`, `DefaultSyncLimits`
and `StructureSyncIssues`, with their types. Why a library: the viewer's events are the surface every host already uses;
a port would duplicate them for one consumer, and a headless host (a server coloring a structure) needs the translation
alone.

- **A target** (`StructureTarget`) is a Location whose parts are on units, or on a molecule unit's atoms: structured
  clone carries it across `postMessage`. A binding (`StructureBinding`: `{ complex, via? }`) names a complex and the
  maps its paths take; absent, every complex with no `via`.
- **Facts** (`structureFacts(model, complex)`): per unit its address, id, `chain`, `auth_chain`, `operator`, label, its
  molecule and start on it, length, class and kind (`polymer`, `branched`, `other`), observed runs (sorted, within it;
  every position when it names none), copy group (EC-65), its `auth` numbering both ways (`label(position)`,
  `position(label)`: the unit's own, else its molecule's shifted onto it, EC-8), and its residue rows (`authSeqId`,
  `insCode` split from the label, `labeled`, `observed`, a row a position). Why: an adapter joins them against its own
  residue list once and per message looks up, never walks the model. A description read by `structureNodes` gives back
  its chains' ids and residues as described (tested on the zinc finger and every mmCIF fixture).
- **To units** (`toUnits(model, loc, complexes?, options?)`): each part (a rect's sides as parts; an atoms part
  collapsed onto its residues) carried onto each bound unit along the model's paths with the binding's `via`; a
  molecule's part reaches every unit of it (EC-48). Through the molecule: a part is carried onto a unit's molecule once
  and steps down to each unit (its start subtracted, clipped to it), so an assembly of `maxUnits` copies searches once a
  part, not once a copy (a search from a sequence reads every unit's step); linear in the parts and the units reached,
  each part read once and a unit no part reaches making nothing. The copy rule as the search keeps it (EC-25): a part on
  a unit of the molecule, or a path climbing out of one onto it, never steps down to another; a unit with maps of its
  own (a step that is no unit step, the spec's or the overlay's: `Overlay.edgesFrom`) is searched on its own (its map
  may be the shorter path). A path ambiguous or past its bound is said once (`seq.resolve.path.*`); one that is none is
  not (most systems reach no structure). Then: unobserved positions cut (`cut` says them; nothing to light, EC-118); a
  molecule unit whole (`ligand: whole`, the default) or, with `ligand: atoms`, its atoms (an atoms part its own atoms, a
  residue its atoms); the size gate.
- **The size gate** (EC-167): past `maxRuns` runs or `maxLength` positions after the cut, whole units (`over: units`:
  each unit lit, an observed position or atom reached, its observed runs: a unit reached only where nothing was modeled
  is never sent) or nothing (`skip`); counting stops once over, said as `seq.structure.gated` and in `gated`. Monotone:
  a bound raised never gates what a lower one passed, and a gated target holds the exact one. `limits: null` gates
  nothing (rows and colors are written once). Why here: a 4,000-unit assembly's selection must not become 4,000 × n runs
  a frame, whatever the host.
- **From units** (`fromUnits(model, target, complexes?, { atoms? })`): a part on a bound unit kept (clipped to it); on a
  unit's atoms collapsed onto its residues unless `atoms` (a section shows the atoms form); on a molecule spread over
  every bound unit of it; the rest counted (`unplaced`, `seq.structure.unplaced` once). Round trips (property tests):
  `fromUnits(toUnits(loc))` is the Location's observed part on the units, as the model's own projection gives it; and
  `toUnits(fromUnits(t))` is `t`'s observed part.
- **Rows** (`unitRows(model, store, annotation, complexes?, { want?, ligand? })`): a ranges or pairs annotation's items
  (a range's rows, a pair's two ends) as targets on units with their fields, items reaching no unit left out. Inline
  rows at once; a source's from the tiles the store holds at its finest level over the whole system, read as the layers
  read them (`readTiles`, Render, Data): a tile not held, failed or past the store's bound on a want's tiles `pending`
  (`tooMany` says the bound), never read as empty (EC-169); `want` asks for them (the library holds no consumer).
  Another kind is `seq.structure.rows`.
- **Colors** (`colorRuns(model, store, layer, complexes, { colors, tokens?, want?, ligand? })`): a letters, values or
  ranges layer's color channel as its track compiles it (`colors`: `Renderer.layerColors()`, Render, Layer colors: the
  renderer's schemes, measures, scales, held domains and the reader's substitution, the type's default and a colourless
  layer's categorical place), with the light and the dark tokens (the built-in ones unless given; `color`, `dark`: a
  scheme's dark variant); evaluated per position (letters, values) or per item (ranges), grouped by color, carried onto
  units ungated, runs in the order of each color's first position. Unobserved positions cut, uncolored ones (a missing
  value, a letter not held) in no run, data not held `pending` as rows' are. Another layer is `seq.structure.layer`. A
  snapshot: a host rebuilds on the layer's changes for live colors. Why the view's compile: a structure colored by a
  track agrees with the track.
- **Performance**: a 2OCJ hover onto its four units is a few µs once the paths are kept, and a selection on 32,768
  copies tens of milliseconds (quick, bounded in the heavy lane); a hover on a `maxUnits` assembly is one search a part
  and molecule and a unit step each (heavy). A unit's observed runs are made once a unit; its other facts when asked.

| Item               | Default                                 | Why the default; the alternatives                                         |
| ------------------ | --------------------------------------- | ------------------------------------------------------------------------- |
| `complexes`        | every complex, no `via`                 | the spec names no host's binding; `{ complex, via? }[]`                   |
| `limits` (hover)   | 2,000 runs, 20,000 positions, `units`   | the precedent's bounds, a frame's worth; any; `skip`                      |
| `limits`           | 20,000 runs, 200,000 positions, `units` | a selection, sent once (`DefaultSyncLimits.select`); `null` gates nothing |
| `ligand`           | `whole`                                 | a ligand's position is its residue; `atoms`                               |
| a map's `priority` | 0 (Maps)                                | equal maps of one pair stay ambiguous until ranked; any integer           |
| a binding's `via`  | none                                    | a per-use choice, read as a section's (Paths, Via); a map or a list       |
| unobserved         | never sent (fixed)                      | nothing to light                                                          |
| the vocabulary     | Locations on units (fixed)              | a structure viewer's types are the host's adapter's                       |

### Issues

`seq.structure.*` (`StructureSyncIssues`): `complex` (a binding naming no complex), `gated`, `unplaced`, `rows`,
`layer`.

## Resolve

`src/resolve/` (`@molstar/seqx-core/resolve`): the static model, tier 1: pure, over a checked spec, nothing from the
store. `resolveSpec(spec, { runtimes, config })`, or `createResolver` to keep entries and paths between specs.

### Configuration

`ResolveConfig`, each choice with its recommended default (`DefaultResolveConfig`):

| Item                | Default     | Why the default; the alternatives                                                                                                                                     |
| ------------------- | ----------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `maxPathSteps`      | 4           | counted steps a path may take (EC-22): bounded search                                                                                                                 |
| `maxPathCandidates` | 16          | candidate paths kept (EC-22): bounded search                                                                                                                          |
| `encoding`          | `cheapest`  | each map in its cheapest form; `runs` or `expanded` hold every bounded map so, for tests and comparisons                                                              |
| `openQueryLength`   | int32's max | what the query forms read on a system of unknown length, so labels resolve before its source opens; layout never reads it (an axis of unknown length is pending)      |
| `maxNestedDepth`    | 4           | the levels of nested sources the overlay requests (EC-179): a spec's sources are untrusted                                                                            |
| `maxUnits`          | 65,536      | the units a complex holds (a large capsid's assembly); a complex past it is refused whole (`seq.resolve.complex.units`), a source's or a spec's; larger at its memory |

### Systems, sources, annotations

- **Systems**: each `sequence`, with its alphabet, cyclic flag, its numberings' ids (the built-ins `index`, `ordinal`,
  then its declared ones), its `name` (the record its source knows it by: the param, else its id), and its length when
  its letters are inline; else the length is unknown, pending until the source is opened, never a guess. Why: nothing in
  tier 1 waits, and nothing lays out on a made-up length. A molecule's residues and atoms and a complex's units are
  systems too (Molecules and units); a data reference to a system (`target`, `a`, `b`, `x`, `y`, a map's `from` and
  `to`) takes any of a sequence, a molecule and a unit; `pairs.a`/`b`, `ranges.target` and `map.to` an alignment's
  columns too (its address). An alignment's systems are made before the maps, pairs and ranges, which read them. A
  ranges source on an alignment is read on the system named `columns` (its columns system's name, as an alignment
  source's coverage says it; a records source names it so with `system`): pending the author's decision (a source shared
  by two alignments cannot tell their columns apart; TODO.md).
- **Numberings**: each `numbering` child of a sequence or a unit (`numberings`, by its address;
  `numberingOf(system, id)`): inline labels, an `offset` (`sequentialNumbering`, so labels snap with no scan) and
  `segments` (`maps-numbering.ts`, arithmetic a run, scanned by the forms) a numbering at once (over the open length
  while the system's is unknown, the length in the entry's key), which the system's query form answers (`numbering(id)`,
  `numberingIds`); a source's holds its source, its field and an access plan on the sequence's name, its labels tier 2's
  (Render, Numberings from a source). Why not in tier 1: nothing there waits for a source.
- **A length from coverage**: once a source is open, the length its coverage gives a system (by name or alias, an
  alignment's columns and rows too) is passed to the next resolve (`sourceLengths(model, coverage)`,
  `ResolveOptions.lengths`, `Resolver.resolve(spec, lengths)`), and the system resolves as one with inline letters of
  that length: its rows wrap the origin, its codon maps and grids are read, its sections and blocks get their pieces.
  `Lengths` is opaque: each length is kept by its system and its source node's hash, so a model whose source changed (an
  edited `uri`) never reads the old file's length; the system is pending until the new source opens. The lengths are
  part of the model's hash and of every entry's key over the system. Why again through resolve: what reads a length is
  the model's, so a cyclic `seq` file works as inline letters do. A sequence whose opened source names systems but none
  after it (its `name`, else its id, nor an alias) is kept by `sourceLengths` as unnamed, and the next resolve says
  `seq.resolve.system.unnamed` at it (an error, its path the `name` param, else `source`): its length stays unknown and
  nothing over it lays out. Why its own code: a spec fault with its own fix (name the sequence), never a wait, where
  `length.unknown` is transient while a source opens.
- **Opening a model** (`openModel(spec, store, options)`, `render/export.ts`): the one implementation of the two phases.
  It resolves (or takes the caller's first resolve), opens the sources that hold systems (`scope: 'view'`: the view's
  only), resolves again with their lengths, then opens the rest, and gives the settled model, its view and every issue
  (the model's, the view's, the openings'). `renderSvg` and the tests use it; the viewer composes the same pieces
  reactively (its store opens on each model change, `sourceLengths` on each store change). Why one helper: four
  hand-written copies had drifted (a source opened twice, a late model swapped in mid-render).
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
  it expected. Pairs and grids read by region wait for the Genomes milestone. **Links** (a `pairs` child, Language,
  Molecules and complexes) are held so too: a molecule's on its residues (`target` and `ySystem` the molecule); a
  complex's with each end on its row's unit (`PairsData.aUnit`, `bUnit`: indices into the entry's `units`, the complex's
  units in order; a row naming no unit of it is left out), so a query reads each end on its own unit and an end on one
  copy never lights another.
- **Wrapped rows**: a row that wraps a cyclic origin is two rows of one item, in the items and the table (EC-9). On a
  sequence of unknown length it is left out with `seq.resolve.length.unknown` until the length is known.
- **The query forms' model**: `system`, `annotation` (ranges; inline items at once, a source's pending over what is
  asked); the molecule forms' reads (`components`, `compoundClass`, `annotationsOn`, `units`, `molecule`, `atoms`:
  `resolve/molecule-query.ts`, one implementation for the model and the overlay); `node` answers none (no overlay in
  tier 1), and the form reports `seq.query.node.missing`.
- **Caching**: an entry is keyed by its address and the hashes of the nodes it reads (a map: itself and its two
  sequences); an equal spec gives the same model, its paths and resolved views with it.

### Molecules and units

`resolve/molecules.ts` (the entries), `StaticModel.molecular` (`molecules`, `complexes`, `units`, `atoms`, `modified`).
Seq\* holds pointers and labels: codes and identifiers are strings.

- **A molecule** is two systems: `<molecule>` (its residues, a code each; `kind: 'molecule'`, no alphabet) and, with
  atoms, `<molecule>/(atoms)` (`atomsSystemOf`; `kind: 'atoms'`). Its class is `ligand` when it names none.
- **A unit** is a system (`kind: 'unit'`) on its molecule's positions over its range, counted from 0: its position i is
  the molecule's `start + i` (`SystemEntry.parent`), its length its range's (its molecule's when it names none, unknown
  while that is), its letters its sequence's over the range: inline, sliced; a source's, read on the sequence's
  positions through the unit's step (the layer's data system is then the sequence), pending while unread (EC-169), and a
  component form reading standard components waits for them. A unit of a molecule with atoms has its own
  `<unit>/(atoms)`: the molecule's atoms of its residues. Why from 0: every system counts from 0, and the step up adds
  the start (decided by the author).
- **Steps**: a unit's `up` to its molecule and `down` from it (`unitEdges`), free in the search and under the copy rule
  (Paths); a carried search is checked against the units it read as against the maps.
- **Copies are derived** (EC-65): units of one complex with equal `of` and range share a copy group
  (`ComplexEntry.copies`, by key, each group's units in the complex's order), never declared; the same groups for any
  order of the units.
- **The atoms group map** (`AtomsEntry.offsets`: residue `r` owns atoms `[offsets[r], offsets[r + 1])`): many-to-one, so
  never a path step (a path is one linear map at one scale). It applies at a path's ends only: a part on atoms is
  collapsed onto the residues it touches before projection (`collapseAtoms`, as nucleotides onto residues), and an image
  on a molecule's residues is expanded onto their atoms after it when the target is an atoms system (`expandAtoms`).
  Focus projects through it (one jump); a follow does not (no path enters an atoms system: a follow needs a fixed scale,
  EC-113).
- **The modified index** (`modified(system)`): a sequence's `modified` rows, sorted (position, component, length);
  `modifiedAt` gives the components covering a position (alternatives in the source's order), `componentsAt` those or
  the letter's standard component (`@molstar/core`'s `standardComponent`). A unit's is its sequence's shifted by the
  start and clipped to the range.
- **Bounded**: a complex past `maxUnits` units is refused whole with `seq.resolve.complex.units` (its hint: read a
  subset of the assembly).
- **Compounds** (`compounds`): the `compound` nodes by code (name, class, parents), what `component-class` reads;
  `compound(code)` falls back to the host's names (`ResolveOptions.componentNames`, the `componentNames` entry's
  lookup). A host's identifier resolver is `ResolveOptions.equivalents` (the `identifierResolver` entry): the
  identifiers it holds the same compound as one, which `(molecule :identifier …)` adds; none by default.

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
- **Segments** (`maps-segments.ts`): the runs as written, clipped to `from` and to a non-cyclic `to` (EC-21), split at a
  cyclic `to`'s origin (EC-11); one run within its target is implicit, else runs.
- **Structure**: a segments map plus `MapEntry.structure` (each run's provenance on `from`, the `conflict` and
  `unobserved` residues on `to`); `structureStatus(map, k, observed?)` is a residue's: `unobserved` by the unit's
  `observed` when given (the unit wins), else the map's; then `conflict`; then `paired` where a run holds it, else
  `unpaired`; `structureProvenance(data, k)` a position's run's word. A structure map onto a unit whose `unobserved`
  rows are not the unit's (on the residues it pairs) is `seq.resolve.map.observed`, once. Why the unit wins: observed is
  the structure's own fact. Status changes no path: a projection through a run pairs its unobserved and conflicting
  residues as any other (a hover on the sequence lights an unobserved residue in Seq\*); cutting what has no atoms is
  the structure translation's (`toUnits`; maps.md, structure-syncing.md: unobserved positions never sent). Why: one path
  for every use, and a residue missing from the model is still the sequence's.
- **Priority**: `MapEntry.priority` (0 when absent), on each of its edges.
- A map that pairs nothing is `seq.resolve.map.empty`; a partial one is no issue.

### Paths

- **A path** is one `Linear` between two systems and a `step`: `1` (residues to nucleotides), `-1`, or `0`. Across a
  codon step the residue side is in fine places (`3k + j`), every map on it scaled by 3. Why: one linear map answers
  every question of the path; a protein through a chromosome is as sparse as its coding sequence (EC-27).
- **Search** (`findPath`): by counted length from 0, at most `maxPathSteps` counted steps, unit steps free, at most one
  codon step (EC-19), each system once. Candidates that pair nothing are dropped; equal compositions count once; of
  distinct ones the one ranked highest is taken: each candidate's maps' priorities sorted weakest first, compared in
  turn (`pathRank`, `compareRanks`; unit steps none; equal compositions count at their best rank), so a map the
  candidates share ties and the maps that differ decide; a full tie is `seq.resolve.path.ambiguous` (EC-23, EC-24:
  ambiguous unless ranked or named), its hint naming `priority` and `via`. Rank orders candidates of one length only
  (they step through as many maps), so a priority never makes a longer path beat a shorter one (EC-22). Why sorted
  priorities (the author, 2026-11-20): weakest first keeps a path as trustworthy as its least trusted step, and an
  unranked map shared by both (an offset before SIFTS and an alignment) no longer hides the rivals' ranks; why no rule
  ranks kinds or specific over general maps: a wrong map taken quietly is worse than the issue, which teaches. More than
  `maxPathCandidates` is `seq.resolve.path.limit`; none within the bound is `seq.resolve.path.none`. A path issue's
  `address` is the system it starts from and its `places` are both systems. The copy rule is in the search: after
  climbing out of a unit to its molecule a path never steps down to a unit of it, so a longer legal path wins (EC-25); a
  path ending at a ranged unit is clipped by its step (EC-26).
- **Edges**: `edgesFrom(system)` lists the spec's maps from a system, both ways: what a search reads, so the overlay
  searches over them and its own.
- **Via** (`maps-via.ts`): `path(from, to, via?)` with map addresses leaves out their rivals (other maps joining the
  same two systems, either way; `viaRivals`, `withoutRivals`) for that search; paths that reach no such pair are
  unchanged. A path that steps through none of the named maps is chosen as if the via were absent (`viaPath`: leaving
  out a rival elsewhere never settles a tie the via did not name), and `PathResult.unused` lists the named maps the
  chosen path does not step through, a name that is no map too (the author, 2026-11-20: "Warn"). A section reports each
  named map on none of its paths (follow, along, layers, grid y, the links onto it and from it) as
  `seq.resolve.via.unused`, a warning at the section with the map as its place (`viaUnusedIssue`; Resolved view, Via).
  One entry a pair and `via` set (sorted, deduplicated: `viaList`); its carried search's deps are every map of the
  systems it read, rivals too. `project(loc, to, { via })` and the overlay's `path(from, to, via?)` take it, and a
  structure binding's `via` (Structures) is read as a section's: by `toUnits`, `unitRows` and `colorRuns` on that
  binding's paths.
- **Kept**: a search is done once per model (the reverse read as the inverse), and carried to the next resolve with the
  map entries of every system it read; it is kept while those are the same objects (entries are kept by hash, a
  `priority` in it, so a changed priority drops the search), so an edit elsewhere keeps the path object. Why: a map
  added or changed where the search looked can make it shorter or ambiguous; nothing else can.
- **Projection** (`projectRanges(path, ranges, out?)`, `StaticModel.project`): ranges along the path, into the caller's
  buffer when given; from nucleotides to residues, the residues they touch (EC-18). A part with no path is left out with
  its issue; rects are not projected yet. `project` applies the atoms group map at the ends (Molecules and units).
- **A path read the other way** (`inversePath`) reads each step the other way: a unit's `up` is `down`.
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
`map.empty`, `map.observed`, `via.unused`, `map.codon_length`, `map.table`, `length.unknown`, `system.unnamed`,
`path.none`, `path.ambiguous`, `path.limit`, `complex.units`, `ref.system` (a reference to a sequence, molecule or unit
that made no system: a unit of a refused complex or of nothing held; what reads it is left out), and the overlay's
`node.invalid`, `node.ref`, `node.depth`. A test checks the catalog against the sources.

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
  search per version, never on hover. `Overlay.edgesFrom(system)` lists what such a search reads from a system (the
  spec's maps, the overlay's units' steps and item maps): a reader telling a unit with a map of its own (the structure
  translation) reads the overlay's too.
- **Nested sources** (EC-179): `nested(needs)` lists the source nodes at a depth some template on their root needs,
  bound to the first runtime that reads their format (`OverlayOptions.runtimes`; none: `seq.resolve.source.runtime`);
  one deeper than `maxNestedDepth` is `seq.resolve.node.depth` and not listed. Why on need: a nested source is a
  request, and a view that shows no domains asks for none.
- **Named systems**: `systemNamed(name, near)` is the system a nested source's coverage names: an item system of a node
  its parent source (`near`) brought, else a spec system of the name.
- **Holds**: `holds(addresses)` is the nodes of the given addresses (a node; an item system, a unit or an atoms system:
  its node) and, transitively, the nodes their `target` and `source` name and a complex's units' `of`: what is shown
  holds what it is drawn through.
- **Molecules and complexes** (`overlay-molecules.ts`, `Overlay.molecular`): a source's `sequence` and `molecule` nodes
  are systems at their node's address (a molecule's atoms at `<node>/(atoms)`), its `compounds` node its compounds, and
  a `complex` node's units systems at `<complex node>/<unit id>` with their steps (Resolve, Molecules and units), their
  inline numberings (`auth`) answering the query forms. Sequences and molecules are taken before complexes, complexes
  before annotations, at each depth. Each unit passes `checkUnit` against its molecule's length and each links node
  `checkLinks` against the complex's units, as a spec's do; one with an error is left out with its issues, and a unit
  whose `of` names nothing held is `seq.resolve.node.ref`. A complex past `maxUnits` units is refused whole. A
  structure's unit names its entity's node by a string key (`(node :key "1")`), as `nodeAddress` prints it: a number
  never names a node.

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
- **Via**: a section's `via` (`ResolvedSection.via`, sorted; absent when none) is what its paths take where two maps
  join one pair of systems (Resolve, Paths): its follow, its `direction: { along }` (where the map does not reach its
  system itself), its layers' data (templates' instances too), its grid's y, its navigators' paths, and the links onto
  it and from it (hover, selection; the viewer's hover fallback and focus hull read it too). A link takes both ends'
  `via` (the author, 2026-11-20), so two sections of one chain, each naming one map, light different residues for one
  hover (a comparison), and a hover on either reaches a section naming none through its map; between two sections naming
  none the pair stays ambiguous and there is no link; two naming rival maps of one pair leave it tied, so the target's
  map carries it and `seq.view.via.conflict` says so at the target (never a link gone silently). Why both: a section's
  via says how its positions meet the others', whichever way the projection runs. A named map on none of those paths is
  `seq.resolve.via.unused` (a warning at the section, the map its place), and those paths were chosen as if the via were
  absent. Why warn: a via that settles nothing is a figure's mistake the reader cannot see.
- **Grids**: a grid section's `y` (`GridAxis`) is its second axis, down the plot: `of_y`'s system whole (`of`'s when
  absent), top down, and what it shows: `follow_y: x` (the default) x's window (`followsX`; through the path from x's
  system when `of_y` is another), a section's region through `follow` (the path from its system), `none` a window of its
  own (`viewport.y`). A y that cannot follow (no path) is `seq.view.follow.unplaced` and keeps its own window. A section
  whose `of` names no system (with blocks: its first block's) is `seq.view.section.system`; its layers are not placed.
- **A grid over blocks** (a grid over a complex): a grid section with no `of` and no blocks of its own that `follow`s a
  section of blocks takes that section's blocks and axis (their display units, on the leader's address), its follow the
  identity on those units, and its y the same blocks (`followsX`). Its tracks' layers are placed per block as a section
  of blocks' are. Why the leader's blocks: the section of blocks already orders the chains and their copies on one axis.
  Refused (`seq.view.section.system`): a grid with blocks of its own, and one following blocks in a stack or on a ring
  (the plot's axes are a row's units).
- **A section of a complex** (`of` a complex; `view/blocks.ts`): a block a unit, at `<section>/<unit id>`, ordered by
  the section's `classes` (`[polymer, branched, ligand, ion]` when absent: water left out, EC-101), polymers in the
  complex's order, then each class in order; a class not listed is left out. Why left out: hundreds of waters say
  nothing on an axis. **An override**: a `block` child with no `of` whose id is a derived unit's (a merged copies
  block's: its first unit's) sets that block's `form` and `label` only; the rest stay derived, from the complex and
  `classes`, units added later included. When every `block` child is an override the section lists every derived unit;
  otherwise the `block` children win (a writer picks units, forms and order). Why: one form switched (`block.form`) must
  not freeze the list of units (the author's decision, "override one unit").
- **A block's form** (`form`): as written, else `whole` for a molecule or a unit of one and `residues` for a sequence or
  a unit of one; `atoms` puts its axis on the atoms system (`<unit>/(atoms)`), its `range` the atoms of the residues in
  range, and falls back to the default where there is none (the language says so). **Its label**: as written, else its
  unit's `label`, else (a molecule) its `name`, else its one component's code (`ATP`), else its id.
- **Blocks**: each its axis (its system, its `range`, the section's direction), its form, its first display unit in the
  section's (`offset`: 0 in a stack), its label, the blocks it stands for (`members`) and their systems (`systems`: a
  merged copy's unit, what an end on a copy is drawn on), its copy group (`copy`, internal: the pairs layer's), and each
  of the section's tracks' layers resolved onto its system (a path a block: placed on each block its data reaches; a
  layer that reaches no block says why once, by its first block's issue). Copies (the same system, or units of one copy
  group, with the same range and form) are one block labeled `×N` with `copies: 'merge'` (the default: a copy group's
  sequence said once, its cost one block's, EC-65, EC-67), a block each with `each`. The section's axis is the blocks'
  display units in order: on one system in a row (or on a ring), the blocks' pieces end to end (their breaks between
  them, a follow along one axis); else units `[0, total)` on no system (its `system` the section's address,
  `pieces: 'blocks'`), the total the sum in a row and the longest block in a stack. Why: a position on them is no
  block's, so a selection, a focus, a bookmark or a cursor on a system never reads them as its own; a block's hits carry
  its system and its index (`block`), and a consumer maps them through the block's axis and `offset`. A block that names
  no system is `seq.view.section.system` and left out.
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
- **Links**: every ordered pair of the view's sections a static path joins, with the path (`fromSystem`, `toSystem`: the
  systems it joins), composed when the view is resolved so hover and selection never compose; the same system is the
  identity. A section of blocks links through each block's system (its own axis is on no system), so what is on a
  sequence lights the units of a complex's section (EC-25) and back.
- **Projection** (`projectOntoView(model, view, loc, { close?, own?, linked? })`, `project.ts`): a Location on the
  view's sections' systems (a section of blocks: its blocks' systems), one Location keyed by system (two sections on one
  system show the same): each part carried along the view's link from its system to each shown system, else, for a
  system no section shows, along the model's static path; rects are not carried. `close` closes it first among the
  view's systems (a hover's: a nucleotide lifts to its codon and gains its residue); `own` keeps the parts already on a
  section's system (a hover's and a highlight's), else only what it adds where it has no part (a selection's);
  `linked: false` keeps the own parts alone. Why keyed by system, not a function per section: the renderer reads a part
  per axis with an indexed loop, nothing called or made per draw.
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
- **A second consumer** (EC-183; `packages/seq/examples/src/_test/consumer/`, tests only): letters, ranges, values and
  variants of a section drawn as plain SVG, compared with core's SVG as positions, text, heights and colours, forward,
  reversed and across a tile edge. What an outside drawer needs is public: the resolved view, the model, the store, the
  theme compile (`theme`), and from `render` the data access (`layerSource`) and the token table (`LightTokens`); a test
  fails on any other import from `render` (layouts, draw lists, backends) and on any path into `src/`. Colours agree
  partly through shared code (the theme compile), positions, text, heights and the dedupe do not. Missing, replicated
  from the code or the spec text: the layer types' default channel forms (a letters scheme by alphabet, `token:value`,
  `token:range`, `token:variant`: render-private), the open values domain (the values in view down to zero: stated only
  in the code), the row filter (`compileRowFilter` is not exported: a consumer cannot apply a layer's query), the
  complement of a nucleotide on a reversed axis and the item key that drops a range seen in an earlier tile (both spec
  text, read again). Not drawn: a layer behind a path (the maps are public, `projectRanges`; not built).
- `seq.view.*` (`ViewIssues`): `section.system`, `layer.type`, `layer.data`, `layer.shape`, `layer.unplaced`,
  `follow.unplaced`, `section.pieces`, `direction.unresolved`, `direction.mixed` (info); a test checks the catalog.

## Render

`src/render/` (`@molstar/seqx-core/render`): tier 3, the resolved view laid out for a size, a viewport, a target and the
store's entries, as draw lists; hit-testing on the same layout; the recording and SVG backends; `renderSvg`.

### Configuration

`RenderConfig`, each choice with its recommended default (`DefaultRenderConfig`, `renderConfig(over)`). The sizes are
the UX foundations' numbers, CSS px, settled by the foundation stories; where no reason is given, none is recorded.

| Item                                                           | Default                                                                                                                | Why the default; the alternatives                                                                                                                                                                                                                                                                                                                                                                      |
| -------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `inset`                                                        | 8                                                                                                                      | the page's inset around everything                                                                                                                                                                                                                                                                                                                                                                     |
| `labelWidth`, `labelGap`, `labelColumn`                        | 112, 8; measured, 80 to 160                                                                                            | one label column a view, so the data columns align; measured to the widest track label between `labelColumn.min` and `.max`, an ellipsis past it, never a hard clip (decided from the defaults audit: a fixed 112 clipped real names); a view with an alignment's rows, or `labelWidth` given alone, keeps `labelWidth`                                                                                |
| `trackGap`, `sectionGap`                                       | 4, 20                                                                                                                  | gaps between a section's tracks and between sections: a section is a unit, grouped by proximity (20, the audit; was 8)                                                                                                                                                                                                                                                                                 |
| `trackPadding`                                                 | top 0, bottom 0                                                                                                        | a track's padding inside its well when its spec gives none (a host's default); never left or right: the horizontal is the axis                                                                                                                                                                                                                                                                         |
| `trackShape`                                                   | `squircle`                                                                                                             | a track's well is a surface, so a squircle, its data clipped to it (the author's look); `square`, a box                                                                                                                                                                                                                                                                                                |
| `heights`                                                      | letters 20, ranges 36, values 32, ruler 22, secondary 16, variants 32, frame 16, rows 240, logo 48, pairs 64, grid 240 | a track's fixed height by its layers' types (the largest): what `heights: fixed` gives, and a lane-less track's (Track heights)                                                                                                                                                                                                                                                                        |
| `laneHeight`, `laneGap`, `lanePadding`, `maxLanes`             | 14, 2, 4, 5                                                                                                            | a lane of ranges (a 12 px mark and its gap), the padding above and below a band's lanes, the most lanes a content height takes (one lane 22 px, the ruler's band; five 78): Track heights                                                                                                                                                                                                              |
| `density`                                                      | letters 10, cells 3, dots 1, bars 1/64                                                                                 | each density step's lowest pixels a position (below)                                                                                                                                                                                                                                                                                                                                                   |
| `densityUnits`                                                 | letters `css`, the rest `device`                                                                                       | legibility is in CSS px; what bounds the lists is a mark a device pixel                                                                                                                                                                                                                                                                                                                                |
| `letterSizes`                                                  | 9 to 14                                                                                                                | a short scale; the largest that fits the cell is used                                                                                                                                                                                                                                                                                                                                                  |
| `monoAdvance`                                                  | 0.6 em                                                                                                                 | a letter's size known from its cell with no measuring                                                                                                                                                                                                                                                                                                                                                  |
| `labelFont`, `itemFont`                                        | `12px sans`, `10px sans`                                                                                               | a family role a backend maps to its stack                                                                                                                                                                                                                                                                                                                                                              |
| `labels.minChars`                                              | 4                                                                                                                      | a label cut inside its mark keeps at least so many letters, else it is dropped (the author: “4 characters as default is fine, should be configurable”)                                                                                                                                                                                                                                                 |
| `letters.modifiedMark`                                         | `shape: dot`, `size: 3`                                                                                                | a modified position's mark (`letters` `modified: mark`): a dot of 3 px at the cell's top right in the letter's ink, never over the letter; a bar over a component of several positions; a tick at a strip's top below the cells (the mark survives a zoom out); `underline` a bar along the cell's bottom                                                                                              |
| `fade`                                                         | 12                                                                                                                     | how far past its max width a text fades out                                                                                                                                                                                                                                                                                                                                                            |
| `squircleRatio`, `squircleMax`, `squircleMin`                  | 1/3, 8, 6                                                                                                              | a radius a third of the height, capped; none below the minimum                                                                                                                                                                                                                                                                                                                                         |
| `outline`                                                      | hover 1.5, selection 1, cursor 1                                                                                       | outline widths                                                                                                                                                                                                                                                                                                                                                                                         |
| `emphasis.selectionPlace`                                      | `over`                                                                                                                 | a tint under the data is hidden by opaque cells and a selection then redraws the data; `under` (to confirm with the author)                                                                                                                                                                                                                                                                            |
| `emphasis.selectionShape`                                      | `column`                                                                                                               | a selection is of positions: one band a section across its tracks and gaps reads as one (the author's look); `tracks`, a box a track                                                                                                                                                                                                                                                                   |
| `emphasis.selectionTint`, `emphasis.selectionOutline`          | 0.18, 0.6                                                                                                              | a soft tint and one thin outline (the author's look)                                                                                                                                                                                                                                                                                                                                                   |
| `emphasis.highlightStripes`                                    | 4                                                                                                                      | stripes a band splits into where channels overlap; the rest share the last                                                                                                                                                                                                                                                                                                                             |
| `emphasis.highlightBand`                                       | 3                                                                                                                      | a `band`-style channel's bar, px                                                                                                                                                                                                                                                                                                                                                                       |
| `emphasis.highlightTint`                                       | 0.3                                                                                                                    | a highlight channel's tint                                                                                                                                                                                                                                                                                                                                                                             |
| `emphasis.projected`                                           | 0.5                                                                                                                    | what is projected onto a linked section, a factor of its own look's opacity: the same look, lighter                                                                                                                                                                                                                                                                                                    |
| `emphasis.hoverOffset`                                         | 1                                                                                                                      | the hover's outline outside its mark, concentric, never touching the selection's                                                                                                                                                                                                                                                                                                                       |
| `emphasis.cursorInk`                                           | `msui-text`                                                                                                            | the cursor neutral, never the selection's accent (the author's look)                                                                                                                                                                                                                                                                                                                                   |
| `axis.directionMark`                                           | `label`                                                                                                                | a reversed axis says so in its first label (an arrow and "reversed", or "along" what it reads along), in the picture and the export; `none` for a host whose own chrome says it                                                                                                                                                                                                                        |
| `header.slots`, `.slotPx`, `.minHeightPx`, `.sectionPx`        | 2, 20 px (an icon of 14), 16 px, 20 px                                                                                 | the eye and the color fit beside a label in 112 px, a target of 20 px on a 24 px row; below 16 px an icon overlaps the next track (its states then marks); a section's line a line of the label font (Track headers)                                                                                                                                                                                   |
| `header.actions.show`                                          | `hover`                                                                                                                | the track hovered (its header or its data), the cursor on its header, the header the viewer names: a figure at rest shows names and states, not buttons; `always`, `never` (the menu only)                                                                                                                                                                                                             |
| `header.states.export`                                         | `include`                                                                                                              | the states are facts about the picture, exported as drawn; `none` leaves them out of a figure (the slots and the "⋯" never are)                                                                                                                                                                                                                                                                        |
| `section.header`                                               | `auto`                                                                                                                 | a section's header line when it has a state or a slot: no line that says nothing, records without hidden tracks unchanged; `always`, `none`                                                                                                                                                                                                                                                            |
| `margin`                                                       | 0.5                                                                                                                    | what a want adds on each side, a fraction of the view (the store's `margin` priority)                                                                                                                                                                                                                                                                                                                  |
| `followArc`                                                    | `shortest`                                                                                                             | a cyclic follower takes the shortest arc holding the image (EC-12); `hull`                                                                                                                                                                                                                                                                                                                             |
| `maxTemplateTracks`                                            | 200                                                                                                                    | the tracks one template makes; past it `seq.render.template.limit`: a source's nodes are untrusted                                                                                                                                                                                                                                                                                                     |
| `maxNumberingPositions`                                        | 2^18                                                                                                                   | the positions a numbering from a source covers at most; past it `seq.render.numbering.limit` and nothing asked for: its labels are held whole                                                                                                                                                                                                                                                          |
| `nodes.holdShown`                                              | true                                                                                                                   | a layout lists the nodes it shows, for its owner to hold in the store; false: nothing held by showing                                                                                                                                                                                                                                                                                                  |
| `ranges.overflow`                                              | `extra-lanes`                                                                                                          | only the items past the lanes collapse into a density row; `layer`, the whole layer                                                                                                                                                                                                                                                                                                                    |
| `values.aggregate`                                             | `minmax`                                                                                                               | below a device pixel a position, each column's extent, so negative extremes show; `max`                                                                                                                                                                                                                                                                                                                |
| `values.coarser`                                               | 1                                                                                                                      | the coarsest level whose bin is no wider than k device pixels' positions (below)                                                                                                                                                                                                                                                                                                                       |
| `values.direct`                                                | true                                                                                                                   | data on the section's own system placed with no mapping (below); false: the general mapping, the same lists                                                                                                                                                                                                                                                                                            |
| `ruler.labelGap`                                               | 48 px                                                                                                                  | the least between two ticks (the foundations' number); labels clear each other by an eighth of it                                                                                                                                                                                                                                                                                                      |
| `text.scale`, `text.groupGap`                                  | 12 px a position, half a cell                                                                                          | the letter size a text opens at (its zoom); a space between words                                                                                                                                                                                                                                                                                                                                      |
| `text.lineGap`, `text.bandGap`                                 | 10, 3                                                                                                                  | a line block reads as one when the line gap is the larger (10 and 3, the audit; were 8 and 2)                                                                                                                                                                                                                                                                                                          |
| `text.collapseEmpty`                                           | true                                                                                                                   | an empty band separates what belongs together (EC-96); false: uniform lines                                                                                                                                                                                                                                                                                                                            |
| `text.maxLines`, `text.maxHeight`                              | 100,000; 600 px                                                                                                        | a figure, not a genome; a text in a host page scrolls inside its box past it                                                                                                                                                                                                                                                                                                                           |
| `text.pieceBreak`, `text.numberGap`                            | `none`; 8                                                                                                              | lines stay aligned, the boundary marks the edge; the gap between a line and its numbers                                                                                                                                                                                                                                                                                                                |
| `text.heights`                                                 | letters 18, ranges 4, values 3, ruler 14, secondary 6, variants 10, frame 14, rows 18, logo 24, pairs 24, grid 120     | a text's tracks thin: a ranges track reads as underlines, a values track as a strip                                                                                                                                                                                                                                                                                                                    |
| `rings.hole`, `rings.maxDiameter`                              | 0.35 of the radius; 720 px                                                                                             | room for the title, chords readable; a map, not a wall (`RenderConfig` preset; 0 to 0.8; px)                                                                                                                                                                                                                                                                                                           |
| `rings.blockGap`                                               | 1°                                                                                                                     | a gap between blocks reads as a boundary (Circos)                                                                                                                                                                                                                                                                                                                                                      |
| `rings.letters`                                                | `upright`                                                                                                              | readable at any angle; `tangent`                                                                                                                                                                                                                                                                                                                                                                       |
| `rings.center`                                                 | `title`                                                                                                                | the name and length in the hole; `none`                                                                                                                                                                                                                                                                                                                                                                |
| `rings.labelGap`, `rings.leaderBend`                           | 8, 12 px                                                                                                               | leaders fan out over the ruler, never through one point                                                                                                                                                                                                                                                                                                                                                |
| `rings.cluster`                                                | `over: 6`, `degrees: 12`, `fanMax: 12`                                                                                 | more than 6 of one layer's labels within 12° (across the origin too) are one "n items" label, its fan of members placed with the layout and opened by a hover (the `tip` pass; nothing laid out on a hover); 0 degrees never groups; a group of up to `fanMax` fans out (past a dozen lines a fan is a hairball: the author's choice), a larger one is listed in its card; `Infinity` fans every group |
| `rings.window`                                                 | `followed`                                                                                                             | a window means only what a follower shows; `always` (idle, quiet, when nothing follows)                                                                                                                                                                                                                                                                                                                |
| `pairs.maxArcs`, `pairs.hitPx`                                 | 4,096; 4                                                                                                               | bounded by the screen: past it arcs binned by their ends' device columns (`seq.render.pairs.limit`); a curve is thin, so a hit within 4 px                                                                                                                                                                                                                                                             |
| `grid.rectMinPx`                                               | 12                                                                                                                     | rects with outlines only where cells are few and large enough to read; image tiles below                                                                                                                                                                                                                                                                                                               |
| `grid.buildBudget`                                             | 65,536 samples a frame                                                                                                 | a frame while moving under 8 ms; a tile past it draws a coarser level meanwhile, never empty                                                                                                                                                                                                                                                                                                           |
| `navigator.window.minPx`, `.hitPx`, `.touchHitPx`, `.handlePx` | 16, 12, 24, 5 px                                                                                                       | a window under a few px cannot be seen or grabbed; a press in the middle still pans; a finger is wider                                                                                                                                                                                                                                                                                                 |
| `navigator.window.floor`                                       | `centred`                                                                                                              | the drawn window says where; `thumb` (a scrollbar's mapping)                                                                                                                                                                                                                                                                                                                                           |
| `navigator.trackScale`, `navigator.trackMinPx`                 | 0.25, 2 px                                                                                                             | a miniature reads as the view; a thinner track is not seen at dpr 1                                                                                                                                                                                                                                                                                                                                    |
| `navigator.marks`                                              | selection, highlights, cursor                                                                                          | where is what I picked; any subset                                                                                                                                                                                                                                                                                                                                                                     |
| `navigator.rail.grownPx`, `.grabShare`                         | 6 px, 0.4                                                                                                              | a target in use; an alignment ruler's labels select columns, its bottom pans                                                                                                                                                                                                                                                                                                                           |
| `navigator.overview.minPlotPx`, `.minPx`                       | 160, 48 px                                                                                                             | a narrow frame gives way by the plot's width: the overview narrows to `minPx` to keep the plot `minPlotPx` wide, hides below it (the rows rail kept), then the label column and the tree, then the frame scrolls (the author's choice); any px, `minPlotPx` 0 never hides                                                                                                                              |
| `internedStrings`                                              | 4,096                                                                                                                  | strings draw lists keep interned across frames, so a frame drawn again interns nothing; past it a full redraw starts the table again                                                                                                                                                                                                                                                                   |
| `imageBytes`                                                   | 64 MB                                                                                                                  | the shared image cache: every `image` tile a layer draws (grid, rows, navigators) in one LRU, since two caches could not lend each other room (The shared reduction)                                                                                                                                                                                                                                   |
| `rows.labelMinPx`                                              | 8 px                                                                                                                   | the least row height a label, a number, a chip or an icon is drawn at; below it the label column is empty (the rows layer, decided)                                                                                                                                                                                                                                                                    |
| `rows.pinInPlace`                                              | true                                                                                                                   | a pinned row is also drawn in its place in the order, marked `·`, so the order reads whole; false: in the pinned band only                                                                                                                                                                                                                                                                             |
| `rows.buildBudget`                                             | 2^18 cells                                                                                                             | a moving frame builds at most so many cells of the rows' new level, the level drawn last drawn meanwhile; a settled frame builds the rest (a 10,000 × 300 level took two frames of about 100 ms)                                                                                                                                                                                                       |
| `rows.columns.floorPx`, `.gapPx`, `.dataMinPx`, `.titlePx`     | 32, 6, 200, 16 px                                                                                                      | a row column yields to a narrow view down to its floor, right ones first (EC-107); the gap before each; the data's least width; a column's title band                                                                                                                                                                                                                                                  |

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
  (`seq-ss-helix`, `-strand`, `-turn`), `seq-variant`, `seq-ruler`, a modified position's outline (`seq-modified`, the
  letters' `modified: tint`), the navigator's window, grip and wash, the centromere, a letter cell's edge), as
  `0xRRGGBBAA`, the logo's chroma floor (`seq-logo-chroma`, a number) and the font stacks of the two family roles. Node
  has a light and a dark table (`LightTokens`, `DarkTokens`); the viewer reads its own from CSS.
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
  `selection`, `hover`, `cursor`, `labels`, `tip` (`Passes`; `window` to `cursor` and `tip` are `EmphasisPasses`).
  `window` is the first of the emphasis, under the highlights: a window is no data emphasis, and a highlight channel
  stays above it. It is the layouts' own (their chrome): no layer draws there. `tip` is the last, over the labels: what
  a hover reads whole over its neighbours (a column cell cut short), on a surface of its own, so no label shows through
  it. The recorder writes an empty pass as its header, but for an empty `tip`.
- **Primitives** as structs of arrays in growable typed arrays kept frame to frame (`reset` zeroes counts): `rect` (x,
  y, w, h, fill, stroke and width, a squircle radius, 0 square, and which ends are round, `Round`; `cutTop` and
  `cutBottom` square a round end's top or bottom corner: a grid's cell inside its plot, a row cut by its band), `path`
  (move, line, cubic, close in shared arrays; fill, stroke, dash), `glyphs` (a run of letters in one font, an x a letter
  at its center, one baseline, a color a letter or one; a **placed** run gives each letter its own y, angle and vertical
  scale: letters along a ring, a logo's stacked letters, a row's letters; one form for the three, so one code a
  backend), `text` (an interned string, x, y, anchor, font, color, its measured width and a max width: faded past it,
  never ellipsized), `image` (RGBA bytes the layer owns), `clip` (a box, or a squircle with a radius and round ends: a
  track's well) and `unclip`. A pass's `ops` keep the order its primitives were added in. Colors are `0xRRGGBBAA`;
  `NoColor` (0) is none.
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

`Emphasis` (`layer-type.ts`), drawn by the view (`tracks-emphasis.ts`) in the emphasis passes, a strip at a time in its
frame through its warp (Section layouts), over the strip's **band column** (its bands and the gaps between them;
`SectionLayout.columns`, in tracks `SectionLayout.band`; the shape `column`, the default) or over each band's box
(`tracks`), `emphasis.selectionShape`, on the strip's transform and system. Why `column` by default: a selection is of
positions, not of a track, and one band reads as one selection; why a strip: the selection is one look in every layout
(a band on a line in text, a sector over the rings on a ring). A layer with a second axis (`axes: 2`) paints the parts
on its y systems itself, in the emphasis passes' `draw` (it alone knows its y transform); the x painting stays the
strip's.

- **The order of emphasis** is fixed, below to above: data, highlight channels, selection, hover, cursor; a higher one
  never hides a lower one entirely (tints are translucent, outlines are thin lines that never touch).
- `highlights`: the channels in order, each a Location, its element keys (`elements`: a layer of keyed items draws those
  it holds, the pairs layer's arcs, chords and dots, whatever the Location covers), a color and a `style` (`tint`, the
  default; `band`; `outline`), on every system its Location has (the `highlight` pass), each style a squircle (Shapes:
  round where a channel's range really ends, square where a stripe change cuts its run; a `band` bar, below the minimum,
  square). **Overlapping channels stay apart** (ux-foundations section 6): a position lit by one channel takes the whole
  band; by several, the band splits into equal stripes, one a channel present there, in channel order, at most
  `emphasis.highlightStripes`, the channels past the last sharing it (a ring's band radially). The style draws in the
  channel's stripe: `tint` fills it at `emphasis.highlightTint`, `band` a bar `emphasis.highlightBand` px along its foot
  in the color as given. A fill's run is one rect while its slot and the stripe count stay: it is cut only where the set
  of channels lighting changes (`highlight-stripes.ts`, a merge of the channels' sorted boundaries: nothing allocated).
  **An `outline` takes no stripe** (it has no fill, so it does not count among the channels lighting): one stroke
  (`outline.selection`) over the channel's own ranges on the whole band, above the fills, never fragmented by another
  channel's boundary. Stripe edges and the bar are on device pixels, one edge between two stripes (no overlap, no seam).
  A range wholly in a hidden stretch lights its break once a channel, at its first run. **The grid** stripes its
  highlighted rows the same way across the plot (a stripe a channel, in channel order; `band` a bar along a stripe's
  left edge, `outline` one stroke over the row); the columns are the strip's, drawn with the tracks'. **The navigators'
  ticks stay one tick a channel** (a tick is a position mark, no band to split). Why stripes: stacked tints of two
  colors blend into a third that names neither.
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
  hover from elsewhere (a structure, no layer) is drawn on every section (the `hover` pass). Every hover a squircle
  (Shapes): a cell's (letters, logo, values) by the cell's size, square on a side the view cuts and below the minimum,
  its radius the cell's plus the offset; a letter's codon parts round on the element's outer sides alone; a chip's
  concentric with its pill; the projected hover round where its range really ends; a secondary structure's mark (its
  radius plus the offset, its round ends: a helix round, a strand or a turn square), a band's (square at the view's
  edge), a rows cell's (square where the view or its band cuts it); a grid's hovered or selected cell by its size, its
  crosshair's column and row square (they end at the plot: a cut).
- `cursor`: a 1 px outline of its cell (`outline.cursor`) over the band, in a neutral ink (`emphasis.cursorInk`,
  `msui-text`), never the accent (the `cursor` pass); a squircle by the cell over the band, inside it by half its width
  (concentric), square on a side the view cuts. Whether it is shown is the caller's (a viewer shows it once the keyboard
  moved it).
- **A box a track** (`selectionShape: 'tracks'`) in a layout that fills its bands: the highlights' and the selection's
  tints drawn inside the track's well clip, so they take its shape; then, after the clip, the outlines (the selection's,
  an `outline`-style highlight's; the hover's in its own pass): an outline around a mark is never clipped, so a cut end
  keeps its square corners. Under the data (`under`) the same, in the background pass, the clip in the page's frame.
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
  outlines. **Every emphasis is a squircle** by the same rule (Emphasis): the hover, the projected hover, the selection,
  the highlights in each style, a cell's and a chip's hover, the keyboard cursor. **A track's well is a surface, so a
  squircle** (`render/wells.ts`, `trackShape`): the rule over the well's height, round at both ends (square below the
  minimum, or `square`), filled with `seq-track` in the background (a stack's later rows in their first strip's data
  pass), and each track's data is clipped to it with the same path, one clip a track and strip (a block's strip in a row
  square where it stops short of the well's ends); the data inside keeps its own shape (cells square: a letters or grid
  track loses a sliver at its corners). **A grid over blocks** has one well, the plot's (its cells together, grown by
  its padding): each cell's clip is square on an edge inside the plot (`cutTop`, `cutBottom`, and a side short of the
  plot's ends), so only the plot's outer corners round; the top row's and the bottom row's clips reach into its padding.
  A ring's well is its ground (a sector, the band grown radially by its padding), not a squircle; text draws no well;
  both clip their data to the well's box (the band grown by its padding). A merged run is round only where its first and
  last items really end. An end cut by the view, a piece of the axis or a piece of the item's own image (an exon) is
  square; a mark cut at both ends is square. The squircle is one curve, `squirclePath` (each corner a cubic on |x|⁴ +
  |y|⁴ = 1), that every backend draws. Outlines are concentric.

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

- **The ranges layer's modules**: `render/ranges.ts` (the layer type: its prepared form, lanes and marks, hit, legend)
  and beside it, a concern each, `ranges-lay.ts` (what it keeps), `ranges-content.ts` (its content: Track heights),
  `ranges-beside.ts` (labels beside), `ranges-draw.ts`, `ranges-density.ts` (the heat strip) and `ranges-theme.ts` (its
  channels).

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
  its protein) and its parts on the section's system as `hull` (its codon); its hover outlines each part. **Molecules**
  (`render/letters-molecule.ts`, `MoleculeLettersType`: the letters wrapped, registered under `letters` in the core
  letters' place by the default compositions, `BuiltInRenderRegistry`; `renderSvg` loads it as a chunk of its own only
  for a spec with molecules, complexes or modified residues; `BuiltInLayerTypes` and a minimal composition keep the core
  letters, so a view of plain letters never reads it): on a molecule's residues (or a unit of one) a cell a residue
  labeled by its component code, on an atoms system a cell an atom labeled by its element symbol, where the label fits
  at the letter sizes (no letters are read: the labels are the model's); the cells' fill `seq-sequence` or a uniform
  `color`. **Modified positions** (`modified`, on the section's own letters or through a one-to-one path; the data
  system's modified index, alternatives one mark): `mark` (the default) a dot at the cell's top right in its letter's
  ink (`letters.modifiedMark`), a bar over a component of several positions (CRO), a tick at the strip's top below the
  cells; `code` the component's code in place of the letter where the cell has room for three letters (the mark
  elsewhere, EC-117); `tint` the cell outlined in `seq-modified`; `none`. No value changes a cell's width; every
  modified position in view has its mark at every step. A modification's color is never the residue scheme's.
- `ranges` **over links** (a complex's or a molecule's `pairs` child, `data` naming it): each end on the layer's system
  (a block's unit) a mark at its site, its label the other end (`→ Glycan 1`), its category the link's kind: the entry's
  `ends` (`LinkEnds`), an inline ranges plan a system (`<links>#<system>`), built when the view's `layerOf` places such
  a layer on that system and not before (the store and the export open the plans built). Why lazy: most links are only
  arcs, and a plan encodes a table. Placed on each section or block with an end (`seq.view.layer.unplaced` where none).
  A data `pairs` node has no ends plan: its ranges reader is not built.
- `ranges`: items in sticky lanes, carried onto the section along the layer's path (a domain on its gene: a piece an
  exon); labels inside marks 24 px or wider (faded past the mark), or inside a narrower one when they fit whole with 3
  px a side (a modified residue's code over its one position; none below 8 px), drawn in the `labels` pass. At `bars` a
  lane's marks touching within a device pixel merge into a run; at `density`, or a density table, one heat strip, a
  device column a pixel. Items past the lanes collapse into a density row in the last lane, the rest drawn as items
  (`ranges.overflow: 'extra-lanes'`), or the whole layer becomes the strip (`'layer'`). In a band under a lane's height
  (a text line's underline) the strip uses only the light `ThinStripShare` (0.3) of the density ramp, so it reads as a
  band, not as a dark rule.
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
  `seq.render.numbering` and `ordinal` is used; a source's draws no tick until its labels are in, Numberings from a
  source). A tick at each label that is a multiple of a step of 1, 2 or 5 × 10ⁿ, the least whose ticks are
  `ruler.labelGap` px apart; a label that does not fit (between the data box's edges and the labels already drawn, 6 px
  apart at least, majors first) is skipped, never shrunk. Majors are the step's round values only: the step never counts
  an insertion, so it holds across them, and a round value the numbering lacks (an IMGT gap) has no tick or label, no
  insertion standing in for it. A label with an insertion code (`100A`, `111.1`) is a minor tick from 4 px a position
  (`MinorScale`), labeled where it clears the labels drawn and the data box (one rule for the ruler: the order decides
  who wins), majors first, then the minors nearest a drawn major first, then the shorter, then by position: a number
  between every pair of majors where one fits, the insertion beside a round value (IMGT's 112.1 beside 112) before its
  neighbour, a label wider than its own position too; decided at layout, never on a hover, so no jitter, and the hit
  card names the rest (the author, 2026-11-20, replacing "within its own position's width", which left the comb of ticks
  unnumbered at a mid zoom). A ruler track is unlabelled unless its spec names it; with more than one ruler over a view
  (two chains' numberings, an against view) the spec names each (its chain, its scheme) and the name stays drawn: it is
  the only sign whose numbers are which (intended, never a label the viewer adds). It labels the positions where they
  are: a reversed axis counts down, a collapsed one labels its pieces. Where a text section prints a group's end number
  (`group_labels` not `none`) the ruler draws the tick without its label: a position carries one label, the group's
  (only where the ruler's numbering is the one the group labels print). Over a piece of more than 4,096 positions, or on
  a built-in numbering, the labels are looked up by number (a tick a step), not a position at a time. **Units**
  (`units`): `none` labels the number as its numbering does; `bases` labels an `ordinal` ruler in bp, kb or Mb chosen by
  the step (`@molstar/core/util/number`'s `formatGenomic`: 7,668,402 at a base step, 7,668 kb at a kilobase one, 7.67 Mb
  from 10 kb), one unit and one count of decimals a step; another numbering keeps its own labels. Why: the writer of a
  chromosome view says so, never inferred from the alphabet. **Where the numbering does not reach** (`unnumbered`;
  decided, the author, 2026-11-19, the Antibody's light chain: Kabat and IMGT number only the variable domain): a
  position past the numbering's labelled extent (before its first label or after its last, found once a numbering, of
  any length of piece) is in an unnumbered run; a position with no label between two labelled ones stays blank. `mark`
  (the default): the ruler's line over the run is muted and dashed, one muted note at its start (`not numbered`,
  ` by <numbering>` added where the ruler's track label does not name the numbering) where it fits clear of the labels,
  no ticks. `fallback`: the run is labelled by `fallbackNumbering` (`ordinal` when absent; a numbering the system
  declares, built in or inline) in the muted style (ticks and labels at a lower opacity), no line and no note. Under
  `fallback`, a `fallbackNumbering` the system lacks is `seq.render.numbering.fallback` and `mark` is used. Why: a blank
  reads as a bug, and a fallback numbering is the writer's choice, never a silent default (its numbers are not the
  scheme's). Canvas and SVG alike (one draw list). A hit is the position under the pointer; the hover card and the live
  region say its position and `not numbered by <id>` there. Fallback: none (omitted).
- `secondary`: the items of a `ranges` layer (the class from the `kind` field, else `category`; an inline table cannot
  name a column `kind`, the builder's key for a node's kind and reserved by the S-expression codec, so it writes
  `category`, while a source table's `kind` column reads: DSSP codes or words: helix, strand, UniProt's `Helix`,
  `Beta strand` and `Turn`, anything else a turn), in one band: a helix a block, a strand a block with an arrowhead on
  its C-terminal end, a turn a thin bar. **The arrowhead** is on the piece holding the C-terminal end of the item, found
  among all its rows (an inline annotation's; a tiled source's rows may lie outside what was read, so its items have no
  head: a head on the wrong piece is worse than none) (the end of its last row in reading order: across the origin of a
  cyclic system, the row that starts at it), at the screen side the residues run to (the axis's direction against the
  layer's path's), and not where the view clips that end (EC-130: reversed, across the origin, through a minus-strand
  codon map). An end cut by the view, a piece of the axis or of the item's image is square. The class is the whole code
  or word (DSSP `H G I`, `E B`, `T S`; mmCIF `HELX_*`, `STRN`, `TURN_*`; helix, strand, turn). **Items in a hidden
  stretch** are the break's tick, as `ranges`'. **Copies** (rows of one extent and class) are drawn once, the first hit
  (EC-67). From the coarse level the source's density bins are one strip. Fallback: `ranges`.
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
  it, so the line stops at it; the hit boxes are the marks' lanes all the same. Whether an item has its arrowhead is its
  own mark's width, at the items and the bars steps alike: 4 px or wider, fixed (narrower, the head would be all of the
  mark); an arrow is never merged into a run of bars, and the density step has none.

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

### Layer colors

`render/layer-color.ts`: a layer's color channel as its track compiles it, one path for the layer types and any other
reader (the structure translation's `colorRuns`). The letters, values and ranges types' defaults (the alphabet's scheme,
`token:value`, a colourless ranges layer's categorical place by its key's first appearance in the view, else
`token:range`) and the view's theme context (`viewTheme`: its scales, the held domains, the registry's schemes and
measures, the reader's substitution) are defined there and read by `layoutView` and the types.
`layerColors(view, theme)` compiles a layer's channel through its kept channels with any tokens;
`Renderer.layerColors()` is it with the renderer's theme, held domains (another renderer's when it reads those) and its
last frame's substitution. Why one path: a structure colored by a track agrees with the track, an open scale's domain
the one the view holds, not the data's extent.

### Numberings from a source

`render/numberings.ts`: a declared numbering whose labels are a source's rows (Language, Numberings).

- **Read whole**: the sequence's every position at level 1 (`readSourceNumbering`, internal), a row a position (its
  `start`: a `position` row is `[p, p + 1)`), the label its `field` (a category's text, a number's spelling; `id` the
  row's id; a row with no label is a position with none); the first row of a position keeps it. Held once every tile is
  in, never in part. Why whole: labels need no order, so a part would snap and print other positions than the whole
  does.
- **Bounded**: a sequence longer than `maxNumberingPositions` (2^18) is `seq.render.numbering.limit`, and its tiles are
  never asked for nor its labels held. Why: the labels are held whole, a string a position.
- **Kept**: held labels stay while the source's coverage is the same, whatever becomes of their tiles: an evicted tile
  is never asked for again for them, so a ruler does not ask for every tile each frame. Before they are held, a read is
  kept while its tiles are the same entries.
- **Blank, said**: a tile that failed, or a source whose coverage has no system of the sequence's name (or no length for
  it), is `seq.render.numbering.source`; held labels with rows past the sequence's end (left out) or rows for another
  count of positions than the sequence has are `seq.render.numbering.rows`, each once a read (the ruler reports them at
  layout). A blank ruler labels by `ordinal`, as for an unknown numbering.
- **Pending is waiting**: the query model (`withSourceNumberings(base, model, store, { maxPositions, want })`) answers a
  held numbering, says one not in yet is pending (`QuerySystem.numberingPending`: a form naming it gives nothing, the
  whole system pending, `seq.query.numbering.pending`, so the address is `waiting`, never unknown), and leaves a blank
  one out of the ids an unknown numbering lists. `want` asks for the tiles of what a form named: whoever resolves asks
  (the viewer's addresses), so a link in a numbering no ruler shows resolves. A unit's molecule is wrapped too, so a
  molecule's source numbering numbers its units. The static model alone says a source's numbering is pending: nothing
  there brings it.
- **Who asks**: the ruler at layout (no ticks while they come, never ordinals in their place) and the viewer's
  addresses; a text's numbers read inline numberings only (a section layout reads no store).

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
  ring; a letter's center warped, 0.35 em above its baseline, the baseline put back below it in the letter's own turned
  frame: a warped baseline puts an upright letter's center off its cell, along the ring), an image refused (its layout
  reports it), a clip applied in the strip's frame before warping (a clip by path would change the lists' form for one
  case). Nothing allocated once warm.
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
  (with blocks) starts a line, a one-position block flows as a chip (Blocks).
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
- Blocks: each starts a line, chips flow like words (Blocks).

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
  `bands` track filling it, `rows` the section's first rows track with no header (no tree, numbers or row columns: a
  narrow picture would hide the columns and draw their chip) filling what the other tracks leave (an overview's whole
  column), its rows fitted to the band (`rowScale` the band's height over the rows shown, the view's folds); a picture
  that fills is made again for another height. The same system with no `system` param: the bound section's own axis.
  Colors are the view's: the picture reads the view's renderer's held domains (`domains`, its `layoutState`, a picture
  key mapped to the view's), never a domain of its own over the whole axis, so a value has one color in both. Why the
  layout interface: a miniature drawn by other code drifts from the view it summarizes.
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
- **The zoom strip** (`zoomStrip`, `DefaultZoomStrip`: `minPx` 48, `heightPx` 16; the viewer's `navigator.zoomStrip`, an
  export's `navigator` option): a positions navigator along the view (`top`, `bottom`: a side part's length is the
  view's height, not its axis; none there) whose reachable span is narrower than `minPx` on the picture is context and a
  zoom strip, one navigator (one renderer, hit, frame and draw list set). The span: the hull of its bound sections'
  pieces, on a cyclic axis the shorter arc covering them (over the origin its end past the total); a bound section on
  the navigator's own axis is the whole axis; known limit: two far-apart regions make a wide hull, so no strip. The
  picture's view gets a second section, the same form over the same system (an ideogram's bands; else the axis, a
  ruler), its viewport the span's units (over a cyclic origin past the total: its pieces twice, the strip runs on over
  the origin), laid out below the box at the box's bottom whatever the context's content (`zoomHeight`, exactly
  `heightPx`: a gap of at most `ZoomGapPx`, 6, and at most half, then the strip's track; `heightPx: 0` is off; the
  frame's `size` is the box, the lists the box and the strip; an export stacks both). The windows, their wash, floor,
  grips, marks, hits and `sectionUnitAt` are the zoom strip's (a drag moves the view at its scale); the context draws
  the span as a mark in the `window` pass (an outline at least `window.minPx` wide, always centred on the truth whatever
  `window.floor`, so the truth is inside it, its hairline at the truth; over a cyclic origin two pieces; no grips, no
  wash) and two connector hairlines (`seq-window`) from its ends to the strip's; a hit on the context is `context`
  (`mark` over it, else `outside`), part `outside`, no target; the gap below the context is the strip's (a finger a few
  px above it lands on it); `markText` names the span (over the origin the arc, `chr:59,001–1,000`). On under `minPx`,
  off past `minPx + ZoomHysteresisPx` (8: a width about the threshold flickers not); made again when it comes or goes;
  the span read each frame (none, a bound section pending: the last kept, mark and connectors with it); the span moving
  resets the strip's viewport. The tokens last set are the renderer's own, so a picture made again keeps them. Why: a
  region on its chromosome is a few kb of tens of Mb, a window there can neither show a move nor take a drag; the
  chromosome's context and a working window are two scales.
- **Draw lists**: the base (`background`, `data`) is the picture's, keyed as any frame; the `window` pass draws the
  washes, the windows (the one the pointer last acted on on top) and a sweep's outline; then the marks: the selection's
  and the highlights' ticks at the picture's foot (`navigator.marks`), the keyboard cursor's hairline, the pointer's
  hairline. A navigator never draws in a section's lists.
- **Hits** (`hit`): the window part under a point, its section and axis, the position on the navigator's system, the
  display row, a band's name (the `bands` layer's hit `state`). A navigator's hover is its own: a readout.
- **The ruler rail** (`Rail`, `FrameInput.rails` and `railsActive`, by section address, `key`; the viewer's): drawn just
  below the first shown ruler's bottom border, in a row of its own (`SlimRailPx`) its well grows by, never over the
  ruler's ticks or labels (decided, the author, 2026-11-19: "Below the border"); the row only where a rail is drawn live
  (the canvas with rails on; with the rail toggle off none), never in an export ("Leave it out of exports": a figure has
  no empty strip, and its included rails lie on the border), else in a slim row of its own under the section's tracks
  (`SlimRailPx`, 12, fixed: the form's size, the least a rail and its grips need; a hidden ruler too; under, not above:
  the tracks keep their place); the whole axis across the data box, the window a segment with a grip each end (over a
  cyclic origin two pieces, a grip at each true end, the axis's ends body), collapsed breaks a notch; its grab zone the
  ruler's height, up from the rail row's foot (on an alignment's columns its lower `grabShare`: above it a drag selects
  columns), the slim row whole; grown (`grownPx`) while hovered or dragged (`Emphasis.drag`); idle (not active) a quiet
  hairline with no grips and no target (at the fit it reads as the line); in a figure (an export) a plain segment, no
  grips, no idle state. Its hit is the `window` target with `navigator: 'rail'` (touch the wider hit: `hitView`'s
  `touch`), its unit `railUnitAt`. Rails are part of a section's key by identity (a slim row takes room); which can move
  is not: the renderer sets it on the rails each frame, an emphasis redraw too (an idle flip redraws the window pass
  alone). `railSections(view, hosted)`: every tracks section no shown navigator binds along its positions: a declared
  navigator and a host's replace its rail, a hidden one brings it back.
- **The rows rail** (`placeRowsRail`, `SectionLayout.rowsRails`; `FrameInput.rails` by the rows track's address): the
  same `Rail` upright (`vertical`: along the page's y, its line and grab zone across on x; one geometry, drawn and hit
  by the same code, the coordinates turned at the rect and at the hit), `SlimRailPx` wide at the right edge (right of
  the last numbers' gutter: the gutter widened by it while any rows rail is on, so the data keeps its width idle or
  live), along the band's scrolled area below its pinned rows, in display rows below the pinned band; its window the
  rows in view. Idle only while every row is shown and the rows cannot zoom in (navigators.md's idle rule: a rail zoomed
  out to the whole band keeps its grips, so a grip drags it smaller again) (`whole`, the layout's: every row shown and
  still at `MaxRowScale`, 40 px a row, the rows' zoom-in limit with `MinRowScale`, 1/64, their least; the renderer sets
  `active` only when the viewer says it can move and not `whole`). Its hit is the `window` target with
  `navigator: 'rows-rail'` (`RowsRail`), its `layer` the track, `position` the display row.
  `rowsRailTracks(view, hosted, mode)`: each tracks section's tracks with a rows layer shown; by `mode` (the viewer's
  `navigator.rowsRail`) drawn beside a shown navigator over the section's rows too (`always`, the default), replaced by
  one (`unlessOverview`: a hidden one brings it back), or never (`none`). Why: the design's "simple slider navigator for
  vertical scrolling", the rail a linear section has, so a tall alignment can be scrolled and zoomed with no declared
  part; kept beside an overview because the overview is a picture and the rail the precise control (the replace rule
  holds on the positions axis alone, where the rail and a strip are two pictures of one axis in one place).
- **A narrow frame** (`render/narrow.ts`: `narrowColumns`, `overviewWidth`; `ViewLayout.narrow`): what gives way when
  the plot (the data's width) runs short, keyed on the plot's width, never the frame's, each rule with a floor (the
  ruler's: never a column squeezed to nothing), cheapest first: (1) a vertical overview at its form's size beside the
  canvas narrows from 96 toward `navigator.overview.minPx` (48) to keep the plot `navigator.overview.minPlotPx` (160)
  wide; (2) below that it is hidden, the rows rail kept (the way to move is never lost); (3) the label column shrinks to
  its floor (`labelColumn.min`, 80, or its own width if less), then the tree collapses to its toggle
  (`LayoutContext.treePx`: 0; `RowsLay.toggleX`: a small tree in the branches' ink and width at the band's top, the
  labels after it, hit as `tree-toggle`; a band's `tree`, the reader's opening, session state, keeps the tree whole:
  `treeKept`, the least width then holding it); (4) the plot keeps its floor and the frame scrolls: `narrow.minWidth` is
  the least frame width with the columns at their floors. `narrow.plotPx` and `.columnsPx` are at the columns' own
  widths (row columns ungranted too, whatever the frame), so a navigator reading them beside the canvas never flips with
  the columns it caused to shrink. A view with a tracks section has it. Why the order: a miniature before names, names
  before data.

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
  sink. **A dense fan is one label** (`rings.cluster`: more than `over` of one layer's labels within `degrees`, found
  sweeping from the layer's largest angular gap, so a run across the origin is one): the grouped label ("n items", at
  the run's middle angle) is placed, and so is the fan of its members as if no group existed (`placement.fan`), both
  with the layout. The grouped label hits as the group (`target: 'group'`, `members`: its card lists them, a click
  selects them all); a hover on it or on a member draws the fan over the labels in the `tip` pass (the grouped label
  covered, each member on a ground of its own) and the fan's labels hit as their items while it is open; an export draws
  the grouped form. Why precomputed: the layout never runs on a hover, and nothing else moves. A fan opens for a group
  of up to `rings.cluster.fanMax` members (12) and pivots about the group's middle angle: its members a line apart at
  the bend (`labelLine / (r + labelGap)` radians, or their own spacing if wider), in their order along the ring, each
  leader from its member's own angle; a larger group opens no fan, its card lists every member (no "more": the card is
  its only opening). Why: a dozen leaders over the group's own few degrees run together; past a dozen lines a fan cannot
  be read.
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
- **A one-row image is a ring band**: a density strip (a ranges layer past its lanes, its overflow row, a variants
  count) is drawn radially, each run of one color an annular sector through the warp, so a ring never refuses its
  density. Why: the strip's columns are the frame's, so the band is exactly what a track shows. **A taller image is
  refused** (a matrix: the warp counts it; `seq.render.ring.image`, said once a section by the frame).
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
- **A chip** (`isChip`): a block of form `whole` (a molecule's default: its every residue in one slot) or of one
  position. A chip is **one unit** of its section however many residues it stands for (`blockUnits`; layouts.md): its
  strip shows its first position's cell, its pill and its hit stand for every position (the hit's `extent`), and a
  pair's end on it is drawn at its middle. Why one unit: the window, the rail and the navigators skip a chip as a gap,
  and 200 residues of a glycan never take the polymers' share.
- **A one-position block is a chip** (`isChip`, `chipWidth`, `chipSlots`, `Strips.chip`): a fixed slot its label's width
  (its label with its `×N` in the label font, `BlockMetrics.chipPad` 5 px either side), never a position's share: zoom
  does not scale it, so it reads at the fit view. Why: a one-position share is a sliver at the fit view (the Complex's
  ATP) and the floor (a 2 px mark, unlabeled) says nothing. **A floor**: the chips take at most half the width
  (`chipShare`; a ring: half its circle) together; past it they shrink by one factor, never below `chipMin` (16 px, less
  only when even that would pass the half), their labels elided; the polymers keep the rest.
- **A row** (`BlockRow`): each block's edges in view on the device grid (`snap`, a chip's kept: its slot is its label's
  width), as a grid's row under it, so the two align at any dpr; a block left none wide has no strip (nothing of it
  drawn or hit; a hit beside it names the block next to it). A chip stands where its block does (a ligand after the
  chains). The chips in view take their px from the width as the gaps do, so the row's x is not linear over a chip
  (EC-102); the window, the rail and the navigators read the section's units, a chip one unit, as across a gap. The row
  is continuous in its window: a chip partly in view takes px by the part of it in view, a gap by its boundary's
  distance to the nearer edge (none at the edge, whole 12 px in), so a pan moves nothing by a jump. **Pan and zoom go
  through the row's inverse** (`rowPan`, `rowZoom`, by bisection on the window's start): a pan brings the unit `dx` px
  from a reference point (the middle less half the move, moved past a gap) to it, so what is under the pointer moves
  with it while nothing enters or leaves; a zoom keeps its anchor's px as chips and gaps enter or leave. Why: a linear
  `dx · span / width` jumps over a chip.
- **A stack**: the chips share one row after the other blocks', side by side from the first position's x, whatever the
  window. A row's ground and its track labels are drawn once, at its first strip (a later chip's would paint over the
  pill before it).
- **The pill**: a letters layer on a chip draws its label in a pill over its cell (`LayoutContext.chip`, by the rule of
  a label inside a mark, `longestCut`: whole, or cut keeping at least `labels.minChars` letters, else the pill without
  text, its label in the card; the label row above follows the same rule): the cell's fill with the letter cells'
  hairline (`seq-cell-edge`, dark), its ink the tokens' when it holds 4.5:1 against the fill, else black or white; one
  pill over a chip's cells. A molecule's cells take the layer's `color` once a block, its theme reading the molecule's
  fields (`class`, `name`, `component`: its first residue's code), so one letters layer colors each chip on its own
  (`(rule (eq component ATP) …)`); a rule that holds for none gives the layer's default. A chip is hit as its whole
  block (the hit's `extent` its every position, made once a layout: a hover over it is one hover), so a click selects
  every residue and, through the atoms map, every atom. The label row above a chip is empty only where a pill was drawn
  (a letter laid out on a band with height); else (no letters track, a pending tile, bars, a collapsed band) the label
  is there. Another layer draws its mark across the chip where it has something at the position, else nothing.
- **Labels**: a block's label above its tracks is elided to its strip's width, dropped when nothing of it fits: never
  over the next block's. Hit, hover, selection and a link read a chip as its one position; `place` at a unit where one
  strip ends and the next starts gives the next (a chip at its own first unit is the chip). Contacts (where a ligand
  binds) are a `pairs` layer a builder or an application adds, never drawn by default.
- **Text**: every block starts a line on its own system (its lines from its first unit, `blockStarts`); a run of chips
  flows like words, as many on a line as fit, each its slot (at most a line) with a word's gap (`text.groupGap`)
  between. The placement reads a line's units from its first section unit (`StripPlacement.lineBlocks`), a line's chips
  by x. A block's label is inline on its first line (EC-95): its width and a word's gap taken from that line
  (`TextBlocks.lead`, the line's units fewer, its letters after the label), so it takes part in the breaking; its
  continuation lines have none. A chip's label is its pill.
- **Rings**: a block an arc (above); a chip a fixed-angle arc its slot on the outer ring (the frame's px), its letters
  ring a pill, the other blocks sharing the rest of the circle by their lengths.
- **Follow across systems** is refused: a section of blocks on more than one system neither follows nor is followed
  (`seq.render.follow.blocks`, it shows its own range). Why: which block's map a follow would carry is undecided. Blocks
  on one system follow as one axis.

### Tracks layout

- **Its modules**: `render/tracks.ts` (the layout pass, `layoutView`, and the public names) and beside it, a concern
  each, `tracks-types.ts` (its input, its layouts, `LayoutState`), `tracks-axes.ts` (a section's axis, a follower's),
  `tracks-cache.ts` (a section's reuse), `tracks-strips.ts`, `tracks-columns.ts` (the label column),
  `tracks-templates.ts`, `tracks-heights.ts` (Track heights), `tracks-headers.ts` (Track headers: the facts and the
  lines), `tracks-draw.ts`, `tracks-emphasis.ts` (Emphasis) and `tracks-hit.ts`; the header line's own modules are
  `headers.ts` (its parts and geometry) and `headers-draw.ts` (drawn and hit).
- **A cyclic axis's window past its total** (a linear view panned across its origin, `wraps`): the transform takes the
  whole system twice, so the window's units go on over the origin (two runs, no break between); a follower of it wraps
  as of a ring's window.
- **On the interface** (`layout/tracks.ts`): one strip a section (a block in view with blocks: Blocks), translated by
  nothing (its frame the page), its bands the tracks stacked with `trackGap`; its chrome the label column (each track's
  label and header line, the section's header line, the direction mark: Track headers). Its draw lists are the same
  before and after the restatement but for the `window` pass's header: the committed records are its guard.
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
  it brings arrive; that level and a nested source's (below) are bounded as a layer's is (coarser while the view would
  want more than `maxTilesPerWant` tiles), so a long nested system at the fit is a coarse tile, never refused. The
  overlay's issues are the layout's. A template below its source also wants each nested source's systems (`systemNamed`)
  in view, carried along the overlay's path from the section, so a nested source that sends its nodes by region is asked
  for them.
- **Shown nodes** (`nodes.holdShown`, true): the layout lists the nodes it shows (`ViewLayout.shown`: the data and the
  path systems of each laid-out layer, through `overlay.holds`), for its owner to hold in the store; with `false`, none.
  Why: a node stays while it is shown, not only while the region that brought it is cached. **Shown means in view**: a
  template's instance is laid out only while its node's region is in view: the node's top node (the one of the
  template's source it came under) is listed up front in that source's coverage, or is carried by a loaded tile of that
  source over the view at any level whose tiles there are within the store's bound; while the template's tiles over the
  view are pending (a zoom's new level), an instance shown last stays. Out of view its track goes and its hold with it;
  the node is evicted as any other when the store needs room. Why: an instance held its node, so every place browsed
  left its track (29 empty after five places on a chromosome); a track whose node's region is on screen never vanishes,
  a zoom included.
- Cell edges snap to the device grid (`snap`). One label column for every track of a view (`labelColumn`, measured over
  the view's own track labels: a template's instances, named once expanded, and a direction mark fade or end in an
  ellipsis at its end, so a reverse or an arriving region never moves the data box); a track's label is centred on its
  well (its padding with it), in tracks and on a text's first line.
- **A track's padding** (its `padding`, else `trackPadding`): px inside its well above and below its data, never left or
  right. The well (`TrackLayout.well`) is the track's height plus its padding; its data's band (`TrackLayout.box`, the
  strips' `bandY`, `bandH`: what every layer lays out, hits and emphasizes in) is inset by it, and `Strips.padT`, `padB`
  say by how much a band (`PlaceContext.pads`, set by each layout). A point in the padding is in no band: no item. Rows
  and blocks stack wells with `trackGap` between them; text lines hold them (a collapsed band drops its padding too); a
  ring's padding is radial, scaled with its rings; a grid's plot is inset by its bands' largest padding.

### Track headers

- **Parts given, geometry laid out.** A header is a track's line in the label column, or a section's line above its
  first track. Its parts come from `LayoutInput.headers` (`HeaderPartsOf`: facts in, parts out, once a layout a header;
  part of the key by identity, `FrameInput.headers`, `RenderSvgOptions.headers`): `slots` (an entry's name, label, icon
  as path data in a 24-unit box, enabled), `states` (name, text, count, the address it is about), `menu` (any entry
  applies) and `hosted` (a host's DOM parts draw them over the laid-out boxes). The facts (`HeaderInputFacts`; the
  viewer's own facts are its entries') are DOM-free: the layers a filter thinned with their counts, the hidden layers,
  the section's hidden tracks read from the spec (each once, in spec order). Absent, `builtInHeaderStates`: a track's
  `seq:filtered` ("n filtered", every thinned layer's rows) and `seq:hidden-layers` ("1 layer hidden", "n layers
  hidden"), a section's `seq:hidden-tracks` ("n hidden"); no slots, no menu. Why given: which entries apply and whether
  they run is the viewer's (its registry and commands); where they go is the layout's, so a canvas and a host's DOM
  parts agree. A navigator's picture lays out no header.
- **Where** (`header.*`): `actionColumns`: the "⋯" at the label column's right end, `header.slots` columns of
  `header.slotPx` left of it, as many slots as leave three columns of label (an 80 px column holds one), none under
  `header.actions.show: never`; one set a view, counted from the "⋯" leftwards (the first slot next to it), so a slot's
  x is never the label's and every track's eye is in one column whatever its line's count. Why the "⋯" last: a kebab
  menu is looked for at the row's end (Mol\*'s order; the author's). The slots are the first given
  (`HeaderSlotBox.part`); the "⋯" when `menu`. A line reserves only the slots it has (`HeaderLine.actions`: its slots
  and the "⋯"): a line with none (a section's, whose built-ins are menu only) has its chips up to the "⋯". Why: at the
  default 112 px a reservation of every column left "1 hidden" a mark. A track under `header.minHeightPx` (`short`)
  takes no slots and its states are marks (6 px accent squircles, their words the card's). States beside the label end
  before the line's own actions when it has any, else at the column's end, the last nearest; a state with no room left
  is a mark, so a slot never overlaps a state: **a narrow column draws its states as marks** (decided: an 80 px column
  with actions holds one slot, its chips marks), and a mark with no room left in the column is left out (never left of
  it; its words the card's). On a track two lines high the states take a line of their own under the label (the label
  whole, the two lines centred together; lower by the slot's half when the header has actions, so the line clears them).
  The label fades before what is beside it (`HeaderLine.label`). **A rows track** draws no states and no slots and gives
  no label hit: its label column is its tree and row labels, hit as before; its header is its title row's (not built,
  the TODO).
- **The section's line** (`section.header`: `auto` when its parts have a state or a slot, `always`, `none`):
  `header.sectionPx` high above its first track, the tracks moved down by it (`SectionLayout.header`); kept when every
  track of the section is hidden. Why `auto`: no line that says nothing, records without hidden tracks unchanged.
- **Drawn**: the states in the labels pass (always: facts about the picture; a stack of blocks draws them once, on a
  track's first row, where they are hit), the actions in the tip pass over a backdrop of the label column's background
  (`msui-bg`: the column has no well) covering the label's tail (the label never moves): each slot's icon stroked in the
  muted ink (an entry not enabled in the line ink, the hovered one in the accent; no icon: its label's first letter),
  the "⋯" three dots. Actions show for `always`, or on `hover` for the track the hover's hit is in (its label line, its
  slots, its states, its data; a section's line only its own hits) and the cursor's header (`cursor.part: 'header'`,
  `cursor.track`): an emphasis redraw, no layout, nothing allocated (`headers.quick.ts`). No emphasis field of its own:
  the hover drives it. A header's hover lights none of its track's items (EC-146: refused by default). A figure
  (`figure`) lays out no slots and no "⋯"; its states as drawn, none with `header.states.export: none`. A `hosted` line
  (a host's DOM parts) draws nothing in its slot, "⋯" and chip boxes (the parts are the DOM's); a figure of it draws its
  states.
- **Hit** (`hitChrome(section, x, y, hover)`, before the bands): a slot or the "⋯" `header-slot` **only while its
  track's actions show** (by the hover `hitView` is given: the label line hovered first, then the slot), a state
  `header-state` (`count` its count), elsewhere on a track's line `header-label` (no `state`), on a section's line
  `section-header`; `state` the entry's name (`HeaderMenuName` the "⋯"), `layer` the track's address (a section's for
  its line), `section` its section. Each hit is laid out with its part (`hit`, a line's `HeaderLine.hit`): a hover makes
  no object.

### Track heights

- **A track with no `height` fits its content, in steps** (`render/sizing.ts`): a layer type's `content(prep, ctx)` says
  what its items need over the whole data, never the view (ranges: the lanes first-fit by start packs, and whether a
  label fits none of its marks at the fit: the section's own whole extent, its pieces once, across the data box, never
  the frame's axis, which a pan over a cyclic origin or a follower's range changes; a zoom only widens marks), so a pan
  never resizes a track. A track's content is one lane block for all its layers: the most lanes any of them needs and
  one label lane when any of them labels beside; its height is those lanes, each `laneHeight`, plus `lanePadding` above
  and below, between one lane and `maxLanes` (22, 36, 50, 64, 78), at least its other shown layers' fixed heights (a
  hidden layer adds none, nor names its kind); past `maxLanes` the rest go to the density row (`ranges.overflow`).
  Values, letters, rulers, pairs and rows keep their fixed heights (no lanes). A source's items (a file, a tiled genome)
  are not known whole: its track keeps its fixed height and its lanes are centred by the most the layer has shown (it
  only grows: its marks never move back on a pan). Why: a default fills its space, and heights decided by the whole data
  stay put as the view moves.
- **The section's `heights` names the sizing**: `content` (the default: tracks of one kind, their shown layers' types in
  order, share the tallest lanes, each with its own label lane, `SizedTrack.label`: a Reverse strand whose labels all
  fit inside reserves no empty lane beside a Forward one's), `uniform` (every content-sized track the section's
  tallest), `fixed` (each its types' fixed height, `heights` above), or an `ext:` sizing (`TrackSizing`, a
  `track-sizing` registry entry; one not registered is `seq.render.heights` and the tracks fit their content;
  `ContentSizing` is public for one to wrap). A track's `height` overrides any. **The host's `trackSizing` port**
  (`TrackSizingHook`: the section, its `SizedTrack`s with their content heights, the named answer) has the last word; an
  answer of another length is ignored. Text keeps its own thin bands. A navigator's picture sizes its tracks with the
  view's sizings and port, and has no label lane.
- **Kept while their key is**: a section's heights and its tracks' contents are kept in the layout state by the section
  (its shown tracks, the px a position of its extent, the config, the layout, the metrics' version, the registry, the
  port), so a pan runs neither the sizing nor the port; a layer's content is kept by its prep (its label check by the
  fit, the font and the metrics' version).
- **Marks centred, never below the band**: the band's lanes below its padding, the label lane above them, the block
  centred in the band by the whole data's lanes (`RangesLay.laneY`, `labelY`), the marks in those lanes only (a band
  taller than its content, by `uniform`, kind or its own height, draws none below them). A sticky lane that would leave
  an item out of lanes enough for the whole data is packed fresh before any item goes to the density row. A band too
  short for its lanes (a track's own height, a text line's underline, a navigator's miniature) drops the label lane,
  then compacts its lanes over the whole band.

### Labels inside if they fit, else beside

- One rule in both layouts (rings: the sink, Rings layout). On a linear track an item's label is drawn inside the first
  of its marks it fits whole in (the mark's part in the box, less its arrowheads and a 3 px inset each side: a mark cut
  at the section's edge keeps its label within the box, never clipped mid-glyph); else in the track's label lane at its
  leftmost mark's x, moved right clear of a label placed before it with a short leader from its own mark (at most its
  own width, or 24 px; an item in a lower lane whose start is under a mark above rises from the first x of it with no
  mark above, `LabelLane.clearX`, else is not moved there); **a leader never crosses another label or a mark above its
  item** (the label lane's labels, the marks of this and the earlier layers: `LabelLane.crosses`), so a label that would
  need one goes on, as does one the lane cannot hold, or one on a track with no label lane: **inside cut** with an
  ellipsis in the first of its marks that holds at least `labels.minChars` letters (4) and the ellipsis, the most
  letters that fit (a label of `minChars` letters or fewer only whole), else dropped (the hover card names it, whole;
  the address too). The canvas and the SVG cut alike (one layout, the text-metrics port); a cut is made once
  (`range-labels.ts`). Why: a label vanishing whole from a box a step of zoom narrower reads as lost data (the author).
  **A run is labeled once**: items of one layer with the same label, each less than the label's width from the one
  before it in x, get one label at the run's start and none of their own (`firstOfRuns`; the label's hit is the run's
  first item; the hover card names each). **One label lane a track**: its layers' labels beside share it, placed by one
  priority, the layer's order, then the longer item. An item is labeled once. The same font and size inside and beside;
  the ink is the contrast rule's inside a mark, the text color beside it (decided). A placed label's hit is its item's.
- **Placed per layout** (the layers', in their order, allocation-free once warm): a pan in the same level lays the
  layers out again and places their labels again in view, as their marks are; never on a hover.

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
- **Size** (EC-109): x is the data box (the view's x, aligned with the sections it shares x with); the height is the
  data box's width over `aspect` (1), capped by the room, the larger of `heights.grid` (a floor, decided) and nine
  tenths of the view's height below the section's top. At least 40 px. `aspect: 'none'` takes the tallest track's
  height. **A plot's box is fixed** (decided): no zoom or pan of x or y changes it, so cells are rectangles when x and y
  show different spans. Why: a box that followed y's window resized the plot under the pointer, moved the page and the
  tracks below it, and let the wheel fall through to the page. Why the floor: a view sized to its content would
  otherwise feed its height back into the room and shrink.
- **Chrome**: the label column (the first track's label) and y's first and last positions at the plot's top and bottom,
  their strings made at place (a draw allocates nothing).
- **A grid over blocks** (its section's blocks a leader's, "Resolved view"): a strip a pair of blocks in view, x's block
  across (`strips.block`) and y's down (`strips.yBlock`, its first unit `strips.yFirst`, in the lays' key), each block
  laid as a row of blocks is (`BlockRow`, gaps between them), y over x's window and the plot's height; the strip's bands
  that pair's cell; `at` finds the strip by x and y. Its follow is the identity on the leader's units (the tracks
  layout's "several systems cannot follow" does not apply; a follower that is no grid still is refused). The y blocks'
  transforms are made once a place and shared by the x blocks' strips. A hit on it carries its x block (`Hit.block`) and
  its y block (`Hit.yBlock`). Each block's edges are on the device grid (`BlockRow.snap`: its strip's x and width, its y
  block's top and height), as the tracks section's row above is, so the columns above and the plot's align and the
  cells' clip and the tiles' edges meet exactly: no cell past its block's edge, no gap's ink inside it; a cell's rect is
  cut at its axes' edges. A chain boundary is the blocks' gap with no line in it (decided: the gap and an off-centre
  hairline together read wrong); a `seq-break` hairline a device pixel wide only where two blocks in view touch with
  less than a device pixel between them (down the plot at an x boundary, across it at a y one). Chrome: the first
  track's label at the plot's top, each y block's label at the label column's right edge at its top (below the track's
  label when it starts there).

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
  clade; `column-title` a rows layer's column title, `state` its field; `columns-hidden` its hidden-columns chip,
  `state` its text; `header-slot`, `header-state`, `header-label` and `section-header` a header's parts, `state` the
  entry's name, `layer` the track's or the section's address: Track headers): `layer` the layer's key (its address, else
  its track's and `:index`), `section` the section's id, `element` an item's `{ owner, key }` with its `extent` (a
  Location on the data's system), `position` on the section's system. Cells by arithmetic on the axis; items by
  `LaneIndex` (a lane's marks by left edge with a running maximum of right edges). In a track, the topmost layer with an
  element under the point wins, else the topmost with a position.
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
- **A complex's links** (`pairs` under a complex): each end on its row's unit (`PairsData.aUnit`, `bUnit` into the
  entry's `units`), placed on the block showing that unit, or on the block a merged copy of it is merged into (its
  `systems`); its hit's extent on the two units.
- **Ends on blocks** (`render/pairs-ends.ts`): an end is on the first block showing its system: the block of its system
  or its row's unit, the block a merged copy of it is merged into, or a block of a unit of it (a sequence or a molecule
  reaches every unit over it, EC-48, by the resolver's step down to the unit: within its range). A pair whose ends are
  on one block of one system is drawn again on each copy of that block (`ResolvedBlock.copy`, the resolver's copy
  groups: one copy group, or one system and one range and form), a caret a copy where its other end is on no block:
  split copies each show a pair within one molecule, a unit named by the data shows it alone, a zoom block (another
  range) is no copy. A pair between two systems with an end on split copies its data names no unit of is drawn on no
  copy, never a guessed one, and `seq.render.pairs.copy` (warning, once a layer) gives the pairs shown so (a filter's
  left out) and each system with its copies; one whose other end is on no block also says `seq.render.pairs.systems`.
  Why: the data says a pair within one system holds on every copy; across two it does not say which copies touch, and a
  contact drawn on the wrong copy states a fact the data does not hold (derive contacts per pair of units).
- **Bounded by the screen**: past `pairs.maxArcs`, arcs binned by their ends' device columns, one a bin (widest first),
  the rest a density strip of ends, said once (`seq.render.pairs.limit`, info). Ends on two systems that are not the
  section's own report `seq.render.pairs.systems`.
- **Feet**: on the band's edge toward the sequence (`above`: its bottom, the arcs rising; `below`: its top), inset by
  the stroke so the whole arc is in its band; with `split`, on the band's middle, each side half the band. A view puts
  the pairs track on the letters' side (the tRNA's: over its letters).
- **On an alignment's columns** (`pairs.a` an alignment, a family's `SS_cons`): the section of the alignment is the
  pairs' own system, so its arcs stand over the rows, a foot over the column each end names (the rows' hit there names
  it, at any zoom and across an image tile's edge); a hover on a rows cell lights the arcs ending in its column, and a
  pair's extent is on the columns, carried to a chain's section through its map onto a seed row (and a residue's hover
  back onto the arcs); a chain's ranges show on the columns along that map, as any layer's data on its path. Pairs are
  not projected: on another system's section they are `seq.render.pairs.systems`.
- **Dots** (with `ctx.axisY`, in a grid): a rect a pair of single positions at every zoom, at least 4 CSS px, outlined
  in the text's ink (it reads over any cell color in both themes); a diagonal segment a stretch, with its partner when
  both ends are on one strand (`a_strand`, `b_strand`), against it otherwise; `mirror` (false: a pair is drawn where its
  data puts it; the tRNA builder sets it) adds (j, i), a self-pair once; `triangle` keeps a side. Pairs at one cell
  share it side by side (EC-63).
- **Themes and filter**: its `color` channel (`channelForms`, held extents, its legend), by default `field type` (a
  categorical palette fixed by value; the builder lists an inline table's kinds, as for ranges' `category`); a pair the
  channel gives no color takes the ranges' ink. The layer's data filter (`<pairs>/(except :of (item :type …))`) is read
  pair by pair at prepare (`id`, `type`, `label`, the strands and ends as fields): a pair left out is not drawn, split
  or hit, counted on its track's "n filtered" (a `section` layer's too); a form not read pair by pair is
  `seq.render.filter`.
- **Split order**: at one place the others' closing ends, each self-pair's two ends, then opening ends; equal intervals
  nest: a self-pair, a duplicate and a shared end cross nothing (held against a brute force).
- **Hit**: candidates by hull over the point's x, then distance to the curve (`pairs.hitPx`), and the carets by the
  distance to their boxes, a tie the later drawn (the carets after the arcs); the pair's element, its extent both ends.
  **Hover**: the hovered pair; as arcs, every pair with an end at a hovered position (a base, or either position of a
  cell) or in the hover projected onto the section; as dots, the pair at the hovered cell. **Highlight**: the pairs a
  highlight channel lists by key (`EmphasisHighlight.elements`) drawn over in the channel's color in the `highlight`
  pass, as the selection's in the accent (a pair two channels list: the later's); the pairs each channel lights are made
  again only when the channels change (by identity), into a reused array. The pairs pass ignores
  `EmphasisHighlight.style`: a pair is a stroke over its arc or dot in any style (a tint or a band has no room on a
  curve).

### The grid layer

`render/grid.ts` (`GridType`, `axes: 2`): a `grid` entry's values over the strip's x and `ctx.axisY`.

- **Only on its own systems**: x and y of the section must be the grid's (`seq.render.grid.systems`); no strip y is
  `seq.render.layer.axes`. Why: a grid is not projected; a generator computes it on the systems it is drawn on. Over
  blocks a grid draws on the strips of its own pair of systems only, silently on the others: a matrix over the complex
  is one grid node a pair of chains. A lay on another pair's strip is `off`: it draws nothing, its emphasis neither.
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
  transform, across x's extent); the hovered cell outlined, its row and column tinted and edged in the hover's ink
  (`seq-hover`, the theme's contrast token: they read on either end of a palette); a selected rect drawn as its rect,
  and on a section showing both its systems never its sides projected back (`rectOver`: the strip's x painting and the
  layer's rows skip them). Linked sections light both positions (the hover's location has a y part). Nothing allocated
  once warm: positions and colors go through typed scratch; a hit reads its point from `at`
  (`AxisTransform.positionOf`). Over blocks the hover's column is painted on every strip of its x system, its row on
  every strip of its y system, the cell outlined on the one strip of both.
- **A summary over a rect** (`gridSummary(view, rect)`, `view` a model and a resolved view): the mean, least and
  greatest value and the count of the cells inside the Location's rects of every grid a grid layer of the view draws
  (each grid once, on a rect of its two systems), each cell once however many rects hold it (an added drag's overlap),
  within the layer's `triangle`, missing values left out; NaN with no cell. A library function, no state: a host's
  readout reads it (a selected rect's mean PAE).

### Tensions (pairs and grids)

- **Images in SVG**: an `image` is runs of rects, so a zoomed-out large grid exports large; a PNG inside SVG needs an
  encoder in core. Rects in this slice; measured with the PNG export, which decides.
- **No navigator for a grid's own y window**: a grid whose y has its own window has no navigator over y; open until a
  Hi-C view needs one.
- **Arcs between sections** (ribbons, EC-209) stay deferred: `follow_y` covers a y following a section of the same view.
- **Emphasis along y** is the two-axis layers' own, through the y path (`axes: 2`).

### The rows layer

`render/rows.ts` (the layer type) and `render/rows-model.ts` (what it computes without drawing: clades, display rows,
measures, consensus, the profile); beside them, a concern each, `rows-lay.ts` (what the layer keeps, its rows'
transform), `rows-read.ts` (held, lazy, joined), `rows-cells.ts`, `rows-header.ts` (the tree, labels, numbers),
`rows-hit.ts`, `rows-draw.ts`, `rows-emphasis.ts` (Emphasis along y), `rows-columns.ts` (Row columns). An alignment's
rows over its columns, its own y inside its band (`axes: 2`).

- **Data**: the alignment's `cells` on its rows system. Held when rows × columns is within `maxBlockValues` (2^24: the
  10,000 × 300 alignment is held): every tile wanted at level 1 and joined once (kept by the entries' versions); larger
  ones are lazy: the rows in view only, in data order (another order is `seq.render.rows.order`, info: ordering needs
  every row's field), a row not held labeled `row 5,124` and drawn pending; below a device pixel a row, a lazy source
  with no row bins draws a neutral fill (`seq.render.rows.reduced`, info; never a pending stripe that never resolves).
  Pins and hidden rows are by key, and a lazy row's key comes with its tile: on a lazy alignment neither applies
  (`seq.render.rows.lazy`, info).
- **Order** (`order`): `data`, `tree` (the tree's leaf order, rows no leaf names after them in data order), or a field
  of the rows table (`label`, `first`, an own field, a joined one, a measure: Row columns; numbers by value, text by
  code points, missing last in either direction, ties in data order; built again when the joined fields come and on
  every edit of the layer). Clades fold only in tree order.
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
  on the page, at least `heights.rows`, 240, at most every display row at `RowHeight` with its title row when the rows
  are counted; a figure takes the default). Why the cap: a page as tall as its content gives the last frame's content as
  the room, so a band taking it grew by the sections below on every frame, an empty well under its rows.
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
  pinned band. While the view moves a level not built yet is built `rows.buildBudget` cells a frame (2^18, about 4 ms;
  `Pyramid.levelWithin`, its finer levels first) and the level pair drawn last stays drawn meanwhile, its tiles placed
  at the new zoom (the rows in view drawn pending while no level of the base is drawn yet, a zoom out of the letters or
  a fold: never an empty band, EC-169); a settled frame builds the rest (`rows-levels.ts`). Why: a 10,000 × 300 level
  from its base took two frames of about 100 ms each zoom (performance.md, Levels). A level whose sample covers under 64
  rows counts its cells' classes exactly, from the base; a larger one is the majority of the majorities below it (shares
  averaged), approximate within the measured tolerance, tested over seeds 1 to 12: the share within 0.05 on average and
  0.3 at worst, the letter exact where the exact winner leads by more than 0.25. Why: a majority is not mergeable across
  levels, and a majority of majorities is wrong for small blocks (measured on an unbiased generator, seeds 1 to 12: the
  letter missed at leads up to 0.41, shares off by up to 0.41); a few rows a sample are cheap to count; and below 3
  device px a cell the shade reads, not the letter. The pinned rows are drawn on their own over a band that does not
  scroll (rects past `pinnedFrom`), the scrolled rows clipped below it. The reducer is the theme's classes, fixed (no
  item, decided by the author 2026-11-20); a reduction by letter is in the package's TODO.
- **The header**, in the label column: the tree (`tree_width`, 120) and the first residue numbers widen the view's label
  column by their width (the widest of its tracks sections, so every section's data starts at one x), the labels keep
  the column's own; the tree by branch lengths, else depth, a folded clade a wedge at its summary row, its geometry in
  display rows kept by the display's version (a pan makes nothing sized by the nodes). **The tree's look**: every branch
  one width, whole device pixels, its edges snapped to them (a branch runs half its width past each end, so it meets its
  node's vertical square); a wedge a triangle from the clade's node to its deepest leaf's x (by branch length, 6 px at
  least, started 6 px before its node at the tree's edge), its far side half the row high, its corners softened by the
  shape rule, filled with a light tint of the branch ink (0.18 of its opacity, 0.36 hovered) and outlined in it at the
  branches' width; its hit the summary row's middle 0.8 from the node to that x. Why: a 1 px rect centred on a
  fractional x smeared over two device pixels, so widths varied, and the solid triangle at the tree's full width read as
  a heavy blade. The tree is drawn down to the clades whose rows take `TreeMinPx` (4, fixed): a node's segments, the
  branch to it and a folded clade's wedge only then (by its shown rows, kept with the geometry; a node's vertical spans
  the node and its drawn children only, none when none is drawn), so zoomed out the deep clades end in their branch's
  tip. Why: below it branches sit under 2 px apart, and 10,000 leaves at the band's whole zoom drew a solid dark slab
  (every leaf's branch). Labels in the items' font (a row's label is data; a Pfam label, `NAME/1-147`, reads whole) from
  `rows.labelMinPx` (8 px a row), the label column empty below it (no label, number, chip or icon a row), cut with an
  ellipsis past their room, a summary row "first +n" (its first label cut, never its count: the count is what it says);
  the first residue numbers right-aligned before the data box and the last after it, in a right gutter the view keeps
  when a rows layer draws numbers (`numbers: gutters`): the numbers of the residues **in view** (the columns the cells
  are drawn for: a pan and a zoom change them; a gap is never counted; a row with no residue in view has none), the
  row's own at the fit view (kept per data, or per lazy row once). The right gutter is clipped to the band vertically,
  as the labels and cells are: a partly shown last row ends at the band's foot on every side.
- **Row columns** (`columns`: a list of `{ field, side, as, width, title }`, `rowColumnsOf`; `render/rows-columns.ts`):
  data columns beside each row's header. `left` (the default) in the header between the label and the first numbers, the
  label column widened by them (a text figure's too); `right` after the last numbers' gutter, the right gutter widened
  by them and by a gap before a rows rail; a gap (`rows.columns.gapPx`, 6) before each. `as`: `text` (the default, cut
  with an ellipsis), `number` (right-aligned; an integer whole, else two decimals; the measure `identity` a percent:
  formatted by the field's origin, never its name), `bar` (the value's share of the column's range from 0, in the
  accent) or `chip` (a swatch of its category's color, `seq-cat-*` by the category's index in code-point order, and its
  text); `width` by `as` (text 96, number 44, bar 56, chip 56), `title` the field. A column whose values all read as
  numbers (numbers, or text that parses as a finite number: a `#=GS` TAXID, an untyped table's field) draws and orders
  as numbers (`numericValues`); the data's own types are unchanged. Why: a Stockholm or TSV field is text, and a bar or
  an order over it means its numbers. **Fields** (`RowsLay.fieldValue`, by origin): `label`, `first`, the rows' own
  (`HeldCells.fields`: `#=GS` tags, an inline table's other columns), the measures (`length` the residues; `gaps` the
  row's `-`; `identity` to the reference row: the same residues over the columns where either row has one, a gap against
  a residue a mismatch, gap against gap and an insert on either side, a lowercase letter or a `.`, left out;
  `ColumnMeasures.row*`), then the fields of a records source whose system is the rows system (its `system` param, or
  its records' `target`: `/data/aln/(rows)`) joined by row key (a record's `id` names a held row by label, else its
  position is the row). Why that identity (PID1-like, over aligned positions): it does not inflate a fragment's identity
  by the columns it lacks, and an insert is no aligned position; the author may choose another. Only the fields a column
  or the order names are joined; one an own field or a measure shadows, or that a second source has, is not read and
  said once (`seq.render.rows.field`, info). The join is kept by the held rows, the prep and the entries' versions,
  compared in place. A missing value is an empty cell, a summary row's cells are empty; a field no held row has is said
  once (`seq.render.rows.column`, info) once no source that could join (its `system` the rows system, a runtime bound)
  is still to attach, or the view has settled (a source that failed holds nothing). Columns read held rows only (a lazy
  alignment's are empty). **The title row**: a rows layer with columns keeps a row of its own above its rows
  (`rows.columns.titlePx`, 16: its track grows by it; a `fill` band's rows take what is left; a text figure's every
  line), the titles in it (a figure's first line only), muted, the ordered field's marked `↑` or `↓`. Why its own: the
  titles never reach the track above, nor leave the page when the rows are a section's first track. A title hits as
  `column-title` by its drawn text (`state` its field: the viewer orders the rows by it). **Narrow** (EC-107,
  `fitRowColumns` in the view's layout, kept while the view, the width and the room it read are): what the view lacks
  for the data's least width (`rows.columns.dataMinPx`, 200) with every column at its own width is taken from the
  columns, the right ones first and then the left, each from the last written, down to `rows.columns.floorPx` (32) and
  then hidden; the label keeps its room. A chip in the title row names those hidden ("Species, Identity hidden";
  "Species +1 hidden" short of room, else its first name cut, its count kept, in the roomier band; no room for a name,
  the count alone, "2 hidden": a column is never hidden with no chip), on the side they live and never over the data
  box: hidden on the right only, at the end of the right columns' band (from the numbers' gutter at the least; the right
  titles cut for it); else, or with no room there, in the label column's title band, ending before its first left title
  shown (else at the label column's end). Why: a count of "columns" beside an alignment read as alignment columns, and
  the chip hung over the data. It hits as `columns-hidden` (`state` its text: the viewer's card lists each hidden column
  and why, "no room at this width"; a click does nothing: showing them is not designed). The widths reach the layer as
  `LayoutContext.rowColumns` (with the right columns' x); a text section's lines leave the right columns their room; an
  export attaches the sources that could join. A cell's text, its cut (kept by the column's room) and its widths are
  made once a row and kept with the resolved values: a layout once warm, a draw and a hover allocate nothing
  (`rows.quick.ts`).
- **Hits** (the header's through the band's gutter: a layer of two axes is asked there), written into one scratch hit (a
  hover allocates no hit object; a clade's element and extent made once a display): a cell is `cell` with its column
  (`position`) and row (`ySystem`, `yPosition`); a header's label is `row` (its `position` the data row on the rows
  system; on a row column's cell, left or right, the same with `state` naming its value, `Species: Homo sapiens`, made
  at layout); a column's title `column-title`, the hidden-columns chip `columns-hidden` (`LayerType.hitOutside`: what a
  layer draws outside its band, asked by `hitView` before a layout's chrome and the bands); a summary row `summary-row`,
  and `tree-node` the inner node whose branch segment is nearest within `TreeHitPx` (8, fixed: decided: a junction's
  size, a target a pointer finds) by distance to the segment, a wedge included (a leaf's branch hits none), each its
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
  band and outlines the cell, and draws a column cell cut short whole in the `tip` pass only while the pointer is on
  that cell, on the copy under it (a pinned row kept in place is shown twice: the layer's last hit says which; the
  cursor on the row's header: all of them, on every copy it outlines), never from a hover on the row's data: on an
  opaque surface (the page's background, then the track's well, outlined in the hover ink, a squircle by the shape
  rule), growing right up to the view's edge where a longer text ends in an ellipsis (a number column's grows left); the
  selection's parts on the rows system as row bands, its rects of rows × columns as their cells, read over the display
  rows in view only (a binary search a row: a selection of any size costs the rows in view); a summary row draws the
  partial mark (the tint over the share of its rows selected, the outline when all are); the cursor with a row
  (`Emphasis.cursor.ySystem`) outlines its cell, or on its `part` `header` the row's label column (`labelX` to
  `numberX`), a squircle as the strip's (`squircleRadius`), on its own section only (`Emphasis.cursor.section`: another
  section of the system marks its column, as the strip draws a cursor without a row, and nothing while it is on a
  header); on the row it stands on (`cursorRow`: its own, a folded row's summary row, a hidden or filtered row's nearest
  shown), a pinned row's in the band and in its place. `displayLabel` is a display row's label whole (a summary row's
  `first +n`, `summaryMore` its count, which the header elides apart); `hitRow` the display row of the last hit (which
  copy of a pinned row a click was on).

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
- **SVG** (`drawListsToSvg`): a string, no DOM; passes as groups, clips as clip paths (a squircle clip the core's path,
  as the Canvas2D backend traces it), squircles as the core's path, images as `<image>`s
  (`image-rendering="optimizeSpeed"`: crisp; an attribute, since resvg ignores a style it cannot read) of a PNG data URL
  (`render/png.ts`: stored deflate blocks, synchronous, no DOM), fades as gradient masks, metadata as `data-*`
  attributes. Every text and attribute from a spec or data is escaped; a character XML forbids (a control but tab,
  newline and return, U+FFFE and U+FFFF, a lone surrogate) becomes U+FFFD, one for one, so a glyph run stays aligned
  with its x list. Ids take a prefix (`idPrefix`, name characters only), a metadata key is written only as a plain name.
  Held by hostile labels, ids, titles, fonts and metadata parsed back as XML (`_test/svg-escaping.test.ts`).
- **Accessibility** (`accessibility.ts`, `accessibleOf`; `SvgOptions.accessible`): an export's root is `role="figure"`,
  `aria-labelledby` its `<title>` (the spec's `title`, else `Sequence view`) and `<desc>` (the sections, their systems
  and track counts, and how many items are not named); over the picture a group `data-seq-accessible` holds a
  `role="group"` a track (`aria-label` and `<title>`: its label, its item count and its section) and in it a
  `role="img"` a named item (`LayerType.items`, which the ranges, secondary, variants, rows (a row an item: its label; a
  folded clade's summary row by its summary label, `HBA +12`) and pairs (a pair an item) layers have: the label, else
  the key, with its extent, 1-based and inclusive: `DNA-binding, p53 94–292`), then a `role="img"` a highlight
  (`Highlight 1: 10 positions, p53 11–20`). A track names at most `MaxNamedItemsPerTrack` (200) items and a figure
  `MaxNamedItems` (2,000), saying "n more not named" on the track and in the description (a figure of 2,000 named items
  adds about 190 bytes each: roughly the size of the picture itself for dense ranges). Each box is
  `fill="none" pointer-events="none"`: the layer never takes a pointer from the picture beneath. Why boxes over the
  picture and not marks with names: the draw lists carry no semantics, and a layer's items are already what the hit test
  names. **Ids**: the default `idPrefix` is a hash of the spec, the title, the size, the token table, the view, the
  highlights and the selection, so exports that differ in them can share a page; two equal exports share their ids, and
  a host inlining both passes each its own `idPrefix`. The viewer's announcements use the same words (viewer `_spec`,
  Announcements).
- **PNG** (`renderPng`, `export-png.ts`): the SVG export rasterized, so a PNG is the SVG's picture (the same options,
  the same figure); the rasterizer is a port (`Rasterizer`: the SVG string and its size and `scale` × `dpr` to PNG
  bytes) so the core carries no native code, and the application composes one (the render CLI's, `cli/seqx-render`, over
  resvg). The SVG paints its own background, so the pixels are opaque. A rasterizer loads no system fonts: the text is
  drawn with fonts bundled for the metrics the layouts measure with (Liberation Sans and Mono for the bundled table's
  Arial and 600-advance mono), so the same SVG draws the same on every machine. Why rasterize the SVG and not draw the
  lists again: one picture, one set of tests; a second renderer would drift. Held by the render CLI's probe pixels per
  track in both themes (`cli/seqx-render/src/_test/png.test.ts`).
- **The render CLI** (`seqx-render`, `cli/seqx-render`): a spec file (JSON, S-expression or MessagePack, by content) to
  an SVG or a PNG, its assets read under `headlessAssetPolicy` (files under the spec's folder, no network unless
  `--allow-host`; a spec's `uri` is untrusted), its issues printed one a line to stderr, exit codes 0 done, 1 spec, 2
  usage, 3 unsettled, 4 file, 5 internal. Cold start, process start to a written file for p53's tracks, under 1 s (about
  0.12 s), held in the heavy lane.

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
- `exportSvg(spec, { size, types, layouts, navigatorExport?, … })` draws with what it is given (none of a kind when
  absent: each layer said); `renderSvg(spec, options)` is `exportSvg` with the built-ins (`BuiltInLayerTypes`,
  `BuiltInLayouts`, `NavigatorSvgExport`) unless the options name their own. Included navigators are drawn by the
  `navigatorExport` part (`render/navigator-export`); without one an export asked to include them leaves them out,
  `seq.render.navigator.export`.
- **Call options** (`RenderSvgOptions`; options of a call, none a viewer setting or a spec param):

  | Option                                                                      | Default                                                         | Why the default; the alternatives                                                                           |
  | --------------------------------------------------------------------------- | --------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------- |
  | `size`                                                                      | required                                                        | CSS px; `height: 'auto'` the sections' content's                                                            |
  | `viewport`, `viewports`                                                     | the fit view                                                    | a window of display units for every leading section, or positions a section; a viewer's sections by address |
  | `view`                                                                      | the spec's first view                                           | an export names one view                                                                                    |
  | `window`                                                                    | the view's `export.window`                                      | a ring's window drawn; false: an exported map is a figure                                                   |
  | `legend`                                                                    | the view's `export.legend`                                      | `none`, `below`, `right`                                                                                    |
  | `navigators`, `navigator`                                                   | the view's `export.navigators`, else `omit`; `DefaultZoomStrip` | a window is session, not a figure; `include` stacks them by `place`                                         |
  | `cvdSafe`                                                                   | false                                                           | an export draws the figure's colors unless asked for the reader's substitution                              |
  | `folds`                                                                     | each rows layer's `collapsed`                                   | a viewer passes its own, as it passes the selection                                                         |
  | `selection`, `emphasis`                                                     | none                                                            | an address resolved by the export once settled; emphasis as the viewer draws it                             |
  | `idPrefix`                                                                  | `s` and a hash of what the export depends on                    | exports that differ can be inlined in one page; equal ones share ids                                        |
  | `theme`, `tokens`, `text`, `dpr`                                            | `light`, the built-in tokens, the bundled metrics, 1            | a reader's theme and metrics                                                                                |
  | `config`, `resolve`, `bounds`                                               | the defaults (`DefaultRenderConfig`, `DefaultSpecBounds`)       | overrides field by field                                                                                    |
  | `runtimes`, `assets`                                                        | inline data through `seq` only; `NoAssets`                      | nothing read unless the caller allows it (security)                                                         |
  | `timeoutMs`, `clock`, `store`                                               | 10 s, real timers, a new store                                  | one deadline; a viewer's store and clock to share its data                                                  |
  | `types`, `layouts`, `sizings`, `trackSizing`, `registry`, `navigatorExport` | the built-ins in `renderSvg`, none in `exportSvg`               | the composition's registry; a minimal composition passes its own                                            |

- `renderSvg(spec, { size, viewport?, runtimes, … })`: checks the spec and opens its model (`openModel`,
  `scope: 'view'`: the sources its view reads, `attachSources`: inline data through `seq`, a source node through its
  runtime, a template's source, the source of a numbering a layer names, which no layer's data does: a ruler's; assets
  refused unless a reader is given, `NoAssets`), wants through a want set of its own (EC-125), each round reads the
  overlay and syncs its sources (let go at the end), lays out until its wants and the store's node version hold, waits
  for them to settle, and draws once more with the view settled. It waits for its view's wants alone, never their
  margins.
- **One deadline** (`timeoutMs`, 10 s, on `clock`: real timers by default, a shared store's clock or a fake one) covers
  opening the sources and settling; past it the opens are aborted and the export is `seq.render.unsettled`. The model
  and view resolved again are taken only when the opening wins the race with the deadline: a source that opens later
  never swaps the model under a render (the first resolve's, its lengths pending, is drawn).
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
  `theme.color`, `scale.log` (info), `legend.unlisted` (info), `filter`, `line-style` (info), `numbering`,
  `numbering.fallback`, `numbering.limit`, `numbering.source`, `numbering.rows`, `layer.type`, `layout`,
  `layer.letters`, `axis.pending` (info), `text.lines`, `follow.leader`, `follow.blocks`, `ring.image`, `follow.cycle`,
  `source`, `template.limit`, `unsettled`, `navigator.of`, `navigator.path`, `navigator.rows` (info), `navigator.layer`
  (info), `navigator.export`; a test checks the catalog.

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

- `composeSeqContext({ registry?, config?, ports? })`: the registry's runtimes, layer types and layouts as listed, the
  configuration and the ports, and through them `check`, `resolve`, `view` (the resolve passes) and `renderSvg`
  (`exportSvg` with the registry's types and layouts, no navigators). `createSeqContext(spec)`, the default: the core's
  layer types and layouts (`BuiltInRenderRegistry`) listed before the spec's own entries, and the navigators in its
  exports (`NavigatorSvgExport`). `retire()` disposes what the context made or was told to dispose. Why two: what a
  composition does not list is not in its bundle (a viewer composes the base), and the headless default draws
  everything.
- **Data adapters, one a format**: a format two registered runtimes read is an error listing every conflict
  (`runtimesOf`), so which one reads a source never depends on the order they are listed in;
  `withAdapter(registry, adapter)` puts a host's adapter in place of the defaults for its formats (a runtime left with
  other formats keeps them under its name, one narrowing however often it is listed). Format runtimes are keyed by their
  formats and their identity, not by name: a narrowed runtime and an adapter may share one. The resolver and the overlay
  bind a source to the one runtime of its format (`runtimeFor`); a list given them with two for a format binds none
  (`seq.resolve.source.runtime`). Why: extensibility's "never resolved by order".
- **Registry entries** are data: `{ kind: 'format-runtime', runtime }`, `{ kind: 'layer-type', type }`,
  `{ kind: 'layout', layout }` and `{ kind: 'track-sizing', sizing }` (an `ext:` sizing a section's `heights` names:
  Render, Track heights), each keyed by its provider's `name` (a format runtime by its formats, above); and the query
  language's `{ kind: 'query-form', … }` and `{ kind: 'query-property', … }` (Queries, Extensions), composed into
  `SeqContext.query` (`queryExtensionsOf`), what a host passes `evaluateQuery` as `extensions` (the viewer's address
  model does not yet: a followup). The same object again is no change; another under one name is an error listing every
  conflict (nothing registered), so a host's own type of a built-in's name goes in `composeSeqContext`. The canvas and
  the SVG export draw the same registry's layer types: a layer of a type not listed is said (`seq.render.layer.type`)
  and not drawn, a section of a layout not listed is said (`seq.render.layout`) and laid out as tracks. Run-time, atomic
  and undoable registration is the viewer's followup.
- **The tracks layout is the one fixed fallback**: a section whose layout is not listed (or every section, with no
  layouts given) is laid out as tracks, and the viewer's navigation reads tracks likewise; tracks needs no entry. Why:
  the tracks layout is the frame's own (the label column, the rails, the bands), a section must be placed somehow, and
  the issue says it was not the one asked for.
- **No built-in table a renderer falls back to**: `Renderer` and `NavigatorRenderer` take `types` and `layouts`;
  `layoutView` and `exportSvg` given none draw none. `BuiltInLayerTypes` and `BuiltInLayouts`
  (`@molstar/seqx-core/render/built-ins`, and `BuiltInRender` to spread into a renderer's options) are the default
  compositions' lists, imported by them alone: `createSeqContext`, `renderSvg` and the viewer's default spec.
- **Configuration** is `[item, value]` pairs (`SeqConfigItem`, `configItem(key)`); the core's items are `SeqCoreConfig`
  (`theme`, `render`, `resolve`). Items of other keys are another package's (the viewer's `group.item` keys).
- **Ports** are factories with ownership, each called once; unset, the default: no reads (`NoAssets`) and the bundled
  text metrics. A result's `disposeWithViewer` is true unless `false`: a host that passes its own keeps it.
  - **assets**: either the application's `AssetManager` behind `policyAssets` (`ManagedAssetsPort`, with a policy), or a
    reader of its own (`ReaderAssetsPort`: a test adapter that refuses the network, the Node fixtures' adapter). Why
    two: a suite per port needs a second adapter, and an `AssetManager` is not what a test or a CLI has.
  - **text metrics**: `measure(font, text)`; the context wraps it in the per-font cache (`createTextMetrics`).
  - **structure**: none in this version (`runtime` unset).
  - **track sizing** (`TrackSizingPort.heights(section, tracks, named)`): the host's last word on each section's tracks'
    heights (Render, Track heights); the context's `trackSizing`, which its renderer and `renderSvg` read. Why a port:
    an application fitting the viewer into its own panel decides heights from what only it knows.
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
- A port factory that throws leaves nothing made: the ports already made are disposed, last first. There is no structure
  port: syncing with a structure viewer is a library over the viewer's events and commands (Structures).
