# Tooltips: the card composed by the application, the spec's words on rendered marks

## Intent

The author, 2026-11-22: "we should be able to associate arbitrary tooltips with data points. MolViewSpec has this too
with .tooltip() ... we should find a good way to do this here too"; then "basically, we should be able to associate
'tooltips' with addresses in general. ... or rather a tooltip by range, etc... the situation here is a bit different. i
think we should design this first."

On the first draft: "i am missing the 'tooltip from data/uri' option here"; "i think we can do markdown? we do have
packages/core/src/util/markdown.ts"; "the placement: maybe configuration on 'viewer spec'? after all, mol* plugin also
composes tooltips via providers. maybe we should adopt similar principle? in general, tooltips should be composed by
'application layer' anyway"; on SVG export, "not sure we need this", then "or well .. maybe <title> in SVG could work
too ..."; on a note outlining its range: "i think we should only be able to associate tooltips with data points/ranges
etc. .. so it would always be tied to something rendered?"

So two things, kept apart: **the spec says which words attach to which data** (portable, read the same by every
reader), and **the application composes the card** (which providers speak, in which order, how links act).

## Precedents

| Precedent | Taken | Not taken |
|---|---|---|
| MolViewSpec (`mvs-tree.ts`; runtime `annotation-tooltips-prop.ts`, `custom-tooltips-prop.ts`) | `tooltip { text }` on a component; `tooltip_from_uri` (`uri`, `format` json, cif, bcif, `schema`, `field_name` default `tooltip`, `text_format` an f-string, `field_remapping`, the CIF block params); `tooltip_from_source` (the same table as a CIF category in the structure's own file: `block_header`, `block_index`, `category_name`); a row with a field undefined gives no text; every matching tooltip shown | texts joined by `" | "` (here a row each); a tooltip a child of each component (here one node) |
| Mol\*'s plugin (`spec.ts` `lociLabels`, `state/manager/loci-label.ts`) | hover labels composed from **providers listed in the plugin spec**, each given the hovered thing, returning a label or nothing | labels as HTML or React elements (here typed rows) |
| Seq\*'s viewer spec (`default-spec.ts`, `header-entries.ts`) | the viewer composed from data; `header-action` entries, a viewer-only registry kind ordered by `order`, then name | — |
| `@molstar/core/util/markdown.ts` | a typed tree, never HTML, linear, nesting capped; `parseMarkdownInline`, `markdownInlineText` | its blocks in a card row; its `isSafeMarkdownHref` as the link rule (a relative `/refs/dbd` passes it: Text) |
| `@molstar/core/util/string-format` (`FormatTemplate`) | MVS's f-strings, one implementation | a second template language |
| `describe.ts` (the card's sections) | `Parts`: a `title`, `rows`, a separate `words` list (the live region's sentence, not one per row: `ItemSection` writes three or four rows and one phrase, `PlaceSection`'s `index` row is never spoken), `whole` (an item, a region, a pair said whole: later place sections skip) | — |
| The SVG export (`core/src/render/svg.ts`) | an accessible item is a named box whose `<title>` is the browser's tooltip and its accessible name | — |

## Decision

### Tied to something rendered

The author's rule: a tooltip "would always be tied to something rendered". So:

- **A tooltip answers only on rendered marks**: an item, a range, a pair, a cell, a letter, a value bar some layer
  draws, or a bin of them. A position no layer marks gets none (a track's empty band, a gap). Why: the words are about
  data the reader sees; a note on nothing visible is a note the reader cannot find.
- **No outline of its own**: the mark's hover emphasis is the tooltip's (one hover, one outline).
- **Two ways to attach, one result** (items with a location and fields; a row's text from its fields):
  - **by data**: a data node's (or a layer's) items speak their own fields (`on: /data/p53/features`, `text:
    "{label}: {note}"`): each item's text on its own marks;
  - **by range from a table**: a table of locations with text (inline rows, or a `source`: a URI, a CIF category):
    each row's text on every mark its location covers (a letter, a feature, a value bar, a cell), in any layer of the
    view. An address is a one-row table (`on: /refs/dbd`, MVS's `tooltip { text }` on a component).

### The spec: a `tooltip` node

Under `views` (every view, beside `scale`) or a `view` (that view only); ids required (a story step adds or removes
one by patch).

| Param | Says |
|---|---|
| `on` | by data: a data node of items or a layer (its marks only, its filter applied); by range: an address or a list (a ref, a range, Keys) |
| `target`, `rows` | by range, the tooltip's own inline table (location columns, text fields) |
| `target`, `source` | by range, a table read from a `source` node under `data` (a reference: data stays in `data`) |
| `schema` | how a table's rows locate: `ranges` (absent: `start`, `end`, a sequence name on an assembly), `address` (an `address` column, forms relative to `target`), or MVS's schemas (`whole_structure`, `entity`, `chain`, `auth_chain`, `residue`, `auth_residue`, `residue_range`, `auth_residue_range`, `atom`, `auth_atom`, `all_atomic`), each read as MVS reads it, as `(elements …)` items |
| `block`, `category` | a CIF table: MVS's `block_header` or `block_index`, `category_name` |
| `field` | what `{}` reads (`tooltip`, MVS's `field_name`) |
| `text` | always an f-string (`"{}"`, MVS's `text_format`) |
| `label` | the card row's label; absent, the text is the row |
| `format` | `markdown` (absent) or `text` |

- **MVS's forms, mapped**: `tooltip { text }` on a component is `on` its selection as an address with `text`;
  `tooltip_from_uri` is a `source` node (`uri`, `format`, `fields` from `field_remapping`) and a tooltip with `source`,
  `schema`, `field`, `text` (and `block`, `category` for CIF); `tooltip_from_source` is the same over the `source` node
  the structure is read from. **Where Seq\* differs**: a sequence spec has no structure file of its own; a CIF category
  is read only through a structure's `source` (a molecule's or a complex's), and until the structure runtime reads
  categories such a tooltip is unavailable (`seq.lang.tooltip.category`, said once).
- **Schema is said, never guessed**: a `label_asym_id` and `beg_label_seq_id` file under `schema: chain` is whole
  chains, as in MVS; a table whose columns do not fit its schema is `seq.lang.tooltip.schema`.
- **One mechanism with a `ranges` column**: `on` a `ranges` node with a `tooltip` column and a tooltip's own table are
  one items table in the store (one index, one region read, one projection). Inline `rows` is the only table that
  lives with its tooltip (MVS has inline tooltips; a never-drawn note list has no other home); a read table is a
  `source` under `data`.
- **GA's assemblies and names**: ranges on an assembly name their sequence by the names rule (aliases, `chr`);
  addresses take `(locus …)` and `(block :name …)`; a name that is not an id is reached by its form (EC-200).
- **Tiers** ([composition.md](composition.md)): tier 1 holds the tooltip node resolved to references (`on`, `target`,
  `source`) and its parsed template, never rows. Every table, inline too, goes through the store as items, as every
  data does ([data-store.md](data-store.md)), and its interval index belongs to the store entry, built when it loads.
- **Read with what it decorates**: a table no layer draws is wanted by no one, so its tiles are requested beside the
  marks it decorates: the viewport's request for each layer whose marks the table can cover (on its system, or
  reaching it by a path) adds the table's tiles over the same window, and they are kept as those are (never evicted in
  view). States per data-store.md: pending gives a `loading` row (EC-169); empty, no rows; failed, no rows and the
  issue on the tooltip, retried by the policy (EC-174). Read whole (addresses, MVS's schemas) only within the asset
  policy's size of one asset and `SpecBounds.tooltipRows` (1,000,000).
- **Assets**: a `source` is read through the assets port's policy ([security.md](security.md)); a refused read is an
  issue on the tooltip, and the rest opens.

### The application: the card composed by providers

- **A `card-provider` registry entry** in the viewer spec, beside `header-action`:
  `{ kind: 'card-provider', name: '<ns>:<name>', order?: number, describe(at, out) }`, ordered as header entries are
  (by `order`, then name: one registry rule). `at` is the hit and its place (the section, the layer, the system and
  position, a cell's y, the element, a bin, the block) and the models; `out` is today's `Parts`: a `title` (a provider
  sets it when empty), `rows`, a separate `words` list (the sentence's phrases, never one per row), and `whole`.
- **The viewer's sections become built-in providers**, `BuiltInCardProviders`, with orders that keep today's order:
  `seq:region`, `seq:pair`, `seq:group`, `seq:item`, `seq:columns-hidden`, `seq:row`, `seq:grid`, `seq:residue`,
  `seq:place`, `seq:tooltips`, `seq:block`. Their rows and words are unchanged, so the live region's sentences and the
  card tests hold as they are.
- **`whole` is the place's, not the writer's**: `seq:tooltips` ignores it, so a note on 94 to 292 shows over a
  feature, a cell or a pair the item and region sections said whole. A host's provider reads it as it chooses.
- **What stays where**: the spec says which words attach to which data; the viewer spec says which providers run,
  their order, the caps (`hover.cardItems`, `hover.tooltipLines`), whether the card shows (`hover.card`), and how a
  link acts. A host adds an entry, filters one by name, or replaces `seq:tooltips` with its own over the core's
  `tooltipsAt(view, place, out)`.
- **The hover contract holds** (HC: every mark answers with a hit, an emphasis and a card): providers change rows,
  never the hit or the emphasis. Its floor is `seq:item` and `seq:place`; a composition without either is
  `seq.viewer.card.floor` (a warning). The HC sweep runs on the default spec.

### What the reader sees of the spec's tooltips (`seq:tooltips`)

- **The mark's place**: the hit's item by key; the position under the pointer on the mark (a cell: both positions; a
  pair: both ends; a matrix cell: its position and state); a bin's stretch. A range matches by overlap of one position
  or more, a cell or a pair by either place, a rect by both; projected through the hover's own projection (a note on
  the protein answers on its gene's codons).
- **By data**: the hit's item's own row. **By range**: every table row covering the place, in position order.
- **Rows and the cap**: a tooltip's single row on the hovered mark is a plain row, never cut. Rows listing several
  matches (over a bin, or several table rows over one mark) are item rows sharing the one `hover.cardItems` cap with
  the place's items, after them (`cutItems` counts both). Past the cap, `HoverCard.more` holds the rows made and a new
  `HoverCard.moreCount` the matches not made (counted by the index), so the `+k more` row says `k = more.length +
  moreCount`; the live region says the rows shown and the count, never a million texts.
- **Stacking**: `views`' tooltips before the view's, each in spec order; none hides another or a built-in row.

### Text: Markdown in the template, links inside the view

- **Inline Markdown** (`format: markdown`): strong, emphasis, code spans, links, line breaks (`parseMarkdownInline`);
  no blocks (a row is a line of a compact card). `format: text` shows the text as written.
- **Markdown applies to the template only; values are escaped**: each `{}` value has Markdown's escapable characters
  backslash-escaped before the parse, so `A*02:01 and B*07:02` stays text, and a data file adds no mark or link.
- **Links act only through commands** (security.md: "links act only through commands on addresses"): only Seq\*'s
  link grammar (`!seq-focus=/refs/dbd`, `!seq-select=…`, the links codec's) is a link; anything else is its text. Why
  not `isSafeMarkdownHref`: a bare `/refs/dbd` is a relative URL to it, and the asset policy is a read policy
  (`https:` reads, no `mailto`), not a navigation rule. External links are an open question.
- **Cleaned after parsing**: `parseMarkdownInline` decodes entities, so C0 and C1 controls (but breaks) and the bidi
  overrides and isolates (U+202A to U+202E, U+2066 to U+2069) are removed from the text nodes it returns.
- **Plain readers**: the live region, a headless card's `value` and an export's `<title>` read `markdownInlineText`;
  `rich` holds the tree for an adapter that draws it.
- **Bounds**: a template at most 1,000 characters (`SpecBounds.tooltipText`); a formatted text cut at 1,000; a row
  wraps to `hover.tooltipLines` (3) lines, then `…`.

### The enterable card (EC-166's bridge)

- **A card with a link is entered**: the pointer may cross from the anchor to the card within a grace (300 ms, the
  safe triangle between them), the card stays while it is over it, the hover frozen on the anchor; leaving hides it. A
  card with no link takes no pointer (`pointer-events: none`), so it never steals a hover from the data.
- **From the keyboard**: the cursor's card's links are listed in the cursor's context menu (Shift+F10, EC-156).

### The SVG export's `<title>`s

- **An option**: the view's `export.titles`, `name` or `card` (renderSvg's `titles` option overriding it, as
  `export.legend`). With `card`, each exported item's named box (the accessible export's: "a box nothing paints and a
  pointer still hits: its `<title>` is the tooltip and the accessible name") takes the composed rows as plain text
  (`markdownInlineText`); with `name`, as today. Items only: bins, image tiles and summaries get none (not items; a
  bin's rows depend on the zoom the figure froze).
- **Cost**: a provider call per exported item, once per export, never per frame; size about the rows' text per item
  (a 60-character note on 10,000 items, 0.6 MB). The viewer's export passes its composed providers; the core's
  headless export passes `tooltipsAt` alone. The default is an open question.

### The builder

```ts
const notes = data.ranges('notes', { target: p53, rows: { start: [93, 324], end: [292, 355], tooltip: ['DNA-binding', 'tetramerisation'] } });
views.tooltip('notes', { on: notes, label: 'Note' });                              // by data
views.tooltip('dbd', { on: refs.ref('dbd', { to: `${p53}/(range :from 94 :to 292)` }), text: 'binds **DNA**' });
views.tooltip('curated', { target: p53, schema: 'address', rows: { address: ['(at :label 175)'], tooltip: ['R175H, [a hotspot](!seq-focus=/refs/hotspots)'] } });
views.tooltip('mvs', { target: tup, source: data.source('ann', { format: 'json', uri: 'ann.json' }), schema: 'auth_residue', field: 'note' });
```

`check()` reports a template field the table lacks. Viewer side: `{ kind: 'card-provider', name: 'lab:assay', order,
describe }` in the registry.

### Performance

- **Per hover step, 0 B warm**: a step is the hit and the emphasis; providers run when the card shows (at rest) and
  again only when the place or the store's version changes; a place re-hit each frame keeps the card it has.
- **The lookup**: the store entry's interval index (the runtimes' `intervalIndex`, moved to a shared internal module
  of the core with its test); matches into a reused `Int32Array`, stopping at `cardItems + 1`, the rest counted:
  O(log n + k), 0 B.
- **Parses cached by text** (64 entries): a cold card parses and allocates; a warm one does not. The card's
  allocation per step, 2,168 B, is the allocation tier's reading of the card show (the viewer's card case); the tooltip lane measures a cold and a warm card show with tooltips in the same tier,
  records both in its report, and holds the warm one at today's figure.
- **Targets**: the lookup over 1,000,000 ranges within 1.5× of 10,000 (heavy, `itWithin`); a hover step at 0 B warm
  with a 10,000-row inline table and a 1,000,000-range source; an export with `card` titles over 10,000 items within
  1.2× of `name` (heavy).

## Defaults

No `tooltip` node: the card as today. With one, a row after the place's rows on the marks it covers:

```
┌──────────────────────────────┐
│ p53 175                      │
│ index 174 · uniprot 175      │
│ Letter R                     │
│ Note  R175H, a hotspot       │
│ Chain A                      │
└──────────────────────────────┘
```

The default viewer spec lists every built-in provider; nothing in an empty band; an export's titles are item names
unless `export.titles: card`.

## A separate change: the matrix's member labels

The author, in passing: a matrix layer should optionally show its member (sample) names, as a categorical grid axis
does. A layer param, `member_labels`: `none` (absent) or `left`, the states' names down the track's label column,
drawn by the categorical axis's `placeMemberLabels`, `drawMemberLabels`, `hitMemberLabel` (thinned past
`rows.labelMinPx`, following `bandTop` and `rowScale`, a label a `target: 'row'` hit selecting the state); off by
default (thousands of states are sub-pixel rows); on a ring not drawn (`seq.render.matrix.labels`, info). It shares
nothing with tooltips but the language bump.

## Rejected alternatives

- **A tooltip on any position, marks or not** (the first draft): the author's rule ties words to what is rendered.
- **An outline of a note's range**: the mark's emphasis is the tooltip's; a second outline is a second hover.
- **Plain text only** (the first draft): the core's Markdown is safe by construction; blocks in a row stay rejected.
- **Markdown over values**: a data value's `*` or `_` would become emphasis (`A*02:01 and B*07:02`), and a file a link.
- **`isSafeMarkdownHref` as the link rule**: it passes relative URLs and is no navigation policy; links act through
  commands only.
- **A `key` column**: `address` with `(item :id …)` says the same.
- **Columns deciding the location kind**: MVS says `schema`; a guess reads a chain file as residues.
- **Indexing tables in tier 1**: the resolved view holds no data; inline data goes through the store.
- **Title boxes over a Location in an SVG, a cap on titles**: titles are on items; their cost is the item count.
- **The card's sections internal** (extensibility.md's "chrome parts ... internal"): the application composing the
  card is the outside user.
- **One word list per row**: today's sentence is phrases, not rows; changing it changes every announcement.
- **Mol\*'s priorities**: card providers order as header entries do, one registry rule.
- **A `tooltip` channel on layers, a param on data nodes, a ref's note, a fifth root level, one tooltip per layer,
  texts joined by `" | "`**: per-layer or per-mark costs, data and presentation mixed, names and words mixed, fixed
  levels, copies; a row each reads, caps, wraps and links.

## In the specs

| Decision | Where |
|---|---|
| A tooltip answers only on rendered marks; no outline of its own | core › Tied to something rendered |
| Two ways to attach (by data, by range from a table), one result | core › Two ways to attach |
| The `tooltip` node and its params (`on`, `target`, `rows`, `source`, `schema`, `block`, `category`, `field`, `text`, `label`, `format`) | core › Language › Tooltips |
| MVS's forms mapped; a CIF category before the structure runtime reads one | core › MVS mapped |
| Schema is said, never guessed | core › `schema` |
| One mechanism with a `ranges` column; every table through the store as items | core › Tooltip tables |
| A tooltip read from a `source` | design-only: not built (the core spec's TODO, `seq.resolve.tooltip.source`) |
| Tables read with what they decorate; pending a `loading` row, failed an issue | core › Tooltip tables wanted with what they decorate |
| Assemblies and names (aliases, `locus`, `block`) | core › By range |
| Card providers: a `card-provider` registry entry, ordered as header entries | viewer › Providers |
| The viewer's sections as built-in providers; `whole` the place's; the floor `seq:item` and `seq:place` | viewer › Providers |
| What stays where: the caps and `hover.card` the viewer's | viewer › Configuration |
| `seq:tooltips`: the mark's place, by data, by range, rows and the shared cap, stacking | viewer › Tooltips |
| The place and rows in the core's terms | core › The place |
| Inline Markdown on the template only; values escaped; cleaned after parsing | core › Texts |
| Links act only through Seq\*'s link grammar | core › Texts |
| Plain readers (`markdownInlineText`) and `rich` for an adapter | viewer › Tooltips |
| Bounds: a template of 1,000 characters, `hover.tooltipLines` | core › Texts |
| The enterable card (a card with a link takes the pointer, a grace) | viewer › Behaviors |
| The cursor's card's links in its context menu | design-only: not built |
| The SVG export's `<title>`s (`export.titles`, `card`) | core › A view's `export.titles` |
| Export titles composed by the viewer's providers | viewer › Export titles |
| The builder's `views.tooltip`; `check()` reports a missing field | core › `field`, `text`, `label`, `format` |
| Performance: a warm card show allocates within a few hundred bytes of the card without tooltips | viewer › Tooltips |
| Performance: the lookup through the store entry's interval index; parses cached by text (64) | core › Render › Tooltips |
| Performance targets (1,000,000 ranges, export cost) | design-only: measurement targets for the heavy lane, not rules |
| Defaults: no `tooltip` node, the card as before | viewer › Providers |
| The matrix's member labels (`member_labels`) | design-only: not built |

## Tensions

- **An enterable card vs the next hover**: only a card with a link takes the pointer.
- **Escaped values vs a notes file's Markdown**: an MVS file's `**` shows as written. Open question.
- **A table read for marks it decorates**: one more request per decorated window; bounded by the same tiles.

## Edge cases

New ids continue the catalog (EC-245 to EC-249 are layouts.md's).

- **EC-250** an `on` or a table address that no longer resolves: `seq.lang.tooltip.on`, that row inert, the rest opens.
- **EC-251** a template field the table lacks: `seq.lang.tooltip.field`; undefined or empty in a row: no row.
- **EC-252** a table not fitting its `schema`: `seq.lang.tooltip.schema`; a CIF category before the structure runtime
  reads one: `seq.lang.tooltip.category`.
- **EC-253** a wrapped range on a cyclic system (EC-9): two runs, one row.
- **EC-254** a note over an empty band, a gap or a collapsed break: none there; on the marks around it, as covered.
- **EC-255** a tooltip `on` a hidden or filtered layer: no marks, no rows.
- **EC-256** bidi overrides, controls (entities too), `<script>`, `[x](javascript:…)`, `[x](/refs/dbd)`: removed,
  text, or the link's text; a 10,000-character template refused at check.
- **EC-257** a table pending under a mark: a `loading` row (EC-169); failed: the issue on the tooltip (EC-174).
- **EC-258** a grid cell in a note's range by x only: the row once; a rect by both.
- **EC-259** a million-range source zoomed out over a bin: `cardItems` rows, `+k more` with `moreCount`, 0 B warm.
- **EC-260** a note on an item said whole (a feature, a region cell, a pair): its row shows (`whole` ignored).
- **EC-261** a card with a link: entered across the gap; one without never takes the pointer; a link's hover on its own
  channel (EC-182); Esc hides it.
- **EC-262** a composition without `seq:item` or `seq:place`: `seq.viewer.card.floor`; a provider that throws: its
  rows dropped, said once.
- Related: EC-141, EC-143 (a re-hit parses nothing), EC-156, EC-166 (the bridge, now built), EC-169 to EC-180 (the
  store's states), EC-200.

## Tests

- Quick, core: the `tooltip` kind, its schemas and issues (EC-250 to EC-252); by data and by range against a
  brute-force overlap over drawn marks (a property); MVS annotation files under each schema giving MVS's runtime text
  (`FormatTemplate`); projection gene to protein; a table's tiles wanted with a decorated layer's window, kept and
  evicted with it; `tooltipsAt` headless.
- Quick, viewer: the built-in providers giving today's cards and sentences unchanged (the existing tests as is); a
  host's entry by `order`, one filtered, one replaced; the floor issue; a throwing provider; `whole` ignored
  (EC-260); escaped values; cleaning after parsing; in-view links only (EC-256); the cap and `moreCount`; the
  enterable card (fake timers, jsdom); the cursor's menu links.
- Medium: the HC sweep on the default spec with a tooltip by range on each system: every card on a mark holds its
  row, none on an empty band; the SVG with `card` titles over every story (items only).
- Perf and heavy: the hover step at 0 B warm; a cold and a warm card show recorded; the lookup's scaling; the export.
- A story: p53 with a curator's inline notes (Markdown, a link to a ref), hotspots by data, an MVS annotation file;
  the Genomes dataset with a BED's names as notes by region; a host's provider.

## Consequences

- **Language 0.9.0**: the `tooltip` kind (parents `views`, `view`; `on`, `target`, `rows`, `source`, `schema`,
  `block`, `category`, `field`, `text`, `label`, `format`) and a view's `export.titles`; with the separate change, the
  matrix layer's `member_labels` (or the next minor after any lane that lands first).
- The viewer: the `card-provider` registry kind and `BuiltInCardProviders` (`describe.ts`'s sections moved, rows and
  words unchanged); `HoverCardRow.rich`; `HoverCard.moreCount`; config `hover.tooltipLines` (3).
- The core: `tooltipsAt`, the shared interval index, a tooltip table's wants beside its decorated layers',
  `renderSvg`'s `titles`, `SpecBounds.tooltipText` and `tooltipRows`.
- Designs to update with the change: [security.md](security.md) (tooltips' Markdown, cleaning after parse, links
  through commands only), [extensibility.md](extensibility.md) and [extension-openings.md](extension-openings.md)
  (card providers opened), [data-store.md](data-store.md) (wants for tables that decorate).
- MolViewSpec converts both ways, as mapped above.

## Answered (the author, 2026-11-22: "Yes, as recommended")

1. An item's own `tooltip` field shows only when a `tooltip` node names it (MVS's rule): a column is not a promise to
   show it.
2. The SVG export's titles default to item names (a figure's size as today); the composed card rows when asked
   (`export.titles`).
3. No external links in version 1: links act only through commands (security.md); in-view `!seq-…` links only.
4. Values in templates are escaped: a value is data (`A*02:01` stays text); Markdown applies to the template.
