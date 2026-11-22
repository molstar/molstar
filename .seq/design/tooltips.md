# Tooltips: the card composed by the application, the spec's words on places

## Intent

The author, 2026-11-22: "we should be able to associate arbitrary tooltips with data points. MolViewSpec has this too
with .tooltip() ... we should find a good way to do this here too"; then "basically, we should be able to associate
'tooltips' with addresses in general. ... or rather a tooltip by range, etc... the situation here is a bit different. i
think we should design this first." In passing: a matrix layer should optionally show its member (sample) names, as a
categorical grid axis does.

On the first draft (2026-11-22): "i am missing the 'tooltip from data/uri' option here"; "i think we can do markdown? we
do have packages/core/src/util/markdown.ts"; "the placement: maybe configuration on 'viewer spec'? after all, mol*
plugin also composes tooltips via providers. maybe we should adopt similar principle? in general, tooltips should be
composed by 'application layer' anyway"; on SVG export, "not sure we need this", then "or well .. maybe <title> in SVG
could work too ...".

So two things, kept apart: **the spec says which words attach where** (data: portable, drawn the same by every reader),
and **the application composes the card** (which providers speak, in which order, how many rows, how links act).

## Precedents

| Precedent | Taken | Not taken |
|---|---|---|
| MolViewSpec (`mvs-tree.ts`; runtime `annotation-tooltips-prop.ts`, `custom-tooltips-prop.ts`) | `tooltip { text }` on a component; `tooltip_from_uri` / `tooltip_from_source`: a separate annotation table (`uri`, `format` json, cif or bcif, `schema`, `field_name` default `tooltip`, `text_format` an f-string, `field_remapping`); a row with a field undefined gives no text; every matching tooltip shown | texts joined into one string by `" | "` (here a row each); a tooltip a child of each component (here one node with an address) |
| Mol\*'s plugin (`plugin/core/src/spec.ts` `lociLabels`, `state/manager/loci-label.ts`) | hover labels composed from **providers listed in the plugin spec**, each given the hovered thing and returning a label or nothing; the application chooses which run (MVS's providers added by its behavior) | a numeric `priority` for order (here the listing's order); labels as HTML or React elements (here typed rows) |
| Seq\*'s viewer spec (`viewer/src/default-spec.ts`, `header-entries.ts`) | the viewer composed from data: registry entries in order, behaviors with params, config items; `header-action` entries as the precedent for a viewer-only registry kind a host adds to or filters | — |
| `@molstar/core/util/markdown.ts` | a typed tree, never HTML, linear in its input, nesting capped; `parseMarkdownInline`, `isSafeMarkdownHref` (http, https, mailto, relative), `markdownInlineText` for plain readers | its blocks in a card row (Text) |
| `@molstar/core/util/string-format` (`FormatTemplate`) | MVS's f-strings, one implementation | a second template language |
| Themes ([themes.md](themes.md)), refs ([addressing.md](addressing.md)) | a rule's `in` (overlap by one position); nodes under `views` shared by every view; an address that no longer resolves reported, the rest opening | a tooltip as a mark's channel (Rejected) |
| The SVG export (`core/src/render/svg.ts`) | an accessible item is a named box whose `<title>` is the browser's tooltip and its accessible name | — |

## Decision

### The spec: which words attach where

A `tooltip` node, under `views` (every view, beside `scale`) or under a `view` (that view only). Ids required, so a
story step adds or removes one by patch. One mechanism: **a tooltip is a set of items, each with a location and its
fields; a text is formatted from a row's fields**. The items come from one of three spellings:

| Spelling | Params | Like MVS's |
|---|---|---|
| on an address | `on`: a ref or a list (a Location, Keys, a data node of items, a layer) | `tooltip { text }` on a component |
| a table of its own, inline | `target`, `rows` (a table: location columns and text fields) | `tooltip_from_uri` over an inline annotation |
| a table of its own, read | `target`, `source` (a `source` node: `uri`, `format`, `fields`, `index`) | `tooltip_from_uri`, `tooltip_from_source` |

Common params: `field` (`tooltip`: what `{}` reads, MVS's `field_name`), `text` (`"{}"`: an f-string, MVS's
`text_format`; on a bare Location a literal, `{{` for a brace), `label` (the card row's label; absent, the text is the
row), `format` (`markdown`, absent, or `text`).

- **Location columns of a tooltip's table**, decided by the columns present (`seq.lang.tooltip.locate` when none or
  two kinds):
  - **ranges**: `start`, `end` on `target` (a sequence, a unit, an alignment's columns), or on an assembly with a
    sequence name column (`target` role through `source.fields`, `chrom`, aliases and `chr` by GA's names rule);
  - **addresses**: an `address` column, each row an address or a form relative to `target`
    (`/refs/dbd`, `(range :from 94 :to 292 :numbering uniprot)`, `(elements :chain A :beg 10 :end 20)`);
  - **keys**: a `key` column naming items of `target`, a data node (a feature's id, an alignment row's id, a matrix
    state);
  - **MVS's schema fields** (`label_asym_id`, `beg_label_seq_id`, `auth_seq_id`, …, by the short and long names of
    [properties.md](properties.md)): each row an `(elements …)` item, so an MVS annotation file is a tooltip table as
    is.
- **Read whole or by region**: a ranges table over an indexed format (BED, GFF, VCF with tabix; bigBed later) is read
  by region as any `ranges` source; addresses, keys and schema rows are read whole (they resolve before they can be
  indexed), bounded by the asset policy's size of one asset and `SpecBounds.tooltipRows` (1,000,000).
- **Assets**: a `source` is read through the assets port's policy ([security.md](security.md)): `https:` and relative
  locations in a browser, no credentials, nothing over the network headless; a refused read is an issue on the
  tooltip, and the rest of the spec opens.
- **One mechanism, not two, with a `ranges` node**: `on` a `ranges` node whose rows carry a `tooltip` column and a
  tooltip's own `rows` or `source` resolve to the same thing (an items table, one interval index, one region read,
  one projection), and a tooltip's own table *is* a `ranges` shape with an extra location kind or two. Why both
  spellings: a table that is also drawn belongs in `data` and is named by `on`; a notes table that is never drawn
  (MVS's case) stays with its tooltip, and an MVS spec converts one node to one node. A tooltip's own table is not
  addressable as data (`/views/notes` is the tooltip) and never drawn.
- **On a layer** (`on: /views/main/p53/features/f`): its marks only, its filter applied: the one scope finer than a
  view, given by the address.
- **Resolved in tier 1** ([composition.md](composition.md)): each view's tooltips with `on` resolved, inline tables
  indexed, templates parsed and their fields checked against what the table declares. Data read by region is never in
  tier 1.

### The application: the card composed by providers

- **A `card-provider` registry entry** in the viewer spec (beside `header-action`, a viewer-only kind):
  `{ kind: 'card-provider', name: '<ns>:<name>', before?: name, describe(at, out) }`. `at` is the hit and its place
  (the section, the layer, the system and position, a cell's y, the element, a bin, the block) and the models;
  `out` takes rows (`label`, `value`, `rich?`, `item?`, `words?`) and may mark the place said whole (an item's own
  card). A provider that has nothing to say writes nothing.
- **The viewer's own rows become built-in providers**, in `BuiltInCardProviders`, in today's order: `seq:region`,
  `seq:pair`, `seq:group`, `seq:item`, `seq:columns-hidden`, `seq:row`, `seq:grid`, `seq:residue`, `seq:place`
  (numberings, letter, chemistry, values, items at the position), `seq:tooltips` (the spec's), `seq:block`. The
  default spec lists them; `describe.ts`'s sections become these entries without a change of words.
- **Order is the listing's**: a host spreads `BuiltInCardProviders` and inserts, filters or replaces by name; an
  extension's entry says `before` a name, else it comes last but `seq:block`. Why not Mol\*'s numbers: two hosts'
  priorities collide, and an order in a list reads in the spec.
- **What stays in the spec, what in the viewer spec**: the spec says what text attaches to which places (portable:
  another reader, a headless card, an export reads it). The viewer spec says which providers run, their order, the
  caps (`hover.cardItems`, `hover.tooltipLines`), whether the card shows (`hover.card`), how links act, and the export's
  titles. A host that wants the spec's tooltips in its own card replaces `seq:tooltips` and calls the core's
  `tooltipsAt(view, place, out)`; one that does not want them filters it out.
- **The hover contract holds** (HC: every mark answers with a hit, an emphasis and a card): providers change the
  card's rows, never the hit or the emphasis. The floor is `seq:item` and `seq:place` (every hit names an item or a
  position): a composition without either is `seq.viewer.card.floor` (a warning at composition); the HC sweep runs on
  the default spec, and on a host's in its own tests.
- **One description**: the card's `text` (the live region's sentence) is every row's words in order; a row without
  `words` is said as `label: markdownInlineText(rich)`, or its value. The keyboard cursor's card is composed the same.

### What the reader sees of the spec's tooltips (`seq:tooltips`)

- **The place a tooltip matches is the hit's**: the item under the pointer by its key; else the position under the
  pointer (a cell: both positions; a pair: both ends; a matrix cell: its position and its state; a row header: its
  row); else a bin's stretch. A Location matches by overlap of one position or more, a cell or a pair by either place,
  a rect only by both; Keys match the hit's element and a position inside a key's extent. **Projected** through the
  hover's own projection (a note on the protein answers on its gene's codons), never a second one.
- **On every mark whose place it covers, in every layer** (a mark-less position is still a hit), except a tooltip
  `on` a layer. When the hit is an item of the tooltip's table, that item's row alone; else every item over the place,
  in position order, as item rows under `hover.cardItems` with `+k more`.
- **Stacking**: `views`' tooltips before the view's, each in spec order; none hides another or a built-in row.
- **Zoomed out**: a bin matching many items gives the first `hover.cardItems` and `+k more` with the index's count;
  the live region says the rows shown and the count. **Pending**: a row `label: loading` (EC-169), replaced while the
  hover rests. **No emphasis of its own**: one hover, one outline.

### Text: Markdown, inline

- **A card row shows inline Markdown** (`format: markdown`, the default): strong, emphasis, code spans, links, line
  breaks (`parseMarkdownInline`). Blocks are not read: a row is a line of the card, so headings, lists, fences and rules
  stay their text. Why inline only: the card is compact (the author, 2026-11-20); a note with structure is a link to
  where it is written. `format: text` shows the text as written.
- **Values are inserted as written**: a notes file's Markdown renders as its writer meant (MVS's file tooltips are the
  file's); the parser is safe by construction (a typed tree, no HTML, nesting capped, linear), so a value can add a
  mark or a link, never markup.
- **Links point at two things**:
  - **inside the view**: a bare address (`/refs/dbd`, `/data/p53/(range :from 94 :to 292)`: focus) or Seq\*'s link
    grammar (`!seq-focus=…`, `!seq-select=…`): run as commands, the hover of a link on its own channel (EC-182);
  - **outside**: `isSafeMarkdownHref`, and allowed by the assets policy's `schemes` and `hosts` (the viewer's one
    list of where a spec may reach); opened by the host (`openLink` on the hover card port, default a new tab with
    `noopener noreferrer`); headless, none. Any other link is its text.
- **Plain readers**: the live region, a headless card's `value` and an export's `<title>` read
  `markdownInlineText`; a link is said as its text. `rich` holds the tree for an adapter that draws it.
- **Parse once**: a literal text is parsed at resolve; a templated row is formatted and parsed when a card first shows
  it, kept in a cache by (tooltip, row) (64 entries), so a resting hover re-hit each frame (EC-143) parses nothing.
- **Cleaned and bounded**: C0 and C1 controls (but `\n`) and the bidi overrides and isolates (U+202A to U+202E, U+2066
  to U+2069) removed; a template at most 1,000 characters (`SpecBounds.tooltipText`); a formatted text cut at 1,000;
  a row wraps to `hover.tooltipLines` (3) lines, then `…`.

### The enterable card (EC-166's bridge)

- **A card with a link is entered**: the pointer may cross from the anchor to the card within a grace (300 ms along
  the path between them, the safe triangle), and the card stays while the pointer is over it; the hover stays on the
  anchor (frozen) meanwhile. Leaving the card hides it. A card with no link takes no pointer (`pointer-events: none`),
  so it never steals a hover from the data. Why only then: a card in the way of the next hover is the cost of a link,
  never of words.
- **From the keyboard**: the cursor's card's links are listed in the cursor's context menu (Shift+F10, EC-156), so a
  link is reached without a pointer.

### The SVG export's `<title>`s

- **An option**: the export may write each exported item's composed rows as its `<title>` (the providers' output as
  plain text, `markdownInlineText`), in the named box the accessible export already writes ("a box nothing paints and
  a pointer still hits: its `<title>` is the tooltip and the accessible name"). A browser then shows the figure's
  notes on hover, and a screen reader says them; the default title stays the item's name.
- **Where it goes**: an item (a feature, a pair, a row) a title; a spec tooltip on a Location with no item, one named
  box over its extent per track; bins, image tiles and summaries no title (they are not items, and a bin's rows depend
  on the zoom the figure froze).
- **Cost, bounded**: a provider call per exported item, once per export, never per frame; at most
  `export.titles.max` (10,000) items titled, past it the item's name only, said once (`seq.render.export.titles`).
  Size: about the rows' text per item (a 60-character note on 10,000 items, 0.6 MB); the cap holds it.
- **Composition**: `renderSvg`'s `titles` option takes a describe function; the viewer's export passes its composed
  providers; the core's headless export passes `tooltipsAt` alone (the spec's words, no viewer rows). The view's
  `export.titles` (`name` or `card`) is the spec's say, as `export.legend`. The default is an open question.

### The builder

Mirrors MVS's `.tooltip()`, `tooltip_from_uri` and `tooltip_from_source`:

```ts
const notes = data.ranges('notes', { target: p53, rows: { start: [93, 324], end: [292, 355], tooltip: ['DNA-binding', 'tetramerisation'] } });
views.tooltip('notes', { on: notes, label: 'Note' });                                  // a drawn table's column
views.tooltip('dbd', { on: refs.ref('dbd', { to: `${p53}/(range :from 94 :to 292)` }), text: 'binds **DNA**' });
views.tooltip('curated', { target: p53, rows: { address: ['(at :label 175)'], tooltip: ['R175H, a [hotspot](/refs/hotspots)'] } });
views.tooltip('mvs', { target: tup, source: data.source('ann', { format: 'json', uri: 'ann.json' }), field: 'note', text: '{}' });
```

`check()` reports a template field the table lacks. Viewer side: `{ ...DefaultSeqViewerSpec(), registry: [...,
{ kind: 'card-provider', name: 'lab:assay', before: 'seq:block', describe } ] }`.

### Performance

- **Per hover step, 0 B warm**: a step is the hit and the emphasis; providers run when the card shows (at rest) and
  again only when the place or the store's version changes; the card for a place re-hit each frame is the one held.
- **The lookup**: a Location's tooltips as one interval index per system (the runtimes' `intervalIndex`, moved to a
  shared internal module of the core with its test); keys as a set by owner and key; a source's tiles index themselves
  when they load and are evicted with them. Matches go into a reused `Int32Array`, stopping at `cardItems + 1`, the
  rest counted: O(log n + k), 0 B.
- **Targets**: the lookup over 1,000,000 ranges within 1.5× of 10,000 (heavy, `itWithin`); a hover step at 0 B warm
  and a card show within today's 2,168 B with a tooltip over 10,000 inline ranges and over a 1,000,000-range source
  (the allocation tier); an export with titles over 10,000 items within 1.2× of one without (heavy).

### The matrix's member labels

- **`member_labels`** on a `matrix` layer: `none` (when absent) or `left`, the states' names down the track's label
  column, a state a row, under the track's name line. Off by default: a cohort of thousands of states is sub-pixel
  rows, thinned to nothing; today's records stay.
- **Drawn as the categorical axis draws them** (`placeMemberLabels`, `drawMemberLabels`, `hitMemberLabel` in
  `render/layout/grid-members.ts`), thinned past `rows.labelMinPx`, following `bandTop` and `rowScale`, cut to the
  column (the card says the name whole); a label hits as `target: 'row'` on the states system (the interaction
  matrix's "a state": click selects it). On a ring not drawn (`seq.render.matrix.labels`, info, once).

## Defaults

No `tooltip` node: the card as today. With one, a row after the place's rows (links drawn as links):

```
┌──────────────────────────────┐
│ p53 175                      │
│ index 174 · uniprot 175      │
│ Letter R                     │
│ Note  R175H, a hotspot       │
│ Chain A                      │
└──────────────────────────────┘
```

The default viewer spec lists every built-in provider; an export writes item names as titles unless asked for the
card's rows; a matrix draws no member labels.

## Rejected alternatives

- **Plain text only** (the first draft): the core's Markdown is safe by construction, and a note's emphasis or link is
  the writer's meaning. Blocks in a row stay rejected (compact cards).
- **The card's sections internal** (the first draft and extensibility.md's "chrome parts ... internal"): the author's
  rule that the application composes tooltips is the outside user that opens them.
- **Mol\*'s numeric priorities**: colliding numbers, an order not readable in the spec.
- **Card composition in the spec** (which rows, their order): a reader's chrome, not the figure's meaning, and another
  reader composes its own card.
- **Two mechanisms for "a tooltip table" and "a `ranges` column"**: two indexes, two region reads, two projections
  for one thing.
- **A `tooltip` channel on layers**: per layer, compiled per mark in the layout's hot loop though read only on a
  hover; silent on other layers and on mark-less positions.
- **A param on data nodes** (`ranges.tooltip`), **a ref's note**, **a fifth root level**, **one tooltip per layer**:
  a drawn table is named by `on`; a ref names, a note speaks; the levels are fixed; `on` a layer scopes it.
- **Joining texts by `" | "`** (MVS): a row each reads, caps, wraps and links.
- **Titles on bins in an export**: a bin's rows depend on the frozen zoom and multiply its size.

## Tensions

- **A tooltip's own table outside `data`**: data is under `data` by the composition's rule; a never-drawn notes table
  stays with its tooltip for MVS's shape. Open question 3.
- **An enterable card vs the next hover**: only a card with a link takes the pointer.
- **Compact cards vs notes**: three lines, `hover.tooltipLines`.

## Edge cases

New ids continue the catalog (EC-245 to EC-249 are layouts.md's).

- **EC-250** an `on` or a table address that no longer resolves: `seq.lang.tooltip.on`, that row inert, the rest
  opens; never re-bound.
- **EC-251** a template field the table lacks: `seq.lang.tooltip.field` at check; undefined in a row: no row; empty:
  no row.
- **EC-252** a table with both `start`/`end` and `address`: `seq.lang.tooltip.locate`.
- **EC-253** a wrapped range on a cyclic system (EC-9): two runs, one row.
- **EC-254** a note over a collapsed break: answers where a position of it is shown.
- **EC-255** a tooltip `on` a hidden or filtered layer: no rows; on its data node, by position.
- **EC-256** bidi overrides, controls, `<script>`, `[x](javascript:…)`: removed, text, or the link's text; a
  10,000-character template refused at check.
- **EC-257** a pending tile: a `loading` row, replaced while the hover rests (EC-169).
- **EC-258** a grid cell in a note's range by x only: the row once; a rect note by both.
- **EC-259** a million-range source zoomed out: `cardItems` rows, `+k more`, 0 B warm.
- **EC-260** a source refused by the policy (a host not listed, headless network): an issue on the tooltip.
- **EC-261** a card with a link: entered across the gap; a card without one never takes the pointer; a link's hover
  on its own channel (EC-182); Esc hides it.
- **EC-262** a composition without `seq:item` or `seq:place`: `seq.viewer.card.floor`; a provider that throws: its
  rows dropped, said once, the card from the others.
- **EC-263** an export past `export.titles.max`: names only after it, said once.
- **EC-264** a matrix with `member_labels: left` on a ring, or more states than room: not drawn, or thinned.
- Related: EC-141 (one hover owner), EC-143 (a re-hit formats nothing), EC-156 (the cursor's menu lists the card's
  links), EC-166 (the bridge, now built), EC-200 (names through forms).

## Tests

- Quick, core: the `tooltip` kind's spellings and issues (EC-250 to EC-252); each `on` type and each location column
  kind against a brute-force overlap (a property); an MVS annotation file's rows giving MVS's runtime text
  (`FormatTemplate`); projection gene to protein; the interval index moved with its tests; `tooltipsAt` headless.
- Quick, viewer: the built-in providers giving today's cards unchanged (the existing card tests as is); a host's
  provider inserted `before`, a built-in filtered and replaced; the floor issue; a throwing provider; the live region's
  sentence equal to the rows' words; Markdown rows, links by kind and policy (EC-256), the enterable card and its
  grace (fake timers, jsdom); the cursor's menu links.
- Medium: the HC sweep on the default spec with a tooltip on `(all)` of each system: every card holds its row; the
  SVG titles over every story (one row per item, none on bins); the matrix labels (canvas and SVG).
- Perf and heavy: the hover step 0 B warm; the lookup's scaling ratio; the export with titles.
- A story: p53 with a curator's notes table (Markdown, a link to a ref), hotspots from a query, an MVS annotation file
  as a tooltip source; the Genomes dataset with a BED's names as notes by region; a host's provider.

## Consequences

- **Language 0.9.0**: the `tooltip` kind (parents `views`, `view`; `on`, `target`, `rows`, `source`, `field`, `text`,
  `label`, `format`), a view's `export.titles`, the matrix layer's `member_labels`; the core `_spec` Versions list (or
  the next minor after any lane that lands first).
- The viewer: the `card-provider` registry kind and `BuiltInCardProviders` (describe.ts's sections moved, words
  unchanged); `HoverCardRow.rich` and `.tooltip`; the hover card port's `openLink`; config `hover.tooltipLines` (3);
  extensibility.md's and extension-openings.md's lists gain card providers.
- The core: `tooltipsAt`, the shared interval index, `renderSvg`'s `titles` option, `SpecBounds.tooltipText` and
  `tooltipRows`.
- MolViewSpec converts both ways: `tooltip { text }` is `on` its selection, `text` escaped; `tooltip_from_uri` is a
  tooltip with its `source`, `field` and `text`; `field_remapping` is the source's `fields`.

## Open questions for the author

1. An item's own `tooltip` field with no `tooltip` node: shown by default, or only when a node names it (MVS's rule,
   recommended)?
2. The SVG export's titles by default: the item's name (recommended: a figure's size and its reader's hover stay as
   today unless asked) or the composed card rows?
3. A tooltip's own table (`rows`, `source` on the node, MVS's shape, recommended) or only `on` a table under `data`
   (one place for data, MVS converting to two nodes)?
4. Order of providers: the listing's with `before` (recommended), or Mol\*'s numeric priorities, so Mol\*'s and Seq\*'s
   providers read alike?
5. Should a note by range outline its range while it answers (a second, quieter emphasis), or stay words only
   (recommended: one hover, one outline)?
