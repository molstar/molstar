# Tooltips: an author's words on places, in the hover card

**Draft for the author's review** (2026-11-22): not yet agreed, not in the index; its open questions are at the end.

## Intent

The author, 2026-11-22: "we should be able to associate arbitrary tooltips with data points. MolViewSpec has this too
with .tooltip() ... we should find a good way to do this here too"; then "basically, we should be able to associate
'tooltips' with addresses in general. ... or rather a tooltip by range, etc... the situation here is a bit different. i
think we should design this first." In passing: a matrix layer should optionally show its member (sample) names, as a
categorical grid axis does.

The hover card says what the viewer knows about a place (its numberings, its letter, its item, its value). A tooltip
says what the writer knows: "the tetramerisation interface", "R175H, a hotspot", "a note from the curator", one text
or a column of a million. It is data in the spec, so the canvas, a headless card, the live region and another reader
say it the same.

## Precedents

| Precedent | Taken | Not taken |
|---|---|---|
| MolViewSpec (`packages/mvs/builder/src/tree/mvs/mvs-tree.ts`; runtime `annotation-tooltips-prop.ts`, `custom-tooltips-prop.ts`) | `tooltip { text }` on a component (a selection) and `tooltip_from_uri` / `tooltip_from_source` over an annotation's rows: `field_name` (default `tooltip`) and `text_format`, an f-string over the row's fields (`{}` the field, `{score:.2f}`); a row with a field undefined gives no text; several tooltips on one place all shown | the texts joined into one string with `" | "` (here a card row each); a tooltip a child of each component (here one node, `on` an address); Mol\*'s loci labels as HTML (here text) |
| `@molstar/core/util/string-format` (`FormatTemplate`) | the f-string grammar MVS uses, reused as is (variables only, format specs `deEfF%cs`, `{{` `}}` literals; any undefined variable: no text) | a second template language |
| Themes ([themes.md](themes.md)) | data, never code; a rule's `in` clause (a mark whose position or extent overlaps a location by one position); a node under `views` shared by every view (`scale`); compiled once in tier 1, looked up per hover | a channel per layer: a tooltip is not a mark's look (Rejected) |
| Refs ([addressing.md](addressing.md)) | `to` as an address or a list of addresses; a node whose address no longer resolves is reported and the rest opens | a ref's `label` doubling as a note (Rejected) |
| The hover card ([ux-foundations.md](ux-foundations.md), viewer `describe.ts`, `hover-card.ts`) | one description, the card's rows and the live region's sentence; sections in order; item rows capped at `hover.cardItems` with `+k more`; the port draws text | — |
| Genome browsers (IGV's popup text, UCSC's `mouseOver` field, JBrowse's `description`) | a text per feature from a column, shown on hover; a bigBed's extra column as the text | HTML popups from data (UCSC's `mouseOver` takes markup: an injection path) |

## Decision

### What a tooltip attaches to: an address

- **A tooltip is a text on what an address names** (`on`): a Location (`(range …)`, `(locus …)`, `(rect …)`, a ref
  to either), Keys (`(item :type domain)`, `(row :species human)`, `(elements …)` through its groups), or Nodes (a
  data node of items: `ranges`, `pairs`, `matrix`, an alignment's rows; a layer). Why one param: an address and a query
  are one language ([locations-and-queries.md](locations-and-queries.md)), and the result type already says whether
  the text is one for a place or one per item.
- **One text for a place, or one per item**: on a Location the text is one string for the whole of it; on Keys or a
  data node each item gets its own text from its fields (`{}` its `field`, `{label}`, `{score:.1f}`). So a table of
  ranges with text is a `ranges` node (inline rows with a `tooltip` column, or a source read by region), a column of a
  source is its `field`, and a query (`(item :type domain)`, `(filter :test …)`) picks which items speak. Why items
  through data nodes: a `ranges` node already has sources, region reads, maps, assemblies and names
  ([genome-coordinates.md](genome-coordinates.md)); a second table of "ranges with text" would duplicate all of it.
- **On a layer** (`on: /views/main/p53/features/f`): its marks only, its items' fields, its filter applied. Why: the
  one place a writer wants a note on one drawing of data and not another; the address is the scope, no extra param.

### Where it lives: a `tooltip` node under `views` or a view

```
{ kind: "tooltip", id: "notes", params: { on: "/data/p53/notes", field: "tooltip", text: "{}", label: "Note" } }
```

| Param | Type | Default | Says |
|---|---|---|---|
| `on` | ref, or a list of refs (one result type) | required | what it answers on |
| `text` | string, an f-string (`FormatTemplate`) | `"{}"` | the text; on a Location a literal (`{{` for a brace) |
| `field` | string | `tooltip` | the field `{}` reads (MVS's `field_name`) |
| `label` | string | none | the card row's label; absent, the text is the row (no label) |

- **Parents `views` and `view`**: under `views` it answers in every view (beside `scale`, shared as a scale is);
  under a view only there. Ids required (`/views/notes`, `/views/main/notes`), so a story step adds or removes one by
  patch. Why under views: composition's levels are fixed (data, views, refs, snapshots); a note is how a figure
  explains data, the views' part, while what it names stays an address into data.
- **GA's assemblies and names**: `on` takes `(locus :text "chr17:7,668,402-7,687,550")` and `(block :name …)`; an
  assembly's `ranges` (BED's name or a bedGraph's column by `source.fields`) is a tooltip's data by region; a name that
  is not an id is reached by its form, never by loosening the id pattern (EC-200).
- **Resolved in tier 1** ([composition.md](composition.md)): the resolved view lists each view's tooltips with `on`
  resolved to a Location, a key set or a data node, the template parsed and its fields checked against what the data
  declares (an inline table's columns, a source's `fields` and `types`). Data read by region is never in tier 1.

### What the reader sees

- **The place a tooltip matches is the hit's** (the hover contract's): the item under the pointer by its key; else the
  position under the pointer (a cell: both positions; a pair: both ends; a matrix cell: its position and its state; a
  row header: its row); else a bin's stretch. A Location matches by overlap of one position or more (the themes' `in`),
  a cell or a pair by either of its two places, a rect only by both. Keys match the hit's element, and a position inside
  a key's extent. **Projected**: the place on the tooltip's system through the hover's own projection (a note on the
  protein answers on its gene's codons), never a second projection.
- **On every mark whose place it covers, in every layer** (letters, ranges, values, a ruler, a track's empty band),
  except a tooltip `on` a layer (its marks). Why: a note on residues 94 to 292 is about those residues wherever they
  are drawn; a mark-less position is still a hit under the contract, so a note answers there too.
- **Its own item first**: when the hit is an item of the tooltip's data, that item's text alone; else every item
  over the place, in position order, each an item row (`hover.cardItems` cap, `+k more`).
- **Stacking**: every tooltip that matches adds its rows, `views`' tooltips before the view's, each in spec order. No
  tooltip hides another or a built-in row.
- **Order in the card**: a new section after the place's (after the item, pair, grid, residue and place rows) and
  before the block's: what the viewer knows first, the writer's words next, the chain last. The section ignores
  `whole` (an item's card gets its notes too). A run of small blocks says its block and the run alone (no tooltips).
- **The live region** says the same words in the same order: `…, Note: tetramerisation domain` (`label: text`, or the
  text alone), one sentence with the card's (the card's `text` stays the announcer's).
- **Zoomed out**: a bin or a band matching many items gives the first `hover.cardItems` in position order and a
  `+k more` row with the count from the index (k exact for inline data and loaded tiles; past a tile, "more"); the live
  region says the rows shown and the count, never a million texts.
- **Pending**: a tooltip's data not loaded at the place gives a row `label: loading` (never nothing, EC-169); the
  card re-describes when its tile lands while the hover rests.
- **No emphasis of its own**: the hover's emphasis is the hit's (one hover, one outline); a tooltip adds words only.
- **The keyboard cursor** is a place like a hover (`cursorPlace`): its card and announcement carry the tooltips.

### Text

- **Plain text**: never markup, never a link, in this version (security.md: "labels from data never become markup").
  The default adapter writes `textContent`; a host's adapter receives plain strings, documented as untrusted.
- **Templates are MVS's f-strings** through core's `FormatTemplate`; any variable undefined gives no row (MVS's
  rule); a number with no format spec is said as the card says numbers (`valueText`: up to three significant
  decimals); a value is inserted as text, never read again as a template.
- **Cleaned**: C0 and C1 controls (but `\n`, a line break in the card and a pause in the sentence) and the bidi
  overrides and isolates (U+202A to U+202E, U+2066 to U+2069) are removed, so a note cannot reorder the card's other
  text. Why: a spec and its data are untrusted ([security.md](security.md)).
- **Bounds**: a template at most 1,000 characters (`SpecBounds.tooltipText`, `seq.lang.tooltip.text`); a formatted
  text cut at 1,000 with `…`; a row wraps to at most `hover.tooltipLines` (3) lines in the card, then `…` (the live
  region says it to the bound). Tooltips in a spec count in its size.

### The builder

Mirrors MVS's `.tooltip()` where it fits: a tooltip is added beside what it names, not under it.

```ts
const notes = data.ranges('notes', { target: p53, rows: { start: [93, 324], end: [292, 355], tooltip: ['DNA-binding', 'tetramerisation'] } });
views.tooltip('notes', { on: notes, label: 'Note' });                     // MVS's tooltip_from_uri: data + field
const dbd = refs.ref('dbd', { to: `${p53}/(range :from 94 :to 292)` });
views.tooltip('dbd', { on: dbd, text: 'binds DNA' });                     // MVS's tooltip { text }
view.tooltip('hotspots', { on: `${variants}/(item :hotspot true)`, text: '{label}: {count:d} tumours' });
```

`views.tooltip` and `view.tooltip` (on the views and view builders, beside `scale`) return the node's address (a
string, as every builder address); `check()` reports a template field the inline table lacks. A literal `text` on a Location is written with braces escaped by the builder (`tooltipText('a {b}')`).

### The hover contract (HC) and the card's placement

- **A tooltip is never a hit and never a mark**: it adds rows to the card of whatever hit the contract gave. The
  sweep stays as it is (every mark: a hit, an emphasis, a non-empty card), and gains one property: with a tooltip `on`
  `(all)` of each system, every swept card holds its row.
- **Rows carry `tooltip: true`** (`HoverCardRow`), so a host's card styles them apart; `cutItems` counts a tooltip's
  item rows with the place's (one cap). The card keeps its place rule (`placeCard`, never over the anchor); a tooltip
  row only makes it taller, and the adapter places by its measured size (already its rule). A drag suspends it.
- **Off by a host**: `hover.tooltips` (true); the spec never says how a reader's card looks.

### Performance

- **Resolve, once**: a Location's tooltips as one interval index per system (the runtimes' `intervalIndex`, moved to
  a shared internal module of the core with its test: short ones by start, long ones apart); Keys as a set by owner and
  key (core's hash); an inline data node's items through the index its store entry builds when it loads.
- **Read by region**: a source's tiles hold their own interval index, built when a tile loads and evicted with it
  (performance.md: summaries belong to their entry); a hover reads the tile under the place. Nothing over a million
  ranges is built at resolve.
- **Per hover step, 0 B warm**: the lookup writes matching row numbers into a reused `Int32Array` scratch, stopping at
  `cardItems + 1` and counting the rest by the index; O(log n + k). Text is formatted only for the rows a card shows,
  when it rests, and kept in a small cache by (tooltip, row) (64 entries), so a resting hover re-hit each frame
  (EC-143) allocates nothing. The card's allocation per step today (2,168 B) is not raised.
- **Targets**: the lookup over 1,000,000 ranges within 1.5× of 10,000 (heavy, `itWithin`); a hover step at 0 B warm
  with a 10,000-range inline tooltip and a 1,000,000-range source (the allocation tier).

### The matrix's member labels

- **`member_labels`** on a `matrix` layer: `none` (when absent) or `left`, the states' names down the track's label
  column, one a state row, under the track's name line. Why off by default: a cohort of thousands of states is
  sub-pixel rows, thinned to nothing; the column is the track's header; today's records stay.
- **Drawn as the categorical axis draws them**: `placeMemberLabels`, `drawMemberLabels`, `hitMemberLabel`
  (`render/layout/grid-members.ts`), thinned past `rows.labelMinPx`, following the band's `bandTop` and `rowScale`,
  cut to the column's width (the card says the name whole). A label hits as `target: 'row'` on the states system: the
  interaction matrix's "a state (a matrix row's header)" (click selects the state).
- **On a ring** not drawn (`seq.render.matrix.labels`, info, once): a ring has no label column.

## Defaults

No `tooltip` node: the card and the live region as today. With one, a row after the place's rows:

```
┌──────────────────────────────┐
│ p53 175                      │
│ index 174 · uniprot 175      │
│ Letter R                     │
│ Note  R175H: a hotspot,      │
│       structural class       │
│ Chain A                      │
└──────────────────────────────┘
```

A matrix draws no member labels unless `member_labels: left`.

## Rejected alternatives

- **A `tooltip` channel on layers** (`tooltip: (field :name note)`, rules by `in`): per layer, compiled per mark in
  the layout's hot loop though read only on a hover; silent on other layers and on positions with no mark; a rule is
  not a template. Themes color and size; a tooltip is words.
- **A param on data nodes** (`ranges.tooltip: "{note}"`): data is defined once and drawn by many views; a note is
  presentation, and a Location's note has no data node.
- **A ref's `label` or a `note` on refs**: a ref names a place; a note is a second thing. `on` takes refs, so a named
  region with a note is two nodes and one name.
- **A fifth root level, `tooltips`**: the levels are fixed; scope by `views` or a view covers what it would.
- **Under each layer or track, as MVS's under each component**: one note copied per drawing; `on` a layer scopes it.
- **Markdown and links now**: the card cannot be entered (EC-166's bridge is untested), so a link in it cannot be
  clicked; markup from data is what security.md refuses. Open question 2.
- **Joining texts into one string** (MVS's `" | "`): a row each reads, caps and wraps; a host styles rows.
- **A Mustache or `${}` template**: core has MVS's f-strings; an MVS spec's `text_format` converts unchanged.
- **An index over a source's ranges at resolve**: tier 1 holds no data; tiles index themselves.

## Tensions

- **Compact cards vs notes**: the author asked for more compact cards (2026-11-20); a note wraps up to three lines.
  Bounded by `hover.tooltipLines`, and a note is the writer's choice.
- **`views` scope vs a figure**: a note under `views` speaks in a published figure's view too; place it under the
  browsing view to keep it out.
- **Keys by overlap**: a key's note answers at positions inside its extent, so a long gene's note answers on every
  base; a layer `on` or a filter narrows it.

## Edge cases

New ids continue the catalog (EC-245 to EC-249 are layouts.md's).

- **EC-250** an `on` that no longer resolves: `seq.lang.tooltip.on` (warning), the tooltip inert, the rest opens
  (as refs, addressing.md); a patch never re-binds it.
- **EC-251** a template naming a field the data does not declare: `seq.lang.tooltip.field` at check; at hover a
  field undefined in a row gives no row (MVS's rule); an empty formatted text gives no row.
- **EC-252** two tooltips on one place, one `views`', one the view's: both rows, `views`' first; the same text twice
  is said twice (two writers' rows).
- **EC-253** a note on a cyclic system's wrapped range (EC-9): two runs in the index, one row when either matches.
- **EC-254** a note over a collapsed break or a hidden piece: answers where a position of it is shown; the break's
  chip says nothing of it.
- **EC-255** a tooltip `on` a hidden or filtered layer: no marks, no rows; on its data node instead, it still answers
  by position.
- **EC-256** a text of bidi overrides, controls or `<script>`: removed or shown as text (jsdom); a 10,000-character
  text refused at check, a 10,000-character value cut at 1,000.
- **EC-257** a pending tile under a note's place: a `loading` row, replaced while the hover rests (EC-169).
- **EC-258** a grid cell whose x is in a note's range and y not: the row once; a rect note only when both are in.
- **EC-259** a million-range note source zoomed out: `cardItems` rows and `+k more`, the hover step 0 B warm.
- **EC-260** a matrix with `member_labels: left` on a ring: not drawn, said once; with more states than room: thinned.
- Related: EC-141 (one hover owner: a stale clear removes no row), EC-143 (a resting re-hit formats nothing), EC-166
  (the card cannot be entered: no links), EC-169 (pending), EC-200 (names through forms), EC-182 (a link's hover
  channel: a link to a noted range shows no card, it is not the pointer's hover).

## Tests

- Quick, core: `tooltip` params and issues (EC-250, EC-251); each `on` result type matched against a brute-force
  overlap over random places (a property), the wrapped and reversed cases; projection gene to protein; the template
  over MVS's fixtures' rows equal to MVS's runtime text (`FormatTemplate`, one implementation); the interval index moved
  with its tests.
- Quick, viewer: the card's rows and order and the live region's sentence equal (one description); stacking; the cap
  and `+k more`; the pending row; cleaning and bounds (EC-256, jsdom adapter); the keyboard cursor's card.
- Medium: the HC sweep with a tooltip on `(all)` of each system: every card holds its row; the matrix labels' record,
  thinning and hit (canvas and SVG), the ring's issue.
- Perf and heavy: the hover step at 0 B warm (the allocation tier); the lookup's scaling ratio, 1,000,000 against
  10,000 ranges (`itWithin`).
- A story: p53 with a curator's notes table, hotspot notes from a query, and a note on the DBD ref; the Genomes
  dataset with a BED's name column as notes by region.

## Consequences

- **Language 0.9.0**: the `tooltip` kind (parents `views`, `view`) and the matrix layer's `member_labels`; the core
  `_spec` Versions list and the language description gain them (or the next minor after any lane that lands first).
- `HoverCardRow.tooltip`; config items `hover.tooltips` (true) and `hover.tooltipLines` (3); `SpecBounds.tooltipText`.
- The runtimes' `intervalIndex` moves to a shared internal module of the core; the resolved view carries tooltips, so
  another reader (Nightingale, JBrowse) can show them.
- MolViewSpec converts both ways: `tooltip { text }` on a component is `on` its selection with `text` escaped;
  `tooltip_from_uri` is a `source`, its items and a tooltip with `field` and `text` (`text_format`).
- Not now: tooltips in an SVG export (open question 5); a structure viewer showing Seq\*'s notes through
  `linkStructure` (it can be written as MVS's tooltips later).

## Open questions for the author

1. An item's own `tooltip` field with no `tooltip` node: shown by default, or only when a node names it (MVS's rule,
   recommended: a column is not a promise to show it)?
2. Markdown and links: plain text now (recommended), or a subset (`**`, `*`, `` ` ``, links to addresses only) with the
   card made enterable (EC-166's bridge) in this slice?
3. The rows' place in the card: after the viewer's own rows (recommended), or first, under the title, the writer's
   words leading?
4. Should a note by range outline its range while it answers (a second, quieter emphasis showing what it covers), or
   stay words only (recommended: one hover, one outline)?
5. In an SVG export, write each mark's tooltips as an SVG `<title>` (a browser shows it over the figure), or leave
   exports without them (recommended for now)?
