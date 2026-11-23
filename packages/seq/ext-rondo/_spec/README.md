# `@molstar/seqx-ext-rondo`: the genome ring

A genome ring like Rondo's, as an extension on the public surface: ribbons between stems wide by count, stubs to blocks
not shown, 24 hues. Two entries: `.` (`rondo`, headless: the layer type, the scheme, params, `rondoRing`) and `./viewer`
(`rondoViewer`: `rondo` and its commands and binding). Why two: the CLI and `createSeqContext` compose the first and pay
nothing for the viewer's parts. It imports `@molstar/seqx-core/extension` and `@molstar/seqx-viewer`'s root only: a gap
is the contract's, brought to the core, never a reach into `src`.

The descriptors carry `version` (0.1.0, what a spec's `extensions: { rondo: '0.1' }` is read against) and `requires`
(`'0.10'`, the language version that has `section.center`). A host lists `rondoViewer` in its spec or registers it late
(`viewer.extend`); never in the minimal or default viewer. The viewer entry loads `menu.js` and `card.js` as one lazy
chunk once composed (the stem menu, the partner keys' walk, the hole's announcements, the card's lines).

## What it ships

| Entry                   | Kind                                                                                                                                                             | Fallback                                                  | Reads                                                                  |
| ----------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------- | ---------------------------------------------------------------------- |
| `rondo:ribbons`         | layer type, `scope: 'section'`, shape `pairs`, on `rings` and a `tracks` row                                                                                     | `pairs`                                                   | core's `PairsReader`: bins at a level, items under its gate            |
| `rondo:hues`            | scheme over categories: n hues round the Oklch wheel, not colour-blind safe, substitute `categories`                                                             | the layer type's default                                  | —                                                                      |
| `rondo:partner`         | command (`.`; `rondo:partner-back` `,`): the cursor to the next partner stem, the stem menu's order, said                                                        | —                                                         | the stem under the cursor                                              |
| `rondo:select-partners` | command: the far ends of the selection's ribbons added                                                                                                           | —                                                         | the selection                                                          |
| `rondo:open`            | command, bound to a double-click scoped to `rondo:ribbons`' hits                                                                                                 | —                                                         | the hit's chromosomes, set as `opens`' `include`                       |
| `rondo:cells`           | layer type, `scope: 'section'`, shape `pairs`, with core's `square` cells (nothing drawn on a strip)                                                             | none (only `rondo:square` names it)                       | core's `PairsReader` at the coarsest bins (`gate` 0) over the assembly |
| `rondo:square`          | navigator picture, `axes: both`: one `rondo:cells` track over the bound section's ribbons' data (`of_layer`, else the first `rondo:ribbons`, else `pairs` layer) | `axis` (the lattice), said (`seq.view.navigator.picture`) | —                                                                      |

**The square's cells** (`cells.ts`): core's square lays the picture's track out over the assembly's blocks (a hidden
section of its own) and draws `rondo:cells`' cells over its lattice: the read's bins summed by block pair (a row's `a`
end along x, its `b` end along y), mirrored for unsided data (one assembly), each pair in its own cell for sided data
(two); a cell's colour its count on a square-root ramp of the five `seq-ramp-*` stops over the largest; the diagonal
(links within a chromosome) included. Whatever the ring shows (`include`), the cells are the whole assembly's: the
square is where the rest is chosen. A cell's words `37 links` (`within` on the diagonal; `links loading` before the
bins), the readout's after the pair (`chr4 × chr7 · 37 links`). Not lazy: a chunk loaded when a spec names
`rondo:square` needs resolve's lazy kinds (`withViewKind`), which the extension entry does not export (a contract
change, asked). The stem menu's logic is lazy (`menu.ts`, loaded once `rondoViewer` is composed, its behavior's `loaded`
waiting; its entries `StemMenu` in `viewer.ts`).

## Stems

- **A stem is one bin of a coarser level**: `2^j` × the read's level (bins), or `2^j` (items), the largest no wider than
  `rondo.stemPx` (16 CSS px) of positions at the ring's outer circumference (a row: the strip's own px), the finest
  strip's; never device px, so the stems are the same at any DPR and in the SVG. Aligned on position 0 of each system: a
  stem is a union of whole bins, its counts exact sums. `stems: none`: the read's own level, no trumpets.
- **Its count** is the pairs leaving it (ribbons, stubs, fan); pairs with both ends in it are counted apart ("within")
  and drawn once a finer level parts them.
- **The trumpet**: its mouth on the layer's band's inner edge (the hole's rim) to its neck `rondo.stemDepth` (0.12) of
  the hole inward, two cubics flaring, in its system's hue at `ribbons.fillOpacity`; on a ring of one chromosome or a
  pair (the section shows two blocks or fewer) at `rondo.focusStems` (0.35) of it. Why: there the information is in the
  ribbons and stubs, and full trumpets dominate them (a look at the SVG of a 23-chromosome toy genome, chr4 alone and
  chr4 × chr7: 0.35 read the stubs over the stems with the block's hue kept, 0.25 nearly washed it out, 0.5 kept the
  trumpets the louder; to confirm in a pane).

## Ribbons

- **A ribbon is a stem pair a strand** (the rows' relative strand), its `a` the stem earlier in the systems' order (its
  hue). Width: `k` = `rondo.maxWidthPx` (24) / the largest stem count in view, a ribbon `k ×` its count; under
  `ribbons.minPx` device px, or past the strongest `rondo.maxRibbons` (1,024), it goes to both stems' fans. So a neck is
  the sum of its ribbons, stubs and fan; at least a device px as a whole. Why linear, not core's square root for arcs:
  arcs are strokes that never add, but a neck is a sum a reader checks; the cap (`rondo.maxWidthPx`) keeps the hole
  readable and the fan keeps the sum exact however thin the ribbons under it.
- **Drawn widest first**: the polygons are trumpets, then ribbons by count descending (ties by index), then stubs, then
  fans, the hole last (skipped in the data pass); hits are tested last drawn first. Why: the narrow lie on top, so they
  stay hittable. Ribbons are filled at `ribbons.fillOpacity` with no outline (an outline is the hover's).
- **Items are ribbons of one count**: past the gate each pair is a ribbon from its ends' own extents, still grouped into
  its stems' necks, so the bundle holds while its members show. An end is placed by core's `placesOf` (every block that
  shows the position, each place counted). Layout makes stems, ribbons and hits; hover and a window's move only draw
  from them.
- **Order in a neck**: ribbons by their partner's place clockwise from the stem (the nearest clockwise last), then stubs
  by partner order, then the fan.
- **The twist**: a slice's edges in display order; `s` = strand × the two blocks' directions; `s` +1 (or a strand of 0)
  joins A's low edge to B's high (the flat ribbon: two facing arcs have opposite tangents), −1 low to low. Held by a
  property test over every combination (the sides cross exactly when `s` is −1).
- **The curve**: each side a cubic, its controls at `(1 − bundle)` of the neck's radius from the centre (`bundle` 0.7);
  on a row an arc down into the band by its span. Flattened into `rondo.flatten` (12) segments: the polygon drawn is the
  polygon hit.
- **Read**: core's pairs reader with the gate `rondo.maxRibbons`: bins until the stem pairs with an end in view sum to
  at most it, then items (held pairs past it binned by the same level rule). Pending (nothing landed, the view not
  settled or the reader `waiting` on a tile it wanted): drawn empty and said ("links loading", EC-169); settled with no
  rows and nothing awaited: "no links"; a tile failed (the reader's `failed`: its source's opening failed, the store
  says) and nothing read: `failed`, no hatch waiting, said "links failed" beside the store's issue. Why the reader's
  word too: the view settles for a frame between a source's open and its tiles' wants, and "no links" there is a
  flicker.

## Stubs

- **An end with no place** (its system not in the section) makes its pairs a stub, one a stem a partner: a taper from
  the neck toward the centre `rondo.stubShare` (0.2) of the hole, in the partner's hue, `k ×` its count wide (a fan
  member under a device px), labelled at its tip ("chr7 · 37", centred), the strongest first, a label only where its box
  (padded 2 px) meets no label taken and its tip is `rondo.stubLabelPx` (24) from every labelled tip. Why: tips crowd at
  a hub, and two labels over each other read as neither. `stubs: none`: not drawn (their counts in the fan).
- **Small stubs fold**: a stem's stubs of fewer than `rondo.stubMin` (2) links, when it has two or more, are one stub of
  partner −1 (`foldStubs`; its count their sum, so the neck is still the sum of its parts), in the muted text's colour,
  unlabelled (every stem has one: the text would repeat round the ring; its count in its card), keyed
  `stub:<system>:<bin>|other`; its card lists them ("Partners: chr3 (1), chr9 (1)", `foldedOf`), its double-click with
  `opens` the stem's chromosome alone. A lone small stub stays. Why: on one chromosome a stem's single links to every
  other chromosome are a comb of hairlines that hides the strong partners. Hovered it unfolds (`unfold`, the hover
  pass): its base split among its partners by count, each a taper in its partner's hue to its own tip, the tips at least
  7 px apart; nothing laid out again.

## Colour

- **The layer's `color` channel**, compiled as other layer types compile theirs (core's kept channels): absent, a field
  of `system` (what `rondoRing` writes) or a scheme over `system` (`rondo:hues`, its field `system` when it names none),
  a hue a block, never compiled: a stem its block's, a ribbon its `a` stem's, a stub its partner's (`rondo:hues`, else
  the categorical tokens, by the systems' order). **A fixed colour**: every polygon in it (stems, ribbons, stubs, fans,
  a folded stub), opaque (`opacity` its alpha, never `ribbons.fillOpacity`): a kind drawn over the hues reads as one
  colour, never a tint of each hue under it (a look at the showcase, light and dark: at 0.35 the accent read brown over
  one chromosome and grey-green over another). **Another field** (`type`, any form over the links' fields) is read link
  by link, and kinds are layers, never mixed inside one ribbon: where every row read gives the field (items; held bins,
  which are split by `type`) a stem pair's rows split into a ribbon a colour (its key the colour's hex after the strand;
  past 64 colours in one pair the last share a ribbon and `rondo.ribbons.colors` says so); where the rows carry none
  (bins by another field, a source's bins, which carry no fields) the layer has one colour, the ranges' ink
  (`seq-range`), `rondo.ribbons.color` (info) says so and the legend names it
  (`<label> (bins in one colour: <field> read link by link)`). The issue follows the level read: at the item level every
  link has its colour, so it is said only while such bins are drawn. Stems, stubs and fans take the layer's one colour
  (the ink when its links' values differ). The scan of a field's colours is kept by the read and the channel: a warm
  layout scans none. Why: a bin's count is every link's in it; a colour a value needs a bin a value, which held bins
  have for `type` and nothing else.
- **Hues repeat on a subset**: the default palette is the theme's `categories` (8 tokens, safe in both themes), taken by
  the systems' order (the assembly's, before any `include`), so chr7 has one hue on the genome's ring, on its own ring
  and on a stub pointing at it. On a subset two chromosomes eight apart (chr4 and chr12) draw in one hue: stated, not
  fixed. Why: re-deriving the domain from the blocks shown would give one chromosome two hues across two sections, which
  is worse; the pair's ribbons are labelled and its stubs named at their tips. `rondo:hues` has no repeat but is not
  colour-blind safe, so the layer never defaults to it: the writer names it.
- **`opacity`** (its channel's uniform value) scales `ribbons.fillOpacity`: a faint layer under a strong one.
- **The legend** (`LayerType.legend`, its code in the card's chunk, `card.ts`, loaded with `rondoViewer`: a viewer's
  legends derive once it lands, the settle waiting for it; the headless entry keys none): a field colour its channel's
  values (core's `channelLegends`); a fixed colour one swatch named by the layer's `label` (else its track's); a
  labelled layer in the hues "`<label>`: by chromosome", its label dropped where the legend's title already says it ("by
  chromosome" under "Contacts"), with the first chromosome's hue as its swatch (look 24: a row with no swatch read as
  broken, and its label doubled); an unlabelled one in the hues none (the hue ring keys the chromosomes).

## Kinds of links

- **A link's geometry says nothing of its meaning** (a 3D contact, a structural variant or translocation, a fusion, a
  synteny block): the kind is in the data (`type`), so the view shows it through layers, never inside one ribbon. **A
  kind is a `rondo:ribbons` layer** over the same pairs node, its data the kind's rows (`<links>/(item :type <kind>)`, a
  base: its rows left out are not "filtered"), in a track of its own (its label the kind's name: its legend's title and
  the card's word for it). Layers draw in their tracks' order, the last on top: a writer lists the rarer kind last. Each
  reads, stems, widths and bounds by itself: its own `k` (its widest stem's), its own fans, its necks the sum of its own
  ribbons, stubs and fan; `layer.hidden` hides a kind. Every layer's mouths are on one rim, the innermost band's inner
  edge (the hole's), whichever track holds it.
- **Widths by a layer's own rule and on its own scale** (the author, 2026-11-23: "per layer, said"): its own `k`, stated
  in its legend entry (`1 px = 1.3 links`, or `n q` for a summed field); `scale: section` (a plain param, beside
  `width`, which stays a word: the measure; `rondoRing`'s `scale`) puts the section's layers of `section` and one
  `width` rule on one `k`, their widest stem's (a layer laid out later read as last laid out: a frame late when its
  widest stem grows), so widths compare across kinds; necks stay sums of their parts. The measure: `width: count` (the
  links), or a numeric field summed (`q`): an item its value, a held bin its field's sum (core's
  `PairsReadOptions.sum`); links that carry none (a source's bins) are counted and `rondo.ribbons.width` (info) says so.
  Over links read by region (a source's, its `type` column mapped), core's reads filter their items by each layer's
  base; their bins carry no fields, so a filter is not read on them: the first layer over the links draws the bins, the
  next draws nothing (`empty: shared`, "links drawn by another layer"), and `rondo.ribbons.filter` (warning) says so on
  each. Never once a layer: kinds that cannot be told apart are drawn once. How bins read by region should split by kind
  is open.
- **The card names the kind** (`describe`): a ribbon at the item level lists the fields its links share (one link: its
  own), its `type`, its ends' strands (`Strands: + −`; none for unstranded ends, `.`: a contact's), then every other
  field its data carries (`q`, gene names where the file has them): a shared value, a number's range where they differ
  (`q: 0.001–0.04`), absent values and differing text left out; a ribbon wide by a field says its sum (`Width: q 0.07`).
  A bin's ribbon lists, when the section has more than one `rondo:ribbons` layer, each layer's links between its two
  stems (`Contacts: 140 links`, `Translocations: 3 links`); a stem each layer's links on its chromosome, under a heading
  `chr22, all links` (a line with no value); read from each layer's own stems and ribbons as laid out (their links
  counted at layout, those meeting the ranges, a ribbon with both ends in the range once): bounded by what the layers
  drew, no row read again on a hover, so a kind under another is never lost under its ribbons; a hidden layer is no peer
  (laid out no more). The layers of a section know each other by its resolved section (a new resolve, a new list). These
  lines are a chunk (`card.ts`) `rondoViewer` loads once composed, as the stem menu's: the headless entry pays nothing
  for them (no card there); a card shown before they land is made again when they do (`cardLanded`: the hover-card
  behavior's `framed`).
- **`rondoRing`'s `kinds`** (`{ type, label?, color?, opacity?, width? }[]`) writes those layers, the colour the
  chromosomes' hues when a kind names none; absent, one layer of every link. Which kinds, their order and colours are
  the application's (the showcase: faint contacts in the hues, translocations in one accent over them).

## Hits, emphasis, words

- **Hits**: a uniform grid over the polygons (`rondo.hitCellPx`, 16 CSS px), candidates topmost first (the last drawn),
  even-odd inside: a ribbon its key `ribbon:<a>:<bin>|<b>:<bin>|<strand>`, a stem or its fan `stem:<system>:<bin>`, a
  stub `stub:<system>:<bin>|<partner>`; each hit made at layout, so a hit allocates nothing. A thin ribbon or stub
  (under `2 × rondo.hitSlopPx` wide) covers `rondo.hitSlopPx` (4 CSS px) either side of its centre line: drawn over the
  wide ones, it is hit along it there; where thin covers overlap, the nearest centre line; the first wide cover below
  ends the search. Why: a one-link ribbon is a hairline a pointer cannot stay on.
- **Hover** (its own): a ribbon or stub at full opacity outlined in `seq-hover`, its stems lit; a stem its ribbons lit
  over `seq-wash` at `rondo.fade` (0.6). Another layer's hover on the section, the same wash over what it touches: an
  item (a cytoband) the ribbons, stubs and stems with a stem in its extent; a place with a position (the hue ring's
  ground, the bare axis) its position's stem; a block's arc or label (no position) every ribbon on the block. Why: on
  one chromosome's ring the block is everything, so lighting the block for a band lit all. The emphasis draw keeps no
  lay past itself (module state cleared after each pass). **Selection**: ribbons with an end in it in the accent.
  **Highlights**: ribbons touching each channel's Location.
- **The stem menu** (`menu.ts`, `StemMenu`, the descriptor's `headerActions`, so `viewer.extend` brings it too): a
  secondary click on a stem opens its element menu (a secondary click on a ribbon opens the track's menu: the ring has
  no header in view, so a ribbon is the pointer's route to the track's actions), the viewer's label menu's shape scoped
  to `rondo:ribbons`' stems (`element.type`, a `stem:` key): "Show with a partner ▸" (its partners by chromosome,
  `partnerChromosomes`: a row a partner chromosome with its links in all, most first, `chr9 · 449 links`; a chromosome
  reached by more than one ribbon its pieces in a submenu, each named by its band range and strand,
  `q31.1–q34.2 · 379 links`, `p13.3 · 2 links, reversed`: the cytobands the section draws there, else its span in bases;
  the stem's own chromosome never a partner but last, `Within chr22 · 229 links`, its pieces the ribbons to its other
  stems and the stem's links within itself, `this stem`; each row and piece shows its chromosome with the stem's as a
  pair in the section the layer `opens`, else its own, within alone. Why by chromosome: bare rows `chr9 · 375 links`,
  `chr9 · 65 links` said nothing of what parted them; a level lists at most 8 rows, past them `n more` with the rest in
  the same order, `capped`, so the pointer reaches every partner the keys do), "Select the stem and its partners"
  (`rondo:select-partners` with the element: one selection), "Show only chr4"; only these: the entries are `alone`, so
  the track's actions (Move up, Reverse, Hide track) stay in its header's menu. The entries find the stem in the lays of
  the viewer the menu is in (`HeaderFacts.viewer`: no state of their own, two viewers of one spec each their own). **Its
  route**: a secondary click on a stem opens it at the pointer; on a ring's empty space (no hit) the browser's menu
  opens (the viewer's rule for rings whose hole a section layer draws); the menu key and Shift+F10 with the keyboard
  cursor on a placed stem's range (the topmost layer's) open it at the stem's mouth (`atCursor` on its first entry,
  `stemAtCursor`). Why: a menu a reader cannot reach is none.
- **The keys**: `.` (`rondo:partner`) and `,` (`rondo:partner-back`), the commands' shortcuts on the canvas: the cursor
  to the next or the previous partner stem of its stem in the stem menu's order (its rows' pieces in turn: one list, two
  readers), said by the announcer (the command moved the cursor: "chr17, 40,000,000 · 140 links, reversed").
- **Words** (core's `LayerType` hooks): `describe` (a ribbon "chr4 ↔ chr7", then both ends, a row each labelled by its
  chromosome, "q13.2–q22.1 · 67.11–100.66 Mb": the cytobands the section draws over its stem, `CardLines.band`, read
  from the viewer whose frame holds the lay, its span alone where none or headless; its count and kind; a bin's "16.8 Mb
  bins: zoom in to see them"; a stem "This stem: 24 links to chr9", its one partner chromosome by name, else "to 23
  chromosomes", then within and fan; a stub its partner, "not shown"), `summary` ("612 ribbons, 23 chromosomes;
  strongest chr1 ↔ chr17, 140"), `announce` ("stem of 312 links").
- **The hole's words**: on a ring whose links are pending the hole says "Loading links…", muted (`msui-text-muted`) and
  centred (the labels pass, measured at layout: the SVG too); on a failure the error its links' source raised replaces
  it (the store's issue at the source's address, `HoleIssue`, of the viewer whose resolved view holds the section laid
  out: two rings of one source address each their own; an export, which resolves its own view, by the section's address;
  headless "Links failed"). Each is said once, politely (`viewer.notify` at a rest, a notice a hover never replaces,
  again only after the hole had other words or none). Why: the pending hole was empty and silent (novice 35). The
  square's pending state is the square's: the pending hatch over its lattice with "Loading…", on a failure the same
  issue's text (`HoleIssue.of` with no section), via the cells' `note` (novice 35).

## What is the core's

- **The zoom is a cut by `section.include`** (core's section param and `patch` command, one undo step): "show only chr4"
  patches `include: ['chr4']`, a pair `['chr4', 'chr7']`, "show all" removes it. The ring is laid out again with no
  tween; the hue domain and the viewports by `{ system, position }` are untouched. Why a param and not new `block`
  nodes: one patch to apply and to undo, the assembly's grammar, and rows get the same switch. In place by default;
  `detail: 'below'` of `rondoRing` writes a second ring section under it, so the genome stays in view.
- **The hole**: ribbons own it, so `rondoRing` writes the section's `center: none` (core's section param, the assembly's
  name and length then in the section's label line).
- **Cytobands and the hue ring are core's `bands` layers**, the hue ring a thin second one over the same data coloured
  by `(field :name system)`; the extension has no layer type for them. Why: the bands cover each chromosome, and a
  reader without the extension draws them too.
- **The square's window, its click that shows a pair, `navigator.switch: include` and the lattice over the assembly's
  blocks are core's navigator**; the extension draws the cells only.
- **Words**: the card, the ARIA description and the keyboard cursor are the viewer's; a layer type supplies its lines
  through the optional `describe`, `summary` and `announce`, so built-in layers need nothing.
- **On a row** (a `tracks` section of blocks) the same stems sit below the band and the ribbons are filled arcs between
  them, the polar warp's absence the only difference. A text or a grid layout is not drawn: `rondo.ribbons.layout`
  (warning, once). Why: its fallback is for readers without the type, not for layouts.

## The fallback

| Part                                   | Drawn as                                                                                                                  |
| -------------------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| `rondo:ribbons` (`fallback: pairs`)    | core's chords: a bin pair one chord stroked by count, items as chords under `pairs.maxArcs`, an end with no place a caret |
| `rondo:hues`                           | the layer type's default colour (the categorical palette over `system`)                                                   |
| `rondo:square`                         | core's `axis` picture: the lattice and its window, its gestures the same                                                  |
| the commands, the `opens` double-click | absent; the label menu's "Show only" (core) remains                                                                       |

One `seq.view.layer.extension` issue a view ("1 layer(s) need `rondo` 0.1: drawn as `pairs`"); `rondoRing` writes the
colour channel explicitly (`(field :name system)` and its palette), so the fallback's colours are the spec's.

## The builder helper

`rondoRing(view, { assembly, links, bands?, detail?, square?, hues?, opens?, height?, kinds? })`: the ring section
(`center: none`), the cytobands and the hue ring (two `bands` layers over `bands`), the ribbons with `fallback: pairs`
and their colour (a layer a kind with `kinds`: Kinds of links), `opens` (the detail ring with `detail: 'below'`, else
the ring), and the square navigator over the ring (`of: [ring, ring]`, `content: rondo:square`, `switch: include`). The
spec's `extensions` is the writer's (`createSeqSpec({ extensions: [rondo] })`).

## Configuration

| Item                                   | Whose            | Default                              | Why                                                                                                                                                                                                                  |
| -------------------------------------- | ---------------- | ------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `width`                                | layer param      | `count`                              | the pairs counted; a numeric field summed (`q`)                                                                                                                                                                      |
| `label`                                | layer param      | absent                               | a kind's name in its card and legend                                                                                                                                                                                 |
| `scale`                                | layer param      | `layer`                              | each layer's widths on its own `k` (the author, 2026-11-23: "per layer, said"): a rare kind readable; `section`: one `k` over the section's layers of `section` and the same `width`, so widths compare across kinds |
| the legend's scale                     | fixed            | —                                    | each layer's entry states its scale, `1 px = n links` (or `n <field>`), `, shared` under `section`                                                                                                                   |
| `stems`                                | layer param      | `level`                              | Rondo's trumpets; `none`                                                                                                                                                                                             |
| `bundle`                               | layer param      | 0.7                                  | ribbons of a stem stay together into the hole                                                                                                                                                                        |
| `stubs`                                | layer param      | `partner`                            | where else a stem's links go; `none`                                                                                                                                                                                 |
| `opens`                                | layer param      | absent                               | the section a double-click sets `include` on                                                                                                                                                                         |
| the colour                             | layer theme slot | `(field :name system)`, `categories` | a hue per chromosome, safe; `rondo:hues` named by the writer; a fixed colour or another field (Colour)                                                                                                               |
| `opacity`                              | layer theme slot | 1                                    | a factor of `ribbons.fillOpacity`: a faint kind under a strong one                                                                                                                                                   |
| `rondo.stemPx`                         | setting          | 16 CSS px                            | a stem a reader can point at                                                                                                                                                                                         |
| `rondo.stemDepth`                      | setting          | 0.12 of the hole                     | room for the flare                                                                                                                                                                                                   |
| `rondo.maxWidthPx`                     | setting          | 24 CSS px                            | one px a count for the view                                                                                                                                                                                          |
| `rondo.maxRibbons`                     | setting          | 1,024                                | the items gate and the drawn bound in one                                                                                                                                                                            |
| `rondo.stubShare`, `rondo.stubLabelPx` | setting          | 0.2; 24 CSS px                       | a stub short of the centre; a label where it reads                                                                                                                                                                   |
| `rondo.hitCellPx`                      | setting          | 16 CSS px                            | few candidates a cell                                                                                                                                                                                                |
| `rondo.hitSlopPx`                      | setting          | 4 CSS px                             | a hairline ribbon or stub a pointer can stay on; 3 to 4 px either side                                                                                                                                               |
| `rondo.fade`                           | setting          | 0.6                                  | the hovered stem's ribbons read over the rest                                                                                                                                                                        |
| `rondo.flatten`                        | setting          | 12                                   | hit within a device px at 720 px                                                                                                                                                                                     |
| `rondo.focusStems`                     | setting          | 0.35                                 | the ribbons read over the stems on one chromosome or a pair                                                                                                                                                          |
| `rondo.stubMin`                        | setting          | 2                                    | single links fold into one "other chromosomes" stub a stem                                                                                                                                                           |

**Settings are the descriptor's and the composition's** (`rondoWith(config)`, `rondo` the defaults; `RondoConfigItems`
name them, listed as the descriptor's `config`): a layout reads the composition's `rondo.*` pairs (`settingOf`, from the
render config it is given) over the descriptor's numbers (`composedConfig`, one object a pair: a warm layout makes
none), so `rondoWith` is a convenience, not the only door. A viewer extended late with `rondoViewer` (which lists no
settings: they are never registered late) reads the descriptor's alone.
