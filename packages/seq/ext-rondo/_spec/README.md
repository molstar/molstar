# `@molstar/seqx-ext-rondo`: the genome ring

A genome ring like Rondo's, as an extension on the public surface: ribbons between stems wide by count, stubs to blocks
not shown, 24 hues. Two entries: `.` (`rondo`, headless: the layer type, the scheme, params, `rondoRing`) and `./viewer`
(`rondoViewer`: `rondo` and its commands and binding). Why two: the CLI and `createSeqContext` compose the first and pay
nothing for the viewer's parts. It imports `@molstar/seqx-core/extension` and `@molstar/seqx-viewer`'s root only: a gap
is the contract's, brought to the core, never a reach into `src`.

## What it ships

| Entry                   | Kind                                                                                                 | Fallback                 | Reads                                                       |
| ----------------------- | ---------------------------------------------------------------------------------------------------- | ------------------------ | ----------------------------------------------------------- |
| `rondo:ribbons`         | layer type, `scope: 'section'`, shape `pairs`, on `rings` and a `tracks` row                         | `pairs`                  | core's `PairsReader`: bins at a level, items under its gate |
| `rondo:hues`            | scheme over categories: n hues round the Oklch wheel, not colour-blind safe, substitute `categories` | the layer type's default | —                                                           |
| `rondo:partner`         | command: the cursor to the next partner stem by count, said                                          | —                        | the stem under the cursor                                   |
| `rondo:select-partners` | command: the far ends of the selection's ribbons added                                               | —                        | the selection                                               |
| `rondo:open`            | command, bound to a double-click scoped to `rondo:ribbons`' hits                                     | —                        | the hit's chromosomes, set as `opens`' `include`            |

Not shipped: `rondo:cells` and the `rondo:square` picture. Why: core's square navigator draws its lattice and never a
picture, so its cells would be a private rebuild; `rondoRing` writes core's square (`switch: include`) and the design's
rule without pictures holds. The stem menu waits for element menus in the viewer.

## Stems

- **A stem is one bin of a coarser level**: `2^j` × the read's level (bins), or `2^j` (items), the largest no wider than
  `rondo.stemPx` (16 CSS px) of positions at the ring's outer circumference (a row: the strip's own px), the finest
  strip's; never device px, so the stems are the same at any DPR and in the SVG. Aligned on position 0 of each system: a
  stem is a union of whole bins, its counts exact sums. `stems: none`: the read's own level, no trumpets.
- **Its count** is the pairs leaving it (ribbons, stubs, fan); pairs with both ends in it are counted apart ("within")
  and drawn once a finer level parts them.
- **The trumpet**: its mouth on the layer's band's inner edge (the hole's rim) to its neck `rondo.stemDepth` (0.12) of
  the hole inward, two cubics flaring, in its system's hue at `ribbons.fillOpacity`.

## Ribbons

- **A ribbon is a stem pair a strand** (the rows' relative strand), its `a` the stem earlier in the systems' order (its
  hue). Width: `k` = `rondo.maxWidthPx` (24) / the largest stem count in view, a ribbon `k ×` its count; under
  `ribbons.minPx` device px, or past the strongest `rondo.maxRibbons` (1,024), it goes to both stems' fans. So a neck is
  the sum of its ribbons, stubs and fan; at least a device px as a whole.
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
  rows and nothing awaited: "no links". Why the reader's word too: the view settles for a frame between a source's open
  and its tiles' wants, and "no links" there is a flicker.

## Stubs

- **An end with no place** (its system not in the section) makes its pairs a stub, one a stem a partner: a taper from
  the neck toward the centre `rondo.stubShare` (0.2) of the hole, in the partner's hue, `k ×` its count wide (a fan
  member under a device px), labelled at its tip ("chr7 · 37", centred), the strongest first, a label only where its box
  (padded 2 px) meets no label taken and its tip is `rondo.stubLabelPx` (24) from every labelled tip. Why: tips crowd at
  a hub, and two labels over each other read as neither. `stubs: none`: not drawn (their counts in the fan).

## Hits, emphasis, words

- **Hits**: a uniform grid over the polygons (`rondo.hitCellPx`, 16 CSS px), candidates topmost first (the last drawn),
  even-odd inside: a ribbon its key `ribbon:<a>:<bin>|<b>:<bin>|<strand>`, a stem or its fan `stem:<system>:<bin>`, a
  stub `stub:<system>:<bin>|<partner>`; each hit made at layout, so a hit allocates nothing.
- **Hover** (its own): a ribbon or stub at full opacity outlined in `seq-hover`, its stems lit; a stem its ribbons lit
  over `seq-wash` at `rondo.fade` (0.6). Another layer's hover on the section, the same wash over what it touches: an
  item (a cytoband) the ribbons, stubs and stems with a stem in its extent; a place with a position (the hue ring's
  ground, the bare axis) its position's stem; a block's arc or label (no position) every ribbon on the block. Why: on
  one chromosome's ring the block is everything, so lighting the block for a band lit all. The emphasis draw keeps no
  lay past itself (module state cleared after each pass). **Selection**: ribbons with an end in it in the accent.
  **Highlights**: ribbons touching each channel's Location.
- **Words** (core's `LayerType` hooks): `describe` (a ribbon "chr1:12–28 Mb ↔ chr17:40–52 Mb", its count and kind, its
  bins; a stem its links, partners, within and fan; a stub its partner, "not shown"), `summary` ("612 ribbons, 23
  chromosomes; strongest chr1 ↔ chr17, 140"), `announce` ("stem of 312 links").

## The fallback

| Part                                   | Drawn as                                                                                                                  |
| -------------------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| `rondo:ribbons` (`fallback: pairs`)    | core's chords: a bin pair one chord stroked by count, items as chords under `pairs.maxArcs`, an end with no place a caret |
| `rondo:hues`                           | the layer type's default colour (the categorical palette over `system`)                                                   |
| the commands, the `opens` double-click | absent; the label menu's "Show only" (core) remains                                                                       |

One `seq.view.layer.extension` issue a view ("1 layer(s) need `rondo` 0.1: drawn as `pairs`"); `rondoRing` writes the
colour channel explicitly (`(field :name system)` and its palette), so the fallback's colours are the spec's.

## The builder helper

`rondoRing(view, { assembly, links, bands?, detail?, square?, hues?, opens?, height? })`: the ring section
(`center: none`), the cytobands and the hue ring (two `bands` layers over `bands`), the ribbons with `fallback: pairs`
and their colour, `opens` (the detail ring with `detail: 'below'`, else the ring), and the square navigator over the
ring (`of: [ring, ring]`, `switch: include`). The spec's `extensions` is the writer's
(`createSeqSpec({ extensions: [rondo] })`).

## Configuration

| Item                                   | Whose            | Default                              | Why                                                          |
| -------------------------------------- | ---------------- | ------------------------------------ | ------------------------------------------------------------ |
| `width`                                | layer param      | `count`                              | the pairs counted                                            |
| `stems`                                | layer param      | `level`                              | Rondo's trumpets; `none`                                     |
| `bundle`                               | layer param      | 0.7                                  | ribbons of a stem stay together into the hole                |
| `stubs`                                | layer param      | `partner`                            | where else a stem's links go; `none`                         |
| `opens`                                | layer param      | absent                               | the section a double-click sets `include` on                 |
| the colour                             | layer theme slot | `(field :name system)`, `categories` | a hue per chromosome, safe; `rondo:hues` named by the writer |
| `rondo.stemPx`                         | setting          | 16 CSS px                            | a stem a reader can point at                                 |
| `rondo.stemDepth`                      | setting          | 0.12 of the hole                     | room for the flare                                           |
| `rondo.maxWidthPx`                     | setting          | 24 CSS px                            | one px a count for the view                                  |
| `rondo.maxRibbons`                     | setting          | 1,024                                | the items gate and the drawn bound in one                    |
| `rondo.stubShare`, `rondo.stubLabelPx` | setting          | 0.2; 24 CSS px                       | a stub short of the centre; a label where it reads           |
| `rondo.hitCellPx`                      | setting          | 16 CSS px                            | few candidates a cell                                        |
| `rondo.fade`                           | setting          | 0.6                                  | the hovered stem's ribbons read over the rest                |
| `rondo.flatten`                        | setting          | 12                                   | hit within a device px at 720 px                             |

**Settings are the descriptor's** (`rondoWith(config)`; `rondo` the defaults; `RondoConfigItems` name them). Why not the
spec's `config`: a layer type reads nothing of its context's settings, and a descriptor listing settings is composed,
never registered late (`viewer.extend`).
