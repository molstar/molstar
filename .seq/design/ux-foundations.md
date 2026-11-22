# UX foundations

## Intent

The author, 2026-11-17: "i would like to have better visual fidelity and nicer UI in general. so we should write "ux
foundations" at start too, similar to the ones the packages/ui has .. but target to showing the sequence data, when to
squircles, etc. … basically i want this to also look really nice. and the default behavior should be good too ..
consistent logic for selections, toggling, etc."

`@molstar/ui`'s foundations (its `styles/_spec`: color, type, spacing, sizes, alignment, shape, motion) govern the
chrome: headers, menus, legends, panels. These foundations extend them to the data: letters, cells, marks, tracks, and
what the pointer and keyboard do to them. The chrome uses `--msui-*` as it is; the data's tokens (`--seq-*`) are
derived from them, so a theme switches both.

## Provisional

This is the first statement. The numbers (sizes, steps, durations) are settled by the foundation stories, which draw
each rule on real data beside its violation, as `@molstar/ui`'s do; the grammar in section 8 is settled by its
conformance tests.

## 0. Configurable, with rules that hold

The author, 2026-11-17: "for "ux foundations" and "Extensibility": should "layouting/spacing/etc" be
extensible/configurable?" Yes, by the rule that what can be configurable is ([extensibility.md](extensibility.md)):

- **Every number here is a setting**: the spacing unit and scale, row heights, gaps, the label column, insets, letter
  sizes, density thresholds, radii, motion durations, and every color through the data tokens (`--seq-*`, derived from
  `--msui-*`). They come in presets (compact, comfortable, figure for export, print) that a viewer configuration or a
  section's style picks and overrides field by field.
- **Layouts share one internal interface**: tracks, text, rings and dot plots take a section of the resolved view, the
  settings and a size, and give geometry, hit-testing and draw lists; the interface opens as an extension point when an
  application needs one ([extensibility.md](extensibility.md)).
- **The rules hold whatever the values**, checked by one set of tests over every layout and preset: positions on the configured unit's grid; one label column a view; columns
  aligned to the device pixel across sections that share an axis; density steps that never move content; the emphasis
  order; square continuous data and squircle discrete marks; the grammar's fixed rules (section 8).

## 1. Color

- **Data and chrome never share a color role.** The accent (`--msui-accent`) is for actions and what is picked:
  selection and the focus cursor. Data colors (palettes, residue schemes, values) never use it, and selection never
  uses a data color. Why: a picked residue must read as picked whatever scheme colors it.
- **Palettes are tokens per theme.** Every categorical palette, ramp and residue scheme has light and dark variants,
  checked for distinctness (neighbors in a palette stay apart for the common color-vision deficiencies) and for
  contrast against the track background.
- **Ink follows its fill.** A letter on a colored cell takes dark or light ink by the cell's luminance (the canvas-ink
  rule), at 4.5:1 or better. Why: a scheme that colors cells must not cost the letters.
- **Missing is not zero.** Missing, unobserved and unmapped each have one neutral treatment (a gray family derived from
  the theme), distinct from the lowest value of any ramp.
- **Conflict is one color**, the theme's warning hue, for a mismatch, a conflicting map element, a failed load.

## 2. Type

- **Letters** are set in one monospaced data face with tabular metrics, centered in their cells on the cap height, in
  sizes on a short scale tied to the cell (section 3), **capped at the chrome's label size** (`letters.maxPx`, the
  label font's size by default): optically matched, so a letter's cap height equals a label's (about 0.9× the label
  size for the mono face). Past the cap a wider cell adds room around the letter, never a bigger letter; the cell's
  color, emphasis and hit area still fill it. Why: letters larger than their labels out-shout the chrome (the author,
  2026-11-21: "using bigger font for the sequence letters than the surrounding UI looks bad"); a host that wants big
  letters (teaching, slides) sets `letters.maxPx` to follow the cell. Measured, not a ratio: the mono face's and the
  label face's cap heights, once per configuration. Not capped: a logo's letters (their height is data) and a text
  section's body (a text figure is read as text, at its own size). Two ink weights only: regular by default, and
  strong for a reference row, the cursor or a conflict.
- **Labels and numbers** in the chrome use `@molstar/ui`'s type scale; rulers and numbering use tabular figures and
  never overlap (a label that does not fit is skipped, not shrunk).
- **One letter, one size per row**: a row shows letters for all its positions or none (EC-116).

## 3. Grid

- **The cell is the unit of data.** A position is a cell whose width is set by the zoom; letters show when the cell is at
  least the smallest legible size, then cells only.
- **Rows are on the 4px grid**: a letters row 20px by default (16 compact, 24 comfortable), feature rows 8 to 16px,
  value tracks in multiples of 4. Track gaps 4px within a group, 8px between groups.
- **A track's padding is vertical** (the author, 2026-11-19): `padding` inside its well, top and bottom (a number is
  both; `(top 8)` one side), default 0 and a render config default for hosts; the well grows by it, the data sits
  inside, hover and hit tests read the data's band. Never left or right: the horizontal is the axis, and columns stay
  aligned across tracks.
- **One label column**: track and row labels sit in a column on `@molstar/ui`'s inset and label width, so the data
  starts at one x for every track of a view, in both layouts (tracks and text).
- **Columns stay aligned across sections** that share an axis (scroll groups), to the device pixel.

## 4. Density

What a cell shows steps down as it narrows, at fixed thresholds, without moving anything:

1. letters (with their colored cells);
2. cells only, colored;
3. dots or ticks per position (sequence present, letter hidden);
4. bars per run (a feature, a gap run);
5. density (a heat strip: coverage, conservation, variant counts).

Each step keeps the previous step's positions and colors, so zooming reads as resolution changing, never as the
picture changing. The thresholds are shared by every layer type. Each step's unit is a config item (`densityUnits`); by default letters are
measured in CSS px (10), the other steps in device px (cells 3, dots 1, bars 1/64) (default; to confirm with the author).
Why: a letter's legibility is in CSS pixels (a size that needs 14 device px at DPR 1 needs only 7 at DPR 2), while what
bounds the draw lists is a mark a device pixel. The letters step never starts below the width the
smallest letter size needs, so cells never show without letters.

## 5. Shape: when to squircle

- **Continuous data is square.** Letter cells, alignment grids, matrices, dot plots, value bars: adjacency carries
  meaning (a run reads as one run), so no corner is rounded.
- **Discrete things are squircles.** A feature, a domain, a ligand chip, a bookmark, a label pill, a selection's
  outline, a hover card: rounded with the superellipse `@molstar/ui` uses for surfaces (`--msui-corner-shape`), its
  radius proportional to its height (a third), capped by the surface radius.
- **Small things are square.** Below 6px in either dimension a discrete mark drops its rounding: a corner that can't be
  seen is noise, and adjacent items must not leave gaps.
- **Cut ends are square.** An item cut by the view's edge, an intron collapse or a block boundary has square ends on the
  cut side and rounded ends where it really ends (the frame rule, EC-100).
- **Outlines are concentric.** A selection or hover outline around a squircle mark takes the mark's radius plus the
  outline's offset, so the gap is even.
- **A track's well is a surface, so a squircle** (the author, 2026-11-19): its radius by the rule above (a third of its
  height, capped by the surface radius), its data clipped to it (a letters or grid track loses a sliver at its corners,
  as a rounded card over a table); a well shorter than 6px stays square, and a config switch (`trackShape: squircle |
  square`) keeps the square look for a host. The data inside follows the rules above: continuous data stays square,
  only the well's corners round.
- **Every emphasis is a squircle** (the author, 2026-11-19, "Yes, round all"): the hover, the hover projected onto
  linked tracks, the selection, the highlight channels (each style) and a chip's or a cell's hover outline take the
  rule above: round where the span really ends, square where the view, a piece or a block cuts it, square below 6px,
  an outline concentric with what it surrounds. Why one rule: the selection already did, and emphasis that differs by
  what it lies on reads as two systems.
- **One path, every backend.** Squircles are drawn as a path in the draw lists, never left to a backend's rounded
  rectangle, so Canvas2D and SVG draw the same curve.

## 6. Marks and emphasis

- **Layers of emphasis, in a fixed order** from below to above: data; highlight channels (named, each a tint of its
  color under the data); selection (an accent tint under the data and an accent outline over it); hover (an outline at
  3:1 or better against the track it sits on, a double outline over a fill it fails, EC-138); the keyboard cursor (an accent outline, thicker). A higher one never
  hides a lower one entirely. **Current pass order** (`selection` after `data`): the selection's tint is drawn
  translucent over the data, with its outline, not under it (default; to confirm with the author: the alternative is a
  separate under-data pass). Why over: a pass under the data is hidden by opaque cells.
- **The look as built** (the author's look at the selection, 2026-11-18, lane P): the selection is a soft tint (0.18)
  with one thin outline (1 px, 0.6), one band a section across its tracks and gaps (`emphasis.selectionShape`:
  `column`, or a box a track, `tracks`); the hover's outline is concentric, offset outward (1 px) and drawn above the
  selection; the keyboard cursor is a 1 px neutral line (`msui-text`), shown only after a key, not the thicker accent
  outline the list above names. A highlight channel is a tint (0.3) of its color; what is projected onto a linked
  section takes its own look at half strength (0.5).
- **Overlapping highlight channels stay apart** (the author, 2026-11-19, on the base's several colored marks: "can we
  have this functionality too?", then "ok as recommended"): where two or more channels light one position of a track,
  the band splits into equal horizontal stripes, one a channel present there, in channel order (the order set); up to
  `emphasis.highlightStripes` (4), the rest merged into the last stripe. Where one channel alone lights, it takes the
  whole band. A channel's `style`: `tint` (the default, as above), `band` (a bar along the track's foot, the base's
  band style), or `outline` (its color's outline, no fill). The legend part lists the active channels with their
  colors, each removable (`highlight.clear`, the themes design). Why stripes: stacked 0.3 tints of two colors blend
  into a third that names neither.
- **Hover is the same everywhere**: one outline style for a cell, a range, a row, a pair; projected hovers on linked
  sections take the same style, lighter.
- **Selection is one look** in tracks, text, rings and dot plots.
- **Boundaries** (an insertion point, an indel, a zero-length item) draw as a hairline with a pinned caret, never as a
  cell (EC-127).

- **One hover contract for every mark** (the author, 2026-11-21: "we should not be fixing this one-by-one ... we need a
  general rule for this and global fix"; built by HC): every drawn mark of every layer type, in every layout and at
  every level of detail, answers a hover with (1) a hit naming what is under the pointer (the item, else the bin's items
  or stretch, else the position), (2) an emphasis in the hovered view and every linked view, reaching every layer the
  covered positions reach (the view draws a position box where a layer draws none), and (3) a card with the hit's
  fields. A layer type without a hit is refused at registration (`seq.render.layer.hover`). A drag suspends the card. A
  rail's hover is its own (no card). A medium-lane sweep over every story holds it. Legends work both ways: an entry
  that names items lights them in view; an entry with no action shows the default cursor.

## 7. Motion

- **Navigation that jumps animates**: focus, reveal, a snapshot's view, zoom by a step (base duration, eased in and
  out); the eye follows the content to its new place.
- **Navigation that follows the hand does not**: drag to pan, wheel and pinch zoom track the input directly.
- **Data does not animate**: arriving regions draw in place (no fade on every tile); a density step is instant.
- **Reduced motion**: every animation is a cut. One viewer item decides it for jumps and stories alike: `view.motion`
  (`platform`, the default: the platform's setting; `cut`; `animate`).

## 8. Interaction grammar

One grammar for every layout, layer type and target. A gesture acts on **the element under the pointer at the
granularity of what was hit**: a cell (a residue, a column of an alignment), an item (a feature, a pair), a row (by its
header), a column (by the ruler), a track (by its header).

**The grammar is configurable** (the author, 2026-11-17: "the controls grammar should be configurable"). It is data: a
table of bindings from a gesture (with its modifiers and what it hits) to an action, and every action is a command
(`@molstar/ui/commands`). The table below is the default. An application composing the viewer replaces bindings, adds
its own, or picks a preset (for example one where a drag on data pans, as genome browsers do); a user's changes are a
settings layer over it. The rules that stay fixed whatever the bindings:

- **One configuration per viewer, applied to every view in it**: a gesture means the same thing in tracks, text, rings
  and dot plots of one viewer. Why: consistency is the point of a grammar; configurability must not bring back
  per-view behavior.
- **Granularity comes from the hit, not the binding**: a binding says "select", the hit says a codon or an item.
- **Modifiers are platform-aware** (Cmd on macOS, Ctrl elsewhere), never browser zoom keys or Tab (EC-155).
- **Every action stays reachable** by a command and from the keyboard, whatever the pointer bindings.
- **Conflicts are reported**: two bindings for one gesture is an issue when the configuration is built, not a race at
  run time.

| Gesture | Effect |
|---|---|
| hover | the element is hovered here and on every linked section; nothing else changes |
| click | the selection becomes the element; clicking the only selected element clears it |
| Shift+click, Shift+drag | extends the range the anchor is in to the element, keeping the other ranges (a range in one system; a run of rows) |
| Cmd/Ctrl+click | toggles the element in or out of the selection, keeping the rest; the anchor moves there |
| Cmd/Ctrl+drag | adds the range to the selection, keeping the rest |
| Alt/Option+drag | removes the range from the selection |
| click on empty space, Esc | clears the selection (Esc first cancels a gesture in progress) |
| drag on data | selects a range, snapped to the granularity hit (whole codons over nucleotides, whole items over items) |
| drag on a ruler (the ruler rail) or with Space held | pans ([navigators.md](navigators.md)); on an alignment's ruler only its lower 40% is the rail's, and a drag above it selects columns; a grid's rulers pan their axis |
| drag on a navigator | on its window pans, on a grip zooms with the other end held, from outside sweeps a new window and zooms to it; a click outside centres there; a wheel over it zooms about the window's centre |
| wheel, two-finger scroll | pans; Cmd/Ctrl+wheel or pinch zooms about the pointer |
| double-click | focuses the element (zooms to it); on the element already focused, goes back; on empty space, fits the view |
| Shift+double-click, -, =, 0 | zooms out about the pointer; zooms out and in about the cursor or the centre; fits |
| Alt+Left, Alt+Right, the mouse's back button | back and forward through the view's jumps (focus, fit, reveal, a bookmark or link): session state, not undo |
| arrows | move the keyboard cursor; Shift+arrows extend the anchor's range; Cmd/Ctrl+Shift+arrows start a new range; Enter focuses; Space toggles, on key-up and only when no drag happened while it was held (so Space+drag pans without toggling) |
| `]`, `[` | zoom an alignment's rows in and out about the cursor's row (`rows.zoom-in`, `rows.zoom-out`); `0` fits rows and columns |
| click on a tree node | folds or unfolds its clade (`clades.toggle`; a preset may bind it to select) |
| drag, wheel or pinch on a follower | moves its leader through the follow, so the follower moves as asked (`follow.followerMoves`: `lead`) |

- **Selection is a location**, held per system and projected: what is selected in a protein shows on its gene, its
  structure and its alignment rows.
- **Toggles are uniform**: hiding a track, collapsing a group, pinning a row, showing a layer: one gesture, one state
  shown in the chrome (an icon, never only a color), one undo step when it changes the spec.
- **View state is not history**: hover, selection, the cursor, the viewport and the folded clades are transient (a fold
  is session state, no undo step; `clades.keep` writes the folds shown into the spec as its opening set, one patch); edits to the spec (a
  bookmark, a hidden track, a theme) are patches, one undo step each, with a preview while a gesture is in progress that
  ends in its patch or a cancel.
- **Every gesture has a keyboard way and a command**, and the context menu lists what applies to the element under the
  pointer, in the same words as the commands.

### The interaction matrix

What a user acts on, against the verbs (the principles, "Interaction is a matrix", from the second slice's fold-back):
each cell is the gesture or key and the command, **planned** with the slice that builds it, or **no** with why. A slice
plan names the cells it fills; its gate checks them; a design that adds an action fills its row and column, and its
self-check's negation (an inverse, a symmetric case, a plural) names the cells to look at. As built after the second
slice, with the cells the third slice's designs fill ([axis.md](axis.md), [themes.md](themes.md),
[layouts.md](layouts.md), [grids-and-rows.md](grids-and-rows.md), [navigators.md](navigators.md); each design's own
table has the cells of its sub-rows: a text or ring section, a residue through a path, a break, a pair, a block):

| | select | extend | add / remove | focus | back | zoom in / out | pan | hide | edit |
|---|---|---|---|---|---|---|---|---|---|
| **an element** (a cell, an item) | click; `selection.set` | Shift+click (the anchor's range) | Cmd/Ctrl+click toggles | double-click; `address.focus` | double-click on it when focused | — | — | no: hide a track or layer | no: data is not edited in the viewer |
| **a range** (positions) | drag on data | Shift+drag | Cmd/Ctrl+drag adds, Alt+drag removes; `selection.add`, `selection.subtract` | Enter at the cursor | Alt+Left | — | — | — | — |
| **the selection** | — | Shift+arrows | Cmd/Ctrl+Shift+arrows start a range | Enter | Alt+Left | — | — | — | a bookmark from it (`selection.bookmark`, Cmd/Ctrl+D); copy its letters as shown (`selection.copyLetters`, Cmd/Ctrl+C); color it on a layer: planned (`theme.color-selection`) |
| **a leading section** | — | — | — | — | — | Cmd/Ctrl+wheel, pinch, -, =, Shift+double-click; in display order (a reversed axis) | wheel, drag on the ruler rail, Space+drag, Shift+Alt+arrows; in display order | no: hide its tracks | its params by patch (height); Shift+R `section.reverse`, Shift+C `section.collapse` (patches) |
| **a follower section** | as a leading one | as a leading one | as a leading one | as a leading one | as a leading one | moves its leader (`lead`) | moves its leader (`lead`) | as a leading one | as a leading one |
| **the view** | — | — | — | 0 fits; double-click on empty space | Alt+Left / Alt+Right, the mouse's buttons; `view.back`, `view.forward` | — | — | — | undo and redo through the host's stack |
| **a track or layer** | a row of the rows layer by its header (below) | — | — | — | — | — | — | `track.hide`, `layer.hidden`; a layer's items by a data filter (`layer.filter`), never a theme | height, color, order by patch; previews while dragged; its color and channels (`layer.color`, `layer.channel`, `scale.set`) |
| **a row** (an alignment's, by its header) | click; `selection.set` | Shift+click: the run in display order | Cmd/Ctrl+click toggles | double-click; `address.focus` | Alt+Left | `]`, `[` (`rows.zoom-in`, `rows.zoom-out`), Cmd/Ctrl+wheel over the header | the wheel over the rows band, drag on the header, Page keys | `rows.hide` (patch); shown again by the hidden-rows chip | `rows.pin` (patch) |
| **a clade** (a tree node, a summary row) | click a summary row: its rows | Shift+click | Cmd/Ctrl+click a node or summary row toggles its rows | double-click or Enter on a summary row opens it; a reveal opens the clades holding a row | double-click again, Alt+Left | — | — | fold, session state: click its node, Space on key-up (`clades.toggle`); `clades.collapseAll`, `clades.expandAll`, `clades.collapseTo` | the opening set: `clades.keep` (one patch) |
| **a legend entry** | click, Enter: its items (`legend.select`) | — | Cmd/Ctrl+click toggles | — | Alt+Left | — | — | its items filtered (`layer.filter`) | its color (`theme.value-color`), a ramp's domain (`scale.set`), the scheme |
| **a navigator** (the ruler rail by default; strip, ideogram, vertical overview, miniature) | no: at its scale a pixel is many positions (on an alignment's ruler, a drag above the rail selects columns) | — | — | a click outside the window centres (`navigator.center`); double-click fits or zooms there | Alt+Left; a far drag pushes a step | drag a grip (the other end held); sweep from outside (`navigator.zoom-to`); Shift+arrows while focused | drag the window's body; arrows, Home, End while focused | `navigator.form: none` (host), a node's `hidden` (patch); a hidden ruler falls back to the slim rail | its params by patch |
| **a window** (the `window` target, parts `body`, `start`, `end`, `outside`: a navigator's or a ring's) | — | — | — | — | Alt+Left | drag `start` or `end` | drag `body` | idle (quiet, not a target) when the view cannot move | — |
| **a bookmark or link** | — | — | `refs.add`; removing one (`refs.remove`) | runs its action; into a hidden stretch, its break and the reason | Alt+Left | — | — | — | its label (`refs.label`); its color: planned (themes) |

From the second slice's use gate, now designed: a section's direction (reversed, following its leader), another
system's letters on a section (one axis for protein and gene), collapsed introns ([axis.md](axis.md)). Still planned:
the history's toolbar.

**The history's one rule** (every design that touches it): a step keeps each section's viewport with its axis and layout
signature; onto a changed axis, layout or row order it is restored by position (the centre and the zoom; a text's anchor
and size; a rows band's top row by its key), clamped to what is drawn (a hidden centre goes to its break, a folded row
to its summary row), never silently dropped.

## 9. Chrome

- Built from `@molstar/ui` components and tokens only: track headers, legends, the hover card, menus, toolbars,
  navigators' frames. Section headers are `@molstar/ui` headers (32px); track header rows are rows (24px).
- A track header names its data and shows its state (hidden layers, filters, a pending load) as icons in fixed slots.
  Its actions are registry entries, drawn on the canvas or placed by a host as DOM parts over the header geometry; a
  section header is a 20 px line drawn when it has something to say (`section.header: auto`), replacing the 32 px
  header above; a secondary click opens Seq\*'s menu as the `headers` behavior's `contextClick` says (`everywhere` in
  the default spec), Shift+secondary click the browser's ([track-header-parts](track-header-parts.md); decided, the
  author, 2026-11-19).
- **A ruler in its sequence's own numbering has no label** (the author, 2026-11-19: "we dont need 'position'
  header on the track"): its numbers say what it is; a ruler in a declared numbering is named by it (Kabat, IMGT).
- Legends sit with their track or in one legend panel, never over data.
- Loading regions show a placeholder in place (`@molstar/ui`'s rule), never as empty data (EC-169).
- **The hover card** (the author, 2026-11-19: "Yes, lane HC now", the novice pass's first finding): a `@molstar/ui`
  surface (a squircle, its tokens) that says what the pointer rests on, in words, after the hover rests (the
  announcer's rest); one card, never stacked. Its text is the announcer's description of the place, never a second
  wording (the viewer `_spec`, Announcements): the system's label and the position in each of its numberings (ordinal,
  index, the declared ones: `heavy 97 · Kabat 93 · IMGT 105`), the letter, the item (its label, its extent, its
  category), a value track's value, a pair's two ends, a row's label, a grid cell's two positions and its value, a
  block's chain and copy. A dragged rectangle on a grid shows its summary while it is dragged (`gridSummary`: count,
  mean, min, max). Placed beside the pointer, never over the hovered cell or item (above it, else below), kept inside
  the view; hidden while a drag pans or zooms, on leave, on Esc; at the keyboard cursor when the cursor moves. Not in
  the draw lists, not in an export: DOM over the canvas. **Through a port** (ports.md, the author, 2026-11-19): the
  viewer computes the card as data (title, rows, the announcer's text) and its place (`placeCard`), headless and
  tested without a DOM; the hover card port draws it, `@molstar/ui`'s card by default, a host's own components and
  styles in its place. **Compact** (the author, 2026-11-20: "the hover cards should be more compact in general"): the
  title one line and said once, each row one text line, the value beside its label; an item's card names its track
  (a Track row, and the announcement), so two tracks' items of one label are told apart (FX4: Kabat's and IMGT's
  CDR-H1). **Capped**: past `hover.cardItems`
  (8) item rows the card shows the first and a "+k more" row; `HoverCard.more` keeps the rest, the live region says
  every item; a host's card through the port draws what it likes (the author: "app level overridable").
  styles in its place. A host turns it off (config) or adds rows behind its internal interface (extensibility.md);
  context menus stay later.
- **The keys panel**: `?` (and a toolbar button) opens a panel listing what the view answers to, generated from the
  keymap's commands and the gesture bindings (never hand-written), grouped (move, zoom, select, edit, view), each with
  its key or gesture and its name; Esc closes it. Its list is data too (`describeKeys(viewer)`), so a host can draw its
  own. A view says once, quietly, where to find it (the toolbar button's
  title: "Keys and gestures (?)").

- **No track drawn as an empty well** (the author, 2026-11-21: "Yes, group them"; FX12): a settled track with no items
  in what it covers is one line (`track.emptyPx` 18) naming it with a quiet "no items"; pending is never empty
  (EC-169); a settled source with no items is empty; an explicit `height` does not keep the well. The host's choice
  (`track.empty: note | collapse | keep`, render config; no spec param); a stacked block leaves out its empty tracks,
  the first block too (`blocks.emptyTracks: collapse | keep`).
