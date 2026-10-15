# Interaction and the default UIs

What a user does and what a host sees: input, viewports, hover, highlight, selection, focus, results, the keyboard, and
the default UIs. The behaviors named here are the defaults, each replaceable (`seqspec:viewer#behaviors`); their params
are in [reference/viewer.md](reference/viewer.md).

## Input events

The engine turns DOM events into reused input records and hands them to `ctx.input`; tests feed the same records. The
dispatcher hit-tests against the headless layout and emits **base events** on `ctx.input.events`, for behaviors (a
behavior author's API; consumers read the streams below): `pointer` (every input), `hover` (only when the hit's element
changes), `click` (a release within the drag slop, seen; counted 2 and 3 for a double and a triple click, each within
the double-click time of the one before), `drag` (start, move, end or cancel, once past the slop, or at the first move
of a press a behavior claimed), `leave`, `wheel`, `key`, `focus`, `menu` (the browser's context-menu request, claimable
as `doubleClick`; from the keyboard the `keyboard` behavior claims it and places it at its cursor), and `doubleClick`, a
query a click's subscriber asks once a click (never per move): a behavior whose double click at that hit means something
(a zoom, an activation) claims it. It decides nothing else.

- **Claims.** Subscribers are called in priority order; on a press or a drag start the first to `claim` owns the
  gesture. An unclaimed drag does nothing. A claimed press released within the slop and the double-click time is also a
  click, after its drag ends. Wheel and key subscribers are all called in priority order; one that sets `handled` stops
  the later ones acting, and the engine prevents the page's default.
- **One pointer at a time**: a second finger never drops a drag.
- **A press follows the pointer anywhere.** A press subscribes the page's moves and release until its release, so they
  reach its surface wherever the pointer goes (over the page, another surface, another element that took the pointer); a
  move with a button down that another surface's press owns is that surface's. During a press a point outside the
  surface is clamped to its bounds, how far past them kept on the record (`overX`, `overY`): every behavior reads a
  point on the surface, and a scroll past an edge reads how far. The engine prevents the browser's drag and drop where
  the surface is its target (a draggable ancestor is the host's).
- **A press always ends.** A release the page never delivered (a move of its pointer with no button down, another press
  from any pointer, as the next touch after a lost one) ends the press as a release at the last point it saw: a drag
  ends there, and it is never a click (nor counts toward a double click). A drag is cancelled only by the platform
  (`pointercancel`), a pause or its surface going away.
- **The hover follows the content** under a resting pointer: when a view moves under it (a wheel, a key, an animation, a
  linked section) or data arrives, the point is hit again, once a frame.
- Every record is reused: subscribers read it during the call and never keep it.

## Streams

A consumer outside the viewer (the default UIs, a host such as Prism) follows input through `ctx.interaction`: rxjs
streams `hover`, `press`, `click` (with `count`), `drag` (`start`, `move`, `end`, `cancel`, with the deltas and the
press's part), `menu`, `wheel`, `key`, `leave` and `activate` (a track header activated: the `headers` behavior, the
keyboard's Enter on a header). The dispatcher feeds them after the behaviors ran, so a hover's payload comes with the
hover already set and a press's with the behavior that claimed it (`owner`). They are observe-only: a subscriber watches
a gesture and never claims it.

- **One payload**, a fresh object a consumer may keep: `part` (where: a region of the surface by name, `data`, `header`,
  `gutter`, `navigator`, …, or a part of the chrome; the `section`, `track` and `group` refs, names not indices; a
  chrome node's `key`), `at` (what, on a surface: the `layer` ref, `lane`, `block`, `index`, display `row`, `split`, an
  `item` with its annotation's key and id), `place` (`x`, `y` in surface px, `clientX`, `clientY`), `buttons`, `button`
  and `modifiers` in Mol\*'s terms (`ButtonsType`, `ModifiersKeys`), `time`, `owner`, `origin` (`pointer`, `keyboard`,
  `chrome`) and `native`, the DOM event during the call only (a consumer may `preventDefault`).
- **No cost unobserved**: a payload is built only for a stream with a subscriber; with none, input allocates nothing
  more than it did without the streams.
- **Scoped**: `ctx.interaction.on(stream, { part, section, track })` is the stream kept to a part (or a list), a section
  and a track.
- **The chrome**: the default UIs report a click, a context-menu request, a hover and a leave on their parts (a view
  tab, a section's header, switches and hidden chips, a track's actions, the legend, the player, search, export, the
  issues, the tooltip, the hover card, the description, the annotations) through the same streams, origin `chrome`, with
  the part's section, track and `data-seqv-key`. A part's root carries its kind in `data-seqv-part`; a host's node
  tagged so reports too (`listenChrome`, or `ctx.interaction.chrome` directly). A focused chrome control's keys are its
  own.
- A behavior's bindings are listed, described, as Mol\*'s `Binding` (`ctx.behaviors.bindings()`), for a help panel.

## Replays and carries

A state gives its value to each new subscriber (a replay), and a reassembled document carries the hover, the selection,
the highlights and the results to its new set, with origin `graph` (a carry). Neither is anyone's move: a subscriber
acts on a `graph` emission only where what it holds or last sent differs; it leads nothing, and never moves a link's
lead, a result set's `current`, the hover's owner or a held click. A section a `graph` change moves (fitted after its
data shrank) is set again from its links' last leader: its `follow` link's, its scroll group's. The one exception is a
section's own focus change from a patch (`seqspec:format#focus`): it moves that section, which then leads its links, as
the same change from a snapshot does; when one patch changes the focus of several linked sections, the last of them in
document order leads. A reassembly that only moves the focus's block to another index (a block inserted before it)
changes no focus.

## Viewports

One viewport per section: its zoom (the cell, px a position), its horizontal scroll and its vertical scroll. `live`
changes during a gesture or an animation; `committed` emits once a change settles (on release, at an animation's end,
after a pause in wheel input, at once for a program's move), never per move. `interacting` is true while more changes
will follow.

- **Across updates**: a section whose axis is equal keeps its viewport; a changed axis keeps the center's location and
  zoom when it is still shown and the view was not the whole axis (`Viewport.OnAxisChange`), else fits; a new section
  opens on its focus (wrapped text: scrolled so its lines are shown), else fits.
- A move made before the section has a width is kept as one pending intent and applied once it has one; so is wrapped
  text's scroll to a focus made before it has lines.
- Default gestures: a drag on data pans (in text, it selects: Selection); Ctrl/⌘ or a pinch with the wheel zooms at the
  cursor; a horizontal wheel (or Shift with a vertical one) pans; a double click zooms in on an element or a position,
  then back out; the navigator jumps, pans and zooms. Reduced motion turns animation off.
- A navigator whose view has nowhere to go (the whole axis in view and zoom-in at its limit: it can neither zoom nor
  pan) is **idle**: drawn quiet, in its space (no grips, ticks or hover, its window the quiet line), its region still
  hit, but its press, jump, drag and double click do nothing and its cursor stays the default. It is live again as soon
  as the view can zoom (a narrower view, more positions), with no relayout.

## Linked sections

- **Scroll groups**: sections with one `scroll_group` (equal axes) share zoom and horizontal scroll, each keeping its
  vertical scroll. A viewport of another context may join a group.
- **Follows**: two sections set each other's window through a coordinate map (`seqspec:format#sections`); the follower's
  moves never lead in turn. A section on rings shows its whole axis: it is never moved by its followers and leads them
  only through its window; between two of them it passes the leader's view on.
- **The ring's window**: a section on rings linked by `follow` to a linear section shows that section's view as a
  window, a sector over its rings with a grip at each end on the rim (an outermost ruler's baseline). It follows the
  linear view however it moves; a plain drag of it pans the linear section, a grip's drag zooms it with the other end
  held (at the linear section's zoom limits and its content's ends the window stops); across a whole cyclic unit's
  origin it wraps. A press elsewhere on the rings (another block's included), or a click on the window, is the rings'
  own. On the ring `+` and `-` zoom the linear section about the window's center; the window's drags are never animated.
- **Hover and selection link every section** through projection, whatever their axes.

## Hover

One hover at a time, across all sources: pointer, keyboard, a host (`hover.setLoci`, a 3D view). A later owner replaces
an earlier; a stale owner's clear does nothing, so a late pointer-leave never clears a 3D hover.

- The hovered element is named in its own coordinate system, never the block's: a translation glyph's residue in its
  protein's, a cell's in its row's.
- A hit is by drawn geometry: the pointer's device pixel hits the cell or glyph drawn on it (a cell's snapped edges, a
  glyph's overhang into a gap on either side of its block), and a drag's end takes the same cell.
- `hover.toLoci()` is what is lit: the element closed by the channel's closure, a molecule residue as that residue in
  every unit. These loci are closed already: passed back, they take `element`. It is a getter: a large set resolves only
  when read, once a version (and again on a new set shown or a host's system registered or removed: a `Coordinates`
  hover lights once its system is there). `hover.size()` says how large before it resolves (`LociSize`: items, and the
  positions covered, an estimate), so a reader can leave a large set alone or summarize it; `lociSize` reads any loci's,
  and the DOM's `seq-hover` and `seq-selection` carry it.
- Its version bumps when the element, its closure or its owner changes, or a relayout refreshes the same element; never
  on a move inside one.
- On rings, an item's label outside the rings is hit as its item: hovering or clicking it hovers or selects the item, as
  on the item.
- **A header's hover** stands for what the `headerHover` behavior says, a viewer-spec behavior like any other: removed,
  a header's hover is the header alone (its row band). By default a feature layer's header (ranges, pairs) lights its
  shown items: placed on a block, its filter and row subset kept, a split row's header that row's, a rows track's row
  header the items placed on its row, a group's header its shown tracks' (a hidden track's not); a continuous layer's
  (letters, values, a matrix) the row band alone; past `maxItems` items the row band alone. Its set is item loci with no
  closure, given to linked views by `hover.toLoci()` when `forward`; it is counted when the hover lands (never per
  move), lit runs built only when it lights, and described again when the view or set changes under it. A streamed
  layer's items are not described yet.
- **A data hover's row band** is the `hoverRow` behavior's: a hover on a track's data (linear sections) lights the
  track's row at half a header hover's strength, so a header's still reads as the whole track (by default half over the
  header and a quarter over the data: `header`, `data`). The row is the track's box, a split row's, a rows track's row,
  or the band its renderer gives (`rowHover`: a matrix its hovered state's band; the ruler and a dot plot none).
  Navigators light none. It draws where the pointer is and nothing else: the hover's kind, its loci, its size and its
  version are the element's, so no reader (a 3D link, a host) takes it for a header's hover. An alignment's hovered row
  and a header's band are drawn as before, whatever its params; a selected row stays distinct (the selection's tint and
  outline).
- **Hover emphasis**: the hovered element alone (an item's own pieces, a cell, a matrix band) is outlined on the overlay
  in an ink derived from the color it is drawn in, with at least 3:1 contrast to the surface and, where the fill's
  lightness allows, 1.5:1 to the fill, in both skins. A section's style may turn it off.
- **A strip lane** in the text view (`seqspec:viewer#the-text-view`): a point on the strip, or in the gap above it, is
  its layer's hit (the hit's `track`, `layer` and `index`), so the hover names the position and the hover card shows the
  layer's value there, as on its linear track; hover emphasis outlines the position's cell of the strip. Beside the
  letters the hover is the letters'. A track drawn in full in the text view is hovered as in tracks: a glyph, an item's
  bar or label, a position under a line or an area, by drawn geometry on its line.
- **A dot plot's cell** (`dot_plot`, `seqspec:format#layer-types`) is two positions: the one under x, as on any track,
  and its **partner** under y (the hit's and the hover's `pairBlock`, `pairUnit`, `pairIndex`). The hover lights both,
  each closed as the hover is, in every section and through links (a 3D view lights both residues); its label adds the
  partner and the cell's value. The track shows every channel's lit positions both ways, a column at x and a row at y,
  so the hover is a crosshair through its cell and a selection of two ranges a rectangle where they cross.

## Highlight

**Highlight channels** are free names chosen by their owner: a host's own, `links:<link>` (a command link),
`set:<section>/<n>` (a section's highlight entry), `results:<name>`, `results:<name>/<k>` and `results:<name>:current`;
the overlay draws their union, in a style (`highlight`, `accent`, `band`), and clearing one never clears another.
`highlight.set(channel, target, opts)` with loci or `Coordinates`; anything else is refused (several loci: their union,
or a channel each). A channel that cannot be carried to a new set is dropped alone. A command's target naming several
kinds (an annotation and an alignment's columns) lights a channel each, `<channel>`, `<channel>/1`, …, set and cleared
together (only the parts it set: a host's own `<channel>/1` stays). A behavior clears its channels when removed.

**Following a highlight.** The `revealHighlight` behavior (not in `DefaultSeqViewerSpec`: a host whose other views light
the sequence adds it) reveals a highlight on the channels it follows (`channels`: names, or prefixes ending in `*`;
`links:*` by default), and an external hover of the owners it follows (`hovers`: a linked 3D view's `3d`; none by
default), once it has rested `dwellMs` (200) unchanged, as `ctx.reveal` does (`pad`, `animate`), in each section that
shows it and does not show it already. Never in the section the pointer is in, never during a press or a drag, never
while this viewer's own hover is held (a link lighting it back); a highlight cleared or replaced during its dwell is not
revealed, a newer one starting the wait again.

## Channels

Every channel lights through one pipeline on reused buffers: the element, its **closure** if the channel closes
(`seqspec:model#closure`), **granularity** (selection only), **projection** onto each block
(`seqspec:model#projection`), then **marks**: per block, sorted disjoint runs, so duplicates collapse and there is no
cap, and **carets**, the boundaries lit (an insertion's place, a boundary item, an insert seen from the other side of
its map), per block sorted and distinct, each drawn as a bar on the edge between two cells in the channel's tint.

| Channel                          | Closure                                                                  | Lights                                                                                                                        |
| -------------------------------- | ------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------- |
| pointer and keyboard hover       | the behavior's `closure`, default `containing`                           | the hover; in a rows section, the column and row band; on a track's data, its row band (`hoverRow`, half a header's strength) |
| external hover (`hover.setLoci`) | its option, default `element`                                            | the hover                                                                                                                     |
| header hover (`headerHover`)     | none: its items' rows                                                    | its shown items within `maxItems`; the row band at full strength over the header, half over the data                          |
| item hover                       | none: the item's rows                                                    | the hover; a pair its own arc                                                                                                 |
| selection                        | the `select` behavior closes, then widens; the verbs store loci as given | the selection; rectangles in a rows section                                                                                   |
| highlight channels               | its option, default `element` (command links pass `containing`)          | highlight, accent or band                                                                                                     |
| focus and scroll-to              | none                                                                     | the hull of the projection per block                                                                                          |
| labels and the tooltip           | —                                                                        | the element in its own numbering, then each containing element                                                                |

- A channel's marks are rebuilt only when that channel or the view changes, never on a layout change.
- The **primary** of a channel on a block is its first lit position the layout shows: the navigator's tick and the
  tooltip's anchor.
- Streamed data: a channel over a pending region keeps its extent and lights it when the data arrives. Loci naming an
  evicted node light nothing meanwhile, and keep the node from eviction while held.

## Selection

`selection.current` holds one loci per kind; `pending` is a drag's live range (`remove`: one that will remove it, drawn
as its outline alone).

- **Verbs** (`set`, `add`, `remove`, `toggle`, `clear`) store loci as given; `Coordinates` are converted through the
  registered systems.
- **The click policy** (the `select` behavior): a click on data or a column selects the element, or removes it when it
  is the whole selection of its kind; Ctrl/⌘ adds or removes it; a click on nothing clears. Shift on a position or a
  column adds or removes it (`shiftPositions`, `toggle` by default); with `extend` it extends the range, in every
  geometry: the selection on the clicked block (as lit: positions, an item, columns, made in any view; its runs as
  shown, a cyclic origin's two parts) to the clicked position, one end kept: outside it its far end, inside it the
  anchor's (the last click's position or the last drag's start, a position of the model) where that is one of its ends,
  else the farther; with several runs there, from the anchor. With nothing selected on the block it extends from the
  anchor on another block of the section, else it is a click. Shift on an item or a row (`shiftItems`) adds or removes
  it. Where a double click at the point means something else (`doubleClick` claimed: a zoom on a linear section's data),
  a single click's change waits for the double-click time and applies once no second click comes (at once at the next
  click elsewhere or a key, never after Escape or another change); a double click then never touches the selection.
  Elsewhere a click selects at once, and a double or triple click is that click again. The element is **closed** (a
  nucleotide's residue), then **widened** within its own system by `granularity` (`position`, `range`, `unit`,
  `copyGroup`, `column`, `row`), so a click selects what the hover lit. "This codon in every copy" is `containing` with
  `copyGroup`.
- **Drags**: on data a plain drag pans; on a column track (a ruler, a logo) it selects columns; Alt selects a range;
  Shift with a drag means what it does with a click (`shiftPositions`, over positions, whole columns on an alignment):
  `toggle` (the default) adds the dragged range, or removes it when the press was on a selected position; `extend`
  selects from the anchor to the drag's end. A drag selection closes and widens on release. Who claims a drag is params
  of the `select` and `dragPan` behaviors, never a rule of the dispatcher. The range runs from the press to the pointer
  at each drag event, the first (past the slop) and the release included, so a fast pointer selects what it was released
  over. Its end follows the point as clamped to the surface: on a linear section the column under its x (on rows, the
  row under its y), so outside the surface above or below it still follows x; past the data area's left or right a
  linear section pans, and past a box section's top or bottom rows scroll, each faster the further past, the end going
  with the view.
- **Text** (a wrapped section) is selected as text is: a plain drag on it selects positions, its ends following the
  pointer line by line to the nearest position (past a line's end, in a gap or on an inline label, the one on the drag's
  side of the point: forward the next, backward the one before; over another block's lines, its own block's first or
  last), by x on the line at the point's y: outside the surface above or below, the line at its edge (above the first
  line the first, below the last the last). Past a box section's top or bottom edge its rows scroll while the drag
  lasts, faster the further past, at a speed in time (never per frame); they stop when the point is back over them, at
  their end, or when the drag's surface loses focus. A section sized to its content does not scroll during a drag, nor
  does the page. A double or triple click is a click (nothing takes it in text); Shift-click extends as above. Alt,
  Shift-drag (adds), Escape and the keyboard's Enter act as on any data; an underline under the letters is its item. The
  text sits inside an inset (`wrappedInset`) and a wider strip on its right; a drag there scrolls the rows (a touch
  reader's grip, the rows' vertical only), as do the wheel, the scrollbar, the keyboard and a selection drag past their
  edges. A limitation, not the design: on touch a drag on the text selects too, so a touch reader scrolls a box text by
  its inset and right strip, and a text sized to its content only outside its canvas.
- **On a dot plot** a plain drag selects a rectangle (the `dotPlotSelect` behavior, ahead of `select` and `dragPan`): on
  release its two ranges, the positions it spans across and down, closed as a click's (`closure`), replace the selection
  (Shift adds them); the view does not pan. While it lasts the pending selection holds both (the plot's `plotLayer`, and
  `yBlock`, `yStartShown`, `yEndShown` beside the x range), drawn as the rectangle and, across the tracks, as both
  ranges. A click on a cell selects its two positions. The cursor there is a crosshair (`cursor`'s `dotPlot`).
- A band over many rows stays one rectangle. Selections follow their set across updates.

## The track menu

A track's menu is the host's; Seq* says when one is asked for: `interaction.menu` is the browser's context-menu request
(a secondary press, Ctrl+click on macOS, the menu key, Shift+F10), hit at its point, with the part, the track's ref and
the place. A subscriber that calls `native.preventDefault()` during the call shows its own menu instead of the browser's
(only that request can suppress it); with none, or where none prevents it, the browser's shows. From the keyboard (the
`keyboard` behavior's menu key, or a keyboard context-menu request) it is asked at the cursor's cell, origin `keyboard`,
and once per ask in every browser: the browser's own request that follows a key's on the same section, with no press,
later key, blur or leave between, is not asked again and has its default prevented as the key's was; the key is handled
only where a subscriber prevented it. The default UIs open their menu at the track's header from the keyboard and leave
a right-click to the browser. A macOS Ctrl+click is also a press: `select`'s `toggle` still sees it.

## Focus

`focus.fromLoci` (and `fromExtents`) projects without closure, links and rows included, and frames every block it
lights, from the first lit position to the last, at the largest cell that shows them all; on a cyclic system, the
shortest arc. A section's `focus` node is its initial focus; a patch or snapshot that changes it moves the viewport,
animated by the `focusViewport` behavior. Wrapped text does not zoom: the behavior scrolls the lines of the focus into
view (`reveal`, on by default), over the same transition as a zoom: every line from the first shown position it lights
to the last, on every block it lights (an end outside what is shown is the nearest shown position inside the range; a
reversed piece's lines in their shown order). Where those lines already are shown nothing moves; else they come
`revealPad` px clear of the edge the scroll brings them to, the first line at the top where they are taller than the
view. `ctx.refs.toLoci(ref, args)` resolves a ref to loci (a unit, molecule, complex, row, annotation, track or layer);
range args outside its domain (index args, or labels not there or resolving outside it; a complex's, outside every
member) give none, never a clipped part. A molecule's args name its positions in each of its units.

**Reveal and focus.** A focus frames its target: tracks zoom to it, and it is the section's state (`focus.current`).
`ctx.reveal(target, opts)` (and the `seqv:reveal` command) only brings a target into view, in every section that shows
it or `opts.section`: it changes no selection, focus or zoom. Text scrolls its lines in as a focus does (`pad` px clear,
12 by default); tracks pan at their cell, the target centered, or its start `pad` px in where it is wider than the view.
Where it is shown already (whole, or, wider than the view, its start) nothing moves. It is animated as a focus is
(`animate`, `durationMs`), origin `api` unless given.

## Results

`ctx.results` shows coordinates as results, with no search logic: an application finds them, the viewer lights them
(`results:<name>`), focuses and accents the current one, and steps through them (`2 / 14`). Several named sets coexist
(a search, a section's bookmarks). Only what the application saves becomes a location set: a `ranges` annotation,
proposed as a patch (`results.save`, `bookmark.add`).

## The keyboard

Each section's surface is one tab stop; Tab never moves inside it. Inside, the `keyboard` behavior's cursor reaches
every position, row and header, and sets the hover (owner `keyboard`), so the label, the live region and the tooltip
follow it.

| Key                       | On data                                                           | On a header                            |
| ------------------------- | ----------------------------------------------------------------- | -------------------------------------- |
| ← →                       | one position, across block gaps; ← from the first onto the header | → onto the first position in view      |
| ↑ ↓                       | the entry above or below, same position                           | the header above or below              |
| Home, End                 | the header; the entry's last position                             |                                        |
| `+` or `=`, `-` or `_`    | zoom at the cursor                                                | zoom at the center                     |
| Alt+= Alt+−               | row zoom                                                          | row zoom                               |
| Enter, Shift+Enter, Space | select, as a click; extend; as Enter                              | activate; activate; expand or collapse |
| Shift+F10, ContextMenu    | the track's menu, the item's actions first                        | the track's menu                       |
| Escape                    | clears the selection, puts the cursor away                        | the same                               |

Every move keeps the cursor in view. A relayout that removes its block, entry or row moves it to the nearest one shown.
A relayout moves the cursor, never a hover another owner holds: while the pointer hovers, the keyboard's hover waits,
and comes back at the next relayout or key once nobody holds one. On rings ← → go round and Home is the first position;
`+` and `-` zoom the section the ring's window shows (`seqspec:interaction#linked-sections`), and nothing without one.

## Color schemes at runtime

Where more than one registered color scheme fits a layer (`seqspec:format#color-schemes`), a reader picks one from the
layer's legend group title or its track menu (Colors ▸). The pick is UI state per layer, never in the document or a
snapshot, kept across patches and cleared by a new document; with `UI.ColorSchemeWrites` (or the command's
`write: true`), it is written into the document as a patch. A scheme that does not fit the layer is refused; a pick kept
across a patch after which it no longer fits draws the layer's own color. A layer's color range may be set the same way,
as UI state.

## The default UI

Two UIs with one slot contract: React (`/react`) and the DOM (`/dom`, with the `<seq-viewer>` element). Both are thin:
parts read the context's states, mount surfaces and send actions and commands back; they hold no state beyond leaf state
(a menu open). Anything a UI could not do without its own state is a missing state of the context. The parts, by slot:
`ViewTabs`, `SectionHeader` (a section's label, its filter, its `HiddenChips` and its `SectionSwitches`, each switch an
action), `TrackActions` (the track menu, `UI.TrackActions`), `RowControl`, `Tooltip`, `HoverCard`, `Legend`,
`SnapshotPlayer`, `Description` (plain text only; markdown is a host's part), `Search` (when given a search function; it
keeps its results as a location set and clears every set it kept, in one patch, so removing one never needs the legend),
`Export` (when the context registers `view.exportSvg`), `Issues`, `Busy`, `EmptyState`, `ErrorState`, and `Hint`,
`ScrollArea`, `Menu`; `ViewSurface` is the canvas. `AnnotationsPanel` (the Annotations part,
`data-seqv-part="annotations"`) is a part a host places (no slot, in no chrome's composition): a section's bookmarks set
as annotations, added from the selection with a label, a color and a line style (one of the format's `line_style`s,
`solid` to start; written as the row's `line_style` column, which a layer shows through `field_line_style`), listed (a
row reveals its range, `ctx.reveal`) and deleted, each a patch through the bookmark actions (one undo step). A host that
names several sets (`places`: each a set's ref and the name a person sees) gets a choice in the add row (the first set
to start with); the list shows the sets one after another, each row marked with its set's name, and a delete acts on its
row's set. A set's own track decides where it is drawn: a `ranges` track ordered before the letters' draws above them,
one after, below. The set is drawn by the document's own `ranges` track over it, colored by the rows' `color`
(`field_color`); without a track the set is a `highlights` entry, as `bookmark.add` makes it.

- **Every part is a component a host places**: inside the viewer's composition, beside it, or not at all, against the
  viewer's context. React's take the context as a prop, under the provider (its slots, primitives and surfaces) or
  without one (the default primitives); the DOM's mount into any element, with a mounted UI (`mountSeqViewer`'s, the
  element's `mountedUI`, or `createDomUI`) or none. What a part shows never depends on the chrome or on which other
  parts are placed. The canvas overlays (`Tooltip`, `HoverCard`, `TrackActions`, `RowControl`) are placed at the canvas,
  inside the root; without a UI (no provider, no mounted UI) they show nothing. A host's own composition starts from the
  root (`SeqViewerRoot`, `mountSeqViewerRoot`: the box the overlays are placed in, and the live region).
- A section's switches offer only what changes it: a numbering that labels every position of the section's units as one
  already offered (`label` and `ordinal` on most structures, often `auth` too) is left out unless the section names it,
  and with one numbering left there is no numbering switch.
- A control that hides while it has focus hands focus to its neighbor (Keep and Clear kept to the search field): focus
  never falls to the page.

- Track labels, trees, rulers, data columns and navigators are drawn on canvas; DOM parts are interactive chrome only.
- Pointer-rate values never re-render: the tooltip and row controls move by custom properties set in a subscription,
  from cached strings.
- React and DOM parts produce the same structure and text for the same context state (tested).
- `Issues` counts errors and warnings as issues, its badge marked by the most severe; notes alone are a quiet "N notes".
  Its list keeps the most severe when it is cut (errors, then warnings, then notes), groups them by ref, and lists the
  notes last under their own heading.

## Slots and chrome

A host replaces any part through its slot with its own (same props), `"none"` (left out wherever it is placed, inside
another part too; `Hint`, `ScrollArea` and `Menu` keep their default) or `"default"`. `RowControl` has no default part.
**Chrome** is a UI option (`full` or `clean`), else `UI.Chrome` (`clean` by default): the default composition the layout
places (`ChromeCompositions`), never a mode of a part. `clean` places the canvas (hover, selection, the tooltip and the
hover card, the navigator), the player, the description, search, busy and the states; `full` adds the tabs, section
headers, track menus, the legend, export and issues (`CleanParts`). A hidden legend and `UI.TrackActions: none` leave
theirs out. A slot naming a part (or `"default"`) places it where `full` would; `HiddenChips` and `SectionSwitches` are
placed by `SectionHeader`, not the layout, so theirs place nothing on their own. A value that is not a part falls back
to the default (`seqv.ui.slot`). The `Menu` contract is controlled: opening focuses the first item, arrows move, Escape,
Tab, choosing and a press outside close it, and focus returns to the trigger.

## The legend

Collapsed by default (`UI.Legend`): a "Legend" toggle; a track header's tooltip carries that track's legend. Open, it
lists a group per shown track with something to list: colors, glyphs, a ramp for a continuous palette (never a list of
values; a discrete palette's bands as hard steps), sample sizes. Tracks whose groups would be identical share one,
titled by what their labels share ("Genes (5)"). A group shows `UI.LegendItems` items and the legend `UI.LegendGroups`
groups before "+n more". Below them, one row per section highlight entry, each with a remove button. A group's title is
the color scheme picker where more than one scheme fits.

## The hover card

Items with actions get a card: a host names item actions in `UI.ItemActions` (registered actions whose target is an
item); after `UI.HoverCardDelayMs` on such an item the card opens with the label and a button per action. It stays while
the pointer crosses into it, and closes on leaving both, Escape, a press outside, a pan, zoom or resize of a section, or
another document (an edit included); a reassembly of the same document (a skin, an arrival) leaves it open. With no item
actions (the default) there is no card. The card follows the pointer only; the keyboard reaches the same actions in the
track's menu (Shift+F10), listed first while its cursor is on such an item.

## The DOM element

`<seq-viewer>` makes a context from its `spec` (or takes a host's `context`), mounts the UI and loads its `document` or
`data`; setting a member changes only what it names (a new document object is set once per context). Its `mountedUI` is
what parts a host mounts outside it share: null while disconnected, replaced on each remount (which `seq-ui` announces),
its surfaces the element's across remounts. A mount that fails (a spec whose registration conflicts) shows its error in
the viewer's place and fires `seq-error`; setting a member tries again. The viewer's own format entries are known by a
brand, so a page holding two copies of the viewer (an update while developing, a host bundling it twice) registers
either's spec. Moving the element keeps its context; a context it made is disposed a task after it leaves the page. It
renders in a shadow root with the skin, and fires composed, bubbling events with shallow-frozen details: `seq-context`
(the live context), `seq-ui` (`{ui}`), `seq-error` (`{error}`), `seq-selection` (the selection value, its `Coordinates`
and system descriptors), `seq-hover` (at most once a frame; loci, `Coordinates` and descriptors computed on first read,
null once the hover has moved on: read them in the listener), `seq-snapshot` (`{ref}`). They are fed from the states
(the hover, the selection, the snapshot shown); a DOM host follows input itself through `context.interaction`
(`seqspec:interaction#streams`), the default UI's chrome included (a section's nodes carry `data-seqv-section`, a
track's actions `data-seqv-track`).

## The skin

The viewer's look: tokens (`--seq-surface`, `--seq-ink`, `--seq-accent`, `--seq-hover`, `--seq-selection`,
`--seq-pending`, `--seq-focus-ring`, the fonts, …) and the parts' styles under `.seqv-*` classes, both public API, in
`@molstar/seq-viewer/skin.css`. Defaults are at zero specificity, so a host sets tokens on any ancestor. Light and dark
define the same tokens, following `data-theme` on the root, else Mol\* Prism's `data-msp-theme` (`blue` is dark), else
`prefers-color-scheme`, with no call from the host. The canvas reads the tokens where it is mounted, again (at most once
a frame) when a theme attribute or class changes on the viewer's root or an element above it, or when the host asks
(`refreshSkin`: an inline color or token, which it does not watch, since hosts write inline properties during drags); it
is dark when `--seq-surface`'s WCAG luminance is below 0.179 (a transparent surface decides nothing: the page's theme
does). Data colors are themes, never skin; ink on a colored cell follows the cell's luminance. "Not loaded" is drawn in
`--seq-pending`, never the ground. A window (a navigator's, a ring's) is filled in `--seq-window` with its edges and
grips in `--seq-grip`, never a selection's inks. A feature (`ranges`, `points`, `secondary`) its layer gives no data
color (a `rules` theme with no rule) is drawn in the faint ink, so the skin decides it.

The viewer is a canvas, and its chrome is canvas chrome. What sits on the canvas or opens over it is flat glass: a
translucent tone of the surface and a hairline edge, no drop shadow; what floats over drawn content is blurred, and what
opens for a moment (a menu, the issue list) adds a soft shadow. The legend has no box: its text sits on the canvas. Text
is 12px, hints 11px; surfaces have an 8px radius, controls 4px; 8px lie between the canvas's edge and a strip, and
between a strip and the canvas's content. No 1px borders: the picked item is an accent wash with accent ink, a control a
tone of the ink, a swatch a hairline. Cut-off text fades rather than ending in an ellipsis; checks and chevrons are
icons. The chrome's own tokens are optional, each derived where it is used from the tokens above when unset:
`--seq-chrome`, `--seq-chrome-panel`, `--seq-chrome-filter`, `--seq-chrome-edge`, `--seq-menu-shadow`, `--seq-control`,
`--seq-selected`, `--seq-selected-ink`, `--seq-font-size`, `--seq-font-size-small`, `--seq-radius`,
`--seq-radius-control`, `--seq-canvas-edge` (between the canvas's edge and a strip, on every side) and
`--seq-chrome-gap` (between a strip and the canvas's content, and between strips).

## Accessibility

- Each surface is focusable, labeled by its section's label, with the role description "sequence view" (the surface sets
  them, so one a host mounts on its own is too); a polite live region announces the keyboard cursor's label, never the
  pointer's.
- The tooltip is placed at its owner's anchor: beside the pointer, or at the keyboard cursor's cell.
- **Keyboard focus is never lost**: a part rebuilding a node keeps focus on the same key; when the focused element goes
  away, focus lands on the next element of its group, else the previous, else the surface, else the root. An unusable
  control is `aria-disabled`, never `disabled`.
- **Hints never block the pointer**: a hint takes no pointer events, hides when the pointer leaves its trigger, on a
  press, scroll, wheel or key, and never while the trigger's own menu is open; it shows after a rest (at once on
  keyboard focus, or right after another hint). Every control has an accessible name; hints go through `Hint`, never
  `title`.
- The default skin keeps ink on surface at WCAG AA in both schemes; `--seq-conflict` is at least 4.5:1.
