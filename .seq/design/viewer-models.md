# The viewer's models

## Intent

The second review found the viewer's models unlisted, session state deferred though the first slice builds it, the
viewport and follow undesigned though "p53's protein on its gene" needs them, hover and selection without a model, the
gesture bindings without a type (`@molstar/ui`'s keymap knows key chords only), hover as a command allocating on every
pointer move, and no field for a command's change kind. This design lists the models and settles those.

## Decision

### The models

All are `ReactiveModel`s on `@molstar/core/reactive`, DOM-free, tested in Node ([sources/molstar-ui](../sources/molstar-ui/README.md)):

| Model | Holds | Change kind ([changes.md](changes.md)) |
|---|---|---|
| `SpecModel` | the spec's root (immutable), its revision (the root's hash), the preview (a patch over the root, applied again when the root changes; [patches.md](patches.md), Previews); makes a commit per applied change (`SpecCommit`, the port's checkpoint) and keeps the recent roots by revision (`revisions.keepRecent`, 64) for a host's rebase; shows a story step over the base and the preview (changes.md); applies `repr` patches (validated with the language and Seq\*'s checks, all or nothing; stale ones refused) and records each result as a revision; a participant of the host's history through Seq\*'s history port (a change outside a run is reported to the host) ([changes.md](changes.md), [patches.md](patches.md)) | spec patch, preview |
| `ViewModel` | the view shown; per section its viewport (the one record of [layouts.md](layouts.md)); scroll groups and follows; the folded clades per rows layer ([grids-and-rows.md](grids-and-rows.md)) | session |
| `HoverModel` | the hover: a Location, an element key, its owner, its system | session |
| `SelectionModel` | the selection: element keys and a Location, the anchor, the cursor | session |
| `HighlightModel` | named highlight channels, each a Location or element keys | session |
| `StoreModel` | the data store ([data-store.md](data-store.md)) | data |
| `ConfigModel` | the configuration, read-only after composition | configuration |

- **Hit type**: the viewer's hit is the render hit ([layout-and-draw-lists.md](layout-and-draw-lists.md)) with a `target` added.
- **Session state** is the `ViewModel`, `HoverModel`, `SelectionModel` and `HighlightModel`'s state, expressed in
  addresses and Locations; a host may save it beside the spec (`viewer.session()` and `viewer.restore(session)`), never in
  it.

### Viewport, scroll groups and follow

- **A viewport per section**: the visible range as offsets along its axis in display order (fractional; 0 at the first
  position shown) and the vertical scroll; zoom is the range's width against the section's width. The record is
  `{ start, end, top, width, pending, scale?, bands?, y? }`, owned by [layouts.md](layouts.md): `bands` by band (`bandTop`,
  a band's own scroll, and `rowScale`, px a row: an alignment's rows scrolling in their band), `y` a grid's second axis
  with its own window; each layout reads it its way (text: `start` an anchor, `scale` the letter size; rings: the window),
  and the `ViewModel` calls the section's `nav` for the arithmetic. Why offsets, not
  positions: a section's `range` and a cyclic range over the origin read like a whole system. A saved session stores
  positions (and a band's top row key and `rowScale`), never offsets, since offsets change meaning when the axis does.
- **The view's history** keeps whole viewport records with their axis and layout signature. A step onto a changed axis,
  layout or row order is restored by position (the centre and the zoom; text its anchor and size; a band its top row),
  clamped to what is drawn (a centre now hidden goes to its break, a folded row to its summary row), never silently
  dropped ([axis.md](axis.md), the one rule in every design that touches history).
- **A scroll group** shares the range only among sections on equal axes (EC-113; `viewport.groupEqualAxes`, true); each
  section keeps its own vertical scroll and width.
- **A follow** is a section whose viewport is derived from another's through a map: the leader's range projected along
  the path between their systems, at scale 1, 3 or 1/3 (a codon step); a move asked of a follower leads: it moves the leader through the inverse of the follow, so the follower lands where it was asked (the author, 2026-11-18: "could we have following in both ways? .. i think we should"; `follow.followerMoves`: `lead`, or `ignore`, or `detach`; EC-114). The ends are the continuous images of the
  leader's (maps.md: continuous within a run), the hull only where an end maps nowhere. A cyclic follower's range past the
  origin wraps (`follow.cyclicArc`: `wrap`, its viewport runs past the axis total and `pieces` splits it; `hull` the
  alternative). For "p53's protein
  on its gene", the gene section follows the protein section through the codon map.
- **Viewports survive updates**: an equal axis keeps its viewport; a changed one keeps center and zoom when the center
  is still shown, else fits (EC-111); a move before the section has a width is kept as an intent (EC-112).
- A viewport is clamped to `viewport.maxPxPerPosition`, and a resize keeps the zoom (`viewport.onResize`: `keep-zoom`, or `keep-range`). Collapsed introns are a section's `pieces` ([axis.md](axis.md)): offsets go through the pieces and the direction, a
  break takes no width.

### Hover and selection

- **Hover** is a method on `HoverModel`, `hover.set(hit, owner)` and `hover.clear(owner)`, never a command: it runs per
  pointer move and allocates nothing (reused records). One owner at a time; a stale owner's clear does nothing (EC-141).
- **Projection**: the hover's Location is projected onto every linked section's system when the resolved view links
  them, along their paths (static or overlay), lazily per frame for the sections drawn, under the copy rule (EC-25); a
  nucleotide lifts to its codon (closure, EC-30, EC-144).
- **Selection** holds element keys (items, rows, chains) and a Location (positions); the grammar's toggle acts on
  element keys, extend on positions from the anchor in one system ([ux-foundations.md](ux-foundations.md)).
  Projections of the selection are drawn on linked sections as the hover's are.

### Behaviors (the first list)

`hover` (pointer to `HoverModel`), `select` (clicks and drags to `SelectionModel`), `navigate` (drag, wheel, pinch,
keys to viewports), `follow` (keeps followers in step), `cursor` (the keyboard cursor), `focus` (focus and reveal), each
with its params as config items; `hoverCard` and the keys panel (ux-foundations, 9. Chrome; the third slice's HC). Later:
`contextMenu`, `snapshots`, `reveal-on-highlight`.

**`linkedPairs`** (decided, the author, 2026-11-19: "Behavior, on by default"; asked on Complex › Contacts: "should
highlighting a range also highlight the arcs that belong to it? at least via an optional behavior of the viewer that
runs the queries … basically an extension"): on a hovered or selected range element (a feature, a CDR), it queries the
pairs with an end in the element's extent and writes them to the highlight channel, on every pairs layer shown (arcs,
chords, dots), cleared when the hover or selection goes. Built on the public API only (the query language, the
highlight channel, the hover and selection models), so it is the proof that an analysis is a behavior a host lists:
listed in `DefaultSeqViewerSpec()`, removed to opt out. Its params: `on` (`hover`, `selection`, both by default), the
pair layers it reads (all by default). Why a behavior and not the pairs layer's own rule: a host chooses it, the same
query serves the selection and later the genome ring's "everything touching this block", and the layer stays a
drawing. The pairs layer's own rule stays: a hovered position (a base, a cell's two positions) lights its partners.
What it needed, decided (the author, 2026-11-19, LP stopped on two gaps in the public API): a pairs annotation answers
the query model, a pair one item with a row per end, each row with its own `system` (`QueryItem.system`, optional; a
ranges row's is its annotation's target), so `(items :overlapping <extent>)` on a pairs node gives the pairs with an end
in the extent ([locations-and-queries.md](locations-and-queries.md)); a highlight channel carries element keys beside
its Location (`EmphasisHighlight.elements`), and the pairs layer draws the highlight pass for the keys a channel lists,
in its colour, arcs, chords and dots alike, the SVG export too. A selection's linked pairs draw in the highlight pass in
the channel's colour, never the selection's: they are related to the selection, not in it (a copy of the selection
leaves them out). The hover's query is cached per hovered element key (a hover allocates nothing a frame past its
first). Why: both are general (any host lights items by key; any host asks which pairs touch a range).

### Gestures, bindings and commands

- **A gesture binding** is `{ gesture, modifiers, target, action }`: gestures are `click`, `double-click`, `drag`,
  `wheel`, `pinch`, `key`; targets are the hit's kind (cell, item, row, column, track, empty, axis); actions name either a
  model method (hover, drag-select, pan: continuous, allocation-free) or a command (discrete). Key gestures go through
  `@molstar/ui`'s `Keymap`; pointer gestures (including `hover`) through this table, a `key` binding in it being reported; a `'*'` modifier matches whatever is held and is tried after exact matches (hover is bound with it, Shift+wheel pans, Alt is unbound for the platform); the selection's clicks and drags are in the default preset (the second slice, lane
P; the `genome-browser` preset keeps the clicks and pans on a drag). The table is filled by the grammar's presets
  ([ux-foundations.md](ux-foundations.md), section 8).
- **Commands** are `@molstar/ui/commands` definitions with Seq\*'s metadata: `changeKind: 'patch' | 'session' |
  'navigate'`, kept beside the definition (`@molstar/ui`'s `CommandDefinition` has no such field; Seq\* keeps a table by
  command id, and proposes the field upstream if it proves general).

### Issues

- **An issue of a pending address is not shown** (decided by the author, 2026-11-20, "As proposed"): core reports the
  static model's issues (`model.issues`); the viewer holds those that name an address whose store entry is `pending`
  ([data-store.md](data-store.md) owns the state) until it settles, so `seq.resolve.length.unknown` does not flash while
  a source opens. Why: the viewer is the only holder of "pending".
- **`viewer.modelIssues` is a state**, the model's issues as the viewer says them (one of a source still opening held
  back), beside `viewer.layoutIssues` (the frame's); a host and the example stories' status line read it.

### Configuration items of the first slice

`foundations.preset` (`comfortable`), `density.letterMinPx`, `density.dotMinPx`, `bindings.preset` (`default`),
`store.cacheBytes` (64 MB), `store.inFlight` (6), `store.tileBins` (1,024), `hover.projectLinked` (true), `hover.closeCodons`, `pointer.rehit` (true: a resting pointer is hit again once per frame in which the layout was rebuilt, so a hover follows a view moving under it, EC-143), `viewport.maxPxPerPosition`, `viewport.groupEqualAxes` (true),
`viewport.onResize` (`keep-zoom`), `follow.followerMoves` (`lead`), `follow.cyclicArc` (`wrap`), `navigate.*` (wheel
zoom step, key pan fraction, wheel coalescing). Each a choice with a default; the slice reads them from a plain default
object of this shape, and the `ParamDefinition`s come with the settings form.

The second slice added the groups `selection`, `highlight`, `cursor`, `view` (the history of jumps, its motion and
keys), `address`, `links`, `copyLink`, `revisions`, `edit`, `preview` and `bounds`, and `navigate.doubleClickZoomStep`;
the viewer's `_spec` (Configuration) lists every item with its default.

## Consequences

- The first slice builds `SpecModel`, `ViewModel` with scroll groups and one follow, `HoverModel` with projection,
  `StoreModel`, and the `hover` and `navigate` behaviors; the second slice `SelectionModel` (several ranges, the
  anchor, the keyboard cursor), `HighlightModel`, the `AddressModel` (links, copy link), the history port and
  previews, story steps, followers that lead and the view's history of jumps.
- Tests in Node: viewport rules EC-111 to EC-114, hover ownership EC-141, projection through the codon map, the
  bindings' table resolving gestures by target.
