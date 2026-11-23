# Navigators: one part over an axis or two, its window moving the view

## Intent

The author, 2026-11-18: "have we designed "navigation tracks" for our views?", then "we forgot to even mention the
"vertical maps"? the thing on the right. it should prolly be built with the "navigation tracks"?". Agreed with the author: one part,
set by the axes it covers (positions, rows, or both), in all four forms of EC-110 (the strip, the vertical overview,
the miniature, the ideogram); drawn from the coarse summaries the views already read; its window dragged and resized to
move the view; idle when there is nowhere to go; a floored window drawn wider; a part a host places, like a legend.

A navigator answers two questions at a glance, "where am I in the whole" and "take me there", for any view whose whole
does not fit: a chromosome region, a 10,000-row alignment, a long protein's tracks, a plasmid's linear detail.

## Precedents

| Precedent | Taken | Not taken |
|---|---|---|
| The current Seq\*'s navigation track | one entry made of parts (an axis, the tracks in miniature, ticks of highlights and the selection); the window's parts (body, two end handles as tabs outside the ends, a hit at least 12 px, a floor of 16 px); idle drawn quiet with no relayout; the window filled in a window token, never the selection's ink; the base drawn on its own pass so a pan redraws only the window | the navigator as a registry entry of a preset with twelve switches (here: one part with a form and a few params); a press outside the window jumping, then the drag panning (here: a sweep, Tensions); a floored window placed as a scrollbar thumb (here: centred on the truth, Tensions) |
| The current Seq\*'s ruler rail | the default linear navigator, drawn just below the pinned ruler's bottom border (a 12 px row of its own on the canvas with rails on; never in an export), the window a thicker segment (drag pans, ends zoom); the whole ruler row its grab zone on a positions axis, only the lower 40% on an alignment's columns (above it a drag selects columns); a slim `rail` row without a ruler | — (taken whole: the author's choice); at rest 1 px and muted, the accent only while hovered or dragged (decided, the author, 2026-11-19, the defaults audit: the loudest mark at rest) |
| The current Seq\*'s vertical overview | the rows in miniature down a column, the window the rows and columns in view, top and bottom handles for the rows, left and right for the columns | a render prop of the section (`verticalNav`) beside a scrollbar mode: here it is the same part as the strip, with both axes |
| The current Seq\*'s ideogram | Giemsa stains as grays derived from the theme's ink (here the themes design's `giemsa` scheme), the centromere's two `acen` bands as triangles, a hairline where a block has no bands | bands as application data found by block name (its own TODO wanted them as document data: here they are an annotation) |
| IGV's ideogram and ruler | a red box over the chromosome; a click centres there; a drag on the ideogram or the ruler sweeps a region to zoom to | the ideogram always drawn p-arm left whatever the view's direction (Tensions) |
| UCSC | the chromosome ideogram above the tracks, the region boxed; a drag selects a region to zoom to | a navigation that reloads the page; switching chromosomes from a genome-wide ideogram (refused here, below) |
| Jalview's overview window | the whole alignment, rows and columns, in a small window with the view's box; drag the box; click moves it there; hidden rows and columns shown differently | a window of its own that the application opens (here a part the host places) |
| VS Code's minimap | the content drawn as blocks at a tiny scale; the slider dragged; a click goes there; decorations (search results, errors) as marks at their places; the slider shown on hover only (an option) | characters drawn (at this scale a row is a pixel: a color is what reads) |
| Genome browsers' chromosome navigators (JBrowse's overview scale bar) | several displayed regions as several boxes on one bar; a drag on the bar sweeps a region | polygons joining the bar to the view (no: the window is enough at this size) |
| This codebase: the view's history of jumps, followers that lead ([viewer-models.md](viewer-models.md)) | a click and a sweep are jumps (Alt+Left back); the window is the view's viewport carried to the navigator's axis like a follower's, and a move of it leads, as a follower's does | — |

### Interactions in the precedents

| What a user does | Here |
|---|---|
| click outside the window: the view centres there (all) | built: `navigator.center`, a jump |
| drag on the ruler to pan (the current Seq\*'s rail, genome browsers' rulers) | built: the ruler rail, the default (an alignment's ruler: its lower 40% only) |
| drag the window: the view pans (all) | built, by hand, no animation |
| drag an end: the view zooms, the other end held (current Seq\*) | built |
| drag from outside the window: sweep a region and zoom to it (IGV, UCSC, JBrowse) | built: `navigator.dragOutside: sweep` (the current Seq\*'s jump-then-pan an option) |
| double-click: the whole axis, or a zoom there (current Seq\*) | built: `navigator.doubleClick: toggle-fit` |
| wheel over the navigator (VS Code) | built: as over the view, fixed |
| back after a jump (browsers) | built: Alt+Left, the view's history; a far drag pushes a step (Self-check) |
| marks of the selection, highlights, search results (current Seq\*, VS Code) | built: ticks in the navigator's `selection` and `highlight` passes |
| band names on hover (IGV, the current Seq\*'s TODO) | built: the navigator's own readout, never the viewer's hover |
| several regions, several boxes (IGV multi-locus, JBrowse) | built: a window per bound section, leaders and followers alike |
| hidden rows and columns marked (Jalview) | built as the view shows them: a collapsed clade is its summary row (from the rows design) |
| the slider shown on hover only (VS Code) | refused: a window shown only on hover cannot be read at a glance, and the idle rule already quiets it when it means nothing |
| switching chromosomes on a genome-wide ideogram (UCSC, JBrowse) | built by [genome-coordinates](genome-coordinates.md) (the fifth slice): `navigator.system` an assembly, windows over blocks, `navigator.switch` (EC-184, EC-210) |
| keyboard: arrows move the box (Jalview's arrow keys, scrollbars) | built: the part is focusable, its windows are sliders |

## Decision

### The author's decisions

The author, 2026-11-18: "navigators.md: ok as recommended too"; in general: "a lot of this is "configurable behavior"
and where appropriate, the code should allow the configuration via spec ports etc."; on sizes: "don't overdo it on the
limits (they were sort of self-imposed and estimate), but obviously try to minimize the amount of code". On the review of the decisions: "ok" to its recommendations.

1. Navigators are nodes of the view in the spec, placed by the host (`place` a hint); the default mount is a
   `SeqViewerSpec.ui` entry, the gestures a `BehaviorProvider`.
2. A drag starting outside the window sweeps a new window and zooms to it, splitting [ux-foundations.md](ux-foundations.md) section 8's
   "drag on an axis, navigator … pans" for navigators.
3. A floored window is drawn centred on the true place, not as a scrollbar thumb. The floor is measured within the picture, not the
   reachable span (decided at NV's fix round: chr17's region window was otherwise 0.05 px, unseen).
4. An ideogram on a reversed section is drawn reversed (the drag follows the view).
5. A hover on a navigator stays in the navigator (a readout).
6. `mountNavigator(element, viewer, id)` and `viewer.navigator(def)` are public; no `renderNavigatorSvg` for now.
7. An SVG export leaves navigators out unless `renderSvg`'s `navigators` says `include`; a far drag pushes a back step
   (a default confirmed by use).
8. Hits are the one `window` target with a `part`; the scroll is `bandTop`; the ideogram's grays are the `giemsa`
   scheme; the reduction and image cache are [grids-and-rows.md](grids-and-rows.md)'s.
9. A drag on a ruler, the author: "i think drag on a rule should pan. or rather .. the original solution solved this by
   showing the rail as well ... maybe that's the way to go?" Adopted: the ruler rail is the default navigator of a
   linear section (the slim rail without a ruler); a drag on a positions ruler pans; on an alignment's columns the
   ruler's lower 40% pans and its upper part selects columns; a grid's axes pan. And: "also the rail had "ticks" on
   the edges to resize the window. and this obviously has edge cases to handle too. and then the rail is
   disabled/hidden, etc": the grips, floor, hit, idle, grown, hidden and disabled rules and the axis cases of "The
   ruler rail" (a press outside centres or sweeps, as decided, not the scrollbar's paging).
10. Export choices, the author: "others as recommended": `navigators` is an optional spec param of a view, the
    `renderSvg` option overriding it.

### One part, set by the axes it covers

- **A navigator is a node of the view** (`navigator`, in `views/<view>/navigators`), beside its sections: `of` (the
  section or sections it binds, refs), `axes` (`positions`, `rows` or `both`), `content` (`axis`, `tracks`, `bands`,
  `rows`, or a registered picture `<id>:<name>` falling back to `axis`; `NavigatorModel` the contract of a host-drawn
  navigator, [extension-openings](extension-openings.md)), `tracks` (the tracks drawn in miniature, refs; all by default), `bands` (an annotation, for `content: bands`),
  `system` (the system it draws, when not the section's: a chromosome under a region), `place` (`top`, `bottom`,
  `left`, `right`: a hint to the host and to an export), `size` (its thickness in px), `hidden`. Why a node: a spec file
  carries its navigators to another reader and an export can draw them; why not a section: it has no axis of its own to
  scroll, its window is another section's viewport.
- **The four forms are presets of those params**, a builder helper each:

  | Form | `axes` | `content` | Default size | For |
  |---|---|---|---|---|
  | rail (the default) | `positions` | `axis` as a hairline, the window a thicker segment | a 12 px row just below the pinned ruler's bottom border on the canvas with rails on, none in an export (`ruler` rail); 12 px (`slim` rail, a section with no ruler) | every linear section |
  | strip | `positions` | `axis`, or `tracks` (a navigation track with the section's tracks small) | 20 px (`axis`), 48 (`tracks`) | any long section |
  | ideogram | `positions` | `bands` over a chromosome | 24 px | a genomic section |
  | vertical overview | `both` (the window's rows and columns, four grips; the fourth slice's FX1, as the previous implementation's) | `rows`, columns compressed across | 96 px wide, narrowing to 48 and hidden in a narrow frame | an alignment's rows, the Globins 10,000-row view |
  | miniature | `both` | `tracks` and `rows`, both axes whole | 160 × 96 px | a section with tracks and rows that both scroll |

  Why one part: the window's geometry, gestures, idle rule, floor, keys and announcements are the same along any axis;
  only the picture differs (the author: built "with the navigation tracks"). The rail is EC-110's strip at its
  thinnest, not a fifth kind: the same `window` target and parts (`body` pans, `start` and `end` zoom, `outside`
  centres or sweeps), drawn on the ruler instead of a row of its own.
- **The rail's look**: The look (the author, 2026-11-22, D1 picked over a capsule, a groove and a quiet variant: "i like the simple look of today's one .. it just needs some polish"): the hairline inset by the well's corner radius; the window a heavier segment in a mid ink (`msui-text` at 55%; `navigator.rail.tint: accent` for the accent tint); a round knob with a 1 px ring in the well's ground at each end where the box zooms, growing about its centre (rest 5, hovered 10, dragged 12 px; it may overhang the row); pan-only rails round caps; at the fit the hairline's weight with small knobs. A rail's hover is its own: no card, no data hover, no header actions.
- **The rail is every linear section's navigator by default** (`navigator.form`: `rail`), with no node in the spec
  (a declared navigator replaces it on the sections it binds, as does a host's `viewer.navigator(def)` part; a hidden
  declared navigator brings the rail back): on a section with a ruler track it is drawn just below that ruler's bottom border, in a 12 px row of its own on the
  canvas with rails on (none with the toggle off) and never in an export (a figure has no row; included rails lie on the
  border) (decided, the author, 2026-11-19: "Below the border", so the rail never covers the ruler's ticks or labels;
  "Leave it out of exports", so a figure has no empty strip);
  with none, a `slim` rail row under the section's tracks (decided at NV's review: under, so the tracks keep the
  coordinates they have with no rail; its size fixed, 12 px, the form's). **The rows rail** (decided, the author,
  2026-11-19: "we should also support a simple slider navigator for vertical scrolling too"): a rows layer's rows axis
  gets the same rail, vertical, 12 px wide at the right edge of its band (beside the right gutter's numbers), its
  window the rows in view (top and bottom grips; a drag scrolls, a grip zooms the rows, `]` and `[` as keys), idle when
  every row is shown and the rows cannot zoom in (rows × the largest row scale fit the band: the ruler's rule); the rows' default navigator, as the slim rail is a linear section's. **A declared vertical overview does not replace it**
  (decided by the author, 2026-11-20, "Both, configurable"): the overview is the picture, the rail the precise
  control, both drawn against one window, so a drag on either moves the other in the same frame; the overview sits
  right of the rail (the decision of 2026-11-20), the rail at the band's edge beside the numbers. `navigator.rowsRail`
  (`ViewerConfig`): `always` (the default), `unlessOverview` (the old rule: a declared overview replaces it; an overview
  a narrow frame hides counts as hidden, so the rail comes back), `none`. **The positions rail beside a declared navigator, and its place** (the author, 2026-11-21, a rows view with its vertical overview: "it's weird this view is lacking the 'horizontal rail'. probably should be at the 'navigator' or at the bottom. in general, configurable"): as on the rows axis, a declared navigator no longer removes the positions rail by default: the picture and the precise control, both drawn against one window. `navigator.positionsRail` (`ViewerConfig`): `always` (the default, as `rowsRail`), `unlessNavigator` (the old replacement rule), `none`. **Where it sits**, `navigator.railPlace` (`ViewerConfig`): `ruler` (the default, just below the ruler's bottom border, as built), `bottom` (a 12 px row under the section's last track, the `slim` rail's place), `navigator` (a 12 px row directly under a declared horizontal navigator of that section, else `ruler`); a section with no ruler and `ruler` takes `bottom`. **A list, and a pinned place** (the author, 2026-11-21: "does this support array with multiple placements too?", "maybe even 'static' so it stays in the view always at the bottom?"): `railPlace` takes one place or a list (`['ruler', 'pinned']`), each a rail on the one window (a drag on any moves the others in the same frame), duplicates dropped after the fallbacks; and **`pinned`**: a 12 px row pinned to the bottom edge of the view's visible frame, outside the scrolled rows (as the ruler is pinned at the top), so a tall view keeps its horizontal control in reach at any scroll; it scrolls with nothing and takes no row in an export. One setting for every layout with a horizontal axis (tracks, rows, the blocks of a genome); a ring has no rail. Why `always`: a 2D overview's window is a box to drag, not a precise horizontal control, and the rows axis already keeps both. Why `ruler` stays the default: the rail's place as built, nothing moves for a reader who set nothing. **A narrow frame** (decided by the author, 2026-11-20, "As proposed"; the order in
  [grids-and-rows.md](grids-and-rows.md)): the vertical overview narrows from 96 toward `navigator.overview.minPx` (48)
  as the plot shrinks, and is hidden below `navigator.overview.minPlotPx` (160) of plot width (`RenderConfig` presets;
  0 never hides), keyed on the plot's width, not the frame's; the rail then is all that is left, so the way to move is
  never lost. The overview's picture draws its rows with no row columns, so the hidden-columns chip never lands on it. **Its grab zone**: on a positions axis (blocks, a genome, a protein) the whole ruler row is
  the rail's, so a drag on the ruler pans (ux-foundations section 8's "drag on an axis pans" holds); on an alignment's
  columns only the ruler's lower part is (`navigator.rail.grabShare`, 0.4 of the ruler's height, a `RenderConfig`
  preset), and a drag on the ruler above it selects columns. Why 0.4: the ruler's labels sit in its upper part, where
  a column pick aims, and the window segment sits on its bottom border, where a pan aims; the lower 40% gives the
  segment a grab of at least a third of a 20 px row while leaving the labels' band to columns (the precedent's split).
  A grid's rulers are its axes: a drag on either pans that axis.

### The ruler rail: its window and its edge cases

The author: "also the rail had "ticks" on the edges to resize the window. and this obviously has edge cases to handle
too. and then the rail is disabled/hidden, etc". The rail's window is the shared `window` ([layouts.md](layouts.md)), so every rule
below holds on every navigator form and on the ring's window unless it says rail.

- **Grips resize**: a tick (grip) on each end of the window; dragging one zooms about the other end, which stays put,
  stopped at the zoom limit and the floor, never flipping. The grips are tabs just outside each end
  (`navigator.window.handlePx`, 5) and at most a quarter of the window inside it, so a press in the window's middle
  always pans. Why: on a 2 px rail an end has no room of its own; a tab outside gives it one without stealing the body.
- **Grips light per part**: over an end its grip lights; over the body both do; outside, none. Why: the light says
  what a press will do before it happens.
- **A floor**: the window is drawn at least `navigator.window.minPx` (16 px; at most half the span, so it still
  travels), or the zoom limit's share of the span if that is longer (the view cannot be narrower). Drawn centred on the
  true window (decided), a drag maps by the pointer's travel at the navigator's scale, so it covers the whole range; an
  end dragged inward on a floored window keeps the view (the drawn end is not the true one: shrinking the floor would
  zoom in by surprise), outward zooms out from the drawn end. Why: a window under a few px cannot be seen or grabbed
  (EC-110).
- **The hit**: at least `navigator.window.hitPx` (12 px) wide, centred on the window and inside the span, with touch
  input at least 24 px (`navigator.window.touchHitPx`). Outside the window is the `outside` part: a click centres the
  view there (a jump) and a drag sweeps (decided). Not the scrollbar's paging: a navigator draws where things are, so
  a press there means "go there", not "one page toward there"; the precedent's scrollbar paged, its navigator did not.
- **Idle** (above, EC-110): when the view shows the whole axis and zoom-in is at its limit the window has nowhere to go:
  drawn quiet, no grips, no ticks, no hover, not a target (a press falls through to the ruler: on an alignment's
  columns, a column pick); it wakes as soon as a zoom becomes possible (a resize, a data length arriving), with no
  relayout. Why: a control that does nothing must not look like one.
- **Grown**: the rail is drawn thicker (2 to 6 px, `navigator.rail.grownPx`) while its drag runs, wherever the
  pointer goes, and while the pointer is over its grab zone; the cursor is `grab` there and `grabbing` during the drag,
  by the same rule. Why: a 2 px line is a hint at rest and a target in use.
- **A drag ends wherever the button is released** (the author, 2026-11-19: "when dragging the slider, if you release
  outside the canvas, it doesn't release the event. this needs to be handled by creating a temporary event handler on
  the window/document. other controls prolly have the same issue"): every drag a mount starts (a navigator's window,
  grip or sweep; the view's rails, a drag-select, a height drag) takes pointer capture and also adds temporary
  `pointerup` and `pointercancel` listeners on the element's window (its `ownerDocument.defaultView`, so an iframe's
  own), removed when the drag ends; a `pointermove` with no button held, a `lostpointercapture`, the window's `blur`
  and the document turning hidden end it too. Each path ends the drag once (the input's `up` is idempotent). Why both:
  capture fails silently where it is refused (an iframe boundary, a pointer already released), and a window listener
  alone misses a release outside the browser window.
- **The cursor says the gesture** (the author, 2026-11-19: "when the action is resize on the sliders, it should get an
  appropriate icon"): over a window's grip (`start`, `end`) and while dragging one, a resize cursor along the
  navigator's axis (`ew-resize` on a horizontal one, `ns-resize` on a vertical one: the rows rail, a vertical
  slider); over the window's body `grab`, `grabbing` while dragged; outside it `pointer` (a click centres); an idle
  navigator the default. One rule for the ruler rail, the slim rail, the rows rail, the strip and the miniature.
- **Hidden or disabled**:
  - `navigator.form: none` (`ViewerConfig`, the host's) draws no rail on any section; a declared `navigator` node's
    `hidden` (spec param) hides that one;
  - the embedded layout (the viewer inside a host's cell) shows a button in the canvas corner that hides and shows the
    rails (a `SeqViewerSpec.ui` entry, `railToggle`; its state the reader's session, never the spec);
  - a section with no ruler gets the slim `rail` row (12 px); a ruler that is itself hidden (`track.hidden`) falls back
    to the slim rail, so hiding the ruler never hides the way to move;
  - disabled: where the view cannot be moved (no `navigate` behavior, a full-state story step, an export) the rail is
    drawn with no targets, as idle. Why: one drawing rule for "nothing to do" and "not allowed".
- **Other sections**:
  - a ring has no rail: its window is the ring's own ([layouts.md](layouts.md));
  - a reversed section's rail runs with the axis's display direction (a drag right moves the view right);
  - collapsed introns: the rail spans the drawn pieces (the section's display units), each break a notch;
  - a follower's rail moves its leader by `follow.followerMoves` (`lead`; `ignore` makes it idle; `detach` detaches);
  - several blocks (`arrange: row`): one rail over the blocks laid end to end, a hairline at each block gap;
  - a narrow data area: the rail is drawn in what width there is and still hits at least 12 px; with no data area at
    all, no rail;
  - keyboard: the window is a slider (one tab stop): arrows pan, Shift+arrows resize about the centre, Home and End go
    to the ends, `aria-valuetext` set on settle ([layouts.md](layouts.md), the `window` behavior);
  - touch: the larger hit above; a touch drag on the rail pans, never selects;
  - the grab-zone split on an alignment's columns: `navigator.rail.grabShare` (0.4, above).

- **A host may also make one without the spec** (`viewer.navigator(def)`, a definition of the same shape, not saved):
  a Mol\* plugin adds a strip to a user's spec without patching it. Same model, same tests. A public entry
  (decided).

### The picture: the whole axis, from what is already summarized

- **A navigator lays out its system whole** (the fit window, in the bound section's display order) into its own box,
  through the layout interface: `tracks` content is the tracks layout at the fit window with each track's height scaled
  (`navigator.trackScale`, 1/4, at least `navigator.trackMinPx`, 2 px), whatever layout the section uses (text, blocks: the miniature is
  linear). Why the layout interface: a miniature drawn by other code drifts from the view it summarizes.
- **It reads what the fit view reads**: its own want set (`store.wanter()`, EC-125) at the level its width needs, at
  `view` priority within that set; these are the tiles a fitted view reads, shared in the store. Why `view`, not
  `margin`: the store evicts wanted `margin` tiles first and an export drops `margin` wants (`render/export.ts`), so a
  navigator's tiles would thrash and an exported one would never wait; they are few (a tile or two a layer at the fit
  level, bounded by `maxTilesPerWant`), so they cost the view's slots little. An export that includes navigators
  waits for their wants as for its view's. A
  values layer reads its pyramid ([data-store.md](data-store.md), Summaries); a letters layer draws at its density step (dots or bars per
  run); a ranges layer its coarse level when the source serves one, else a level-1 read within `maxTilesPerWant`, else
  it is left out with `seq.render.navigator.layer` (info, once). Why: a navigator must never cost what a view zoomed
  out costs twice, nor ask for a chromosome's features at level 1.
- **Colors are the view's**: a miniature colors by the view's resolved domains (the themes design), never by a domain
  of its own over the whole axis, so a value has one color in both.
- **Pending is drawn pending** (EC-169), as in the view. A bound section with no viewport yet (no width, or an axis
  whose length is pending: EC-112) draws no window, its sliders `aria-disabled`; the picture draws the axis pending
  until its length is known. A window appears when the viewport is placed, with no relayout of the base.
- **The rows content** reads what the rows layer gives ([grids-and-rows.md](grids-and-rows.md), "What a navigator reads"): the rows in
  display order, the rows band's scroll (`bandTop`, the layouts design's viewport record) and visible
  height, and its row summaries: a coarse grid of row bins × column bins, one color a bin
  (an `image` primitive, reduced by the shared 2D reduction and cached in the shared image LRU of [grids-and-rows.md](grids-and-rows.md));
  falling back to one value a row (identity to the reference, else the non-gap share), then to the rows'
  structure only (clade bands, pinned rows). It shows the rows as the view orders them: a collapsed clade is its summary
  row, a filtered row is absent. Why: the window maps the view's rows; a navigator with other rows would put it in the
  wrong place.
- **The ideogram** is a `bands` layer (a new layer type, fallback `ranges`) over an annotation with a `stain` field
  (UCSC's cytoBandIdeo: `gneg`, `gpos25`…`gpos100`, `acen`, `gvar`, `stalk`), colored by its default theme, the
  `giemsa` scheme (the themes design: grays derived from `msui-text`, so dark flips them; `acen` in `seq-centromere`),
  the two `acen` bands as triangles meeting at the centromere, a hairline where no band is.
  Why a layer type: the bands are data like any annotation, and the same layer may sit on a track of a figure.
- **Marks over the picture**: the selection and the highlight channels as ticks at their places (`navigator.marks`),
  the keyboard cursor as a hairline; each at least 1 device px. Why: "where is what I picked" is a navigator's question
  too (VS Code's decorations).
- **Draw lists** are the navigator's own, in the usual passes plus one, `window`, the first emphasis pass, under the highlights (the layouts design's order;
  the window, its wash outside, its handles), so the selection's and highlights' ticks stay readable over it; the base (`background`, `data`) is keyed by the navigator's axis, size, tokens and the store
  versions it read, so a pan redraws only the overlay passes. A navigator never draws in a section's lists.

### The window: the view's viewport, carried to the navigator's axis

- **A window per bound section**: the section's viewport on its axis (and its rows band's `bandTop` and `rowScale`), carried to the
  navigator's system along the path between them (identity when the navigator draws the section's own system; a region
  into its chromosome; a protein into its gene). A scroll group shares one range, so its sections draw one window; a
  follower's window is its leader's carried along the follow, so a strip may draw a leader's window and its followers'
  (one leading section moves them all), each labelled by its section. **Several leading sections on one navigator**
  (`of` a list naming sections that lead independently: IGV's multi-locus) draw one window each, each labelled by its
  section, each drag moving its own section; the renderer's frame carries a viewport per leading section
  (`FrameInput.viewports`, built by the layouts design, before the navigators). Overlapping windows: the one whose
  section the pointer last acted on is on top and takes the press. Why a path: an ideogram of chromosome 17 over a section on the gene's region, and a protein's strip showing
  where its gene's follower is, are the same rule.
- **The reachable span**: where the bound section's axis covers only part of the navigator's system (the region within
  its chromosome), the rest is washed (`seq-wash`), fixed, so a window cannot be dragged where the view cannot go. Why
  fixed: no second use for a window dragged where nothing can follow.
- **The zoom strip** (the author, 2026-11-19: "the window on the gene strip doesn't work… isn't synced", then "Yes,
  the zoom strip"; UCSC's and IGV's): where the reachable span is narrower than `navigator.zoomStrip.minPx` (48 CSS px)
  on a navigator, a window there could neither show a move nor take a drag (p53's 6 kb region on chr17's 83 Mb is
  1/13,800 of the strip; a pixel there is 78 kb). Then the navigator is **context**: the reachable span drawn as a
  mark (an outline at least `navigator.window.minPx` wide, centred, no grips, the default cursor, not a target for a
  drag; a click on it, as a hover, names the region), and a **zoom strip** under it, the same navigator kind's form
  over the reachable span alone (an ideogram's bands there, a strip's axis), with the real window, its grips and
  gestures. Two connector lines run from the mark's ends to the zoom strip's ends (`seq-window`, hairlines). The strip
  takes `navigator.zoomStrip.heightPx` (16) more; the pair is one navigator (one mount, one hit target, one record).
  Why: the chromosome's context and a working window are two scales; one strip cannot be both.
- **Moving a window moves its section**, through the inverse of the path, as a follower leads ([viewer-models.md](viewer-models.md)): a
  window bound to a follower moves its leader (`follow.followerMoves`: `lead`; `ignore` makes that window idle; `detach`
  detaches the follower). Why: one rule for "a move asked of something that follows".
- **Parts**: the body (a drag pans), the two ends on each covered axis (a drag zooms that axis, the other end held,
  stopped at the zoom limit and at the floor, never flipping), outside (a click centres there; a drag sweeps a new window,
  `navigator.dragOutside`). Handles are tabs just outside the ends and at most a quarter of the window inside it
  (`navigator.window.handlePx`, 5), so a press in the middle always pans; the window is hit over at least
  `navigator.window.hitPx` (12). The rows axis zooms only where the rows layer has a row zoom; else its ends only pan.
- **A floored window is drawn wider** (EC-110): at least `navigator.window.minPx` (16 CSS px, at most half the span),
  centred on the true window and clamped inside the span, with a hairline at the true extent when that is at least a
  device pixel (`navigator.window.floor: centred`; `thumb`, the scrollbar's mapping, the alternative). A drag moves the
  view by the pointer's travel at the navigator's scale, so it reaches either end whatever the drawn width. Why centred:
  on an ideogram or a miniature the window says where; a thumb's position is off by up to half its width (a band).
- **A cyclic axis** (the plasmid's strip): a window over the origin is drawn in two pieces and a drag carries it across
  (follow.cyclicArc's `wrap`; EC-9, EC-12).
- **Display order**: the navigator draws in its bound section's display order (a reversed section's strip runs right to
  left in positions, marked as the axis design marks a reversed axis), so a drag right moves the view right. Collapsed
  introns collapse in a strip over the section's own axis; an ideogram's system is the chromosome, uncollapsed, so the
  window is the hull of the view's image there.

### Gestures, by the grammar's table

The navigator's hits are the one `window` target the layouts design defines for the ring's window, its `part`
(`body`, `start`, `end`, or `outside` on a navigator), with the navigator's id, the position and row under the point at
the navigator's scale, and the band. Its gestures are bindings like any other ([ux-foundations.md](ux-foundations.md), section 8), run by the
`navigator` behavior (a `BehaviorProvider` in the viewer's spec, extending the layouts design's `window` behavior), so a
preset or a host rebinds or leaves them out.

| Gesture on a navigator | Effect | Command or method |
|---|---|---|
| click outside the window | the view centres there, animated (section 7), reduced motion a cut | `navigator.center` (navigate, a jump) |
| click on the window | nothing (it is a press for a drag) | — |
| drag the body | pans with the hand, no animation | `navigator-drag` (a method: per move, allocation-free) |
| drag an end | zooms that axis, the other end held | `navigator-drag` |
| drag from outside | sweeps a new window; on release the view zooms to it | `navigator-drag`, then `navigator.zoom-to` (a jump) |
| double-click | not fitted: fit; fitted: zoom so the window is `navigator.zoomPx` (48) navigator px wide, centred there | `view.fit`, `navigator.zoom-to` (one jump with its first click's: below) |
| wheel, Cmd/Ctrl+wheel, pinch | as over the bound section (the same deltas in that section's px; one leading section, above), zoom about the window's centre, not the pointer: the position under the pointer at the navigator's scale may be far outside the view, and a zoom about it would leap. Fixed: no second use for a navigator without the wheel | the view's `pan`, `zoom` |
| Esc during a drag or sweep | cancels it, the view back where the press found it | — |
| hover | a hairline and a readout in the navigator (the position, the row, an ideogram's band name); the viewer's hover is not set | — |
| a rows window: drag its top or bottom grip | zooms the rows about the other end; never scrolls (the author, 2026-11-20) | the rows' zoom (`rowScale`) |
| a rows window: double-click (the rows rail's, the vertical overview's) | resets the rows' zoom (`rowScale` to its default, the top row kept), one jump; on the overview with the rows not zoomed, the positions' double-click (the author, 2026-11-20: "double click on the window (and perhaps the 'thin slider') should prolly reset the vertical zoom too") | the rows' zoom, a jump |

- **History**: a click, a sweep and a double-click are jumps (Alt+Left goes back); a step a navigator pushes is the
  view's, under its one rule (a step over a changed axis or layout restored by position, clamped to what is drawn,
  never silently dropped). One gesture is one step:
  - a double-click's first click defers its step the way input's first click defers a clear (`_cleared`): the click
    centres at once without a step and keeps the viewports it left; a double-click in time pushes them as its own
    jump's step; no second click in time, they are pushed then;
  - a drag is by hand and pushes nothing, unless on release the view moved farther than one window width from the press
    (`navigator.dragStep: far`, a default to be confirmed by use): the navigator keeps the viewports at the press
    (`ViewModel.snapshot()`, every group's range, as `jump` reads them) and pushes them on release
    (`ViewModel.record(snapshot)`, new: `jump` records only before a move it makes itself). The deferred click uses the
    same two methods. Why: back must undo a trip across a chromosome,
  not a nudge.
- **Idle** (EC-110): an axis of a window is idle when its section can neither pan nor zoom along it (the whole axis
  shown and zoom-in at its limit; all rows shown with no row zoom; a follower window under `ignore`). Idle is drawn quiet
  (no handles, no wash, the window a quiet line), its gestures do nothing, its cursor is the default and its slider is
  `aria-disabled`; it is live again when the view can move, with no relayout. A navigator whose every window is idle is
  still drawn, quiet (fixed, decided by the author, 2026-11-20; `hidden` in `packages/seq/TODO.md` until a host asks): a part that comes and goes
  moves the host's layout.
- **Keys**, while the part has focus (one tab stop; a slider per covered axis and window, roving): the slider's own
  keys, fixed (ARIA's slider pattern): arrows pan by `navigate.keyPanFraction` (Up and Down on rows), Shift+arrows
  resize the window about its centre (a zoom step), Home and End go to the ends (jumps), Enter centres on the readout's position (a jump); the view's keys from
  `view.keys` (`=`, `-` about the window's centre, `0`, Alt+Left and Alt+Right) act on the focused window's section. No
  key of its own outside the part. Platform-aware, never Tab or the browser's zoom keys (EC-155).

### A part a host places

- **The viewer gives parts, not frames**: `mountNavigator(element, viewer, id)` draws one navigator into an element the
  host owns (its own base and overlay canvases, sized by the element, its listeners on it alone; it returns its
  disposer). The default mount, a `SeqViewerSpec.ui` entry (`ui.navigators`, as `ui.legends` is the legend's), places
  the view's declared navigators by their `place` through the same call, as a host would. Why: Seq\* inside a bigger app (a Mol\* plugin, a notebook) must not own the page's layout; the default is a
  convenience over the public call, not a privilege.
- **Headless first**: a `NavigatorModel` (DOM-free, tested in Node) holds the bound windows' frames, the lit part, idle
  per axis, the drag in progress and the readout; the canvas view and the SVG writer read it.

### Light and dark, export, accessibility

- **Tokens** (the themes design's list): `seq-window` (the window's fill and edge), `seq-grip`, `seq-wash` (outside
  and unreachable), `seq-centromere` (a hue, 3:1 against the surface in both themes); the Giemsa grays are the `giemsa`
  scheme's, derived from `msui-text`. The
  window is chrome, never the accent: the accent is what is picked ([ux-foundations.md](ux-foundations.md), section 1), and the selection's
  ticks use it.
- **SVG export**: `renderSvg` leaves navigators out (its `navigators` option: `omit`), settled by EC-138's rule:
  chrome state is not a figure, and the window is chrome state; with `include` it stacks the view's
  navigators by `place` around the view, each with its picture and its windows at the export's viewports as outlines
  and wash, no handles, no hover, no idle state (chrome state is not a figure; EC-138's rule for hover). Each navigator
  is a group with a `<title>` ("Window: chr17:7,668,402–7,687,550 of chr17") and `data-navigator` its id; the export
  waits for the included navigators' wants. A navigator alone (`renderNavigatorSvg`): not now (decided).
- **Announced**: each window slider has `aria-valuetext` in the section's labels ("positions 1,200–1,450 of 2,686";
  "rows 401–440 of 10,000"), set when the viewport settles (`committed`), never per frame; a screen reader reads it on
  focus and on each key. A move made elsewhere announces nothing (noise). Forced colors: the window and handles in
  `CanvasText`, no wash (the accessibility work's rule).

### The window's look: one control row, one stroke, one readout place (amended 2026-11-22, the author: "Yes, as recommended")

Lane DL, the fifth slice, 2026-11-22. A look design before a fix lane. It amends [navigators.md](navigators.md)
("The window", "The zoom strip", "Light and dark", the matrix and Configuration) and
[genome-coordinates.md](genome-coordinates.md) ("The zoom strip over several systems"); F25 (the
fourth slice's novice finding) is decided here with it.

#### Intent

The author, 2026-11-21, on a ring that switches chromosomes: "i don't understand why chr6 gets a 'box' and the other
'sub-section'"; "the 'subsection visual' is really ugly, with different line thicknesses". And 2026-11-22, the BED
cytobands story: "it's better if the label stays outside in top right instead of jumping around all the time. maybe we
can do the same thing for the 'vertical maps' as well, where it appears 'above' the map instead."

What a reader must be able to say after one look: where the control is, which mark is "where I am" and which is "what
I drag", and where to read the numbers. Today none of the three is stable.

#### What is drawn today, and why it reads as two designs

- **The strip comes and goes with the window's width.** Over blocks, `zoomFlips` takes the window's own width
  (`core/src/render/navigator.ts:657`, `px = this.overBlocks ? this.windowPx : …`, `windowPx` set at `:967`) and turns
  the zoom strip on under `navigator.zoomStrip.minPx` (48, `window.ts:28`), off past it plus 8 px of hysteresis
  (`navigator.ts:94`, `:659`). chr6 whole is about 85 px on the genome row (a box with grips); chr13 whole is under 48
  (a strip). The trigger is a pixel count no reader can see. On one system it reads the reachable span instead (a fixed
  property of the section), which is the rule navigators.md states.
- **Four stroke weights.** The window's outline 1 CSS px, 2 when lit (`window.ts:247`); the context mark an outline of
  1 CSS px with no fill (`navigator.ts:1142–1168`, `drawNavigatorWindow` with `FlagOutline`); the connectors one
  device pixel (`g[6] = 1 / dpr`, `navigator.ts:1126`), so half the window's line at DPR 2; the idle outline 1 px at
  0.4 (`window.ts:36`). The lit 2 px also makes the box jump a pixel under the pointer.
- **The mark and the window share an ink but not a role**: the mark is drawn as a window without its fill
  (`seq-window`, full strength), so the context's "where" mark is as loud as the control.
- **The readout moves** (`viewer/src/readout-place.ts:82–122`, `placeReadout`): top right, else top left, else beside
  the windows, else above or below the part, chosen per frame from the window boxes and the pointer
  (`viewer/src/navigator.ts:376–378`, `mount-navigator.ts:89–95`). It jumps as the window or pointer crosses it.

#### Precedents

| Precedent | Taken | Not taken |
|---|---|---|
| IGV: the genome or chromosome ideogram with its red box, the locus box above the tracks | the context row as a quiet "where" mark; the readout in one place, read at rest | a red, loud box (here the window is chrome, never louder than the data) |
| UCSC: the ideogram above the browser image, the position box at the top, unmoving | the numbers in a fixed line above the picture | — |
| JBrowse 2: the overview scale bar with the region highlighted and polygon lines to the view | the funnel as a thin, muted tie between two scales | its filled polygon (it covers the ideogram's bands) |
| This codebase: the zoom strip (navigators.md), the floored window, idle drawn quiet (EC-110) | the strip's parts, gestures and record; the floor | the window-width trigger over blocks |

#### Decision

##### 1. One rule for the second row: by what the section reaches, never by the window's width

- **A navigator has a control row and, sometimes, a context row above it.** The control row holds the window with its
  grips; the context row holds only a quiet mark of where the control row lies in the whole. A reader learns one thing:
  *the grips are always on the bottom row*.
- **When there are two rows** (one rule, readable from the spec, never from a pixel count):
  1. **Over blocks** (`navigator.system` an assembly): **always two rows**. The context row is the genome (each
     chromosome's bands, click to jump, sweep across blocks: unchanged); the control row is the hull of the blocks the
     window touches (one chromosome when the view is within one; the genome's ruler at the fit). So chr6 and chr13
     look the same: both the chromosome below, a quiet mark above.
  2. **One system whose reachable span is narrower than `navigator.zoomStrip.minPx`** (a region on its chromosome, p53's
     6 kb on chr17): two rows, as built. The span is the section's, so this never flips with a zoom.
  3. **Otherwise one row**: the window on the whole axis, its grips there; a window narrower than the floor is floored
     (16 px, centred, its hairline at the truth), never handed to a strip.
- **The height never changes with a gesture**: in cases 1 and 2 the part's height is the context plus the strip, fixed
  for the navigator's life (no hysteresis needed: `ZoomHysteresisPx` goes). Why always over blocks: the genome's and the
  chromosome's scales are two (genome-coordinates.md's reason), and a toggle by width is a rule the reader cannot see.
- **The control row at the fit is the genome's ruler** (case 1 with the whole genome in view): Mb ticks under the
  ideogram, not a blank well (the empty-well rule's "a strip with no content draws its system's ticks"), its window at
  the fit (2. below).

##### 2. A window at the fit is one quiet frame (F25)

- **At the fit** (the window covers the whole control row, whether or not zoom-in is possible) the window is drawn as
  the idle frame: an outline in `seq-window` at the quiet alpha, no fill, no grips, no wash. It is not a drag target;
  **the whole row is a sweep**: a drag draws the region to zoom to (`navigator.zoom-to`), a double-click zooms there
  (`navigator.zoomPx`), the cursor `pointer` (the `outside` part's). Why: at the fit a box with grips at the row's own
  edges says "drag me" and does nothing visible (F25: "nothing"); a sweep is the one gesture that moves the view from
  there.
- **Idle** (cannot move at all) stays as built: the same quiet frame, no sweep, the default cursor.

##### 3. One stroke, two inks, one fill

| Mark | Stroke | Ink | Fill | Grips |
|---|---|---|---|---|
| the window (control row), at rest | `navigator.window.strokePx` (1 CSS px), snapped to the device grid | `seq-window` | `seq-window` at 0.12 | `seq-grip` tabs |
| the window, lit (hover, drag) | the same 1 px (no 2 px step) | `seq-window` | `seq-window` at 0.2 | the lit grip(s) in `seq-window`, one px wider (as built) |
| the window at the fit, idle | 1 px | `seq-window` at 0.4 (quiet) | none | none |
| the context mark (where the control row lies) | 1 px | `seq-window` at 0.4 (quiet) | none | none |
| the funnel (two lines, mark's ends to the control row's ends) | 1 px | `seq-window` at 0.4 (quiet) | — | — |
| the wash outside the reach | — | — | `seq-wash` | — |

- **One stroke weight** for every line of the part: the window, the mark, the funnel, an idle frame, a floored window's
  hairline at the truth stays a device pixel (it marks a sub-pixel truth, not an outline). Lit is said by the fill and
  the grips, never by a thicker line (the box no longer jumps a pixel).
- **Two inks from one token**: the control at full `seq-window`, every "where" mark at the quiet alpha (`IdleAlpha`,
  0.4, made a `RenderConfig` value). So the loudest mark of the part at rest is the data (bands, a miniature's marks),
  then the window, then the context: the defaults audit's "the loudest mark at rest", kept.
- **The funnel kept, thin and muted** (it ties a mark 16 px wide to a row 600 px wide; without it the two rows read as
  two navigators), from the mark's lower corners to the control row's upper corners, never filled.
- **The control row has no frame of its own**: its edges are the picture's; the gap (`ZoomGapPx`, 6) holds the funnel.

##### 4. The readout in one place, outside the part

- **Above a horizontal navigator, right-aligned at its top-right corner**, in a line of its own
  (`navigator.readout.linePx`, 16, the `msui-text-xs` leading). **Above a vertical navigator** (the vertical overview),
  right-aligned to its right edge, in the line above its picture. It never moves with the window or the pointer, so it
  never covers the control in use.
- **The line is the part's own**: the mount lays it out inside the part's box, the canvas under it, so a story's or a
  page's top edge cannot clip it and nothing needs room above the part. A toolbar line (the viewer's, when the host
  places the navigator under it) is not borrowed: two owners of one line would fight on width.
- **What it says**: at rest, the window's region (`chr17:20,145,469–27,700,018 of chr17`; several leading windows:
  the one last acted on); while hovering, what is under the pointer (position, band); while a gesture runs, the
  window's region as it moves, every frame. The value never truncates; the system's name gives way first (as built).
  Always shown at rest, so the line is never a blank band and the part never changes height (Gaps, 1).
- **A vertical overview is narrow (96 px, 48 in a narrow frame)**: its readout takes the short form, the rows only
  (`401–440 / 10,000`), right-aligned; past the part's width it may extend left over the label-free gap above the rows
  band, never over data (Gaps, 3). The rows rail and the ruler rail drawn inside the view have no readout (as built).
- **A ring**: a ring's window is not a navigator part and has no readout (the hover card names a position; the hole
  holds the title). A navigator part bound to a ring section (a strip over the plasmid's system; Rondo's square beside a
  genome ring) follows the horizontal rule: above, right-aligned.
- **Removed**: `placeReadout`'s candidate search, `windowBoxes`, `pointerBox`, `GripReach`, `PointerReach`
  (`readout-place.ts`, its test `viewer/src/_test/readout-place.test.ts` replaced by one on the fixed place), the
  mount's per-frame `readoutAt` (`mount-navigator.ts:91–95`) and its "no `overflow: hidden`" exception (`:36`): the
  readout is inside the part now. FX5's "beside the window, else above or below" goes with them.
- **A bug for the fix lane, not a design question**: the readout lags a gesture (the window at 21 Mb while it reads
  `1–83,257,441`). The text is set only from a hover's hit: a move with a press held returns before `say`
  (`viewer/src/navigator.ts:468–472`), and wheel, keys and moves from elsewhere never call it. Fix: while a window moves
  (any source), the readout reads the moving window's text from the frame each draw; on release the hover's text.

#### Defaults: what a reader sees

Over blocks, a chromosome in view (chr6 or chr13 alike):

```
                                              chr13:1–114,364,328 of chr13   <- readout line, fixed
 ▐▌▐▐▌▌▐▌ ▐▌▐▐▌▐ ▐▌[▒]▐▌▐ ▐▌▐▐ ▐▌ … ▐▌▐▐▌         <- genome: quiet mark (1 px, 0.4), no grips
                   ╱  ╲                           <- funnel: 1 px, 0.4
 ╱                       ╲
▕▐▌▐▐▌▌▐▌▐▐▌▐▌▌▐▐▌▐▌▐▐▌▌▐▌▐▐▌▐▌▌▐▌▐▐▌▐▌▐▐▌▐▌▌▐▏   <- chr13 whole: at its fit, one quiet frame
```

Zoomed into chr13 (the control row's window has its grips; the genome mark floored to 16 px):

```
                                        chr13:48,300,000–49,100,000 of chr13
 ▐▌▐▐▌▌▐▌ ▐▌▐▐▌▐ ▐▌[▒]▐▌▐ ▐▌▐▐ ▐▌ … ▐▌▐▐▌
 ╱                 ╱  ╲                       ╲
▐▌▐▐▌▌▐▌▐▐▌▐▌▌▐▐▌▐▌▐▐▌▐▐ ▌▐[███]▌▐▐▌▐▌▐▐▌▐▌▌▐▌▐▐   <- window: 1 px seq-window, 0.12 fill, ▐ ▌ grips
░░░░░░░░░░░░░░░░░░░░░░░░░░░    ░░░░░░░░░░░░░░░░   (wash outside the window)
```

One system, the whole axis reachable (a protein's strip): one row, a box with grips; zoomed in deep, a floored box:

```
                                                    p53:120–180 of 393
▕───────────────[▒▒▒]──────────────────────────────▏  floored 16 px, a hairline at the true 60 residues
```

A cyclic axis, the window over the origin (the plasmid's strip, pUC19c, 2,600–60 of 2,686): two pieces, a grip at
each true end only (`noStart`, `noEnd`: as built), the wash between them:

```
                                              pUC19c:2,601–60 of 2,686 (across the origin)
▕██▌░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░▐██▏
  ^ end grip                          start grip ^
```

The vertical overview: the readout above it, short form, right-aligned:

```
            401–440 / 10,000
           ┌────────────┐
           │▓▓▓▓▓▓▓▓▓▓▓▓│
           │▓▓[      ]▓▓│  <- window: rows and columns, four grips
           │▓▓▓▓▓▓▓▓▓▓▓▓│
```

#### The interaction matrix's cells

Changed cells only; the rest stand as navigators.md gives them.

| Row | Cell | Was | Now |
|---|---|---|---|
| a navigator | zoom in / out | drag a grip; sweep from outside | the same; **at the fit, a drag anywhere on the control row sweeps** (no grips there) |
| a navigator | pan | drag the body | the same; **at the fit: none** (the frame is not a target) |
| a navigator | focus | click outside centres; double-click | the same; at the fit, double-click zooms to `zoomPx` there |
| a navigator over blocks | focus | click on a chromosome jumps (context row when the strip shows, else the row) | click on the context row's chromosome jumps, **always the context row** |
| a window | zoom, pan | drag `start`, `end`, `body` | the same off the fit; at the fit quiet, not a target |
| the readout | (no row) | — | not a target (`pointer-events: none`, as built); Enter at the readout (navigators.md's focus cell) still centres on the readout's position |

Negations checked: a window narrower than the floor (floored, no strip: case 3); a window across two blocks (the
control row their hull); the fit on a single system with nowhere to zoom (idle, no sweep); a reversed section (the
control row in display order, as built); a follower's window (moves its leader, as built).

#### Configuration

| Item | Whose | Default | Why | Alternatives |
|---|---|---|---|---|
| `navigator.zoomStrip.minPx` | `ViewerConfig` | 48 | **now read on the reachable span only**, never the window: a span too small to drag gets the strip | any px |
| `navigator.zoomStrip.overBlocks` (new) | `ViewerConfig` | `always` | the visible rule: over a genome the chromosome is always the control row | `narrow` (today's: by the window's width, with hysteresis), `never` (one row, floored window) |
| `navigator.zoomStrip.heightPx` | `ViewerConfig` | 16 (over blocks at least the ruler's height, as built) | — | px |
| `navigator.zoomStrip.funnel` (new) | `RenderConfig` preset | `muted` | ties two scales; quiet | `none` |
| `navigator.window.strokePx` (new) | `RenderConfig` preset | 1 | one weight for every line of the part | px |
| `navigator.window.quietAlpha` (new; `IdleAlpha` today, fixed) | `RenderConfig` preset | 0.4 | the context mark, the funnel, an idle or fitted frame | 0–1 |
| `navigator.readout.place` (new) | `ViewerConfig` | `above` | one place, never over the control (the author) | `none` (a host draws `readoutParts` itself) |
| `navigator.readout.linePx` (new) | `RenderConfig` preset | 16 | the `msui-text-xs` leading | px |
| `navigator.readout.rest` (new) | `ViewerConfig` | `window` | where am I, read without a hover | `none` (shown only while hovering or moving; the line stays reserved) |

#### Tensions

- **Two rows over a genome at all times cost 16–22 px** that the old rule saved on a large chromosome. Taken: a fixed
  height and one place for the grips beat the saving; `overBlocks: narrow` keeps the old for a host short of room.
- **The fit's sweep** removes "drag the box" at the fit; a reader who drags at the fit gets a zoom, not a pan (there is
  nothing to pan to). Taken.
- **A readout line always reserved** costs 16 px a navigator even with `rest: none`. Taken: a line that comes and goes
  moves everything below it.

#### Self-check

(a) The trigger is now visible: over a genome, always; on one system, by the section's reach. (b) Inverse: leaving a
chromosome for the whole genome keeps both rows (the control row becomes the genome's ruler at its fit). (c) Symmetric:
the vertical overview's readout above it as the horizontal's is above. (d) Plural: several windows on one control row
(multi-locus) each drawn by (3), the readout the last acted on. (e) Cyclic: two pieces, grips at true ends. (f) Export:
navigators omitted by default (EC-138, unchanged); with `include` the window an outline at `strokePx`, no readout line
(chrome state). (g) Forced colors: `CanvasText` for every line at one weight, no wash, no alpha (as built). (h) DPR 1
and 2: one CSS px snapped, so the window and the funnel match at both.

#### Consequences

- `core/src/render/navigator.ts`: `zoomFlips` keyed on `overBlocks` and the reach only (`:646–670`), the hysteresis
  removed; `drawContext` and `drawMark` at the quiet alpha and `strokePx` (`:1096–1168`); the fit drawn as a quiet frame
  with the sweep as its part (`:1071`'s flags gain a `fit` case).
- `core/src/render/window.ts`: lit by fill and grips only (`:247`); `IdleAlpha` to config.
- `viewer/src/readout-place.ts` removed; `viewer/src/navigator.ts` readout from the frame while a window moves;
  `mount-navigator.ts` lays out the line above the canvas.
- navigators.md: "The zoom strip" and its Defaults picture, "Light and dark" (the stroke table), the matrix and
  Configuration rows above; genome-coordinates.md: "The zoom strip over several systems" (always), its Defaults 8.
- Tests: a twin test (canvas, SVG with `include`) that every line primitive of a navigator has one width; a frame test
  over chr6 and chr13 (two rows, grips on the bottom row only, equal part height); F25: at the fit no grip primitive and
  a drag sweeps; the readout's box fixed across a drag and its text equal to the window's every frame.

#### Decided by the author

2026-11-22, every gap as recommended (the author: "Yes, as recommended").

1. Show the window's region in the readout at rest (not only on hover)? **decided yes**: the line is reserved anyway,
   and "where am I" without a hover is IGV's and UCSC's locus box.
2. Always two rows over a genome (chr6 and chr13 alike), rather than a strip only for small windows? **decided yes**
   (`overBlocks: always`); `narrow` stays for a host short of height.
3. On the narrow vertical overview, a short readout (`401–440 / 10,000`) allowed to run left past the part over the
   gap above the rows band? **decided yes**; the fix lane measures that the gap is data-free in the rows stories.
4. Keep the funnel (thin, 0.4) rather than drop it? **decided: keep**: without it a 16 px mark and a full-width row read
   as two navigators; `funnel: none` for a host.
5. At the fit, make the whole control row a sweep (no grips, drag zooms)? **decided yes** (F25).
6. Drop the 2 px lit stroke (lit by fill and grips only)? **decided yes**: one weight, no pixel jump under the pointer.

## The interaction matrix's cells

The row "a navigator" ([ux-foundations.md](ux-foundations.md), section 8), filled:

| select | extend | add / remove | focus | back | zoom in / out | pan | hide | edit |
|---|---|---|---|---|---|---|---|---|
| no: at a navigator's scale a pixel is many positions; selecting is the view's (on an alignment's ruler, a drag above the rail's grab zone selects columns: the ruler's, not the rail's) | no, as select | no, as select | click outside the window (`window`, part `outside`) centres (`navigator.center`); double-click fits or zooms there; Enter at the readout | Alt+Left, Alt+Right, the mouse's buttons (`view.back`, `view.forward`); a far drag pushes a step | drag a grip (`start`, `end`: the other end held; a floored window's end dragged in keeps the view); sweep from outside (`navigator.zoom-to`); Cmd/Ctrl+wheel, pinch; Shift+arrows, `=`, `-`, `0` while focused; idle: none | drag the `body` part (the rail grown, `grabbing`); wheel; arrows, Home, End while focused; idle: none | `navigator.form: none` (host), a node's `hidden` by patch (one undo step), the embedded corner button (reader); a hidden ruler falls back to the slim rail; a host unmounts its part | its params by patch (form, tracks, size); no gesture |

Cells of other rows it touches: **a follower section**, pan and zoom: its navigator window leads (`lead`); **a leading
section**, pan and zoom: through any window bound to it; **the view**, back: a navigator's jumps are its steps;
**a track or layer**: a `bands` layer is a layer like any (hide, color).

## Configuration

Whose, by the rule the axis design states (spec param, `RenderConfig` preset, `SpecBounds`, `ViewerConfig`, a port, a
registry entry, an option of a call, or fixed with its reason).

| Item | Whose | Default | Why | Alternatives |
|---|---|---|---|---|
| a navigator's `of`, `axes`, `content`, `tracks`, `bands`, `system`, `place`, `size`, `hidden` | spec param (the `navigator` node) | the form's builder defaults (Decision) | what the figure offers its reader | — |
| `navigator.window.minPx`, `.hitPx`, `.handlePx` | `RenderConfig` preset | 16, 12, 5 | a window under a few px cannot be seen or grabbed (EC-110); a press in the middle still pans | any px |
| `navigator.window.floor` | `RenderConfig` preset | `centred` | the drawn window says where (Tensions) | `thumb` |
| `navigator.trackScale`, `.trackMinPx` | `RenderConfig` preset | 0.25, 2 | a miniature reads as the view; a thinner track is not seen at dpr 1 | any |
| `navigator.marks` | `RenderConfig` preset | selection, highlights, cursor | where is what I picked | any subset |
| `navigator.rail.grownPx`, `navigator.window.touchHitPx` | `RenderConfig` preset | 6, 24 | a target in use; a finger is wider than a pointer | px |
| the embedded corner button | `SeqViewerSpec.ui` entry (`railToggle`), its state the reader's | on in the embedded layout | a cell's room is the reader's | off |
| `navigator.form` | `ViewerConfig` | `rail` | a navigator on every linear section at no row's cost (the author's choice, the precedent) | `none` |
| `navigator.positionsRail`, `navigator.railPlace` | `ViewerConfig` | `always`, `ruler` | the precise horizontal control beside any picture (the author, 2026-11-21); the place as built | `unlessNavigator`, `none`; `bottom`, `navigator`, `pinned`, or a list of them |
| `navigator.rail.grabShare` | `RenderConfig` preset | 0.4 | an alignment ruler's labels select columns, its bottom pans (above) | any share |
| `navigator.dragOutside` | `ViewerConfig` | `sweep` | IGV, UCSC and JBrowse (Tensions) | `pan` (jump there, then pan: the current Seq\*) |
| `navigator.dragStep` | `ViewerConfig` | `far` | back undoes a trip, not a nudge | `never`, `always` |
| `navigator.doubleClick`, `navigator.zoomPx` | `ViewerConfig` | `toggle-fit`, 48 | the whole, then back in, scale-free | `fit`, `none`; any px |
| `navigator.zoomStrip.minPx`, `navigator.zoomStrip.heightPx` | `ViewerConfig` | 48, 16 | a reachable span too small to drag gets a zoom strip | any px |
| `navigator.direction` | `ViewerConfig` | `display` | a drag right moves the view right (decided) | `genomic` |
| the gestures | `BehaviorProvider` (`navigator`, on the layouts design's `window`) | in the default spec | a host may leave them out | — |
| the default mount | `SeqViewerSpec.ui` entry (`navigators`) | on in the standalone viewer | the host places parts | `mountNavigator` by hand |
| navigators in an export | spec param of the view (`export.navigators`, optional), `renderSvg`'s `navigators` overriding it (decided) | `omit` | a figure shows data; the window is session | `include` |
| the reachable span washed; the wheel as over the view; the slider's keys | fixed | — | no second use; ARIA's slider pattern | — |

## In the specs

| Decision | Where |
|---|---|
| A navigator is a node of the view (`of`, `axes`, `content`, `tracks`, `bands`, `system`, `place`, `size`, `hidden`); the four forms are presets | core › Navigators |
| Navigators placed by the host; the default mount a `ui` entry, the gestures a behavior | viewer › Navigators |
| `mountNavigator` and `viewer.navigator(def)` public; a host-drawn navigator's contract | viewer › NavigatorModel |
| The rail's look (hairline, heavier segment, knobs) | core › The look |
| The rail is every linear section's navigator; `navigator.positionsRail`, `navigator.railPlace` (a list, `pinned`); its grab zone (`grabShare`) | core › The ruler rail |
| The rows rail beside a declared overview (`navigator.rowsRail`) | core › The rows rail |
| A narrow frame narrows and hides the vertical overview | core › A narrow frame |
| Grips, floor (centred on the truth, measured in the picture), hit, idle, grown | core › The windows |
| A drag ends wherever the button is released | viewer › A drag ends wherever the button is released |
| The cursor says the gesture | viewer › A window's part under the pointer |
| Hidden or disabled: `navigator.form: none`, `hidden`, the rail toggle, a hidden ruler's slim rail | viewer › Navigators |
| Other sections: a ring's own window, reversed, collapsed, follower, blocks, narrow data area | core › The windows |
| The picture: the whole axis from summaries, the view's colors, pending drawn pending, rows, the ideogram, marks | core › The picture |
| Navigator pictures as registry entries | core › Pictures |
| Draw lists and the `window` pass; hits one `window` target with a `part` | core › Draw lists |
| A window per bound section, the reachable span, moving it leads, a cyclic window, display order | core › The windows |
| The zoom strip; two rows by what the section reaches; one stroke, two inks; the funnel | core › The zoom strip |
| Gestures: click centres, sweep from outside (`navigator.dragOutside`), double-click, wheel, a far drag a step | viewer › Behaviors |
| Idle per axis, keys while focused, announcements on settle | viewer › NavigatorModel |
| The readout: one line above the part, right-aligned, at rest the window's region | viewer › Navigators |
| Over blocks and the square: switching blocks | viewer › NavigatorModel |
| SVG export leaves navigators out unless `include` | core › Export |
| Tokens `seq-window`, `seq-grip`, `seq-wash`, `seq-centromere`; the `giemsa` scheme | core › Navigators |
| The configuration items | viewer › Configuration |
| The readout lagging a gesture | design-only: a bug handed to a fix lane, not a decision |
| Removing `placeReadout`'s candidate search | design-only: done, history (the module is gone) |
| A navigator for a grid's own y window | design-only: not built (Tensions) |

## Tensions

- **Spec node or host part.** A node travels with the spec and exports; a host part does not touch a user's spec.
  Decided: a node, and `viewer.navigator(def)` for a host; `place` a hint the default mount honours.
- **Display order against the ideogram's convention.** IGV draws the p-arm left always; here a reversed section's
  ideogram runs reversed so the drag follows the view. Recommendation: `display`, the direction marked; `genomic` keeps
  the convention and inverts the drag. Decided: `display`.
- **A floored window: where it is, or how far it travels.** Centred, it says where but its position is no longer
  proportional to the scroll; as a thumb, it travels the whole span but lies about the place. Recommendation: centred,
  the drag mapped by pointer travel.
- **Sweep against jump-then-pan.** A drag from outside that sweeps cannot also be "go there and keep dragging".
  Recommendation: sweep (three precedents), `pan` an option.
- **What is not summarized is not drawn.** A chromosome's genes with no coarse level, and a 10,000-row alignment whose
  source serves no row summary, draw less in the navigator than in the view. Kept: the bound on work wins
  (§bounded-work); the issue says what was left out. Standing (Self-check).
- **No navigator for a grid's own y window**: the forms cover positions and rows; a grid whose y has a window of its
  own (not following a section) has none over y. Open until a Hi-C view (the Genomes milestone) needs one.
- **The navigator's hover is its own.** A hover there does not light the view, so "what is under my pointer here" is
  answered by the readout only. Why: a navigator pixel is thousands of positions; lighting one residue would be false
  precision. Open to a later use gate.

## Self-check

Claim: "this design covers what a user does with navigators".

1. **Assume it holds.** Every navigator is one part; a user sees the whole, centres, drags, resizes, sweeps, goes back;
   keyboard and screen reader reach the same; idle when nothing moves; exports leave it out.
2. **What it rests on, negated.** (a) one navigator moves one view, negated: one navigator over several views, and
   several navigators over one; (b) the user goes somewhere, negated: the user comes back (the inverse of "centre");
   (c) the navigator's axis is the view's, negated: it is larger or another system (a chromosome under a region, a gene
   under a protein); (d) the navigator reads what the view read, negated: it needs what no view reads at that scale;
   (e) it is chrome, negated: it is part of a figure.
3. **The negated worlds.** (a) Several windows on one strip: IGV's multi-locus and a ring with its linear detail are
   the same case; several navigators on one section (a strip and an overview) must agree on the window. (b) Back after a
   drag across a chromosome returns to before the last click, not before the trip; a double-click's first click and its
   fit are two steps. (c) A window can be dragged where the section cannot go; the window of a follower is moved by
   someone who expects the follower to move. (d) A miniature of genes over 83 Mb would ask for the chromosome at level 1;
   an overview of lazy rows would load every row. (e) An ideogram on a published figure is expected, and a window there
   is meaningless.
4. **Held both.** Integrated into the Decision: a window per bound section (a leader's and its followers'), labelled, and a scroll group's one window; independent leaders each their own window
   (a); `navigator.dragStep: far`, a double-click as one step, and its steps under the view's one history rule (b); the path-bound window, the reachable span, a
   follower's window leading (c); the fallback chains, the shared reduction and image cache, and `seq.render.navigator.layer` (d); `renderSvg`'s `navigators` with
   `include` drawing outlines only, and the ideogram as a `bands` layer usable on a figure's track (e). Standing: (d)'s
   residue, a navigator that draws less than its view where nothing is summarized (Tensions); (e)'s residue, whether a
   navigator belongs in the spec at all (Tensions); and a switch between chromosomes (several systems on one axis),
   refused until the Genomes milestone (EC-184, EC-210). A second negation, (f) "the bound section stays": a spec
   edit removes the `of` section or changes its axis. Integrated: a removed section is `seq.render.navigator.of`
   (warning), the picture drawn with no window, idle; a changed axis is EC-111's (the viewport keeps its centre and
   zoom, or fits) and the window follows it, the base relaid by its key (the axis); an undo restores both. A third,
   (g) "the rail can always be used" (the author: "then the rail is disabled/hidden, etc"), negated: nowhere to go
   (the whole axis shown at the zoom limit), not allowed (a read-only step, an export), no room (a narrow data area, a
   floored window), no ruler (none, or hidden), not wanted (a host's `none`, a reader's corner button), not linear (a
   ring, a reversed or collapsed axis, a follower, several blocks), not a mouse (keyboard, touch). Integrated: the ruler
   rail's rules above (idle, disabled, the floor and its inward drag, the 12 px hit and the touch hit, the slim-rail
   fallback, `navigator.form`, the corner button, each axis case, the slider keys).

## Consequences

- **`seqx-core`**: the `navigator` kind (params above, issues for an unknown `of`, a `bands` content with no bands, a
  `rows` content on a section without rows), the four builder helpers; the resolved view's navigators with their bound
  sections and paths (static, or through the overlay); `render/navigator.ts` (the fit layout through the layout
  interface, the miniature's scale, the rows image, the window pass), `render/window.ts` (the headless window geometry:
  frames, parts, ends, the floor, idle per axis: written by the layouts design for the ring's window, extended here for
  linear and two-axis windows), the `bands` layer type; the `window`
  pass in draw lists, the recording and SVG; `renderSvg`'s `navigators`; issues
  `seq.render.navigator.*` (`layer`, `path`, `rows`, `of`). Tests: property tests (a viewport carried to a window and back
  within a device pixel, at every form and through the codon path; ends never flip; the floored window covers the true
  one; idle exactly when no pan or zoom can move the view); drawing against the data (every band of the chr17 fixture
  has its mark and hit at the fit; every window's outline is its viewport's image); records for a strip, an ideogram
  and an overview; reused buffers compared against fresh ones; quick representatives under `itWithin`: a window move
  redraws in 0.5 ms and allocates nothing warm; the bases: a strip of p53's tracks in 4 ms, chr17's ideogram (about 60
  bands) in 1 ms, the 10,000-row overview from its summary grid in 8 ms. **Size**: the kind, the navigator layout, the
  window pass's linear and two-axis parts and the `bands` layer, an estimated 8 kB minified (the window geometry is
  the layouts design's, the reduction and image cache the rows': nothing copied); measured when it lands, the code kept minimal.
- **`seqx-viewer`**: `NavigatorModel`; the `navigator` behavior provider (click, drag by part, sweep, double-click,
  wheel, Esc) on the `window` behavior, and the `window` target's `outside` part in the bindings; `SeqViewerSpec.ui.navigators`; commands `navigator.center`, `navigator.zoom-to` (navigate) in the change
  kinds table; `ViewModel.snapshot()` and `record(snapshot)`, the double-click's deferred step; `mountNavigator` and the default mount's placement; the sliders,
  keys and value text; Configuration rows. Tests: a headless twin per story on a fake clock (a strip and the view: drag
  the window, the view moves; click, Alt+Left returns; a follower's window leads; a scroll group's one window; two
  windows of independent leaders on one strip, each dragged moving its own section; idle on and off with no relayout; keys only; a removed bound section); the allocation tier for the drag.
- **`seqx-examples`**: p53's gene region with chr17's ideogram (a cytoBandIdeo fixture for chr17, public, a few kB),
  p53's tracks with a strip of the tracks, the protein's strip showing its gene follower's window; the Globins
  10,000-row view with its vertical overview and a miniature; the plasmid's linear detail with a strip wrapping the
  origin. Each view's test checks drawing against the data for its navigators.
- **Storybook**: a story per form in light and dark, canvas and SVG side by side, with its headless twin; a Design page
  "Navigators" (folded from this design).
- **Specs that change**: core `_spec` Language (the kind), Resolved view, Render (Navigators, Draw lists' `window`
  pass, Hits' `window` target and its `outside` part, Export), Configuration; viewer `_spec` Models (NavigatorModel; the history's
  `snapshot` and `record`), Input (the part), Behaviors, Bindings, Commands, Configuration; [ux-foundations.md](ux-foundations.md)
  section 8 (the navigator row as above, and its gesture table's line "drag on an axis, navigator or with Space held |
  pans" split: on a navigator a drag on the window pans, a drag from outside sweeps; a wheel over a navigator zooms about
  the window's centre); [features.md](features.md) (the Navigators row); [edge-cases.md](edge-cases.md)
  (EC-110 kept, its tests).

### What this design needs from the others

- **The rows design** (`[grids-and-rows.md](grids-and-rows.md), "What a navigator reads", aligned): the rows in display order (the count,
  a uniform height, collapsed clades as one row, pinned rows outside the scroll, an order version), the band's
  `bandTop`, visible height and `rowScale`; row summaries as a coarse grid through the shared reduction and image LRU;
  one value a row as the fallback (identity, else the non-gap share).
- **The layouts design**: the layout interface callable at the fit window into a given box with a height scale (the
  miniature); for text, the positions on the lines in view and a move to a position range (a strip over a text section);
  the ring's window written as the shared window geometry (`render/window.ts`), which this design extends.
- **The axis design**: display order with reversed pieces and their mark, collapsed introns, the cyclic window's pieces;
  the path from a region to its chromosome.
- **The layouts design**, again: `FrameInput.viewports` (one viewport per leading section), which several
  leaders' windows on one navigator read.
- **Themes**: the `giemsa` scheme and the tokens `seq-window`, `seq-grip`, `seq-wash`, `seq-centromere`, light and
  dark. **Stories**: a centre's
  animation through the `AnimationPort` once it exists (until then the view's jump animation).

- **A strip with no `content`** draws its axis ruler at its `size`, never a blank well (FX12, the plasmid's strip).
