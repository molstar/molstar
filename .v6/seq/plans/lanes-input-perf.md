# Plan: tracks as lanes in the text view, input as streams, the performance pass

Status: built 2026-10-14 (lanes §2.6–2.8, input and events §3.6, the performance pass §4.3), one PR on
`v6-seq-improvements-261014`; planned 2026-10-14 (David). Takes three items from
[../plan.md §4](../plan.md) and `packages/sequences/todo.md` (Deferred features: "Tracks in the text view, as lanes",
"Input as streams, reconciled with Mol\*'s"; Follow-ups: "The host options' performance pass"). Each leaves those lists
in the commit that builds it. Facts from the code at `318412a94`.

## 1. Order and scope

One PR, in this order. Each step ends with its tests, its spec change (§5) and its todo entry removed:

1. **Lanes** (§2), in Seq\*'s renderer, then Prism's predicted story using them in Text.
2. **Input and events** (§3): the viewer's own input as now, public rxjs streams over it for consumers (UIs, Prism,
   hosts), one events surface.
3. **The performance pass** (§4), last: its guards then cover the new lane kinds and streams as well as the round-2
   host options.

Who: the renderer and Seq\* input changes go to Opus (design and implementation, high); the Prism side and mechanical
consumer updates to Sonnet (medium). Each agent runs only the tests its change touches; the heap tiers run once, in §4,
on a quiet machine.

## 2. Tracks in the text view, as lanes

### 2.1 Today

- `layout.ts` `bands()` (about line 1672) gives each track a band per line: the renderer's `lineHeights` if it has one
  (ranges, secondary and points: `runs.ts` `laneLineHeights`, lanes compacted per line by `lineLanes`), else a fixed
  height: letters `ceil(wrappedFontPx × 1.6)`, **any other non-ruler type `rowHeight`**.
- That fallback is why a `values` track (Prism's pLDDT) takes a lane in text and draws nothing there: `values.ts` draws
  only in `linear`. The same holds for `color` (draws no cells in text), rows, logo, matrix, pairs, map, variants and
  dot plots.
- Underline lanes already draw and hit-test in text (`draw/features.ts` `drawLineLanes`, its hit and hover outline).
- Prism works around it: `sequenceDocument({ tracks: false })` builds a letters-only document for Text, so a style
  switch rebuilds the document (a moment of no text; prism-sequence.md §9, round 3).

### 2.2 What changes

| Layer kind                                   | In the text view                                                                                          |
| -------------------------------------------- | --------------------------------------------------------------------------------------------------------- |
| letters                                      | as now                                                                                                    |
| ranges, features, secondary, points          | underline lanes, as now; heights `wrappedLaneHeight`, the first after `wrappedTrackGap`                   |
| `values`, `color` alone on its track         | **new: a strip** under each line, one color per position from the layer's color theme, `wrappedStripHeight` high |
| dot plot, pairs, trees, rows (multi-row), logo, matrix, map, variants | **no lane**: they need tracks, and reserve nothing in text                       |

- **No lane without a drawing.** `bands()` loses the `rowHeight` fallback: a renderer whose `geometries` leave out
  `wrapped` (`renderer-helpers.ts` defaults them to `['linear']`) reserves nothing; one that includes it gives its band
  through `lineHeights` (or the letters' height).
- **The strip.** `values` declares `['linear', 'wrapped']`; `color` already declares `wrapped` and draws no cells there
  (§2.5). A strip lane is one band per line under the
  letters (after the underline lanes, in track order), split per position like the letters' cells, colored from the
  layer's theme (a `values` layer's `color` bands or ramp, pLDDT's AlphaFold bands); a missing value draws nothing.
  New render prop `wrappedStripHeight` (default 3 px; `compact` 2). Drawn as the letters' cells are: by run of one
  color, `setTransform` placement, no allocation per frame (the performance contract, `seqspec:viewer#performance`).
- **Hit and hover.** A position under a strip is that layer's hit (the hit's `track`, `layer`, `index`), so the hover
  and the hover card show the value, as linear does. `wrappedPositionNear` and the selection need nothing new: lanes
  only add height to a line, which `lines.top` already holds. The hover outline of a strip is the position's cell.
- **Consistent heights**: every lane kind takes its height from a render prop, so a line's spacing reads the same
  whatever is under it; the `wrapped` and `compact` presets set them.

### 2.3 Prism: the predicted story with lanes

- **One document for both styles.** `sequenceDocument` drops its `tracks` option: Text and Tracks show the same
  document, and a style switch only changes the preset (no rebuild, no blank moment). The pLDDT track keeps its
  `encoding: 'line'` (Tracks) and is a strip in Text.
- **`Viewer/With extensions › Predicted`** in Text: the letters in the residue palette, the secondary structure
  as underlines (rose, amber, teal), pLDDT as a strip in `PrismPlddtBands`. The legend is the same in both styles.
- The track gutter stays Tracks-only (Text has no track rows to put icons beside); the right-click menus work on a
  lane as on its track (the hit names the track).
- prism-sequence.md §9 gets a "Round 4 (as built)" entry; the round-3 note "Text has only the letters" is superseded.

### 2.4 Tests and guards

- `wrapped-spacing`: band heights per lane kind; no band for the kinds that do not draw in text; a document with only
  letters spaces as before.
- `draws` (recorded draws): a strip lane, a strip with gaps (missing values), strips under underline lanes.
- `text-select`: a drag that crosses strip lanes; hits on a strip.
- A heap guard for the strip (`wrapped-lanes.heap.ts`, the guard template with its positive control), beside
  `text-select.heap.ts`.
- Prism: the model's tests (one document for both styles; the style switch keeps the document object).

### 2.5 To settle while building

- `color` in text: a `color` layer beside letters is what `wrappedColors` tints them with; a `color` layer alone on its
  track (no letters) is the strip. Check how `color.ts` and `letters.ts` share it today before choosing.
- Where a strip sits relative to underline lanes when one track has both (track order, or strips last).

### 2.6 As built

What differs from §2.2, and the §2.5 decisions (2026-10-14, Seq\*'s renderer and layout; Prism's side is §2.3, next):

- **What text did before** (a correction to §2.1): `values` was `tracks`-only in the format's layer types, so in text it
  fell back to `color` (a `svs.render.fallback` warning) and drew 22 px cell backgrounds; the "no lane" kinds already
  resolved to their fallbacks or were omitted (`pairs` and `variants` as points, `rows` as letters, a ruler and a dot
  plot omitted). Now `values` is `tracks text` (`core/src/format/tables/kernel.ts`, the format reference regenerated) and
  its renderer declares `wrapped`.
- **`bands()`**: letters and `emphasis` the font's band (`emphasis` drew a 22 px row beside the letters: a track of
  letters with emphasis, case 1's, is now 20 px a line, as letters alone); the underline and strip lanes by
  `lineHeights`; the `rowHeight` row is kept only for a renderer that declares `wrapped` with no `lineHeights` (an
  extension's: the testing kit's stars), so declaring `wrapped` never hides a layer. A track whose `lineHeights` give 0
  on every line takes no band (before: zero-height entries).
- **The strip** (`draw/strips.ts`: `stripBand`, `drawStrip`, `outlineStrip`; `runs.ts` `stripLineHeights`, in the
  renderer API): `wrappedTrackGap` + `wrappedStripHeight` on each line of a block the layer has data on (values, rows
  or a stream), none elsewhere; a rect per run of one color and alpha within a group. A streamed `values` block is read
  per position through its regions.
- **§2.5, `color` in text**: `wrappedColors` tints the letters with the letters layer's own theme, not a sibling `color`
  layer. So: on the letters' track a `color` (or `values`) layer colors the letters' cells behind them, as `color` did
  before; alone on its track it is the strip. The renderer tells them by its band (no taller than a strip lane: the
  strip). A `values` layer with no color theme beside the letters draws nothing (its ink would hide them).
- **§2.5, strips and underlines**: track order. Each track is its own band; the layers of one track share it, as in
  linear, so there is no "strips last" reordering.
- **Hit and hover**: no hit code changed (a line's band per track already names the track, its first layer and the
  position; the hover card's value provider reads the track's `values` and `color` layers). `values` and `color` gained
  a wrapped `overlay` outlining the hovered position's strip cell.
- **Presets**: `compact` sets `wrappedStripHeight: 2`; `wrapped` keeps the defaults as its own (setting them there would
  override a host's render props).
- **Guard**: `viewer/_test/wrapped-lanes.heap.ts` (draw frames scrolling, and pointer moves over strips and underlines
  with an overlay frame each, its positive control). Its hover case runs first: after the draw case in the same worker
  the hover's overlay kept allocating past the settle (about 30 bytes a step, with or without the strips' points).

### 2.7 As built (batch 2)

Three items beside the lanes (2026-10-14, Seq\*'s renderer, layout and behaviors; Prism untouched):

- **A. A ligand's pill never vanishes** (`renderers/letters.ts` `pill`): narrower than its box (the cell less its 1 device
  px insets under 2 device px), a filled cell in the faint ink, at least 2 device px wide, centered on the position,
  unlabeled. Tested in `draws` at the fit of narrow views (DPR 1 and 2). At the narrowest fits the last block can fall
  outside the view (the blocks' 28 px gaps): that is the fit's, not the pill's.
- **B. Chips: short units in the text.** A ligand block in the `position` form with one position is a chip: its run's
  one position is as wide as its chemical component's label (`chipFontPx`, half the letters' band, as the pill draws
  it; `CHIP_PAD` 5 px either side; two cells at least), measured once per prepared run by the layout's text measurer
  (no CSS; the SVG export's metrics measurer alike). Representation: each chip stays its own line (one block a line
  holds everywhere: draws, marks, the selection, the keyboard), and the line table gains `left` (a chip's x on its
  row), `follows` (it shares the previous line's row) and `bottom` (its row's). `layoutRuns` packs consecutive chips
  on rows, the group gap apart, wrapping as words do. The layout lays out by row: every line of a row at one top, each
  track's band the row's tallest, chips' rows after the chains at `wrappedGap` (not `stackGap`), no group labels'
  band on a chip row. The views carry the line's `lineLeft` (`AxisView.x`, `DrawView.xd`/`xe`; set where the line is:
  the composer's `place`, the hit, the conformance kit), so every draw, mark and outline places chips unchanged. The
  hit picks, among a track's entries at one top, the chip whose `left` is at or before x (`lineAtX`), as does
  `wrappedPositionNear`; the keyboard's cursor and `Layouts.reveal` (now `lines.bottom`) read the offset too. Not
  done: a cap on how many chips a strip shows (the todo entry's last question): every short unit flows.
- **B, fixed after Prism's first look** (David, AlphaFill CYP3A4 with HEM, ERY, NI, ZN): (1) only the last chip of a
  row was drawn after a scroll: a scroll's redraw of a row its edge cut cleared the whole row per entry, and the first
  entry in view was the last of those sharing the row's top (`entriesIn`'s search), so the chips before it were
  skipped. Now `entriesIn` starts at the first entry of a shared top, a chip after another on its row clears nothing
  of its own in that redraw, and a chip's entry is clipped (and cleared, for a stream's redraw) to its chip. (2) The
  "highlight at the line's start" was the first chip's (HEM at the row's left) with its pill erased by (1); every
  wrapped path places through `place` and the line's `lineLeft` (audited: the composer's draws, marks, the pending
  drag, overlays, the hit, the keyboard's cursor). (3) A chips' row is spaced as a new line: `wrappedGap` under the
  last line's bottom band, then the group labels' band's room, empty (with `compact`'s touching lines, the band alone).
  Tests: `draws` (every chip's pill through a scroll's redraws across a wrap at 220 px, with and without the fix; the
  hover and the selection at a chip's x), `wrapped-spacing` (the row's spacing, touching lines too).
- **C. Reveal without selecting.** `ctx.reveal(target, opts)` on the context, not a manager: it holds no state (a focus
  is state, a reveal only a move) and needs both the projector (the viewports) and the text's lines (the layout), which
  only the context holds together (`reveal.ts`). Text: `Layouts.reveal` (the focus's path, `REVEAL_PAD` 12 px); tracks:
  the new `Viewport.revealBlocks` pans at the cell (centered; wider than the view, its start `pad` px in), nothing
  when shown (whole, or its start where wider); `seqv:reveal` the command. `revealHighlight` (optional,
  `SeqViewerOptionalBehaviors`): params `channels` (names or `*` prefixes; `links:*`), `dwellMs` 200, `pad` 12,
  `animate`; a stable timer callback and an allocation-free scan (`Map.forEach` with callbacks made once) per
  highlight change; it skips the pointer's section, a press or a drag (dropping the wait), and this viewer's own hover
  (owner set, not external).

### 2.8 As built (batch 3): tracks drawn in full in the text

For figures (David, 2026-10-14):

- **Opt-in** (David, during the batch): the render prop `wrappedTracks: 'lanes' | 'full'` (style `wrapped_tracks`),
  `lanes` the default and the `wrapped` and `compact` presets': a track's `height` alone changes nothing in text (Prism's
  pLDDT keeps its strip); `full` is the figure's.
- **Which tracks** (under `full`): a track with a `height` and no `letters` layer (the letters are the text: their track keeps its band
  and lanes; pins over letters stay a linear thing). On a chip's line it keeps its lane.
- **Unit**: rows of the text, a text row the letters' band (1.6 × `wrappedFontPx`, `textRow`): the format's `height`
  is in rows already, and a row of the text scales the figure with its font (`rowHeight` would not follow it).
- **Drawing**: the line stands in for a linear view (`AxisView.trackForm`, `setTrackForm`: `visLo`/`visHi` the line's
  positions, `cell` the run's advance; set where the line is placed: the composer's `place`, the hit, the conformance
  kit, and put back after), so the 1D renderers' linear paths draw over each line through `xd`/`xe` and `placeRun`,
  which already follow the line and its group gaps: `values` line, area, bars (one scale, the layer's domain; a line
  joins across group gaps and breaks at a line's end, no continuation stub: the next line starts at its first
  position, so no value is drawn twice), `points` glyphs, `ranges` items and labels, `secondary` shapes (the coil to
  the line's last cell). Hits: `ranges`' linear hit with `hitWindow` over the line in track form; points as in tracks.
  `variants` stay tracks-only (their fallback, points, draws in full).
- **Axes and labels**: `wrappedTrackAxes`, `wrappedTrackLabels` (style `track_axes`, `track_labels`), off by default;
  the axis is the `values` renderer's `gutter` (in the composer's wrapped gutter pass, now running the layers' gutters
  for a full track), its labels formatted once per prepared state; 30 px left of the lines (`WRAPPED_AXIS_ROOM`); the
  label is the composer's, on line 0.
- **Guard**: `viewer/_test/wrapped-tracks.heap.ts` (its own file), with its positive control.
- **Story** (its style `wrapped_tracks: full`): `Seq*/Text view › Publication tracks` (`storybook/src/_stories/text-view.stories.tsx`,
  `support/seqstar-publication.ts`): p53 from the committed fixture, Kyte–Doolittle as a line, an illustrative
  pLDDT-like area in AlphaFold DB's bands, Zn²⁺ sites as glyphs, UniProt regions labeled, illustrative secondary
  structure, axes and labels on, with the SVG and PNG export card.
- **Annotations UI** (next, another agent): a ranges track colored by an item `color` column needs nothing new
  (`field_color`), drawn as underlines in text or in full with a height. The composer of text tracks is a todo entry.
- **Prism** needs no change: under the default `lanes` its pLDDT track (`height: 2`) stays a strip in Text.

### 2.9 As built (batch 4): the text drawn in full, laid out as a figure

From David's look at `Publication tracks` (2026-10-15):

- **Range labels per line**: in track form a range's label is its item's part on the line (`extent` clamped to the line,
  the glyph and label edges the line's own extent): centered there, ellipsized, dropped under three graphemes, never
  past the line's end. Lanes draw no labels.
- **A line's block**: under `full` each band is `wrappedTrackGap` below the one before (a lane carries it already), the
  next block `wrappedGap` below; the new `figure` preset (`wrappedTracks: full`, axes and labels on, gaps 3 and 18 px)
  is the story's style (`style: 'figure'`). The group labels' band stays topmost; a track ordered before the letters
  sits between it and the letters (as built before; a test now holds it, and the larger gap between blocks makes it
  read as its line's).
- **Empty bands collapse** (`wrappedCollapseEmpty`, default true, under `full`): a features-only track (its layers'
  `lineHeights` are `laneLineHeights`) with no lane on a line takes no band there; `values`/`color`/others keep theirs.
- **Track labels in a row of their own**: `EntryKind.trackLabel`, above the band on the first line the track has one
  (`trackLabelRow`), drawn by the composer's entry pass (`trackLabel`); the left gutter was the alternative, rejected:
  it holds the axes and is 30 px wide.
- Checked in the preview: the story (an Above annotation added), and Prism's predicted story in Text (lanes unchanged).

## 3. Input and events

### 3.1 Today

- **Input**: `canvas/dom-input.ts` (`DomInput`, one per surface) fills reused records (`PointerInput`, `WheelInput`,
  `KeyInput`) and hands them to `InputDispatcher` (`input.ts`), which hit-tests and emits the base events on
  `ctx.interaction` (`pointer`, `hover`, `click`, `drag`, `doubleClick`, `wheel`, `key`, `focus`, `leave`) through
  `SeqEvent`, an allocation-free emitter with priorities: behaviors claim a press in priority order.
- **What a host sees** is spread over five surfaces: `ctx.interaction` (internal-shaped: a raw `Hit` with section
  index, claims); `ctx.events` (`activate`, `menu`, `issue`, `log`, `patch`); states (`hover.live`,
  `selection.current`, `pointer`); the element's DOM events (`seq-hover`, `seq-selection`, …); and nothing for the
  chrome (a click on a legend entry, a track action, a section switch is visible only through what it changes).
  Prism's track gutter re-derives the row under the pointer from the layout for want of a hover stream with its track.
- **Mol\*'s** `InputObserver` (`packages/core/src/util/input/input-observer.ts`) follows mouse and touch events,
  listens for keys on the window, and allocates a payload per event; its `ButtonsType`, `ModifiersKeys` and `Binding`
  (`@molstar/core/util/binding`) are the vocabulary Canvas3D's `interaction` speaks. Seq\*'s packages may import
  Mol\*'s (the restriction to `seq-molstar-poc` came from the separate repository's vendored Mol\* and was lifted
  2026-10-14, in this PR).

### 3.2 Decision: the viewer's input stays, public streams over it

Two layers (David, 2026-10-14):

- **Inside the viewer, nothing changes in kind**: `DomInput` on each surface's DOM events, the reused records, the
  dispatcher, the hit test and the behaviors' priority-ordered claims (`SeqEvent`), under the performance contract and
  its guards as today. Seq\*'s `Binding` stays (allocation-free `match`) and gains a `description`.
- **For consumers outside it** (the DOM and React UIs, Prism, a host): rxjs streams on `ctx.interaction` (Canvas3D's
  name for the same thing), fed by the dispatcher after the behaviors. They may allocate: a fresh payload per event,
  in Mol\*'s vocabulary. The contract they keep is **no cost unobserved**: the dispatcher builds a payload only when its
  stream has a subscriber (`Subject.observed`), so the input path's heap guards hold with no one subscribed.

Mol\*'s `InputObserver` is not used and not changed here; moving it to pointer events (Canvas3D and Seq\* on one
observer) is a separate decision, if it is ever wanted.

**The streams**

- `hover` (the hit changed), `click` (with `count`), `drag` (`start`, `move`, `end`, `cancel`), `press`, `menu`,
  `wheel`, `key`, `leave`.
- **`menu` comes from the DOM's `contextmenu` event** (the secondary press, Ctrl+click on a Mac, the menu key and
  Shift+F10 alike), hit-tested at its point, so a subscriber decides the browser's menu: `e.native.preventDefault()`
  during the call shows its own instead. Only that event can suppress the browser's menu (a `preventDefault` on the
  press cannot), and it fires once per ask, so it is outside the contract. With no subscriber the browser's menu
  shows, as now; the `contextMenu` behavior (`events.menu` today) becomes this stream. Observe-only: a subscriber
  watches a gesture, never claims it; the behaviors' base events move to `ctx.input.events` beside the sink.
- **One payload** (`SeqInputEvent`, a fresh object a consumer may keep):
  - `part`: where: `kind` (`data`, `header`, `gutter`, `dataColumn`, `dataColumnHeader`, `navigator`, `vnav`,
    `scrollbar` on a surface; `legend`, `trackAction`, `sectionSwitch`, `viewTab`, `search`, `snapshotPlayer`, … for
    the chrome), `section` and `track` refs (names, not indices), `key` (a chrome node's `data-seqv-key`);
  - `at`: what: `layer`, `lane`, `item`, `block`, `index`, and `coordinates` computed on first read (as `seq-hover`'s);
  - `place`: `x`, `y` (surface px), `clientX`, `clientY`;
  - `buttons`, `button`, `modifiers`: Mol\*'s `ButtonsType` and `ModifiersKeys` (types from `@molstar/core`), so a host
    linking Seq\* and Canvas3D reads one vocabulary; `time`;
  - `native`: the DOM event, during the call only (a consumer may `preventDefault`).
- **Scoping**: `ctx.interaction.on(stream, scope?)` returns the stream filtered by `{ part?, section?, track? }`.
  `ctx.interaction.hover` and the rest are the unscoped streams.
- **The chrome**: the DOM parts (`dom/parts.ts`, and the React UI's) tag their roots with `data-seqv-part`; one
  delegated listener per mounted UI (at the element's root, in its shadow root) feeds the same streams with `part` from
  the nearest tagged node. It lives in the UI layer (`dom/`), outside the contract; no listener per node.
- **Bindings in Mol\*'s terms for consumers**: a behavior's bindings are listed with their descriptions and their
  triggers as Mol\*'s `Binding` (for the default UI's and Prism's help); the behaviors match Seq\*'s own.

**The test case** that motivated this (David): Prism's track menus (prism-sequence.md §9, round 2) needed the
`contextMenu` behavior added to the spec, a `contextmenu` listener of their own, a second hit test through
`seq.layout.hitTest` with a region set, the track's ref read through `graph.assembled`, and the place converted from
surface px. Done when that is one `ctx.interaction.on('menu', { section })` subscription that reads the track's ref and
the place from the payload and calls `preventDefault` where the track offers actions.

### 3.3 One events surface

| Today                                             | After                                                                                   |
| ------------------------------------------------- | --------------------------------------------------------------------------------------- |
| `ctx.interaction` (base events, claims)           | `ctx.input.events`, for behaviors (documented for behavior authors only)                |
| `ctx.events.activate`, `ctx.events.menu`          | `ctx.interaction` streams (`click` on a header, `menu`)                                 |
| `ctx.events.issue`, `log`, `patch`                | stay on `ctx.events` (not input)                                                        |
| `hover.live`, `selection.current`, `pointer`      | stay (states: what the viewer holds, replayed on subscribe; the streams are what happened) |
| DOM `seq-*` events                                | stay, fed from the same states and streams; `seq-input` added only if a DOM host asks   |

Consumers updated in the same step: Prism's `sequence-control.tsx` (`events.menu` → `interaction.menu`; its own
`contextmenu` listener and second hit test go, `preventDefault` only where the track offers actions) and
`track-gutter.tsx` (`interaction.on('hover', { part: 'header' })` with the hit's track, instead of re-deriving the row),
the storybook's `seqstar-link.ts` and `seqstar-viewer.tsx`, `linkViews`, nightingale's mount if it reads
its old `ctx.interaction`. Seq\* is unreleased in v6, so the rename is a change, not a break; it goes into
`.v6/plans/changelog.md`.

### 3.4 Tests and guards

- Each stream's payload per part (surface regions and chrome parts), scoping, observe-only (a subscriber cannot claim;
  a behavior's claim still wins), coordinates on first read only, `native` cleared after the call.
- No cost unobserved: `input.heap.ts` and `canvas/dom-input.heap.ts` run as today with the streams created and no
  subscriber (each case with its positive control); one case subscribes and unsubscribes and checks the path returns to
  flat.

### 3.5 To settle while building (Opus, a short design note in this file before the code)

- The chrome parts' `kind` list and their keys, from `dom/parts.ts` and `react/`; whether a focused chrome node's keys
  go through `key` with `part`.
- Whether `ctx.interaction.on` filters per subscription (an rxjs `filter`) or the dispatcher keeps scoped subjects.

### 3.6 As built

Built from the design note (Opus, 2026-10-14) with the coordinator's simplifications. Paths under
`packages/sequences/viewer/`.

- **Two layers.** `ctx.input` is `SeqInput` (the sink and `events`, the base events behaviors claim from, as before, plus
  `menu`); `ctx.interaction` is `SeqInteraction` (`src/interaction.ts`, headless): rxjs `Subject`s `hover`, `press`,
  `click`, `drag`, `menu`, `wheel`, `key`, `leave`, `activate`, `on(stream, { part, section, track })` (an rxjs
  `filter` per subscription, no scoped subjects), `observed(stream)` and `chrome(...)`. The dispatcher (`src/input.ts`)
  emits each public event after the base one, only when `subject.observed`; the records carry `native` (the DOM
  event, set by `DomInput` around the sink call). After `next` the payload's `native` becomes a stand-in; headless, a
  call's stand-in records `preventDefault` (what `menu` and the menu key report as handled).
- **The payload** (`SeqInputEvent`, fresh): `part` (`kind`: a `HitRegion` name or a chrome kind; `section`, `track`,
  `group` refs; `key`), `at` (`layer` ref, `lane`, `block`, `index`, `row`, `split`, `item` with annotation key and
  id; null on chrome and keys), `place` (`x`, `y`, `clientX`, `clientY`: the native's, or the surface's rect plus the
  point from the keyboard; NaN unknown), `buttons`, `button`, `modifiers` (Mol\*'s `ButtonsType`, `ModifiersKeys`,
  type-only from `@molstar/core`, now a dependency of seq-viewer; modifiers are 16 frozen objects), `time` (the
  native's `timeStamp`, else the scheduler's), `owner`, `origin`, `native`. **Deferred**: `at.loci` and
  `at.coordinates` (no `Hit.copy`, no `HoverManager.lociAt`, `describeLoci` stays in `dom/element.ts`): todo.md "Loci
  and coordinates on input events".
- **The menu.** `DomInput` listens for `contextmenu` into its own record (`pointerId` -1 when `pointerType` is "" or
  the id is -1) and calls the sink's `menu`, never preventing the default itself. The dispatcher hit-tests it, emits the
  base `menu` (`MenuEvt`, claimable; `keyboard`), then the public one. The `keyboard` behavior's menu key calls
  `input.keyMenu(k)` (handled only where a subscriber prevented it) and claims every keyboard `menu`, placing it at its
  cursor's cell (started as a focus would when there is none). Dedupe, a one-shot flag (no timer, no config): after a
  key's menu the dispatcher remembers the section and whether it was prevented; the next `contextmenu` on that section
  before any press is not emitted and gets `preventDefault` if the key's did; a press, a later key, the section's blur
  or a leave clears it (a browser that fires no request after a prevented key swallows no later one). `contextMenu`
  (`ContextMenuBehavior`), `events.menu`, `events.activate`, `MenuEvent`, `ActivateEvent` and the context's
  `openMenu`/`activate` are gone; `activate` is a stream fed by `input.activate(hit, input, owner)` (the `headers`
  behavior, the keyboard's Enter).
- **The default UIs**: `TrackActions` (DOM and React) subscribe `interaction.on('menu', { section, track })`, open on
  origin `keyboard` and prevent the default; a right-click stays the browser's.
- **The chrome**: `dom/chrome-input.ts` (`listenChrome`, `chromePartOf`), one delegated listener per mounted UI root
  (`mountSeqViewerRoot`, React's `SeqViewerRoot`): `click`, `contextmenu` → `menu`, `pointerover` → `hover` on a part
  or key change, `pointerleave` → `leave`; each checks `observed` first. Every part root of both UIs carries
  `data-seqv-part` (its kind; the shared classes are the fallback), sections `data-seqv-section`, the track box
  `data-seqv-track`.
- **Bindings**: `SeqBehaviorProvider.bindings` (descriptions by param path, `defineBehavior`'s last argument) and
  `ctx.behaviors.bindings()` (`src/bindings.ts`) as Mol\*'s `Binding`; described: `select`, `dragPan`, `wheel`,
  `rowZoom`, `dotPlotSelect`, `keyboard`.
- **Prism**: `sequence-control.tsx`'s track menu is one `seq.interaction.on('menu', { section: SectionRef, part:
  ['data', 'header', 'gutter'] })` reading `e.part.track` and `e.place.clientX/Y`, preventing the default only where
  `viewer.trackActions(track)` is non-empty; `Hit`, `HitRegion`, its hit test and its `contextmenu` listener are gone,
  and `model.ts` no longer adds a behavior. `track-gutter.tsx` follows `on('hover', { section })` (the track of any part
  of the row, as before; not `part: 'header'`, which would have lit the icon only over the header) instead of its own
  `pointermove` and row lookup; its leave listeners stay (the icons are outside the surface).
- **Tests**: `interaction.test.ts`, `context-menu.test.ts` (rewritten), `bindings.test.ts`, `dom/chrome-input.test.ts`,
  `react/chrome-input.test.tsx`, `canvas/dom-input.test.ts` (the `contextmenu` listener, `native`). Heap:
  `input.heap.ts` (a context-menu storm unobserved, flat; element changes with a hover subscriber caught, flat once it
  unsubscribed), `keyboard.heap.ts` (the menu key unobserved, in the staying presses), `canvas/dom-input.heap.ts` (a
  `contextmenu` storm into its record, flat: the record's `time` is not read, a fractional `timeStamp` boxed 12 bytes a
  request). §4's guard (c) is the input storm.
- **Open**: Safari's keyboard `contextmenu` has no `pointerType` and reads as the pointer's; with no subscriber the menu
  key's default is the browser's (the default UIs prevent it); a macOS Ctrl+click is also `select`'s `toggle`; rxjs
  reports a throwing stream subscriber asynchronously, not as `seqv.behavior.error`; `_spec/reference/*` regenerated (the `contextMenu` row left, the `keyboard` and `headers` descriptions changed).

## 4. The performance pass (last)

### 4.1 First: a quick allocation tier (David, 2026-10-14)

The heap guards take minutes (about 4 for the viewer, half an hour for every tier) and need a quiet machine: their
warm-ups are timed (`FRAME_WARM_MS`), so tier-up follows the clock, and they infer allocation from heap growth across
a short and a long run. They become the deep check, run at checkpoints and before a release; the edit loop gets a
cheaper guide.

- **`perf:quick`**: each case measures bytes allocated per step directly, through `v8.GCProfiler` (the heap before and
  after each collection in the window: allocated = Σ(before − previous after) + the final delta), in a worker with a
  small young generation (`--max-semi-space-size=1`), so 32 bytes a step over 20,000 steps shows as several
  scavenges with their sizes.
- **Warm-up by count**, not time: a fixed number of steps, or until two windows agree on bytes a step. No dependence
  on load, so files run in parallel and beside other work.
- **A guide, not a gate**: each case prints bytes a step against a committed baseline (`perf-baseline.json`); a clear
  regression fails (over the baseline by a threshold, about 16 bytes a step), a smaller move is listed, `--update`
  records an intended change. The positive control (`leak`) still shows each case can see allocation.
- **The same cases**: a harness takes the step function each `*.heap.ts` already has, so every guard gets a quick twin
  without new setup; new hot paths get a quick case first.
- **A spike before the rest**: port `wrapped-lanes.heap.ts` and `hover.heap.ts`; done when the quick tier catches the
  order-sensitive hover allocation (about 32 bytes a step, below) and every positive control, in seconds and with a
  load running. Then the remaining guards get twins; the target is the viewer's quick tier well under a minute.
- **Docs**: Seq\*'s `AGENTS.md` ("every hot path has a quick case; the deep guards run at checkpoints"), the tier
  table there, and the spec check's guard-template rule extended to the quick cases.

### 4.2 The pass

From todo.md's follow-up and plan.md §4, measured with the quick tier and confirmed once with the deep one:

- **The cost**: `unobservedAt` (`renderers/letters.ts:206`) scans the unit's `observed` ranges per position per frame.
  Build a mask per block once per prepared state (or walk the ranges with a forward cursor, since positions come in
  shown order); the zoomed-out cell class reads the same mask.
- **Guards**, each with its positive control: (a) letters with `unobservedGray` on a linear section zoomed in over a
  unit with several `observed` ranges; (b) wrapped text with `wrappedColors` and `unobservedGray` (the tinted path), now
  with strip and underline lanes; (c) a storm of secondary presses with `contextMenu` (the dispatch, not the stream);
  (d) text selection and hits on a centered `wrappedMaxWidth` layout.
- **Runs**: the quick tier throughout; once at the end, on a quiet machine, `test:heap`, `test:browser-heap`,
  `test:slow`; numbers that move are fixed, not loosened.
- **Root-cause, not order around**: `viewer/_test/wrapped-lanes.heap.ts`'s hover case (pointer moves over a text
  section, an overlay frame each) allocates about 32 bytes a step when it runs after the file's draw case (full
  frames), with or without the strips' points, so not the strip code; it passes run first, which is how the file runs
  it now (§2.6). Find what the draw case leaves (a shared helper's feedback turned megamorphic, a module the reset
  does not reload) and remove the order.

### 4.3 As built

The quick tier, every guard's twin, count warm-ups for the deep tier, and the pass (2026-10-15; Node 24, M4). Paths
under `packages/sequences/`.

- **The quick tier** (`test:perf` in core, viewer, nightingale and molstar-poc; `perfTests` in `vitest.shared.ts`):
  `core/_test/model/alloc.ts` measures the bytes a step allocates exactly through `v8.GCProfiler` (Σ heap before a
  collection − the previous after, plus the final delta; the code and trusted spaces apart; the measurement's own
  bytes calibrated once a process). Flags: a 1 MB young generation, `--single-threaded`, lower tier-up thresholds. The
  warm-up ends only on two quiet windows that compiled nothing, re-measuring when the measured windows compile; a
  steady non-zero read runs to the case's cap (the "two windows agree" rule latched onto lower-tier plateaus and lulls).
  `quickCase` (`core/_test/model/quick.ts`) measures the step and then its control `leaking(step)` on the same state
  (opt out: `control: false` with a `why`); a case fails over its baseline (`_test/perf-baseline.json`) by 4 bytes a
  step or a tenth of a larger baseline; the report (`perf-reporter.ts`) lists smaller moves; `--update` rewrites it.
- **Twins**: every Node `x.heap.ts` has `x.quick.ts` beside it, sharing its cases (`x-cases.ts` or the guard's own
  fixture, no copies); kit guards cycle a case's frames in one step (`kit-frames.ts`) or measure each frame. The spec
  check enforces the twin, the `quickCase` call and each opt-out's `why`. The browser's guards stay deep-only.
- **The deep tier**: timed warm-ups are gone. `settleFrame` reads exact windows until they are quiet and compile
  nothing over 6,000 steps (capped by count); `frameScaling` warms, settles and measures; `cheapScaling` first runs
  60,000 steps for a step of tens of microseconds (the rings, overlays, small draws), whose climb has quiet spans of
  thousands of steps. `wrapped-lanes.heap.ts` lost its hover-first order: "cases share V8's feedback" was the climb to
  TurboFan (16 to 130 bytes a step for 13,000 to 60,000 steps under V8's defaults, longer after other draws), which a
  timed warm-up ended inside under load.
- **§4.2's cost**: `LettersState.observed`, a mask per unit with `observed` ranges built once per prepared state; the
  letters' cells, their zoomed-out classes and the gray text read it (`unobservedAt`); a test with several ranges.
- **Found by the quick tier and fixed**: a wrapped track's hit boxed `lineAtX`'s x (16 bytes a hit on text; also
  a select drag's follow); ranges' and pairs' hits passed the point to `hitWindow` (now `hitReach` and `hitWindowAt`,
  exported for renderers); `boundaryAt` took the lane height as an argument; the ring window's drag passed the
  circumference to `dragTo` (17 bytes a move); `runXInto` inlines `endOf`. Found by the deep tier: the dot plot's
  drag passed the pointer's x and y to `follow`, 18 bytes a move under V8's defaults and none under the quick flags
  (the quick tier's thresholds can change what V8 optimizes, so the deep tier stays the check on the defaults); it reads
  the input's record now. Test artifacts, fixed in the cases: hits
  called `hitTest(section, x, y)` with fractional coordinates from test code (fingerprint 8 bytes a hit, layout 0.5);
  cases call `hitTestAt`.
- **Not fixed** (`todo.md`, Costs): a frame at a fractional DPR allocates 16 bytes (the overlay pass; most likely its
  `dpr` argument, a backend API change); a frame allocates in V8's middle tier for its first 15,000 to 60,000 frames.
  By design and baselined: the hover label's strings (575 to 1,032 bytes a hover change), a press-drag-release cycle's
  records (48, select 1,952), a cursor key's 80, a SeqQL run's result (3 to 42 KB), a select drag's autoscroll restart
  (about 4 bytes a step averaged).
- **New quick cases** (§4.2 and this PR's paths; `viewer/_test/perf-pass.quick.ts`): (a) unobserved cells on a linear
  section over 40 observed ranges; (b) unobserved text, tinted and plain, over underline and strip lanes; (c) the
  unobserved context-menu storm, the streams unobserved and after an unsubscribe (`input.quick.ts`); (d) a centered
  `wrappedMaxWidth` column, a select drag and hits; chips (hover, scroll); full text tracks (line, area, bars, points,
  ranges); strip lanes (`wrapped-lanes.quick.ts`). All read 0 bytes a step with their controls at +88.

| Tier (2026-10-15, Node 24, M4)            | Cases                                   | Wall time                                              | Result                         |
| ----------------------------------------- | --------------------------------------- | ------------------------------------------------------ | ------------------------------ |
| `test:perf` viewer                         | 451 (26 above 0.5 bytes a step, by design or `todo.md`) | 53 s (load 5 to 10); 64 to 68 s (load 10 to 41)        | pass, no move over two runs    |
| `test:perf` core                           | 97 (52 above 0.5: builds, reads, SeqQL runs) | 33 to 36 s                                        | pass under load 41             |
| `test:perf` nightingale, molstar-poc       | 10, 1                                   | 1 to 2 s, 2 s                                          | pass                           |
| `test:heap` viewer (every guard file)      | 70 files, 323 tests                     | 25 min (was about 4 min with timed warm-ups)           | pass after the dot plot's fix  |
| `test:heap` core, nightingale, molstar-poc | 14, 4, 1 files                          | 92 s, 5 s, 4 s                                         | pass                           |
| `test:slow` viewer                         | 27 tests                                | 14 s                                                   | pass                           |
| fast tiers, Seq\* and prism/viewer         | 337 files, 3,687 tests                  |                                                        | pass                           |

The deep viewer tier is slower: its warm-ups now wait for the steady state by count (`cheapScaling`'s 60,000 steps,
the settle's 6,000-step quiet span, caps of 100,000 for a step that never settles), where the timed ones measured
inside the climb under load. The quick tier is the edit loop's; the deep one runs at checkpoints.

## 5. The specs

Seq\*'s spec says what, never how; each change lands with the step that makes it.

| Spec                                   | Change                                                                                                   |
| -------------------------------------- | -------------------------------------------------------------------------------------------------------- |
| `_spec/viewer.md`, Render props        | the wrapped geometry gets its own heading ("The text view"): lanes per layer kind (§2.2), no lane without a drawing, `wrappedStripHeight` |
| `_spec/viewer.md`, Extension points    | a renderer's `geometries` decides whether it takes a lane in text; `lineHeights`                          |
| `_spec/viewer.md`, Events              | `ctx.interaction` (streams, payload, scoping) and what stays on `ctx.events`                                   |
| `_spec/viewer.md`, Behaviors           | behaviors read `ctx.input.events` (claims); bindings with descriptions                                           |
| `_spec/viewer.md`, Performance         | the streams' rule (no cost unobserved); the new guards (strip lanes, the round-2 options)                                      |
| `_spec/interaction.md`, Input events   | the public streams vs the behaviors' events (observe-only, after the claims), the chrome's parts                              |
| `_spec/interaction.md`, The track menu | `interaction.menu` replaces `events.menu`                                                                     |
| `_spec/interaction.md`, The DOM element | `data-seqv-part`, what the `seq-*` events are fed from                                                   |
| `_spec/interaction.md`, Hover          | a strip lane's hit                                                                                       |
| `_spec/reference/viewer.md`            | regenerated (`viewer/_test/reference.test.ts`): config, behaviors' bindings                               |

Outside Seq\*'s spec: Seq\*'s `AGENTS.md` (public streams may allocate, never when unobserved), prism-sequence.md §9
(round 4), plan.md §4 (the three rows leave), todo.md (the three entries leave), `.v6/plans/changelog.md`. Done already: the Mol\* import restriction
lifted (`biome.json`, `core/_test/layers.test.ts`, Seq\*'s `AGENTS.md`, `_spec/README.md`, `libraries.md`,
`model.md`).
