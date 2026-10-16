# The viewer: the host contract

What a host can rely on when it embeds `@molstar/seq-viewer` and extends it. Commands, behaviors, config items and the
built-in registry entries are in [reference/viewer.md](reference/viewer.md), generated from the code. What a user does
and sees is [interaction.md](interaction.md).

## The viewer spec

A viewer is made from a **spec**, in the shape of Mol\* v6's `PluginSpec`: registry entries, behaviors with params,
config, and initial render props. `createSeqViewer(spec, ports)` reads it once and never mutates it: it builds the
managers, registers `spec.registry` in order, then creates `spec.behaviors` in order. `DefaultSeqViewerSpec()` returns a
new spec each call (`DefaultRegistry` and the default behaviors); a host composes its own from it:

```ts
const spec = { ...DefaultSeqViewerSpec(), registry: [DefaultRegistry, InteractionsViewerRegistry],
  behaviors: [...DefaultSeqViewerSpec().behaviors, SeqViewerSpec.Behavior(MyLink, { select: false })] };
```

The spec is required: nothing in the viewer imports the default composition (`default-spec.ts`, as Mol\*'s
`DefaultPluginSpec`), and its behaviors and entries are built in pure calls, so a host that lists its own parts bundles
only those (a viewer made from an empty spec is about 720 kB minified, the default about 1,040 kB; `todo.md` "No slim
build yet" has what the floor still holds).

## The context

`SeqViewerContext` holds all state of one viewer: the document and its graph, the managers, the layout, the input
events. It needs no DOM and runs in Node; every behavior is testable headless.

- One context shows one document at a time; contexts share nothing that outlives a call (module scratch serves one
  synchronous draw, and a module cache keys on all its value depends on, never on one context being alone), so a page
  holds several and scroll groups may link them.
- **Every timer goes through the scheduler port**, so tests drive time (`advance(ms)`).
- `setActive(false)` pauses a covered or hidden viewer: no events, no frames, no keys; an animation jumps to its end and
  the player pauses.
- `dispose()` releases everything in order (behaviors, joins, timers, graph work, streams, pending commands, states, the
  registry). A disposed context throws on use; a mounted surface sees its states complete and disposes itself.

## Ports

A host gives the context its I/O as ports: `loaders` (bytes and streams, `seqspec:sources#the-loader`), `document`
(where proposed patches go: The document), `scheduler` (timers and frames), `text` (a text measurer) and `log`.
Defaults: loaders refuse every request; the document port commits to the context's own base; the scheduler is a manual
one that moves only on `advance(ms)` (a host passes `rafScheduler` for the browser's frames); the measurer is a
monospace estimate (0.6 em a character) until a surface sets the canvas one. `renderSvg` measures with its own font
metrics.

## State

Every piece of state is subscribable:

- **UI-rate state** is a `SeqState<T>`: an rxjs `Observable` with a synchronous `value` that replays its current value
  to a new subscriber and emits only on change. A set made while a state's subscribers are being called is queued and
  delivered after, in order. Sources never error: errors are values.
- **Pointer- and frame-rate values** are `Live<T>`: one record mutated in place by its manager, with a `version` that
  bumps on change. Read `current` when `version` emits; never keep it.
- **Events** (`SeqEvent<T>`) are happenings, not state: nothing replays; subscribers are called in priority order.
  Pointer-rate events reuse one record; edge-rate ones (`activate`, `patch`, `issue`, `log`) emit a new one.
- A host's subscriber that throws: on a state, reported as rxjs reports it; on an event, logged, and the others still
  run. The context guards its own fan-out (`seqv.state.error`; input handlers, `seqv.behavior.error`), so the viewer's
  own code never crashes a Node host.
- **Origins** travel with changes (`pointer`, `keyboard`, `api`, `graph`, `snapshot`, `behavior`, `host:<tag>`), so a
  host skips its own echoes. A size change (a resize, a DPR, the data area's inset) is `graph`: what its clamp moves is
  no one's move, so a linked section is set again from the one that led last.
- Managers and hosts name a section by its ref; hot records carry its index in the current view.

The context's states, by manager: `document`, `graph`, `views`, `snapshots`, `viewports`, `layout`, `pointer`, `hover`,
`highlight`, `selection`, `focus`, `systems`, `results`, `labels`, `ui`, `tasks`, `streams`, `issues`, `history`,
`behaviors`, `registry`, `renderProps`, `skin`.

## Events

Input is `ctx.interaction`: streams of fresh payloads (`hover`, `press`, `click`, `drag`, `menu`, `wheel`, `key`,
`leave`, `activate`), fed after the behaviors, observe-only, scoped by `on`, built only while observed
(`seqspec:interaction#streams`); the track menu is its `menu` (`seqspec:interaction#the-track-menu`). The behaviors'
base events, which they claim, are `ctx.input.events` (`seqspec:interaction#input-events`). `ctx.events` keeps what is
not input: `patch` (every proposal), `issue` (every new issue), `log`. The DOM element adds its own
(`seqspec:interaction#the-dom-element`).

## The document

The shown document is layered, each layer applied by patch so untouched subtrees keep their identity: the **base**, then
the current **snapshot**, then **local** patches, then a **preview**.

- `setDocument(next)` replaces the base by rebasing: the patch from the old base to `next` is applied, so a document
  that changed one node reruns one node. Snapshots, metadata and revision come from `next`.
- `applyPatch(patch, { layer })` applies a patch on a layer, all or nothing; a preview that does not validate is not
  drawn.
- **Patches out.** A UI action makes a patch and proposes it; the **document port** routes it: `commit` (the default
  port applies it to the base; a host's forwards it and later passes the new document back), `local` (kept in this
  viewer only; the owner's document is untouched), or `refuse` (`seqv.port.refused`).
- The default port routes a patch `local` **under a snapshot**, tagged with that snapshot (it applies only while that
  snapshot is shown), and when it is **masked**: it targets a node only a snapshot or a local patch adds, or a field a
  local patch sets. A base or local change drops a local patch made over the snapshot shown that no longer applies
  (`seqv.patch.dropped`); one made over another snapshot is kept and reported (`seqv.patch.skipped`). A snapshot change
  never drops one.
- A committed patch carries the base's revision as its `base`; a stale one is rejected.

## Data mode

A viewer is usable without a document: `setData` takes a `SeqSet` or the model's build inputs, and views written with
the builder against the set's keys, or a default view (a ruler, the letters, a track per annotation by its geometry).
Refs are the set's keys. Validation runs against the incoming data. Everything downstream (patches to views, snapshots,
managers) is as with a document.

## Updates and status

A document becomes a graph of nodes keyed by ref, whose edges are exactly the format's references
(`seqspec:format#what-a-change-affects`): a node's refs as core reads them (`seqspec:format#nodes-and-refs`), a query's
names, a theme's `column` and `profile` and a table's targets among them. An update reruns only the nodes whose own
content, or an input's result, or what they read in the set, changed; results that come out equal stop there. The set it
assembles is the one `compile` gives the same document, and so are its issues, at the document's paths. A view switch
resolves only that view; switching back is free. Registrations and the skin are inputs too: a light/dark switch keeps
every prepared state.

- **Status.** A node that fails is `error` with its issues; its dependents are `unavailable` with a cause; an
  unavailable layer is omitted (`svs.render.omitted`) and its section still draws. A missing provider, format or service
  is not a failure: the node is unavailable with `svs.provider.unavailable`. Nothing cascades.
- A node keeps at most `Issues.PerNode` issues, the most severe first; the rest are counted in one
  `seqv.issues.truncated`, by severity. The issues from no node (registration, behaviors, the engine) are kept up to a
  bound the same way: past it the least severe goes first, so a flood of warnings never hides an error.
- One update runs at a time; inputs arriving meanwhile are folded (the newest document wins, patches compose). Async
  nodes are aborted when superseded. An update always settles, and `graph.assembled` emits once, from which the managers
  remap.
- Held loci (selection, highlights, an external hover) follow their set: within a lineage they are rebound, across
  lineages remapped by key (`seqspec:model#identity-and-lineage`).

## Snapshots and playback

`snapshots.apply(ref)` sets the snapshot layer; `next`, `prev`, `play`, `stop`. Entering one follows the format's rule
(`seqspec:format#snapshots`): its focus changes are animated over its `transition_ms` by the `focusViewport` behavior.
`play` stays on each snapshot `duration_ms` after its transition and update end. An invalid snapshot is listed and
skipped. On load and on entering a snapshot, the next and previous snapshots' focus regions are requested as prefetch
hints.

## History

The context records `{patch, inverse}` for every patch it commits itself (with a port whose `history` is true) and for
local patches. `undo` applies the inverse; a new patch clears redo; an inverse that no longer applies is refused
(`seqv.history.refused`). When a host echoes the context's own commits back by `setDocument`, history is kept; any other
document clears it. A host that owns undo (Stories) sets `history: false`. Interaction state (viewports, selection,
expansion, params) is never in history.

## Commands

A **command** is the core's only request mechanism from outside: an id (`seqv:focus`, an extension's `<ns>:<name>`), a
params table, and its subscribers, as Mol\*'s `PluginCommand` with stable ids.

- `ctx.commands.dispatch(cmd, params)` from code; `ctx.commands.run(invocation, origin)` for plain data
  (`seqspec:format#descriptions-and-invocations`): the id looked up in the registry, params checked by its table.
  Unknown id or invalid params resolve `false` and run nothing.
- A command with no subscriber does nothing. The built-ins are handled by the `commands` behavior; a host's own
  subscriber runs in addition, after it; to replace what a built-in does, leave the behavior out.
- Subscribers run in order and are awaited; the first rejection rejects the dispatch (`seqv.command.error`).
- Targets are refs (with range args), loci (from code only) or `Coordinates` (runs in named systems). Refs of several
  kinds act on each kind (selected, lit, framed together), never on the first alone. Range args outside a ref's domain
  (index args, or labels a numbering resolves; a complex's: outside every member) are refused as a document's focus is,
  never clipped: nothing is done, `seqv.target.range` names the domain, `run` resolves false and `dispatch` rejects.
- Commands are edge-rate: a dispatch allocates its promise. Pointer-rate paths call managers directly.
- `seqv:focus` frames a target (a focus, the section's state); `seqv:reveal` only brings it into view, changing no
  selection, focus or zoom (`ctx.reveal`; `seqspec:interaction#focus`).

## Behaviors

**Behaviors** are policy over the context's mechanism: what a hover, click, drag, wheel or key does, how linked sections
follow, what a command does. Each is a provider with a name, a category, a priority and params (plain data), added and
removed at runtime (`ctx.behaviors.add`, `get(name).setParams`).

- Removing a behavior releases everything it set (its channels, its hover, a pending drag, cursors), so nothing stays
  lit or stuck. It removes the policy only: clicks still emit `click`.
- A preset may turn a behavior off per section, or override its params (`forSection`).
- A behavior that throws while being created, registered or updated is disposed and reported (`seqv.behavior.error`);
  the viewer goes on without it.
- Bindings match modifiers exactly, so Ctrl/⌘ with `+` stays the browser's zoom. A behavior describes its binding
  params, and `ctx.behaviors.bindings()` lists them from the current params as Mol\*'s `Binding` (its description with
  `${triggers}`), for a help panel.
- Behaviors read the base events (`ctx.input.events`) and claim gestures there; the public streams (`ctx.interaction`)
  see what they did.
- The interaction rules are behaviors' params, in the viewer spec (`SeqViewerSpec`, as Mol\*'s PluginSpec): a host
  composes, overrides, replaces or removes them there, or at runtime (`setParams`), and a preset per section. Input
  delivery (the engine, the dispatcher) only delivers: a press followed anywhere, its point clamped, a release never
  seen synthesized; without `select` a click selects nothing. The selection's rules are `select`'s: what Shift does on
  positions with a click or a drag (`shiftPositions`: `extend`, `toggle`) and on items (`shiftItems`), Ctrl/⌘
  (`toggle`), whether a click waits where a double click is claimed (`doubleClickWait`), what a double click is
  (`repeatClicks`), the drag regions, and the scroll past an edge (`autoscroll`, `autoscrollMin`, `autoscrollRate`,
  `autoscrollMax`). A behavior whose double click means something at a point claims the `doubleClick` query there
  (`doubleClickZoom`, `navigator`, `headers`), so `select` asks the behaviors, never a list.
- Some built-in behaviors are not in `DefaultSeqViewerSpec` (`SeqViewerOptionalBehaviors`; the generated reference lists
  them apart): a host adds them. `revealHighlight`, a linked view's highlight brought into view once it rests, is one
  (`seqspec:interaction#highlight`).

## Actions

**Actions** are what the UI shows as menu items and switches and what a host maps into its own commands: a name, a
title, a `when` and a `run`, mostly proposing a patch (`track.hide`, `section.numbering`, `bookmark.add`) or changing UI
state (`view.switch`, `group.toggle`). `ctx.actions.run`, or the `seqv:runAction` command. The origin reaches the patch.
`bookmark.add` keeps loci of the shown set as rows of a section's inline `ranges` set (`<section>.bookmarks`, or the
target's `data`), with an optional `label`, `color` (a `color` column on the rows, shown by a layer's `field_color`) and
`line_style` (one of the line styles, a `line_style` column, shown by a layer's `field_line_style`; no column is written
for either without it); on a set the document already holds it is one `set` of the rows, so one undo step.

## Config

Typed items with defaults (`SeqViewerConfig.Viewport.SettleMs`), read by their managers at use: `ctx.config.get(item)`,
`set(item, value)`, or `spec.config`. A UI reads its items when it mounts; the one write is the first surface's:
`Viewport.ReducedMotion`, when unset, follows the page's reduced-motion setting.

## Registration

`ctx.register(entry)` adds registry entries: pure data, grouped by kind, as Mol\* v6's.

- One call registers every entry's format half (`svs`) first, then the viewer's kinds in a fixed order (node kinds,
  renderers, columns, themes, navigators, presets, labels, commands, actions), entries in order within each.
- Conflicts are checked before any change (a taken name, a renderer whose layer type is registered nowhere, a theme
  serving a kind its format entry does not); on a conflict nothing is registered.
- It returns an idempotent undo. There is no replace: a host that replaces a built-in registers a copy of
  `DefaultRegistry` without it.
- A registration is a graph input: only nodes whose entry changed rerun. A lookup that misses returns nothing; its user
  takes the declared fallback (an unknown preset, theme or renderer geometry is reported, an unknown data-column type is
  drawn as `text`).

## Extension points

What the core promises each kind of entry, and what the entry must do. The built-ins of each kind are in the reference.

- **Node kinds** wrap the format package's stage functions, one to one. A host may replace one (a caching loader); it
  never adds a document kind.
- **Renderers** draw one layer type each, as a Mol\* representation draws one geometry; every built-in is one, with no
  private access. A renderer prepares headlessly and purely once per input change, hits headlessly, and draws through
  the renderer API (`/draw`), which is DOM-free, so the same entries serve a canvas, SVG in Node and tests. It declares
  the geometries it draws (`linear`, `wrapped`, `circular`); without `circular`, a section on rings draws its extents as
  bars and reports it. Its `wrapped` decides whether it takes a lane in the text view (The text view): with it, its band
  per line is its `lineHeights` (from each line's runs: the underline and strip lanes), else a row (`rowHeight`);
  without it, its layers take their fallback in text and no band. It must honor the kernel type's semantics (placement,
  split rows, channels, fallbacks). It declares `content: false` when it never reads residues, so a residue edit keeps
  its prepared state. It declares `zooms` when its free height follows the row zoom (its height then takes the row
  scale), so the zoom moves its track with the rows without a new layout. Its draws allocate nothing per frame. It may
  give the band a data hover lights on its track (`rowHover`: its rows, or none), else the track's box. It declares
  `stable` when its draw of a strip of the view (a pan's uncovered part) equals its draw of the whole view there, so a
  pan copies the rest: a decision over the view (items or density) reads the whole view's range, not the strip's. Its
  other optional parts: `tiny`, its draw below its `minHeight` (else bars); `across`, what it draws across the tracks
  its layer spans, after the entries (a frame's outline: Frames); `miniature: false`, left out of a navigator's
  miniature; `overview`, its rows in the vertical overview; `ringLabels`, its items' labels outside the rings; `cache`,
  a store kept per prepared state and skin (made and disposed by the engine); `bands`, a multi-band layer's band labels;
  `header` and `gutter`, its rows' header and gutter (a rows layer's labels and tree, its line numbers). Renderers run
  the conformance kit (`@molstar/seq-viewer/testing`) as the built-ins do, its edge geometries (`EdgeKinds`: reversed, a
  cyclic origin, stacked, filtered, copies) included: drawn and hit back in each.
- **Data columns** draw a `rows` annotation's field per alignment row; an unknown column type is `text`.
- **Theme providers** are the drawing half of a format theme entry of the same name; one provider serves every channel
  kind its entry declares. A renderer fills a reused theme location (layer, block, unit, index, item, letter, value,
  state); a theme reads it and returns a packed color, a size or a glyph.
- **Navigators** draw the whole axis and the view's window; their gestures are the `navigator` behavior's. A navigator
  whose view can neither zoom nor pan is drawn idle and its gestures do nothing (interaction.md, Viewports). A
  navigator's parts (`NavigatorParts`) are an axis (the window over a strip, or a ruler's rail) and the tracks in
  miniature (each layer drawn small over the whole axis, `seqv:map`); the built-ins are in the reference, and a host
  registers its own, such as `ideogramNavigator(bands)`, a chromosome's cytobands under the window. A section with rows
  gets a vertical overview (the rows in miniature down a column, its window the rows in view) or a scrollbar, by its
  render props (`verticalNav`).
- **Presets** are a section's geometry, render props, navigator and per-behavior overrides, picked by the UI's runtime
  style, else the section's `style` (through `Presets.Alias`; an unmapped key is a preset name), else `Presets.Default`;
  an unknown one is the default (`svs.render.preset`).
- **Label providers** write structured lines for the hovered element into a reused sink, by priority.
- **Commands** and **actions** as above; **behaviors** are not registry entries but are added the same way by a spec or
  at runtime.

An extension ships one registry entry (format and viewer halves) and, when it has policy, behaviors.

## Labels

`ctx.labels` computes the hovered element's label by calling every applicable label provider with the hover and a reused
sink, grouping equal lines (a merged copy group's units become `× n`). One label is built per hovered element; when to
compute and hide it is the `hoverLabel` behavior's.

- **The codon.** On a block a streamed layer reaches through a path whose map takes each of the block's elements to
  three positions of its system (a codon step from the nucleotides, `seqspec:model#codon-maps`), the hovered element's
  label gives those positions' letters in reading order (`codon`; an unknown one `·`, a partial codon's known ones),
  split by `|` where they are not adjacent (an intron; on a cyclic system the origin is no gap), and complemented by the
  `complement` letter map (`complement_rna` for RNA) on the codon step's minus strand; then the positions in the
  system's numbering, the unit's own first, in reading order (high to low on the minus strand) with the strand
  (`codon at`: `chr2 1001–1003 (+)`). A lazy molecule's letters come from its stream, a view's letters layer of it
  first: those not held are asked after the view's regions, each label's replacing the last's, and the line says
  `loading…` until they arrive; `failed` (Retry asks again) or `not served` when they cannot.

## Host coordinate systems

The viewer talks to its host in coordinates, and a host may own coordinate systems (Mol\*'s chains).

- **`ctx.systems`** is core's system registry over the shown set (`seqspec:libraries#the-system-registry`: its own keys
  `seqstar:<ref>` and their descriptors) with a reactive `version`, bumped by every registration.
- **Events describe their systems.** What the viewer reports outward gives extents as `Coordinates` with a plain-data
  descriptor per system (key, kind, domain, refs, a unit's `mvs` fields, numberings, a complex's members), so a host
  needs no Seq\* types.
- **Requests in any registered system.** `hover.setLoci`, `highlight.set`, the selection verbs and `focus.fromLoci` take
  `Coordinates` (`seqspec:libraries#coordinates`). A host registers an external system with its maps
  (`ctx.systems.register`); it joins the view only through them. A structure registers as one complex system whose
  members are its chains'.
- Coordinates naming a key not registered survive as keys in a hover or a highlight, light nothing, and light again when
  it returns. The selection holds loci: its verbs convert Coordinates at once through the systems registered then, so a
  key not registered selects nothing (`set` clears the selection; register the system first).

## The render backend

The composer and every renderer draw through a `RenderBackend`, an interface at the level of what they mean: batches of
rectangles, sectors, paths, glyphs and text, images, passes, clips and copies. Backends: Canvas2D (live), SVG (`/svg`,
platform-free: the same code exports from a page and renders figures in Node), and recording (tests).

- Batches are struct-of-arrays and grow-only; colors are packed integers or skin tokens, alpha its own byte. A frame
  allocates nothing. A path's stroke may be dashed: on and off lengths and a phase, device px, so a pattern is fixed to
  the content (a range's `dashed` and `dotted` line styles); Canvas2D sets its line dash, SVG `stroke-dasharray`.
- **Surfaces** (`mountSection`, one per section) size their canvases inside a frame only, so the page never paints a
  blank canvas, and follow the host's size and device pixel ratio. Each pass canvas carries `data-pass` (`main`,
  `overlay`, `nav`, `nav-overlay`, `vnav`, `vnav-overlay`), marks for a host's CSS: the vertical overview and scrollbar
  are the `vnav` ones.
- **Per-position text** (a block's letters, another polymer's residues, three-letter codes, atom symbols) is set per
  row: one font size and one form for the whole row, from its cells' width and never larger than the row allows nor than
  the cap (the render prop `maxLetterPx`; by default the header labels' size, `labelFontPx`), so zoomed in past it the
  text keeps its size; where the row's widest text does not fit, none is drawn and the cells stay. A row never shows
  some columns' text and not others'. Text never exceeds the part of its residue it is centered on: a residue's
  three-letter code needs a few px of margin inside the narrowest such part in view (a codon an intron splits has a
  smaller one), else the row shows one letter. The text layout's letters have a size of their own (`wrappedFontPx`), the
  cap's by default. A row that draws no cell colors (text only) marks where text does not fit which positions hold a
  residue: a dot in each cell while a cell holds one (about 40% of its width, at most 9 px), then a thin bar over each
  run of residues, so a deletion stays a break; both in the ink at half strength. A row with cell colors shows its
  cells, and marks nothing.
- **A boundary over letters.** Where a `ranges` layer draws over its own track's letters at their size (a design's row
  with its insertions), a boundary is a small triangle pinned to the row's top edge, pointing down at the edge between
  the two cells, and a hairline down that edge, in the item's color and clear of the letters' ink at every zoom; with no
  room above the text for a triangle 3 px tall, the hairline alone (broken around the ink where the letters leave it no
  room beside them). Where no text is drawn, the boundary is its glyph on its lane. At a block's first or last boundary
  a boundary's mark sits inside the data area, never cut. It hits within a few px of its edge over the whole row; beside
  it, the residues.
- `renderSvg` draws a whole view or section headless; an export never commits a viewport. A document's streamed regions
  are awaited at the size the export draws (each section its content's height), so a tall streamed alignment exports
  whole; `timeout` bounds the wait (a static source's load too), drawing what stands and saying so (`unsettled`, with
  the stream requests in flight; also when heights still move after a few rounds); `issues` gives the document's issues
  as it was drawn.
- **The `seqv-render` CLI** renders a document (JSON or MessagePack) as `renderSvg` does:
  `seqv-render <document> [--view <ref>] [--sections a,b] [--width 1200] [--scale 1] [--overlay] [--nav] [--background] [--title <text>] [--gap 8] [--root <dir>] [--timeout 120000] [--png] [-o <file>]`,
  SVG to stdout or `-o`, PNG with `--png` (or `-o` ending in `.png`, through a local Chrome). Sources are files under
  `--root` (default the document's folder), never outside it (by path, links resolved) and never the network; it streams
  nothing, so a streamed source is drawn without its data (`svs.provider.unavailable`). Exit codes: 0 written; 1 nothing
  written (the document does not load, the render or the PNG failed); 2 usage or an unreadable file; 3 written short of
  data (a source that could not be read, a warning or an error of the document's, or sources still loading after
  `--timeout` ms), each said on stderr.
- Rings (`circular`) are the same frame unrolled and mapped onto the circle, so every renderer, mark and hit runs as on
  a line; a track's layers overlay on its ring as they do linearly. Items' labels are set outside the rings, upright,
  each with a leader, and hit as their items; the ring's window (`seqspec:interaction#linked-sections`) is drawn on the
  overlay.

## Text measuring

Labels, ticks, data column text and wrapped lines need widths, and the layout is headless. Text is **prepared** once per
font into a bounded cache (`Text.CacheEntries`); layout is arithmetic over the prepared arrays, with no measuring and no
allocation (pretext's approach). `ctx.setTextMeasurer` replaces the port; the first mounted surface sets the canvas
measurer unless the host opts out.

## Render props

The layout's and engine's sizes and switches (`ctx.renderProps`, `spec.render`, `setRenderProps`), the analog of Mol\*'s
canvas3d props: header width, row height, gaps, fonts, ring geometry, navigation, LOD thresholds. A section's props are
its preset's over the context's, its inline style's keys over both, the host's `Render.Clamp` last. All are JSON base
types. A dot plot's height (its `aspect_ratio`, `seqspec:format#layer-types`) combines three rules in order: the data
area's width ÷ the ratio, capped by the box's rows' room (its visible height under what stays pinned, the navigator and
a pinned ruler) less a margin of a tenth of it, at least two rows, so a plot that fits shows whole at more than one
scroll position (720 px in a section sized to its content), and raised to 120 px; then at least the width ÷ (ratio ×
`dotPlotMaxStretch`, 2.5 by default), which wins over the cap. It follows that target, so a resize lands exactly, except
a bounce: a shrink of less than 32 px back to the target it just grew from (within the cap) keeps the taller height
while the target stays there, so a section sized to its content whose plot's shrink hides the page's scrollbar settles
after one round instead of flickering. Rows have no borders by default (`rowBorders`): a document turns them on per
section (`style.row_borders`), a host in its render props; an SVG export draws as the screen does, and
`renderSvg(doc, { render: { rowBorders: true } })` turns them on for a figure. An expanded split track's rows have quiet
guides regardless: a hairline in the skin's grid tone at the bottom edge of each row but the last, under the data, drawn
only where a row is at least 4 device px tall and shifting nothing (a row border, where there is one, takes its place).
The text view's sizes and spacing are render props too (The text view).

## The text view

A `text` layout's sections draw the `wrapped` geometry: the positions as lines of text, each line a band per track under
it. The text sits `wrappedInset` px (8, the skin's `--seq-canvas-edge` default) inside the canvas on every side, its
line numbers' gutter within; hits, the selection, the drag's nearest positions and its scroll all read that geometry.
The layout reads no CSS (an SVG export draws it the same): a host that changes `--seq-canvas-edge` sets
`render.wrappedInset` to match. Its spacing is render props too: `wrappedGroup` positions a group, `wrappedGroupGap` px
between groups, `wrappedGap` between lines, `wrappedTrackGap` between a line's letters and a track's lanes,
`wrappedLaneHeight` an underline lane, `wrappedStripHeight` a strip lane (Lanes, below), `wrappedMaxWidth` the longest
line in px (0, the default: the whole data area), `wrappedColors` the letters themselves in the layer's color theme
(false, the default: the ink), `wrappedCells` the letters and an alignment's rows on cells in their layer's color theme,
each letter in the ink that reads on its cell (false, the default; the style's `text_cells`; it wins over
`wrappedColors`: a figure's shading is legible only on its cells), `unobservedGray` a unit's unobserved positions
(`observed` ranges that leave them out) in the skin's `unobserved` color, opaque, in the letters' cells at every zoom
(zoomed out, a column's share of them as their own class) and in text, plain or colored (false, the default),
`wrappedAlign` (`start`, the default, or `center`) where the lines sit: centered, the block is as wide as its widest
line, centered in the data area the section ends with (the vertical overview's room only where the text scrolls), and
the margins either side belong to the gutters, hits and the selection reading the same geometry. A line holds whole
groups, and a short last group goes on a line of its own where a full group would not fit. The `wrapped` preset spaces
the text to read; `compact` is the dense one (half the group gap, lines and lanes touching, lower lanes), still in
groups of 10.

**Lanes.** Under each line, every track that draws there takes a band, in track order (under the group labels' band,
`aboveEnd`); the layers of one track share its band, the tallest of theirs:

| Layer kind                                        | Its band on a line                                                                       |
| ------------------------------------------------- | ---------------------------------------------------------------------------------------- |
| `letters`, `emphasis`                             | the letters: 1.6 × `wrappedFontPx` (`emphasis` draws its track's letters, emphasized)    |
| `rows`                                            | a text row (the letters') per row: the rows as text (Rows in text, below)                |
| `frame`                                           | its titles' row on a line where a titled item starts (Frames, below), else none          |
| `ranges`, `points`, `secondary`                   | underline lanes: `wrappedTrackGap`, then `wrappedLaneHeight` per lane with an item there |
| `values`, `color`                                 | a strip lane: `wrappedTrackGap`, then `wrappedStripHeight` (3 px; `compact` 2)           |
| `ruler`, `dot_plot`, the types text does not draw | none                                                                                     |
| but `letters`, with a `height`, under `full`      | the track drawn in full: `height` text rows (below)                                      |

- **Tracks drawn in full** (opt-in: the render prop `wrappedTracks`, the style's `wrapped_tracks`; `lanes`, the default
  and the `wrapped` and `compact` presets', keeps every track a thin lane whatever its height; `full` draws them in
  full): a track with a `height` (and no `letters`: the letters are the text) takes that many of the text's rows on
  every line (a text row is the letters' band, 1.6 × `wrappedFontPx`, so a figure's tracks scale with its text; rows, as
  a track's height is in tracks), and its layers draw their track forms over each line, as over a linear view of the
  line's positions: `values` its line, area or bars, on the layer's one scale for every line, a line broken at each
  line's end (no stub: the next line starts at its first position), from a line's (or a run's after a missing value)
  first cell's left edge to its last's right, the end values held over their half cells as the area's cells are;
  `points` its glyphs, `ranges` its items with their labels, each on its item's part on the line (centered in it, cut
  with an ellipsis, left out under three graphemes, never past the line's end), `secondary` its shapes, `color` its
  cells; `variants` take their fallback, points. An extension renderer that declares `wrapped` draws its linear form
  there (`AxisView.trackForm`, the line's positions as `visLo`/`visHi` and its advance as `cell`). Two-dimensional kinds
  (dot plots, pairs, trees) are not drawn in text; an alignment's rows are drawn as text (Rows in text). On a chip's
  line a track keeps its lane. Hits there are by drawn geometry, as in tracks.
- **A line's block** is one unit: its group labels' band topmost, then its tracks in track order (one ordered before the
  letters between the group labels and the letters), each band `wrappedTrackGap` below the one before; the next line's
  block `wrappedGap` below its last band. A block reads as one where `wrappedGap` is clearly the larger: the `figure`
  preset (text with its tracks in full, the axes and the labels on) takes 3 and 18 px. Under `full`, a track of features
  (`ranges`, `points`, `secondary`) with nothing on a line takes no band on that line (`wrappedCollapseEmpty`, true by
  default: an empty band is room that separates what belongs together); `values` and `color` keep theirs (a gap in their
  data is data, and one scale reads across lines).
- **Axes and labels** (off by default): `wrappedTrackAxes` (the style's `track_axes`: a small y axis left of each band,
  its domain's ends, in 30 px the text leaves left of its lines) and `wrappedTrackLabels` (`track_labels`: the track's
  label in a row of its own above its band, `EntryKind.trackLabel`, on the first line the track has a band: never over
  its data; the left of the lines holds the axes, too narrow for a label); an SVG export draws them alike.
- **Underline lanes** are compacted per line to the lanes with an item on it, each underline 2 px less than its lane. An
  item's line style (`seqspec:format#themes`) draws its underline (a wave over the whole lane); a `highlight` or `box`
  item also draws behind the line's letters (the line's first `letters` track): a translucent fill or an outline over
  its cells, drawn before the letters (`LayerRenderer.behindText`), its lane keeping its own mark; where the letters sit
  on cells (`wrappedCells`) it is drawn over them, the fill translucent, so the cells do not hide it.
- **A strip lane** is on each line of a block the layer has data on: a color per position from the layer's color theme
  (a `values` layer's at its value; with none, the ink), cut at the gaps between groups; a missing value draws nothing.
- **No lane without a drawing.** A type whose renderer leaves out `wrapped` takes its fallback in text
  (`seqspec:format#layer-types`: a streamed alignment's `rows` its letters, `variants` and `pairs` points) or is omitted
  (a ruler, a dot plot), and takes no band; a track whose layers draw nothing on any line takes none.
- **Beside the letters**, on their track, a `color` or `values` layer colors the letters' cells behind them (a `values`
  layer with no color theme there draws nothing: its ink would hide them); alone on its track it is the strip. Their
  cells never change the letters' ink: a figure shades letters by their own theme (`wrappedCells`), where each letter
  reads on its cell.
- **Consistent heights**: every band's height is a render prop's (the letters' the font's), so a line's spacing reads
  the same whatever is under it. The `compact` preset sets them tighter (`wrappedTrackGap` 0, `wrappedLaneHeight` 4,
  `wrappedStripHeight` 2).
- **Chips**: a short unit (a ligand in the `position` form: one position) is a chip, as wide as its chemical component's
  label (measured in the chips' font, the letters' band's half, with 5 px either side) and not a block of its own line:
  consecutive chips flow on shared rows, the group gap apart, wrapping where the next does not fit, as words do; their
  rows follow the chains spaced as a new line is (`wrappedGap`, not the blocks' `stackGap`, then the group labels'
  band's room, empty: chips have no group labels). Each track's band on a row is the tallest of its chips'. A chip stays
  its unit's one position: its hit, hover, hover card, selection, a drag's nearest position and the links read it as any
  position (a drag still stays on the block it started in). Tracks draw a ligand as a pill over its cells, its label
  fitted; narrower than its box, a filled cell at least 2 device px wide, unlabeled, so it never vanishes.
- **Hits** (`seqspec:interaction#hover`): a point on a track's band is that track's; on an underline its item, by lane;
  on a strip, and in its gap, its layer's position, so the hover and the hover card show the value, as on a linear
  track. Lanes only add height to a line: the selection, a drag's nearest positions and its scroll read the lines alone.

### Rows in text

An alignment's `rows` layer in a text layout is a block of text rows on every line of its block, so a line reads as a
figure of the alignment (ESPript's): each row is named and numbered on every line, and its residues line up in columns.

- **The block.** A text row (the letters' band, 1.6 × `wrappedFontPx`) per display row: `order`, `filter` and
  `row_subset` apply; no tree, every clade unfolded; under `lanes` and `full` alike (a `height` on its track is ignored:
  the rows are the text). Only the rows in view are drawn. A streamed alignment's rows take their fallback in text
  (`seqspec:format#layer-types`): text is a figure of rows held.
- **Letters and cells.** Under `wrappedCells` each residue sits on a cell in the layer's color theme (a conservation
  shading: `espript`, `seqspec:format#themes`), its letter in the ink that reads on the cell; with `wrappedColors` the
  letters themselves in it; else the ink. A gap is the quiet ink and no cell (`letter_map` `{"-": "."}` draws dots).
- **Beside the lines.** Each row's label left of every line, in a column as wide as the widest label (at most
  `wrappedRowLabelMax`, 160 px; a longer one cut with an ellipsis), and right of the line the number of the row's last
  residue on it in the row's numbering (`NAME/27-137` starts at 27; none on a line where the row has no residue). The
  layout keeps the room for both.
- **Group labels** (`aboveEnd`, `inlineStart`) number the alignment's reference row (row 0 unless `reference_row` says):
  the group's last (inline: first) residue of that row, none where the group holds none of its residues, as ESPript
  numbers its first sequence.
- **The rows track's other layers** draw on the text row of the row their item targets: a range's `box` or `highlight`
  (`LayerRenderer.behindText`, with `v.row` its lane offset) over the row's cells; nothing of theirs is drawn as a lane.
  Tracks before the rows track are drawn above the block, those after it below; no track sits between its rows.
- **Hits.** On the block, y names the display row and x the column: the row's residue there (`aligned`, the row's own
  position), or on a gap the row at that column. Items of the track's other layers are hit on their row by the position,
  as on a rows track's row. The hover and the selection light columns over the line's bands, as on letters.

### Frames

A `frame` layer (`seqspec:format#layer-types`) outlines named regions across tracks: a figure's domains and indels over
every row of an alignment, in text and in tracks.

- **The outline** is drawn after the entries of each region the composer draws, in its clip (`LayerRenderer.across`,
  `v.top` and `v.height` its span), so a pan's strip, a scroll's band and an export draw what a full draw does. It runs
  from the top of `span.from`'s band to the bottom of `span.to`'s (wrapped, on each line, half `wrappedTrackGap` outside
  them; linear, over the span's scrolled entries, so it follows a row zoom without a rebuild). Where a rows band is that
  edge, the edge is the topmost or bottommost row the item targets (an item on the columns spans every row).
- **Open sides.** A side is drawn only at its item's own end: where a line break (or, in tracks, the view's edge) cuts
  the item, that side stays open and the top and bottom edges run to the line's end. The edges sit on the cells' edges,
  1.5 px, in the item's color; `dashed` and `dotted` line styles are drawn so, the pattern fixed to the content, the
  others solid.
- **Titles.** Each item's label is drawn in the frame's own track, at its left edge, on the line (or in the view) where
  the item starts, up to the next title or the line's end, cut with an ellipsis; in text the track takes a title's row
  only on such lines. The title's band hits the frame's item over its columns; the outline takes no hits from the tracks
  under it.

## Performance

No allocation per frame, per pointer event or per lookup on these paths: pointer moves and hit tests, state and event
emission, frames, hover across sections (closure and projection on cached routes), labels, marks, drag-select, theme
lookups, text layout, and SeqQL runs per position or run. Each such path is guarded, with a positive control (a guard
that is shown to fail on a planted allocation): its bytes allocated per step, read exactly against a recorded baseline,
and its heap growth over a long run once optimized. Commands, the graph and edge-rate events may allocate. The public
input streams (`ctx.interaction`) cost nothing on these paths while a stream has no subscriber (their guards run with
the streams created and none subscribed, and after a subscriber left); a subscribed stream may allocate its payload per
event.

A draw's edge-rate work (an index, a pyramid, a tick list) is bounded by what its block shows and its data holds, never
by the length of the block's coordinate system: a protein's letters on its chromosome index the protein's codons, not
the chromosome.
