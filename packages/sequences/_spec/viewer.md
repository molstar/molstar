# The viewer: the host contract

What a host can rely on when it embeds `@molstar/seq-viewer` and extends it. Commands, behaviors, config items
and the built-in registry entries are in [reference/viewer.md](reference/viewer.md), generated from the
code. What a user does and sees is [interaction.md](interaction.md).

## The viewer spec

A viewer is made from a **spec**, in the shape of Mol\* v6's `PluginSpec`: registry entries, behaviors with
params, config, and initial render props. `createSeqViewer(spec, ports)` reads it once and never mutates
it: it builds the managers, registers `spec.registry` in order, then creates `spec.behaviors` in order.
`DefaultSeqViewerSpec()` returns a new spec each call (`DefaultRegistry` and the default behaviors); a
host composes its own from it:

```ts
const spec = { ...DefaultSeqViewerSpec(), registry: [DefaultRegistry, InteractionsViewerRegistry],
  behaviors: [...DefaultSeqViewerSpec().behaviors, SeqViewerSpec.Behavior(MyLink, { select: false })] };
```

## The context

`SeqViewerContext` holds all state of one viewer: the document and its graph, the managers, the layout,
the input events. It needs no DOM and runs in Node; every behavior is testable headless.

- One context shows one document at a time; contexts share nothing that outlives a call (module scratch serves one
  synchronous draw, and a module cache keys on all its value depends on, never on one context being alone), so a
  page holds several and scroll groups may link them.
- **Every timer goes through the scheduler port**, so tests drive time (`advance(ms)`).
- `setActive(false)` pauses a covered or hidden viewer: no events, no frames, no keys; an animation jumps
  to its end and the player pauses.
- `dispose()` releases everything in order (behaviors, joins, timers, graph work, streams, pending
  commands, states, the registry). A disposed context throws on use; a mounted surface sees its states
  complete and disposes itself.

## Ports

A host gives the context its I/O as ports: `loaders` (bytes and streams, `seqspec:sources#the-loader`),
`document` (where proposed patches go: The document), `scheduler` (timers and frames), `text` (a text
measurer) and `log`. Defaults: loaders refuse every request; the document port commits to the context's
own base; the scheduler is a manual one that moves only on `advance(ms)` (a host passes
`rafScheduler` for the browser's frames); the measurer is a monospace estimate (0.6 em a character) until
a surface sets the canvas one. `renderSvg` measures with its own font metrics.

## State

Every piece of state is subscribable:

- **UI-rate state** is a `SeqState<T>`: an rxjs `Observable` with a synchronous `value` that replays its
  current value to a new subscriber and emits only on change. A set made while a state's subscribers are
  being called is queued and delivered after, in order. Sources never error: errors are values.
- **Pointer- and frame-rate values** are `Live<T>`: one record mutated in place by its manager, with a
  `version` that bumps on change. Read `current` when `version` emits; never keep it.
- **Events** (`SeqEvent<T>`) are happenings, not state: nothing replays; subscribers are called in
  priority order. Pointer-rate events reuse one record; edge-rate ones (`activate`, `patch`, `issue`,
  `log`) emit a new one.
- A host's subscriber that throws: on a state, reported as rxjs reports it; on an event, logged, and the
  others still run. The context guards its own fan-out (`seqv.state.error`; input handlers,
  `seqv.behavior.error`), so the viewer's own code never crashes a Node host.
- **Origins** travel with changes (`pointer`, `keyboard`, `api`, `graph`, `snapshot`, `behavior`,
  `host:<tag>`), so a host skips its own echoes. A size change (a resize, a DPR, the data area's inset) is `graph`:
  what its clamp moves is no one's move, so a linked section is set again from the one that led last.
- Managers and hosts name a section by its ref; hot records carry its index in the current view.

The context's states, by manager: `document`, `graph`, `views`, `snapshots`, `viewports`, `layout`,
`pointer`, `hover`, `highlight`, `selection`, `focus`, `systems`, `results`, `labels`, `ui`, `tasks`,
`streams`, `issues`, `history`, `behaviors`, `registry`, `renderProps`, `skin`.

## Events

`ctx.events`: `activate` (a track header activated), `menu` (a track's menu asked for from the keyboard),
`patch` (every proposal), `issue` (every new issue), `log`. The DOM element adds its own
(`seqspec:interaction#the-dom-element`).

## The document

The shown document is layered, each layer applied by patch so untouched subtrees keep their identity:
the **base**, then the current **snapshot**, then **local** patches, then a **preview**.

- `setDocument(next)` replaces the base by rebasing: the patch from the old base to `next` is applied, so
  a document that changed one node reruns one node. Snapshots, metadata and revision come from `next`.
- `applyPatch(patch, { layer })` applies a patch on a layer, all or nothing; a preview that does not
  validate is not drawn.
- **Patches out.** A UI action makes a patch and proposes it; the **document port** routes it: `commit`
  (the default port applies it to the base; a host's forwards it and later passes the new document back),
  `local` (kept in this viewer only; the owner's document is untouched), or `refuse` (`seqv.port.refused`).
- The default port routes a patch `local` **under a snapshot**, tagged with that snapshot (it applies
  only while that snapshot is shown), and when it is **masked**: it targets a node only a snapshot or a
  local patch adds, or a field a local patch sets. A base or local change drops a local patch made over
  the snapshot shown that no longer applies (`seqv.patch.dropped`); one made over another snapshot is
  kept and reported (`seqv.patch.skipped`). A snapshot change never drops one.
- A committed patch carries the base's revision as its `base`; a stale one is rejected.

## Data mode

A viewer is usable without a document: `setData` takes a `SeqSet` or the model's build inputs, and views
written with the builder against the set's keys, or a default view (a ruler, the letters, a track per
annotation by its geometry). Refs are the set's keys. Validation runs against the incoming data. Everything
downstream (patches to views, snapshots, managers) is as with a document.

## Updates and status

A document becomes a graph of nodes keyed by ref, whose edges are exactly the format's references
(`seqspec:format#what-a-change-affects`): a node's refs as core reads them (`seqspec:format#nodes-and-refs`), a
query's names, a theme's `column` and `profile` and a table's targets among them. An update reruns only the nodes whose
own content, or an input's result, or what they read in the set, changed; results that come out equal stop there. The
set it assembles is the one `compile` gives the same document, and so are its issues, at the document's paths. A view
switch resolves only that view; switching back is free. Registrations and the skin are inputs too: a light/dark switch
keeps every prepared state.

- **Status.** A node that fails is `error` with its issues; its dependents are `unavailable` with a cause;
  an unavailable layer is omitted (`svs.render.omitted`) and its section still draws. A missing provider,
  format or service is not a failure: the node is unavailable with `svs.provider.unavailable`. Nothing
  cascades.
- A node keeps at most `Issues.PerNode` issues, the most severe first; the rest are counted in one
  `seqv.issues.truncated`, by severity. The issues from no node (registration, behaviors, the engine) are kept up to a
  bound the same way: past it the least severe goes first, so a flood of warnings never hides an error.
- One update runs at a time; inputs arriving meanwhile are folded (the newest document wins, patches
  compose). Async nodes are aborted when superseded. An update always settles, and `graph.assembled`
  emits once, from which the managers remap.
- Held loci (selection, highlights, an external hover) follow their set: within a lineage they are
  rebound, across lineages remapped by key (`seqspec:model#identity-and-lineage`).

## Snapshots and playback

`snapshots.apply(ref)` sets the snapshot layer; `next`, `prev`, `play`, `stop`. Entering one follows the
format's rule (`seqspec:format#snapshots`): its focus changes are animated over its `transition_ms` by the
`focusViewport` behavior. `play` stays on each snapshot `duration_ms` after its transition and update end.
An invalid snapshot is listed and skipped. On load and on entering a snapshot, the next and previous
snapshots' focus regions are requested as prefetch hints.

## History

The context records `{patch, inverse}` for every patch it commits itself (with a port whose `history` is
true) and for local patches. `undo` applies the inverse; a new patch clears redo; an inverse that no longer
applies is refused (`seqv.history.refused`). When a host echoes the context's own commits back by
`setDocument`, history is kept; any other document clears it. A host that owns undo (Stories) sets
`history: false`. Interaction state (viewports, selection, expansion, params) is never in history.

## Commands

A **command** is the core's only request mechanism from outside: an id (`seqv:focus`, an extension's
`<ns>:<name>`), a params table, and its subscribers, as Mol\*'s `PluginCommand` with stable ids.

- `ctx.commands.dispatch(cmd, params)` from code; `ctx.commands.run(invocation, origin)` for plain data
  (`seqspec:format#descriptions-and-invocations`): the id looked up in the registry, params checked by its
  table. Unknown id or invalid params resolve `false` and run nothing.
- A command with no subscriber does nothing. The built-ins are handled by the `commands` behavior; a
  host's own subscriber runs in addition, after it; to replace what a built-in does, leave the behavior
  out.
- Subscribers run in order and are awaited; the first rejection rejects the dispatch
  (`seqv.command.error`).
- Targets are refs (with range args), loci (from code only) or `Coordinates` (runs in named systems). Refs of
  several kinds act on each kind (selected, lit, framed together), never on the first alone. Range args
  outside a ref's domain (index args, or labels a numbering resolves; a complex's: outside every member) are refused
  as a document's focus is, never clipped: nothing is done, `seqv.target.range` names the domain, `run` resolves
  false and `dispatch` rejects.
- Commands are edge-rate: a dispatch allocates its promise. Pointer-rate paths call managers directly.

## Behaviors

**Behaviors** are policy over the context's mechanism: what a hover, click, drag, wheel or key does, how
linked sections follow, what a command does. Each is a provider with a name, a category, a priority and
params (plain data), added and removed at runtime (`ctx.behaviors.add`, `get(name).setParams`).

- Removing a behavior releases everything it set (its channels, its hover, a pending drag, cursors), so
  nothing stays lit or stuck. It removes the policy only: clicks still emit `click`.
- A preset may turn a behavior off per section, or override its params (`forSection`).
- A behavior that throws while being created, registered or updated is disposed and reported
  (`seqv.behavior.error`); the viewer goes on without it.
- Bindings match modifiers exactly, so Ctrl/⌘ with `+` stays the browser's zoom.
- The interaction rules are behaviors' params, in the viewer spec (`SeqViewerSpec`, as Mol\*'s PluginSpec): a host
  composes, overrides, replaces or removes them there, or at runtime (`setParams`), and a preset per section. Input
  delivery (the engine, the dispatcher) only delivers: a press followed anywhere, its point clamped, a release never
  seen synthesized; without `select` a click selects nothing. The selection's rules are `select`'s: what Shift does on
  positions with a click or a drag (`shiftPositions`: `extend`, `toggle`) and on items (`shiftItems`), Ctrl/⌘ (`toggle`), whether a click waits
  where a double click is claimed (`doubleClickWait`), what a double click is (`repeatClicks`), the drag regions, and
  the scroll past an edge (`autoscroll`, `autoscrollMin`, `autoscrollRate`, `autoscrollMax`). A behavior whose double
  click means something at a point claims the `doubleClick` query there (`doubleClickZoom`, `navigator`, `headers`),
  so `select` asks the behaviors, never a list.

## Actions

**Actions** are what the UI shows as menu items and switches and what a host maps into its own commands:
a name, a title, a `when` and a `run`, mostly proposing a patch (`track.hide`, `section.numbering`,
`bookmark.add`) or changing UI state (`view.switch`, `group.toggle`). `ctx.actions.run`, or the
`seqv:runAction` command. The origin reaches the patch.

## Config

Typed items with defaults (`SeqViewerConfig.Viewport.SettleMs`), read by their managers at use:
`ctx.config.get(item)`, `set(item, value)`, or `spec.config`. A UI reads its items when it mounts; the one write
is the first surface's: `Viewport.ReducedMotion`, when unset, follows the page's reduced-motion setting.

## Registration

`ctx.register(entry)` adds registry entries: pure data, grouped by kind, as Mol\* v6's.

- One call registers every entry's format half (`svs`) first, then the viewer's kinds in a fixed order
  (node kinds, renderers, columns, themes, navigators, presets, labels, commands, actions), entries in
  order within each.
- Conflicts are checked before any change (a taken name, a renderer whose layer type is registered
  nowhere, a theme serving a kind its format entry does not); on a conflict nothing is registered.
- It returns an idempotent undo. There is no replace: a host that replaces a built-in registers a copy of
  `DefaultRegistry` without it.
- A registration is a graph input: only nodes whose entry changed rerun. A lookup that misses returns
  nothing; its user takes the declared fallback (an unknown preset, theme or renderer geometry is
  reported, an unknown data-column type is drawn as `text`).

## Extension points

What the core promises each kind of entry, and what the entry must do. The built-ins of each kind are in
the reference.

- **Node kinds** wrap the format package's stage functions, one to one. A host may replace one (a caching
  loader); it never adds a document kind.
- **Renderers** draw one layer type each, as a Mol\* representation draws one geometry; every built-in is
  one, with no private access. A renderer prepares headlessly and purely once per input change, hits
  headlessly, and draws through the renderer API (`/draw`), which is DOM-free, so the same entries serve
  a canvas, SVG in Node and tests. It declares the geometries it draws (`linear`, `wrapped`, `circular`);
  without `circular`, a section on rings draws its extents as bars and reports it. It must honor the kernel type's
  semantics (placement, split rows, channels, fallbacks). It declares `content: false` when it never
  reads residues, so a residue edit keeps its prepared state. It declares `zooms` when its free height follows the
  row zoom (its height then takes the row scale), so the zoom moves its track with the rows without a new layout.
  Its draws allocate nothing per frame. It may give the band a data hover lights on its track (`rowHover`: its rows,
  or none), else the track's box. It declares `stable` when its draw of a strip of the view (a pan's uncovered part)
  equals its draw of the whole view there, so a pan copies the rest: a decision over the view (items or density)
  reads the whole view's range, not the strip's. Its other optional parts: `tiny`, its draw below its `minHeight`
  (else bars); `miniature: false`, left out of a navigator's miniature; `overview`, its rows in the vertical
  overview; `ringLabels`, its items' labels outside the rings; `cache`, a store kept per prepared state and skin
  (made and disposed by the engine); `bands`, a multi-band layer's band labels; `header` and `gutter`, its rows'
  header and gutter (a rows layer's labels and tree, its line numbers).
  Renderers run the conformance kit (`@molstar/seq-viewer/testing`) as the built-ins do, its edge geometries
  (`EdgeKinds`: reversed, a cyclic origin, stacked, filtered, copies) included: drawn and hit back in each.
- **Data columns** draw a `rows` annotation's field per alignment row; an unknown column type is `text`.
- **Theme providers** are the drawing half of a format theme entry of the same name; one provider serves
  every channel kind its entry declares. A renderer fills a reused theme location (layer, block, unit,
  index, item, letter, value, state); a theme reads it and returns a packed color, a size or a glyph.
- **Navigators** draw the whole axis and the view's window; their gestures are the `navigator` behavior's. A
  navigator whose view can neither zoom nor pan is drawn idle and its gestures do nothing (interaction.md, Viewports).
  A navigator's parts (`NavigatorParts`) are an axis (the window over a strip, or a ruler's rail) and the tracks in
  miniature (each layer drawn small over the whole axis, `seqv:map`); the built-ins are in the reference, and a host
  registers its own, such as `ideogramNavigator(bands)`, a chromosome's cytobands under the window. A section with rows
  gets a vertical overview (the rows in miniature down a column, its window the rows in view) or a scrollbar, by its
  render props (`verticalNav`).
- **Presets** are a section's geometry, render props, navigator and per-behavior overrides, picked by the
  UI's runtime style, else the section's `style` (through `Presets.Alias`; an unmapped key is a preset
  name), else `Presets.Default`; an unknown one is the default (`svs.render.preset`).
- **Label providers** write structured lines for the hovered element into a reused sink, by priority.
- **Commands** and **actions** as above; **behaviors** are not registry entries but are added the same way
  by a spec or at runtime.

An extension ships one registry entry (format and viewer halves) and, when it has policy, behaviors.

## Labels

`ctx.labels` computes the hovered element's label by calling every applicable label provider with the
hover and a reused sink, grouping equal lines (a merged copy group's units become `× n`). One label is
built per hovered element; when to compute and hide it is the `hoverLabel` behavior's.

- **The codon.** On a block a streamed layer reaches through a path whose map takes each of the block's elements
  to three positions of its system (a codon step from the nucleotides, `seqspec:model#codon-maps`), the hovered
  element's label gives those positions' letters in reading order (`codon`; an unknown one `·`, a partial codon's
  known ones), split by `|` where they are not adjacent (an intron; on a cyclic system the origin is no gap), and
  complemented by the `complement` letter map (`complement_rna` for RNA) on the codon step's minus strand; then the
  positions in the system's numbering, the unit's own first, in reading order (high to low on the minus strand) with
  the strand (`codon at`: `chr2 1001–1003 (+)`). A lazy molecule's letters come from its stream, a view's letters
  layer of it first: those not held are asked after the view's regions, each label's replacing the last's, and the
  line says `loading…` until they arrive; `failed` (Retry asks again) or `not served` when they cannot.

## Host coordinate systems

The viewer talks to its host in coordinates, and a host may own coordinate systems (Mol\*'s chains).

- **Own keys.** The shown set's systems are `seqstar:<ref>` (a unit, a molecule, a complex, an alignment's
  columns by the alignment's ref, a row `<alignment>/<k>`, `k` in digits); a descriptor lists a numbering's name
  once.
- **Events describe their systems.** What the viewer reports outward gives extents as `Coordinates` with a
  plain-data descriptor per system (key, kind, domain, refs, a unit's `mvs` fields, numberings, a complex's
  members), so a host needs no Seq\* types.
- **Requests in any registered system.** `hover.setLoci`, `highlight.set`, the selection verbs and
  `focus.fromLoci` take `Coordinates`. A host registers an external system with its maps
  (`ctx.systems.register`, a `Translator`: `seqspec:libraries#translate`); it joins the view only through
  them. A structure registers as one complex system whose members are its chains'.
- Coordinates naming a key not registered survive as keys in a hover or a highlight, light nothing, and light again
  when it returns. The selection holds loci: its verbs convert Coordinates at once through the systems registered
  then, so a key not registered selects nothing (`set` clears the selection; register the system first).

## The render backend

The composer and every renderer draw through a `RenderBackend`, an interface at the level of what they
mean: batches of rectangles, sectors, paths, glyphs and text, images, passes, clips and copies. Backends:
Canvas2D (live), SVG (`/svg`, platform-free: the same code exports from a page and renders figures in
Node), and recording (tests).

- Batches are struct-of-arrays and grow-only; colors are packed integers or skin tokens, alpha its own
  byte. A frame allocates nothing.
- **Surfaces** (`mountSection`, one per section) size their canvases inside a frame only, so the page
  never paints a blank canvas, and follow the host's size and device pixel ratio.
- **Per-position text** (a block's letters, another polymer's residues, three-letter codes, atom symbols) is set per
  row: one font size and one form for the whole row, from its cells' width and never larger than the row allows nor
  than the cap (the render prop `maxLetterPx`; by default the header labels' size, `labelFontPx`), so zoomed in past
  it the text keeps its size; where the row's widest text does not fit, none is drawn and the cells stay. A row never
  shows some columns' text and not others'. Text never exceeds the part of its residue it is centered on: a residue's
  three-letter code needs a few px of margin inside the narrowest such part in view (a codon an intron splits has a
  smaller one), else the row shows one letter. The text layout's letters have a size of their own (`wrappedFontPx`),
  the cap's by default.
  A row that draws no cell colors (text only) marks where text does not fit which positions hold a residue: a dot in
  each cell while a cell holds one (about 40% of its width, at most 9 px), then a thin bar over each run of residues,
  so a deletion stays a break; both in the ink at half strength. A row with cell colors shows its cells, and marks nothing.
- **A boundary over letters.** Where a `ranges` layer draws over its own track's letters at their size (a design's
  row with its insertions), a boundary is a small triangle pinned to the row's top edge, pointing down at the edge
  between the two cells, and a hairline down that edge, in the item's color and clear of the letters' ink at every
  zoom; with no room above the text for a triangle 3 px tall, the hairline alone (broken around the ink where the
  letters leave it no room beside them). Where no text is drawn, the boundary is its glyph on its lane. At a block's
  first or last boundary a boundary's mark sits inside the data area, never cut. It hits within a few px of its edge
  over the whole row; beside it, the residues.
- `renderSvg` draws a whole view or section headless; an export never commits a viewport. A document's streamed
  regions are awaited at the size the export draws (each section its content's height), so a tall streamed alignment
  exports whole; `timeout` bounds the wait (a static source's load too), drawing what stands and saying so
  (`unsettled`, with the stream requests in flight; also when heights still move after a few rounds); `issues` gives
  the document's issues as it was drawn.
- **The `seqv-render` CLI** renders a document (JSON or MessagePack) as `renderSvg` does: `seqv-render <document>
  [--view <ref>] [--sections a,b] [--width 1200] [--scale 1] [--overlay] [--nav] [--background] [--title <text>]
  [--gap 8] [--root <dir>] [--timeout 120000] [--png] [-o <file>]`, SVG to stdout or `-o`, PNG with `--png` (or `-o`
  ending in `.png`, through a local Chrome). Sources are files under `--root` (default the document's folder), never
  outside it (by path, links resolved) and never the network; it streams nothing, so a streamed source is drawn
  without its data (`svs.provider.unavailable`). Exit codes: 0 written; 1 nothing written (the document does not load,
  the render or the PNG failed); 2 usage or an unreadable file; 3 written short of data (a source that could not be
  read, a warning or an error of the document's, or sources still loading after `--timeout` ms), each said on stderr.
- Rings (`circular`) are the same frame unrolled and mapped onto the circle, so every renderer, mark and
  hit runs as on a line; a track's layers overlay on its ring as they do linearly. Items' labels are set
  outside the rings, upright, each with a leader, and hit as their items; the ring's window
  (`seqspec:interaction#linked-sections`) is drawn on the overlay.

## Text measuring

Labels, ticks, data column text and wrapped lines need widths, and the layout is headless. Text is
**prepared** once per font into a bounded cache (`Text.CacheEntries`); layout is arithmetic over the
prepared arrays, with no measuring and no allocation (pretext's approach). `ctx.setTextMeasurer`
replaces the port; the first mounted surface sets the canvas measurer unless the host opts out.

## Render props

The layout's and engine's sizes and switches (`ctx.renderProps`, `spec.render`, `setRenderProps`), the
analog of Mol\*'s canvas3d props: header width, row height, gaps, fonts, ring geometry, navigation, LOD
thresholds. A section's props are its preset's over the context's, its inline style's keys over both, the
host's `Render.Clamp` last. All are JSON base types. A dot plot's height (its `aspect_ratio`,
`seqspec:format#layer-types`) combines three rules in order: the data area's width ÷ the ratio, capped by the box's
rows' room (its visible height under what stays pinned, the navigator and a pinned ruler) less a margin of a tenth of
it, at least two rows, so a plot that fits shows whole at more than one scroll position (720 px in a section sized to
its content), and raised to 120 px; then at least the width ÷ (ratio × `dotPlotMaxStretch`, 2.5 by default), which
wins over the cap. It follows that target, so a resize lands exactly, except a bounce: a shrink of less than 32 px back
to the target it just grew from (within the cap) keeps the taller height while the target stays there, so a section
sized to its content whose plot's shrink hides the page's scrollbar settles after one round instead of flickering.
Rows have no borders by default (`rowBorders`): a document
turns them on per section (`style.row_borders`), a host in its render props; an SVG export draws as the screen
does, and `renderSvg(doc, { render: { rowBorders: true } })` turns them on for a figure. An expanded split
track's rows have quiet guides regardless: a hairline in the skin's grid tone at the bottom edge of each row but the
last, under the data, drawn only where a row is at least 4 device px tall and shifting nothing (a row border, where
there is one, takes its place). Wrapped text sits `wrappedInset` px (8, the skin's
`--seq-canvas-edge` default) inside the canvas on every side, its line numbers' gutter within; hits, the selection,
the drag's nearest positions and its scroll all read that geometry. The layout reads no CSS (an SVG export draws it
the same): a host that changes `--seq-canvas-edge` sets `render.wrappedInset` to match. Its spacing is render props
too: `wrappedGroup` positions a group, `wrappedGroupGap` px between groups, `wrappedGap` between lines,
`wrappedTrackGap` between a line's letters and its first annotation lane, `wrappedLaneHeight` a lane under the text.
A line holds whole groups, and a short last group goes on a line of its own where a full group would not fit. The
`wrapped` preset spaces the text to read; `compact` is the dense one (half the group gap, lines and lanes touching,
lower lanes), still in groups of 10.

## Performance

No allocation per frame, per pointer event or per lookup on these paths: pointer moves and hit tests,
state and event emission, frames, hover across sections (closure and projection on cached routes), labels,
marks, drag-select, theme lookups, text layout, and SeqQL runs per position or run. Each such path has a
heap guard with a positive control (a guard that is shown to fail on a planted allocation). Commands, the
graph and edge-rate events may allocate.

A draw's edge-rate work (an index, a pyramid, a tick list) is bounded by what its block shows and its data holds,
never by the length of the block's coordinate system: a protein's letters on its chromosome index the protein's
codons, not the chromosome.
