# Seq* viewer reference

Generated from the viewer's tables; do not edit.

## Commands

Run from code (`dispatch`) or as data (`run`, an invocation `{command, params}`).

### `seqv:highlight`

Lights a target on one highlight channel, projected onto every section that shows it; null clears the channel.

| Param | Type | Default | |
|---|---|---|---|
| `channel` | `str` | required | the owner's channel; others are left as they are |
| `target` | `{refs: list<ref<unit \| sequence \| ligand \| complex \| alignment \| row \| annotation \| track \| layer \| block>>, args?: {start?: int \| null, stop?: int \| null, beg?: str \| null, end?: str \| null, numbering?: str \| null, items?: list<str> \| null} \| null} \| {coordinates: {systems: list<str>, runs: list<int>}} \| null` | required | null clears the channel |
| `style` | `"highlight" \| "accent" \| "band"` | `"highlight"` |  |
| `closure` | `"element" \| "containing"` | `"element"` | closed before it lights or acts |

### `seqv:select`

Sets, adds to, removes from or toggles the selection; null clears it.

| Param | Type | Default | |
|---|---|---|---|
| `target` | `{refs: list<ref<unit \| sequence \| ligand \| complex \| alignment \| row \| annotation \| track \| layer \| block>>, args?: {start?: int \| null, stop?: int \| null, beg?: str \| null, end?: str \| null, numbering?: str \| null, items?: list<str> \| null} \| null} \| {coordinates: {systems: list<str>, runs: list<int>}} \| null` | required | null clears the selection |
| `mode` | `"set" \| "add" \| "remove" \| "toggle"` | `"set"` |  |
| `closure` | `"element" \| "containing"` | `"element"` | closed before it lights or acts |
| `origin` | `str` | `"api"` | the caller's: `run` sets it |

### `seqv:focus`

Brings a target into view in one section or in every section that shows it.

| Param | Type | Default | |
|---|---|---|---|
| `target` | `{refs: list<ref<unit \| sequence \| ligand \| complex \| alignment \| row \| annotation \| track \| layer \| block>>, args?: {start?: int \| null, stop?: int \| null, beg?: str \| null, end?: str \| null, numbering?: str \| null, items?: list<str> \| null} \| null} \| {coordinates: {systems: list<str>, runs: list<int>}}` | required |  |
| `section` | `str \| null` | `null` | null: every section that shows it |
| `animate` | `bool` | `true` | false: at once |
| `origin` | `str` | `"api"` | the caller's: `run` sets it |

### `seqv:setSnapshot`

Shows a snapshot of the document, the next or the previous one, or the base.

| Param | Type | Default | |
|---|---|---|---|
| `ref` | `str \| null` | required | a snapshot ref, "next" or "prev"; null: the base |
| `animate` | `bool` | `true` | false: at once |

### `seqv:setView`

Switches to another view of the document; nothing is reloaded.

| Param | Type | Default | |
|---|---|---|---|
| `ref` | `str` | required | a view ref |

### `seqv:colorScheme`

Colors a layer by a registered color scheme, or back to its own color; refused (nothing done) for a scheme that does not fit.

| Param | Type | Default | |
|---|---|---|---|
| `layer` | `ref<layer>` | required | the layer whose color it sets |
| `scheme` | `str \| null` | required | a registered color scheme that fits the layer; null: its own color |
| `write` | `bool` | `false` | true: into the document (a patch on the layer's color); false: UI state |
| `origin` | `str` | `"api"` | the caller's: `run` sets it |

### `seqv:colorRange`

Sets a layer's color scale over [lo, hi], wider or narrower than its data; null: its own range.

| Param | Type | Default | |
|---|---|---|---|
| `layer` | `ref<layer>` | required | the layer whose color range it sets |
| `range` | `json` | required | [lo, hi] with lo < hi; null: the layer's own range |
| `origin` | `str` | `"api"` | the caller's: `run` sets it |

### `seqv:legend`

Opens, closes or toggles the legend; any other value is refused.

| Param | Type | Default | |
|---|---|---|---|
| `open` | `json` | required | true, false or "toggle" |
| `origin` | `str` | `"api"` | the caller's: `run` sets it |

### `seqv:retry`

Requests a track's failed stream regions again now, or every track's.

| Param | Type | Default | |
|---|---|---|---|
| `track` | `ref<track> \| null` | `null` | null: every track with a failed region |
| `origin` | `str` | `"api"` | the caller's: `run` sets it |

### `seqv:runAction`

Runs a registered action on a target, as a menu or a hover card would.

| Param | Type | Default | |
|---|---|---|---|
| `name` | `str` | required | a registered action |
| `target` | `json` | `null` |  |
| `origin` | `str` | `"api"` | the caller's: `run` sets it |

## Behaviors

In `DefaultSeqViewerSpec` order; params are defaults, set per behavior or per preset.

| Name | Category | Priority | What it does | Params |
|---|---|---|---|---|
| `scrollbar` | interaction | 100 | Dragging the thumb scrolls the rows; a press beside it pages. | `{}` |
| `navigator` | interaction | 100 | A press jumps; dragging the window pans; dragging an edge zooms; a double click toggles the full view and a zoom there. | `{"doubleClick":"toggle","zoomPositions":60,"zoomCell":0}` |
| `ringWindow` | interaction | 100 | On rings, the window of the linear section linked by `follow`: dragging it pans that section, dragging a grip zooms it with the other end held (a plain press of the primary button; a click still selects). | `{}` |
| `verticalOverview` | interaction | 100 | The window on both axes: top and bottom handles zoom rows, left and right columns; a press jumps. | `{}` |
| `headers` | interaction | 80 | A group's header, a split track's chevron or a clade toggles on a single click; activating a track header emits `events.activate`. | `{"activate":"doubleClick"}` |
| `select` | interaction | 50 | A click selects (again: clears a whole selection), held for the double-click time where a double click means something else (a behavior claims `doubleClick`: a zoom), which then never selects; Shift (a click or a drag) extends a range of positions or columns, and adds or removes an item or a row; Ctrl/⌘ adds or removes; a drag it claims selects a range, its end following the point (clamped to the surface) and the view scrolling past an edge. In text, a drag follows the pointer line by line. | `{"click":[{"button":"primary"}],"extend":[{"button":"primary","shift":true}],"toggle":[{"button":"primary","ctrl":true},{"button":"primary","meta":true}],"shiftPositions":"extend","shiftItems":"toggle","doubleClickWait":true,"repeatClicks":"ignore","clearOnEmpty":true,"toggleWhole":true,"closure":"containing","granularity":"position","dragRegions":["column"],"dragOn":[{"button":"primary","alt":true}],"dragColumns":[{"button":"primary","shift":true}],"autoscroll":true,"autoscrollMin":120,"autoscrollRate":15,"autoscrollMax":1440}` |
| `dotPlotSelect` | interaction | 60 | On a dot plot a drag selects a rectangle: the positions it spans across and down, as two ranges (Shift adds them). | `{"dragOn":[{"button":"primary"}],"add":[{"button":"primary","shift":true}],"closure":"containing"}` |
| `dragPan` | interaction | 10 | Dragging pans both axes, each with a lock; with `dominantAxis` held (none by default), the dominant axis only. | `{"dragRegions":["data","header","gutter"],"dragOn":[],"x":true,"y":true,"dominantAxis":null}` |
| `wheel` | interaction | 10 | Zoom at the cursor with Ctrl/⌘ and pinch; a horizontal wheel pans, and so does a vertical one with Shift; a vertical wheel scrolls a box section's rows. | `{"zoom":[{"ctrl":true},{"meta":true}],"zoomSpeed":1,"panSpeed":1}` |
| `doubleClickZoom` | interaction | 40 | A double click on an element shows its whole extent, then the focus window there, then the full view; on the data elsewhere, the focus window there, or zoomed in to it or closer, the full view. | `{"doubleClick":"toggle","zoomPositions":60,"zoomCell":0}` |
| `rowZoom` | interaction | 20 | Zooms rows vertically, every row at once. | `{"zoom":[{"alt":true}],"speed":1}` |
| `keyboard` | interaction | 0 | A cursor over positions, rows and headers (the hover, owner `keyboard`): arrows, Home and End; +/- zoom, Alt+=/Alt+- row zoom; Enter and Space select or activate; Shift+F10 or ContextMenu asks for the track's menu (`events.menu`); Escape clears the selection and puts the cursor away; a press puts it away too. Tab leaves the surface. | `{"bindings":{"left":[{"key":"ArrowLeft"}],"right":[{"key":"ArrowRight"}],"up":[{"key":"ArrowUp"}],"down":[{"key":"ArrowDown"}],"home":[{"key":"Home"}],"end":[{"key":"End"}],"zoomIn":[{"key":"+"},{"key":"="}],"zoomOut":[{"key":"-"},{"key":"_"}],"rowZoomIn":[{"code":"Equal","alt":true}],"rowZoomOut":[{"code":"Minus","alt":true}],"select":[{"key":"Enter"}],"extend":[{"key":"Enter","shift":true}],"toggle":[{"key":" "}],"menu":[{"key":"F10","shift":true},{"key":"ContextMenu"}],"clear":[{"key":"Escape"}]},"closure":"containing"}` |
| `hover` | interaction | 0 | The element under the pointer is the hover; cleared on leave and during a drag another behavior owns. | `{"regions":["data","column","header","gutter","dataColumn","dataColumnHeader","navigator","vnav"],"whileDragging":"clear","closure":"containing"}` |
| `headerHover` | interaction | 0 | A header's hover lights what it labels: a feature track's shown items (a split row's, a group's tracks'), a continuous track's row band; the row band alone past maxItems; forward gives linked views the set. | `{"features":"items","continuous":"rowBand","types":{},"maxItems":5000,"forward":true}` |
| `hoverRow` | interaction | 0 | A hover on a track's data (linear) lights its row band: the track's box, a split row's, a rows track's row, or its renderer's own (a matrix's hovered band); none on the ruler, a dot plot and navigators. A drawing only: the hover, its loci and what linked views get stay the element's. | `{"band":"entry","data":0.25,"header":0.5}` |
| `hoverLabel` | interaction | 0 | Computes the hovered element's label (labels.current); none over chrome and while dragging. | `{"regions":["data","column","header","dataColumn"],"whileDragging":"hide"}` |
| `cursor` | interaction | 0 | The surface's cursor by region (a window's body, edges and outside on the navigator, overview and scrollbar), and for a whole drag the press's (`grabbing` where it was `grab`); `grab` on data only where dragPan is on; on a rows layer's header, `rowHeader` off its tree's nodes; on a dot plot's data, `dotPlot`. With the cursor, the window part lit (`ui.window`: a grip, the body), the press's for a whole window drag. | `{"byRegion":{"data":"grab","column":"crosshair","header":"pointer","navigator":"pointer","vnav":"pointer"},"rowHeader":"default","dotPlot":"crosshair","dragging":"grabbing"}` |
| `focusViewport` | view | 0 | Moves the viewports of the sections that show a changed focus, over its transition. Wrapped text does not zoom: with `reveal` it scrolls the lines of the focus into view, over the same transition (`revealPad` px clear of the edge it scrolls to), and not at all when they already are shown; a focus before the text has lines waits for them. | `{"animate":true,"pad":0.1,"reveal":true,"revealPad":12}` |
| `scrollGroups` | view | 0 | Copies the zoom and horizontal scroll among a scroll group's members; each keeps its vertical scroll. | `{}` |
| `follow` | view | 0 | Sections that follow each other (`follow`) set their windows through a coordinate map: the other centers the image of the leader's center element at the leader's cell times the route's scale; each keeps its vertical scroll. | `{}` |
| `snapshotControls` | view | 0 | A press or a key on a surface stops the player. | `{"stopOnInput":true}` |
| `commands` | commands | 0 | Handles the built-in commands (SeqCommands) with the managers. | `{}` |

Built in, not in `DefaultSeqViewerSpec` (`SeqViewerOptionalBehaviors`): a host adds them.

| Name | Category | Priority | What it does | Params |
|---|---|---|---|---|
| `contextMenu` | interaction | 85 | A press of `on` (the secondary button) on a track (its data, its header or its gutter: `regions`) asks for the track's menu (`events.menu`, origin `pointer`, with the pointer's place), as the keyboard's menu key does; where no track is under the pointer nothing is asked. It asks on the press. The menu is the host's; not in `DefaultSeqViewerSpec`. | `{"on":[{"button":"secondary"}],"regions":["data","header","gutter"]}` |

## Config

Items with defaults, read by their managers at use (`ctx.config.get`, `set`).

| Item | Key | Default | What it sets |
|---|---|---|---|
| `Viewport.SettleMs` | `viewport.settleMs` | `150` | ms after the last interaction before the viewport settles and commits its view |
| `Viewport.AnimationMs` | `viewport.animationMs` | `250` | an animated focus or zoom, ms |
| `Viewport.OnAxisChange` | `viewport.onAxisChange` | `"keepCenter"` | when a section's axis changes: keep the view's center, or fit the whole axis |
| `Viewport.ReducedMotion` | `viewport.reducedMotion` | `null` | true: no animation; null: the engine reads `prefers-reduced-motion` |
| `Input.DragSlopPx` | `input.dragSlopPx` | `3` | CSS px a press moves before it is a drag |
| `Input.DoubleClickMs` | `input.doubleClickMs` | `300` | the longest gap between the clicks of a double click, ms |
| `Tasks.BusyDelayMs` | `tasks.busyDelayMs` | `750` | ms a task runs before the UI shows it as busy |
| `Issues.PerNode` | `issues.perNode` | `100` | issues a graph node keeps, the most severe first; the rest are counted in one `seqv.issues.truncated` |
| `Stream.SampleBudget` | `stream.sampleBudget` | `null` | samples per view at most, where the source lists a level that coarse; null: one per device pixel of the data area |
| `Stream.MaxInFlight` | `stream.maxInFlight` | `6` | stream requests in flight at once |
| `Stream.CancelGraceMs` | `stream.cancelGraceMs` | `250` | ms a request no view needs any more stays in flight before it is aborted (sooner when its slot is wanted); 0: at once |
| `Stream.CacheBytes` | `stream.cacheBytes` | `134217728` | bytes of streamed regions held before eviction |
| `Stream.Margin` | `stream.margin` | `1` | view widths requested on each side of the view |
| `Stream.LettersPerPx` | `stream.lettersPerPx` | `16` | a letters stream is asked only where a device px holds at most this many positions |
| `Stream.RegionBytes` | `stream.regionBytes` | `1048576` | a data request's bound on what its provider reads (`max_bytes`); past it the provider answers with counts |
| `Stream.ShimmerFps` | `stream.shimmerFps` | `30` | frames a second of a pending region's shimmer; 0: a still band |
| `Stream.Prefetch` | `stream.prefetch` | `true` | request ahead: a view width past the margin on the side last panned toward, a zoom's next level, and the regions of the snapshots next to the one shown |
| `Text.CacheEntries` | `text.cacheEntries` | `20000` | measured strings the text measurer keeps |
| `Presets.Default` | `presets.default` | `{"blocks":"tracks","columns":"msa","text":"wrapped"}` | the preset by section kind (blocks, alignment columns, text layout) when neither the UI nor the section's `style` names one |
| `Presets.Alias` | `presets.alias` | `{}` | a section `style` key → a preset name |
| `Render.Clamp` | `render.clamp` | `null` | bounds on the numeric render props, over presets and styles; null: none |
| `History.Limit` | `history.limit` | `200` | undo steps kept |
| `UI.Chrome` | `ui.chrome` | `"clean"` | the default layout's composition: `clean`, the canvas (hover, selection, the tooltip and hover card, the navigator), the player, the description, search, busy and the states; `full` adds the tabs, section headers, track menus, legend, export and issues |
| `UI.Legend` | `ui.legend` | `"collapsed"` | the legend `collapsed` (its toggle only), `expanded` or `hidden`; a UI's own option, the person's toggle and `view.legend` go first |
| `UI.LegendItems` | `ui.legendItems` | `24` | items a legend group shows before its "+n more" |
| `UI.LegendGroups` | `ui.legendGroups` | `8` | groups the legend shows before its "+n more" |
| `UI.TrackActions` | `ui.trackActions` | `["rows","filter","focus","colors","hide","retry"]` | the track menu's items in order; `none`: no menu (`remove` is built in but not listed by default) |
| `UI.ColorSchemes` | `ui.colorSchemes` | `true` | the color scheme pickers (a legend group's title, a track's Colors) |
| `UI.ColorSchemeWrites` | `ui.colorSchemeWrites` | `false` | true: a picked scheme is written into the document; false: UI state only |
| `UI.ItemActions` | `ui.itemActions` | `[]` | the registered item actions a hover card offers, in order; empty: no card |
| `UI.HoverCardDelayMs` | `ui.hoverCardDelay` | `400` | ms the pointer rests on an item with actions before its card opens |
| `UI.OwnRowHover` | `ui.ownRowHover` | `"residue"` | a bar on a split row on a unit of its own: `item` hovers the whole item, `residue` the residue under the pointer |

## Built-in registry entries

What `DefaultRegistry` registers, by kind.

### Node kinds

| Name | Stage |
|---|---|
| `source` | load |
| `sequence` | load |
| `ligand` | load |
| `unit` | load |
| `alignment` | load |
| `mapping` | load |
| `annotation` | load |
| `complex` | complex |
| `computed` | compute |
| `#assemble.base` | assemble |
| `#assemble.mapped` | assemble |
| `#assemble` | assemble |
| `#assemble.stream` | assemble |
| `#entry` | load |
| `#set` | assemble |
| `section` | section |
| `layer` | layer |
| `#prepare` | prepare |
| `#theme` | theme |
| `#template` | layer |
| `#view` | view |

### Renderers

One per layer type; a geometry it lacks draws the layer's extents as bars and is reported.

| Layer type | Geometries |
|---|---|
| `letters` | linear, wrapped, circular |
| `ranges` | linear, wrapped, circular |
| `points` | linear, wrapped, circular |
| `emphasis` | linear, wrapped, circular |
| `variants` | linear, circular |
| `secondary` | linear, wrapped, circular |
| `values` | linear, circular |
| `color` | linear, wrapped, circular |
| `ruler` | linear, circular |
| `pairs` | linear, circular |
| `rows` | linear |
| `logo` | linear, circular |
| `matrix` | linear |
| `dot_plot` | linear |
| `seqv:map` | linear |

### Data columns

An unknown column type is drawn as `text`.

| Type | Width (CSS px) |
|---|---|
| `number` | 64 |
| `fraction` | 84 |
| `text` | 110 |
| `category` | 100 |

### Themes

| Name | Label | Category | Channels |
|---|---|---|---|
| `uniform` | uniform | kernel | color, glyph, size |
| `field_scale` | field_scale | kernel | size |
| `field_palette` | field_palette | kernel | color |
| `field_color` | field_color | kernel | color |
| `letter` | letter | kernel | color |
| `by_unit` | by_unit | kernel | color |
| `by_copy_group` | by_copy_group | kernel | color |
| `field_map` | field_map | kernel | glyph |
| `field_glyph` | field_glyph | kernel | glyph |
| `rules` | rules | kernel | color, glyph, size |
| `scheme` | scheme | kernel | color |

### Navigators

| Name |
|---|
| `rulerRail` |
| `band` |
| `strip` |
| `rail` |
| `none` |

### Presets

| Name | Geometry | Navigator | Behaviors changed |
|---|---|---|---|
| `tracks` | linear | `rulerRail` | none |
| `msa` | linear | `rulerRail` | none |
| `strip` | linear | `strip` | `select` |
| `wrapped` | wrapped | `none` | `select`, `dragPan`, `cursor`, `wheel` |
| `compact` | wrapped | `none` | `select`, `dragPan`, `cursor`, `wheel` |
| `circular` | circular | `none` | `wheel`, `dragPan` |

### Label providers

| Name | Priority |
|---|---|
| `position` | 100 |
| `copies` | 95 |
| `row` | 90 |
| `item` | 80 |
| `layer` | 75 |
| `value` | 70 |
| `pair` | 70 |
| `column` | 60 |
| `stream` | 50 |
| `codon` | 90 |

### Actions

| Name | Title | Key |
|---|---|---|
| `track.hide` | Hide track |  |
| `track.show` | Show track |  |
| `track.remove` | Remove track |  |
| `track.filterView` | Show only what this track covers |  |
| `section.clearFilter` | Show all positions |  |
| `section.numbering` | Numbering |  |
| `section.copies` | Copies |  |
| `section.ligandForm` | Ligand form |  |
| `section.rowSizing` | Row sizing |  |
| `section.style` | Style |  |
| `layer.colorScheme` | Color scheme |  |
| `layer.colorRange` | Color range |  |
| `view.legend` | Legend |  |
| `view.switch` | Switch view |  |
| `group.toggle` | Expand or collapse |  |
| `clade.toggle` | Expand or collapse clade |  |
| `selection.clear` | Clear selection |  |
| `viewport.fit` | Fit |  |
| `stream.retry` | Retry |  |
| `history.undo` | Undo |  |
| `history.redo` | Redo |  |
| `bookmark.add` | Bookmark |  |
| `bookmark.remove` | Remove bookmark |  |
| `bookmark.step` | Next bookmark |  |
| `bookmark.go` | Go to bookmark |  |
| `highlight.remove` | Remove highlight |  |
| `results.save` | Keep results |  |

## The viewer's layer types

### `seqv:map`

other tracks in miniature over the whole axis, the view as a window

| Param | Type | Default | |
|---|---|---|---|
| `tracks` | `list<ref<track>>` | `[]` | the tracks; empty: all others |

## The viewer's issue codes

| Code | Severity |
|---|---|
| `seqv.patch.dropped` | warning |
| `seqv.patch.skipped` | warning |
| `seqv.history.refused` | warning |
| `seqv.port.refused` | warning |
| `seqv.issues.truncated` | warning |
| `seqv.ui.slot` | warning |
| `seqv.render.labels` | info |
| `seqv.renderer.error` | error |
| `seqv.behavior.error` | error |
| `seqv.command.error` | error |
| `seqv.state.error` | error |
| `seqv.target.range` | error |
