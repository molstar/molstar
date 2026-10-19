# The sequence view and the 3D view

How the viewer keeps Seq\* and the Mol\* plugin in step, how Undo works across both, and the sequence band in the shell.
Seq\*'s host contract is `packages/sequences/_spec/` (`viewer.md`, `sources.md`, `libraries.md` § The Mol\* adapter).
The code is the truth for how.

## 1. Two engines, one viewer

- **Neither engine knows the other.** Mol\* (`PluginContext`) owns data, parses, structures, representations, camera,
  hover, selection and focus as loci; Seq\* (`SeqViewerContext`) owns the sequence views, their layout, viewports and
  its own hover, selection and focus. The viewer joins them three ways: data through Seq\*'s `Loaders` port (§2),
  interaction through the link (§4), content through sequence behaviors (§3). Actions on the 3D view go through the
  viewer's API (§5).
- **The viewer owns the sequence as its own feature** (`src/sequence/`), not an extension: no new package, no prism-ui
  abstraction; prism-ui stays free of Seq\*. Dependencies: `@molstar/seq-viewer`, `@molstar/seq-core`.
- **The Mol\* adapter is the viewer's** (`src/sequence/molstar/`, exported as `@molstar/prism-viewer/sequence/molstar`;
  its contract is Seq\*'s `libraries.md` § The Mol\* adapter): structures read headless as Seq\* complexes, the
  structure formats, services and providers (`MolstarSvsEntry`; `afdb_pae` among the formats), the link
  (`StructureLink`, registered on a Seq\* system registry: the Seq\* viewer's `ctx.systems`, or one with no viewer;
  `seqspec:libraries#the-system-registry`; and `linkViews`) and SeqQL over Mol\*. The Seq\* Storybook's Mol\* demo
  imports it through the export; its MVS composition and bare plugin are the demo's own.
- **The sequence model owns every Seq\* context** (`createSequenceViewers`, `sequence/model.ts`), keyed by id with
  `'default'` the only one: one sequence viewer, switched between Text and Tracks. It takes what it follows as plain
  arguments (the ink, the host structure, the actions, a scheduler): no DOM, no plugin, testable headless.
- **One composition**: `PrismSequenceSpec()` (`sequence/spec.ts`) is the one place Seq\* is composed and trimmed:
  Seq\*'s default spec (every layer type and behavior, so a loaded SeqViewSpec draws in full), the Mol\* structure as a
  source (`MolstarSvsEntry`), Prism's palette (`PrismRegistry`), `revealHighlight` on `links:*` and the `3d` hover, and
  the text version's render props.

## 2. Data: read once, through the plugin

- **A runtime chunk.** Seq\* and everything over it is `sequence/runtime.ts`, loaded (`spec.sequence.runtime`,
  viewer-spec.md §2.3) when a structure is first shown, or an operation or a saved state needs it. The core holds only
  the plugin side (`sources/sequence.ts`); a spec without a sequence bundles none of it.
- **One parse.** The document names its structure by host file `structure:<key>`; the loaders' `parsed` answers it with
  a `StructureSource` over the trajectory the 3D view drew, so `StructureLink` maps the very `Structure` the 3D view
  draws. One source per trajectory, with the CIF block it was parsed from (`cifFrameOf`: SIFTS numberings, a predicted
  model's origin, the scheme's auth labels). A new key rebuilds the view.
- **A frame keeps the document.** A trajectory step or the model-index animation makes a new `Structure` of the same
  atoms (unit ids, element indices and model hierarchies equal, same trajectory): it keeps the key, so no document is
  made and the band does not empty. Only `SequenceViewer.frame` changes; the link rebinds `StructureLink` to it (its
  systems, maps by unit id and chains stay) and starts its views again, as the 3D view's loci are of the frame drawn.
  Tracks read from the structure (B-factors, secondary structure) stay those of the frame the document was made for. A
  structure of other atoms (another assembly, a model with another hierarchy) is a new key. On 20 steps of a 10-model
  11.5k-atom ensemble: 40 structure changes and 168 system registrations → 0 (measured, headless).
- **Two caches, both read-through, neither fetching for the sequence:**
  - **Bytes, in core**: Mol\*'s asset manager reads URLs through a URL-keyed, ref-counted byte cache
    (`packages/core/src/util/assets.ts`): the plugin's downloads, MVS's and the viewer's loaders share one fetch.
  - **Parses, in the viewer**: `dataByUrl` (`sources/data-by-url.ts`, `SequenceSource.dataByUrl`) indexes the plugin's
    data cells by address (a Download's URL, a read file's name) with the trajectory parsed under each. A data-tree
    change only marks it stale; the next lookup rebuilds it (edge rate), so a large scene pays nothing while no document
    loads.
- **`load` reads what the plugin holds**: a URI (resolved against a loaded document's own address) or host file the
  index has is answered with its bytes; a `molstar:structure` source it parsed is shared through `parsed` (one
  `StructureSource` per trajectory), never parsed again; a loaded `.svsx`'s relative URIs come from the archive (§6).
  Anything else is refused, saying to load it in the viewer first (`svs.source.load`): the sequence fetches nothing.

## 3. Content: the document and sequence behaviors

- **What the sequence follows** (`createSequenceSource`, `sources/sequence.ts`): a structure `show` names (§6), else the
  one `sequence.follow` set, else the focused one, else the one whose selection grew last, else the first loaded. No
  picker in the band: switching is rare, and the follow rule and `sequence.follow` choose.
- **The built document** (`sequenceDocument`, `sequence/document.ts`): one view, one section (`main`) over a
  `molstar:structure` complex (ref `structure`); with no blocks the section takes every polymer unit in order, then the
  ligands the 3D view draws (`ligand_comps`). Tracks: ruler, letters, secondary structure
  (`molstar:secondary_structure`); the data node has ref `data`. One document serves Text and Tracks (§8).
- **Behaviors add to it by patch.** A `SequenceBehavior` (`sequence/behaviors.ts`) attaches with a context: `shown`
  (what the document shows, a new value per document: structure, complex, data node, section, the last track),
  `patch(shown, ops)` (Seq\*'s `applyPatch` on the base layer, origin `host:<id>`; dropped for a replaced document) and
  `offer` (its tracks' menus by track ref and color theme, and its legend blocks). A new document has none of a
  behavior's nodes: `shown` changes and it adds them again.
- **Refs start with the behavior's id**; ids are `<owner>.sequence.<name>`, unique (a taken id throws). The spec's own
  behaviors come from `spec.sequence.behaviors`; an extension brings its own (`ViewerExtension.sequence`, viewer-spec.md
  §3; `extensions/host.ts` checks the prefix), overridable by id as any part: model-archive's pLDDT track (README
  decision 27, `extensions/model-archive/src/prism-viewer.ts`).

## 4. Interaction: the link

- **One link per shown document** (`linkSequence`, `sequence/link.ts`): `StructureLink.attach` over the shown structure
  and complex, `linkViews` for hover and selection both ways and focus; made again for each document, subscribed while
  the band is rendered and dropped with it. The replay rule makes remounting safe (a subscribe sends the selection,
  never a clear; `seqspec:interaction#replays-and-carries`).
- **The 3D side is `view3d`** (`SequenceSource.view3d`): reads hover and selection from the plugin; writes highlight,
  selection and the camera's focus as operations with origin `link` (`highlight.set`, `select.set`, `camera.frame`),
  never on the plugin directly. View changes, not steps (api.md §1, `history: false`). One channel per kind, in order; a
  write superseded while it waited is skipped.
- **State lives in the engines**: selection and focus are Mol\*'s, mirrored into Seq\*; the sequence's viewport and
  style are Seq\*'s. No third copy.
- **A ligand** position lights and selects its whole residue in 3D (`StructureLink.toLoci`).
- **Emptying the selection from the sequence** (origin `pointer` or `keyboard`: a click, Shift+drag removing the last
  stretch, Escape) resets the camera once (`SequenceView3D.resetCamera`, the person's `camera.reset`); an emptying from
  3D, the selection bar, the API or a carry leaves it.
- **Marks** ([marks.md](marks.md)) show in the sequence as highlight channels `marks:<owner>:<id>` in their tone, and
  the view's hover as the `3d` hover.
- **Size gates** are `linkViews`' defaults (hover 2,000 items, selection 20,000; past them, whole chains).

## 5. From the sequence back to the viewer

- **Every action on the 3D view is an API call** (api.md §11), attributed and recorded, as the command and the tool that
  do the same. Selection needs none: the link mirrors it (§4).
- **A track's coloring** (its menu, its gutter icon) is `HostActions.colorStructure` → `representation.color` on the
  followed structure's polymer (`default` is `chain-id`); `HostActions.coloring()` reads the polymer's current theme
  from the plugin (one theme across its non-ball-and-stick representations).

## 6. Both documents: MolViewSpec and SeqViewSpec

- **An MVS needs no sequence part**: it builds the plugin state and the sequence follows it (§3).
- **The viewer loads an MVS or an SVS, not both at once.** A document holding both is Stories'.
- **`load.svs`** (the `loaders` feature, `features/loaders/operations.ts`): a URL (`.svsj`/`.svs`, JSON or MessagePack;
  `.svsx`), a file (also dropped, after the files dropped with it) or data; `svs-url` in the app's URL. A document
  without a `molstar:structure` complex is refused.
- **Its first `molstar:structure` complex is matched** to a loaded structure (`matchStructure`): the Download cell's
  address (resolved against the document's), the entry id (`service` `pdbe`/`afdb` and `id`; an AlphaFold model by
  accession) or the file name (`host_file`). None matches: it is loaded (`load.url`, `load.entry`, `load.data`) in the
  same step.
- **It replaces the followed structure's document** (`SequenceOverride`, `setOverride`): `svsDocument` makes a copy
  whose complex reads the host's parse, gives refs to the data node, sections and last track (so behaviors attach by
  `shown`) and puts the person's letter colors over its own; the person's style applies to every section; the band shows
  its first section with tracks. The source shows that structure (`SequenceSource.show`) **until the person focuses
  another, a structure is loaded, `sequence.follow` changes, or it goes**; then the built document is back.
  `sequence.get` names it.
- **An `.svsx`** is read with Seq\*'s `readSvsx`: a structure it holds opens as a file of that name; its relative URIs
  resolve inside it (`addArchive`). The files stay in memory for the page's life, out of the saved state.

## 7. Snapshots, history and undo

### 7.1 The sequence's state in snapshots

- **Saved as `prism-viewer.sequence`**: `{ version: 3, style, colors, rail, structure?, override?, edits? }`
  (`SavedSequence`, `sequence/runtime.ts`). `override` is the loaded SeqViewSpec (the document itself, its URL, its
  structure's ref): content, so Undo keeps it. Version 1 reads as 2 without `override`, version 2 as 3 without `edits`;
  another version or none resets to the defaults.
- **The edits overlay** (`SequenceEdits`): the person's changes over the document shown, by track ref
  (`{ tracks: { [ref]: { hidden } } }`), each a change from the document's own value (a loaded SeqViewSpec may hide a
  track itself), so showing a track the person hid drops its entry. Written into every document the model builds and
  into a behavior's inserted tracks (hidden from its first patch, no flash), and brought to the document shown by one
  base patch (`set` of `hidden`, origin `host:prism-viewer.sequence.edits`) when they change; a ref the document lacks
  is kept for the next that has it. So a hidden track stays hidden across another structure and back, new colors, a
  reload and Undo. Seq\* hides by its document (`track.hidden`): no renderer change.
- **A registry, not a field**: `PluginState.registerCustomState(key, { get, set })`
  (`packages/plugin/core/src/state.ts`). `getSnapshot` writes every key into `snapshot.custom`; `setSnapshot` calls each
  `set` after the data tree, with `undefined` when the snapshot has none; transition frames do not apply it. Snapshots,
  their playback and `.molx` / `.molj` sessions carry the sequence with no other hook. Before the runtime is in, the
  value is held as given.

### 7.2 Checkpoints

- **The history uses the plugin's state directly, never its snapshot manager** (which adds entries and clears on
  restore; README decision 24): `createCheckpointPort` (`sources/checkpoint.ts`) takes
  `getSnapshot({ data, behavior, componentManager, canvas3d })`, leaving out interactivity, selection, camera, canvas
  context and animation, minus the followed structure; restore applies it without the focus and keeps each trajectory's
  current frame.
- **Undo restores**: the data tree, behavior settings, canvas props, the sequence's style, colors, rail, override and
  edits (its hidden tracks), the hidden contents. **Never**: the camera, the selection, the selection mode and picking
  level, the focus, which structure the sequence follows, hover, a trajectory's frame, the Snapshots list, open panels.
- **Headless**: `camera: false` keeps the checkpoint from reading the canvas.
- **A lazy extension's load** stays across a restore to before it (its registrations stay; its cells go).

### 7.3 One history manager

`managers.history` (prism-ui, `ui/src/managers/history.ts`) holds every undoable step, the person's and the agent's.

- **A step** is one undoable run (`api.run`, api.md §1) or one group: atomic (`group`: a script's `viewer.run`, a
  command's several calls) or open (`beginGroup(origin)`: an agent's turn, keyed by `origin.turn`). One checkpoint is
  taken when it opens.
- **Not steps**: `history: false` operations (`select.*`, `focus.*`, `camera.*`, `highlight.*`, `view.*`,
  `sequence.follow`, `trajectory.frame`, playing animations, showing a snapshot: `snapshot.apply`, `.next`, `.previous`,
  reads); a result that changed nothing ("it is already so"); `undoable: false` (an upload: no checkpoint, the steps
  before it stay undoable).
- **Inside a step**: each top-level run or group gets a `HistoryContext`; work that calls back in with it (an
  operation's `ctx.api`, a group's API, `viewer.run`'s `api`) is part of that step, with no lock. Calls without it queue
  and are steps of their own: a click during an agent's snippet is never rolled back with it.
- **While a turn runs**, the person's own calls still run; the first closes the turn's entry (its later calls open a new
  one under the same turn). Undo turn undoes the turn's entries while they are the latest, else asks before discarding
  the later steps.
- **Undo while a run is in flight or queued is disabled**, never queued; the button and the command say why.
- **A failed or aborted run restores its checkpoint**: nothing half-built, and the result says so.
- **Redo** puts back what Undo took (the state each Undo left is kept); a new step or a change outside the API empties
  it.
- **Changes outside the API** (`watchOutsideChanges`) make every earlier checkpoint unsafe, and History says so. Not
  outside: runs, the history's own restores, cells behaviors own (the focus representation, volume streaming, custom
  properties), snapshots shown or playing and frames or animations viewed.
- **Showing a snapshot is viewing a story** (David, 2026-10-20): no step and no checkpoint, so stepping through an MVS
  story leaves the steps before undoable. A checkpoint notes the snapshot current with it, and Undo makes it current
  again with the scene it restores (when the list still has it; the list itself is not undone). Like any
  `history: false` run, a failed apply is not rolled back (`setSnapshot` checks the transformers before it changes
  anything).
- **One undo**: Mol\*'s state history is off (`PluginConfig.State.HistoryCapacity` 0, `sources/compose-plugin.ts`).
- Every run and step is also in the session record ([journal.md](journal.md)).

### 7.4 Large scenes

- **Each checkpoint is measured**: time, and the JSON size without data payloads; the time counts the size pass.
- **A step costs what it changed, not the scene** (`SerializedTrees`, `sources/checkpoint.ts`): each transform is
  serialized and measured once per transform object (a new version is a new object) and shared by every checkpoint that
  holds it, its state flags compared on each read (the state's transient tree changes them in place). At 1401 cells: the
  size pass 16 → 0.2 ms a step, 100 checkpoints 14.8 → 1.6 MB. Shared parts are read-only.
- **Over the hard limit, no checkpoint is even taken for 30 s** (or until a clear): each try is a full snapshot thrown
  away, and an agent's turn tried on every call.
- **Per-step checkpoints turn off** for the scene when one is over the size budget, or three in a row over the time
  budget (50 ms, 5 MB: a first guess); the status says once "This scene is large: undo is off. Settings turns it back
  on." (Settings › Undo each step).
- **A group's one checkpoint** (an agent's turn) is still taken under the hard limit (500 ms, 50 MB). When a turn cannot
  be undone, its first changing call asks the person, once for the turn (`agent/policies.ts`).
- **Clearing the scene** turns per-step undo back on (`session.clear`).
- **Undo goes back the latest 10 steps** (`options.ui.undoSteps`, `HistoryOptions.checkpoints`): only they keep a
  checkpoint; an older step stays in History (100 at most, `capacity`) as unsafe, "Too old to undo", its checkpoint let
  go. The record still lists every run (journal.md). Ten because a step's checkpoint holds its whole scene: memory is
  bounded by a few scenes, not a hundred (David, 2026-10-20).

## 8. The band and its menu

- **The band is a frameless top tab** (`id: 'sequence'`, `frame: 'none'`, `sequence/shell-ui.tsx`), registered once the
  first document is in and kept until the host shows no structure (another structure's document being made keeps it);
  the first one shows it. Overridable as `core.tab.sequence`.
- **It lies over the canvas' top**, in the canvas' color behind its lines, fading out below them (`BandFadePx`,
  `BandDissolvePx`; a gradient, no `backdrop-filter`); the canvas keeps its full height and ends with it
  (`--msp-canvas-end`). The lines take the pointer, the fade lets it through.
- **The control has no controls of its own** (`SequenceControl`): Seq\*'s `mountSection` for the shown section (`main`,
  or a loaded SeqViewSpec's), none of Seq\*'s DOM or React UI. Seq\* draws letters, ruler, hover, selection and focus
  and handles wheel, drag, click and keys.
- **The canvas stack's Sequence button** (shown while a structure is) has the menu (`sequenceMenu`): one choice of No
  sequence, Text, Tracks, then the letters' colors (one choice) and Legend…, then the Visible tracks submenu (below).
  Each item runs a command (`sequence.none`, `sequence.text`, `sequence.tracks`, `sequence.colors.<scheme>`,
  `sequence.legend`), so the command menu does the same, the Visible tracks rows the operations their commands run; the
  rail is commands only (`sequence.rail.<ruler|strip|band|none>`). Style, colors, rail and hidden tracks go through the
  API as the person: steps.
- **Tracks is the default style** (`spec.sequence.style`); the switch only sets the preset (`ctx.ui.setStyle`): the
  document and the link survive.
- **Text** is Seq\*'s `compact` with the composition's render props: whole groups per line, centered, letters at most
  the panel text size (`--msp-font-size-sm`, read once with the fonts), in the layer's colors, unobserved positions in
  the skin's `unobserved`; tracks show as lanes (secondary structure underlines, a `values` track a strip;
  `seqspec:viewer#the-text-view`). The edges fade where there is more to scroll and the scrollbar shows on hover, scroll
  or drag, both written when they change, never per step.
- **Tracks** is `prism_tracks` (`sequence/presets.ts`): no header column, the rail a host preset (`TracksRail`, the
  ruler's own by default), a left gutter of `TrackGutterPx` when the host offers colorings.
- **Height**: Tracks follow `seq.layout.contentHeight(section)`, Text is fixed (`TextHeight`); plus the air above and
  the fade. Written only when it changes, kept while a new document is made, restored when the band goes.
- **Track menus**: Seq\*'s `interaction.on('menu', …)` on a track of the section shown opens prism-ui's
  `openContextMenu` at the pointer and prevents the browser's menu: the track's actions (the model's `TrackMenus` and
  the behaviors' offers), then Hide track; elsewhere the browser's menu stays.
- **Hiding and showing tracks** (`viewer.tracks`: the section's tracks, hidden ones included, each named by its document
  label, the viewer's own name, a behavior's offered `labels` or its title): `sequence.track.hide` and
  `sequence.track.show` (a track by ref or name; show with none: every hidden one), steps through the API as the person,
  so Undo takes each back (§7; show with none is one edit, `setTracksHidden`); the commands `sequence.track.hide`,
  `.show` (a choice of the tracks) and `.show-all`. The Sequence menu's **Visible tracks** is a submenu (prism-ui's
  `CanvasMenuItem.items`): a checkbox per track (checked: shown), then Show all while one is hidden. `sequence.get`
  lists the tracks.
- **The hidden tracks' hint**: only while tracks are hidden, a faint "N hidden" at the foot of the track gutter (Tracks
  only); a click opens a menu at it: Show all tracks, then Show each hidden one (a context menu: the canvas stack's
  submenu belongs to its button). Nothing shows while every track does.
- **The track gutter** (`sequence/track-gutter.tsx`, Tracks only): a coloring icon per track row that has one (its
  menu's first action), rows read from Seq\*'s layout when it is rebuilt, the row under the pointer from
  `interaction.on('hover')`; the 3D view's current coloring stays lit in the accent; clicking it again does nothing.
- **The keyboard cursor is said** in a polite live region once the arrows rest.
- **Framing below the band**: `paddedFocus` (`sources/framing.ts`) frames a structure in the canvas below the band, as a
  snapshot provider for `requestCameraReset`.

## 9. The skin

- **The ink is prism-ui's** (`CanvasInk`, `canvasInk`, `ui/src/model/canvas-theme.ts`); the viewer holds it and hands it
  to the sequence model, which follows it while it exists (so the skin is right before the band mounts).
- **An explicit skin**: `skinFor(ink, fonts)` (`sequence/skin.ts`): surface from the ink, text and accent from Prism's
  palette for the ink's theme, selection the accent at 25%; Seq\*'s skin reader (`--seq-*`, `data-msp-theme`) is not
  used. Nothing reads the DOM for colors; the control reads `--msp-font-*` once into the skin (`setFonts`).

## 10. Colors, legend and the hover label

- **The model owns the letters' scheme** (`SequenceViewer.theme`, written into the document; Seq\* draws it). Default
  `prism:residues` (`sequence/palette.ts`: Zappo's classes in Prism's hues, light and dark), registered as a Seq\* color
  scheme from the viewer's side; it lists every polymer alphabet (Seq\* judges fit by the first block). The menu offers
  it, Zappo, Taylor, Hydrophobicity and none (`LetterSchemes`; none: plain cells, `wrappedColors` off).
- **The legend** (`legendOf`, `withTracks`, `sequence/legend.ts`) is read from the scheme's public definition
  (`KernelColorSchemes`), plus the secondary structure's classes and the behaviors' offered blocks; shown in prism-ui's
  `Legend` in an overlay card (`LegendCard`, `sequence.legend`). Seq\*'s own legend model stays inside its UIs.
- **One hover label, Mol\*'s where it can be**: a sequence hover the link maps into the structure lights its residues in
  3D, and Mol\*'s label providers label that highlight, so the label is the same as hovering those atoms in 3D
  (`pickHoverLabel`, `sources/hover-lines.ts`; `SequenceSource.hoverMapped` says the hover mapped; dropping the link
  clears the 3D highlight, and with it `hoverMapped`). **What the sequence knows is added to it**: the link registers a
  Mol\* `LociLabelProvider` (`SequenceView3D.annotateHover`, priority −1000, so after every other) that, while the hover
  maps, adds the track's own lines (an annotation and its fields, a value; under the track's name, `HoverLabel.track`,
  unless a line names it), escaped and on one line each. Never Seq\*'s title or the residue's lines (providers
  `position` and `copies`: numbers, copies, the stored letter; but `copies: varies`, which Mol\*'s label of one residue
  cannot say): Mol\*'s label says those, and a hover with nothing else adds nothing; a new Seq\* label over the same
  highlight labels it again (`LociLabelManager.refresh`), unless a highlight write is on its way (it labels itself with
  the new lines when it lands: no second label, no new lines on the last residue). Any reader of `labels.highlight` gets
  them, and Mol\*'s provider API is unchanged (`refresh` is new). The band's own Seq\* label (`seq.labels.current`, at
  most once a frame and only when its lines change) shows only for a hover that does not map (a gap, an annotation track
  with no counterpart in the structure). The shell's `HoverLabel` listens on the viewport (`data-msp-hover-source`), in
  the corner or at the cursor (`options.ui.hoverLabel`).
