# The sequence view and the 3D view: one state, shared data

How the Prism viewer keeps Seq\* and the Mol\* plugin in step: what each owns, how data is read once, how what the
sequence shows is decided (by behaviors the viewer and extensions add), how a sequence action reaches the 3D view, and
how both document languages, MolViewSpec and SeqViewSpec, fit. Builds on the first slice
([../../seq/prism-sequence.md](../../seq/prism-sequence.md)) and on Seq\*'s host contract
(`packages/sequences/_spec/viewer.md`, `sources.md`, `libraries.md` § The Mol\* adapter).

Status: **accepted** (2026-10-15).

**As built (plan 2.4, 2026-10-15).** `sources/sequence.ts` follows a structure (`sequence.follow` overrides the rule),
hands the sequence the plugin's own parse (`StructureSource.of` on the trajectory the 3D view drew), and is the link's
3D side: hover and selection read from the plugin, highlight, selection and focus written as operations with origin
`link`, one channel at a time in the order sent. The track menus' coloring is `representation.color`. The `sequence`
area (`style`, `colors`, `follow`, `get`) exists when the spec has a sequence; its state is `prism-viewer.sequence`
(§7.1). Still the first slice's: the document is `sequenceDocument` (the follow and secondary-structure rules in one
builder); the behaviors of §3 come with the first extension that needs one (model-archive's pLDDT track, plan phase 5),
and the loaders' `load` (§2) with `load.url` (the `loaders` feature, phase 3).

## 1. Two engines, one viewer

| Engine                        | Owns                                                                                                    | Its document                     |
| ----------------------------- | ------------------------------------------------------------------------------------------------------- | -------------------------------- |
| Mol\* (`PluginContext`)       | data cells, parses, structures, representations, camera, canvas; hover, selection, focus as loci        | the state tree; MVS builds it    |
| Seq\* (`SeqViewerContext`)    | the sequence views (letters, tracks, alignments), their layout and viewports; hover, selection, focus   | a SeqViewSpec document (layered) |

Neither engine knows the other. The viewer joins them in three places, each one way of doing it:

1. **Data** goes through Seq\*'s `Loaders` port, which the viewer implements over the plugin (§2).
2. **Interaction** (hover, selection, highlight, focus) goes through the link: `StructureLink` maps chain instances
   onto Seq\* units, `linkViews` carries the channels both ways, origins stop echoes (§4).
3. **Content** (what the sequence shows) comes from **sequence behaviors**: the viewer's and the extensions' rules
   that turn the viewer's state into patches of the SeqViewSpec document (§3).

Actions that change the 3D view from the sequence go through the viewer's API, never straight to the plugin (§5).

## 2. Data: read once, through the plugin

Seq\* already has the ports: `Loaders.load` (bytes for a URI or a host file), `Loaders.parsed` (a parse the host shares:
`StructureSource.of(trajectory, format)`) and `Loaders.stream` (streamed sources). The first slice uses `parsed` for the
structure the 3D view shows and refuses every `load`. The viewer's loaders, in its `sources/` folder:

- **`parsed(request, format)`**: when the plugin holds a parse of that file (a trajectory cell under a data cell whose
  URL, or host file, is the request's), hand it over. One parse, one `Structure`, so the link maps the very objects the
  3D view draws.
- **`load(request)`**: when the plugin holds the bytes (a data cell from that URL), hand them over. Otherwise download
  **through the API** (`load.url` with origin `link`, so the cell it makes is a recorded step, not a change outside the
  API, §7), so the 3D view finds it when it loads the same file next, and a session saves it once.
- **Two caches, both read-through** (David, 2026-10-15):
  - **Bytes, in core**: a URL-keyed, ref-counted cache in Mol\*'s asset manager (its "TODO: add URL based ref-counted
    cache?"): a miss downloads, a hit returns the held bytes, an entry lives while something holds it. MVS's downloads,
    the plugin's `Download` and the viewer's loaders all read through it, so a URL is fetched once whoever asks first.
  - **Parses, in the viewer**: an index (`dataByUrl`, a source) of the plugin's data cells by URL and host file with
    their parse cells, kept from the state tree's cell events. A miss parses through the plugin (a data cell and its
    parse in the state tree), so the 3D view and a saved session find it.
- **As built (`load.svs`, 2026-10-15)**: `parsed` answers the followed structure's host file, which a loaded SVS's
  first `molstar:structure` complex is rewritten to read (§6); `load` still refuses (an SVS's other URIs are reported as
  issues), and the caches are not built.
- **Seq\* side**: no change needed for files and parses; services (`uniprot`, `pdbe`, `afdb`) fetch through `load`
  too, so they share the same cache. Streamed sources (`stream`) are the host's provider and do not touch the plugin.

## 3. Content: sequence behaviors

What the sequence shows is a SeqViewSpec document. In the viewer it is made, and kept up to date, by **sequence
behaviors**: small rules that read the viewer's state and propose patches.

```ts
interface SequenceBehavior {
  id: string;                                         // 'sequence.follow-focus', 'rcsb-validation.track'
  attach(ctx: SequenceBehaviorContext): () => void;   // returns its dispose
}
interface SequenceBehaviorContext {
  viewer: ViewerSources;                              // structures, focus, selection, custom properties (read-only states)
  seq: SeqViewerContext;                              // its document (applyPatch on a layer), its states
  origin: `host:${string}`;                           // every patch it makes carries it
}
```

- **The viewer's own**: follow the structure the person works with (focus, then latest selection, then in view, then
  first loaded; prism-sequence.md §3); secondary structure.
- **Extensions add theirs** (`extension.sequence?: SequenceBehavior[]`, README §3.1): model-archive's pLDDT track and
  PAE dot plot for predicted models (README decision 27), a validation report as a track (rcsb, pdbe), membrane
  orientation, symmetry, interactions around the focus. The composition overrides them by id, as
  any part.
- **As built** (`sequence/behaviors.ts`, plan 5.1): `attach(ctx)` with `ctx.shown` (what the document shows, a new value
  per document: the structure, the complex's ref, the `data` node's, the section and the viewer's last track), `ctx.patch`
  (onto that document only: a patch for a replaced one is dropped) and `ctx.offer` (its tracks' menus, by track ref and
  color theme, and its legend blocks). The viewer's own content (ruler, letters, secondary structure) is still the
  document the model builds; the follow behavior and the person's overlay are not split out yet.
- **They write the base layer by patch** (Seq\*'s `applyPatch`), so untouched nodes keep their identity and only the
  changed ones rerun. SeqViewSpec patches are all or nothing and nothing cascades (`packages/sequences/_spec/format.md`,
  review F23), so:
  - **Stable refs**: the follow behavior owns the data and the complex under fixed refs (`structure`, `complex`) and
    switches structures by setting their fields, never by removing and inserting them, so other behaviors' tracks that
    refer to the complex stay valid.
  - **A behavior removes what refers to what it removes**, in the same patch; each behavior's own nodes carry its id as
    a ref prefix.
  - **A composition test** runs every built-in behavior together over a structure switch, a load and a clear.
- **The person's edits** (a track hidden, a color chosen, a track's order) are a viewer-owned overlay keyed by track
  ref, applied after the behaviors on Seq\*'s `local` layer and saved with the sequence's state (§7); an edit to a
  track that no longer exists is dropped and said once.

## 4. Interaction: the link

As built in the first slice ([../../seq/prism-sequence.md](../../seq/prism-sequence.md) §2.3), generalized:

- **One link per shown complex**: `StructureLink.attach(structure)` over the structure the 3D view draws (§2's shared
  parse), `linkViews(seq, view3d, link, options)` for hover and selection both ways, highlight and focus from Mol\*'s
  loci; `pluginView3D` lives in `sources/`. The link's writes into Mol\* (selection, highlight, the camera's focus)
  are view changes, not steps (operations.md §1, `history: false`): the viewer's `pluginView3D` makes them through the
  same operations with origin `link` (`select.set`, `highlight.set`, `camera.frame`), never straight on the plugin
  (review F21).
- **Origins**: Seq\*'s changes carry `pointer`, `keyboard`, `api`, `graph`, `host:<tag>`; Mol\*'s side tags what the
  link sends, so neither re-sends its own echo. Replays and carries (`origin: graph`) lead nothing.
- **State lives in the engines.** Selection is Mol\*'s (`managers.structure.selection`), mirrored into Seq\*; focus is
  Mol\*'s (`managers.structure.focus`); the sequence's own viewport and style are Seq\*'s. The viewer's view models read
  both through sources; nothing keeps a third copy.
- **Size gates** stay `linkViews`' defaults (hover 2,000 items, selection 20,000; past them, whole chains).

## 5. From the sequence back to the viewer

Everything the sequence can do to the 3D view is a call to the viewer's API (README §3.3), so it is attributed,
optimistic where it can be, part of the turn's checkpoint, and the same as the command and the tool that do it:

| In the sequence                                           | API call                                                  |
| --------------------------------------------------------- | --------------------------------------------------------- |
| click, drag, Shift (selection)                            | none: the link mirrors it into Mol\*'s selection (§4)      |
| a track's menu: Show in 3D, Color by this track           | `camera.frame(target)`, `representation.color(target, theme)` |
| a structure's name (the strip's one icon)                 | `view.open({ tab: 'structure', reveal: structure })`       |
| a SeqQL result                                            | `select.set(target)`, `focus.set(target)`                  |
| hover card actions (later)                                | the same calls                                            |

`HostActions` (molstar-poc's `structureTrackActions`) take these calls instead of touching the plugin. A **target**
is the shared address (README §3.3): Seq\* positions become residue addresses through the link's own mapping
(`StructureLink`, Seq\* → Mol\* at edge rate; `coordinatesOf` is the other way, Mol\* runs onto Seq\* coordinates), so
a sequence action, a command and a tool name the same residues the same way.

## 6. Both documents: MolViewSpec and SeqViewSpec

- **SeqQL and MolQL meet in the link**: `molstar.select` runs a MolQL (PyMOL, VMD, Jmol, MolScript) selection on the
  linked structure and returns residues on Seq\* units; a SeqQL result goes back to 3D through the link's mapping. The
  agent's tools build either language and get one target.
- **Loading an MVS** builds the plugin state; the sequence follows it through the behaviors (§3), so an MVS needs no
  sequence part to get one.
- **Loading an SVS** (`.svsj`, `.svsx`) sets the sequence document; its `molstar:structure` complexes name their sources
  by URI, which the loaders (§2) resolve against the plugin first, so a structure the 3D view has is not fetched again.
  Whether an SVS also loads its structures into 3D is the follow behavior's choice (yes, by default, with the default
  representation); it does so through `load.svs`, an operation, so those cells are recorded steps (§7).
- **As built (`load.svs`, the loaders feature; David's decisions, 2026-10-15)**: a URL (`.svsj`/`.svs`, JSON or
  MessagePack) or data; `.svs`/`.svsj` files go through `load.files` after the files dropped with them; `svs-url` in the
  app's URL. Its **first `molstar:structure` complex** is matched to a loaded structure by its source: the address
  (resolved against the document's) of a Download cell, the entry id (`service` `pdbe`/`afdb` and `id`), or the file
  name (`host_file`, a data cell's or the model's entry id); none matches: it is loaded (`load.url`, `load.entry`,
  `load.data`) in the same step. The document then **replaces the followed structure's document** (the model's
  override, `SequenceOverride`: `build()` uses a copy whose complex reads the host's parse, with refs given to the data
  node, sections and last track so behaviors such as the pLDDT track attach by `shown`), and the source shows its
  structure (`SequenceSource.show`) **until the person focuses another structure, a structure is loaded, the sequence
  follows another (`sequence.follow`), or the structure goes**; then the built document is back. The person's style
  (Text/Tracks, on every section) and letter colors apply over its presets. `sequence.get` names it. An `.svsx` (a URL, a
  file, dropped) is read with Seq*'s `readSvsx`: a structure it holds at a relative URI opens as a file of that name, and
  its relative URIs resolve inside it (`archiveLoaders`, the override naming the archive by id; the files stay in memory
  for the page's life, out of the saved state). A document without a `molstar:structure` complex is refused.
- **Not both at once in the viewer** (David, 2026-10-15): the viewer loads an MVS or an SVS. A document holding both is
  Stories' (its schema-based document), or a new format later (a "Mol states" zip, `.msx`), not part of the viewer's
  first scope.
- **One sequence viewer** (David, 2026-10-15), switched between Text and Tracks as in the mocks; the keyed model
  (prism-sequence.md §2.1) stays, with `'default'` its only entry.

## 7. Snapshots, history and undo

### 7.1 The sequence's state in snapshots

- **What is saved**: the plugin's state (`PluginState.Snapshot`) and the sequence's own state,
  `{ version: 2, style, colors, structure, rail, override?, edits }` (`edits`: the overlay of §3; not built yet: the
  saved state has no `edits` so far). `override` (version 2, `load.svs`, §6) is the loaded SeqViewSpec: the document
  itself (data has no address to fetch again), its URL when it had one, and the structure's state ref; Undo's
  checkpoints keep it (it is content, not view state) and a session restores it. Version 1 is read as version 2 without
  `override`. The document itself is rebuilt by the behaviors from the plugin's state.
- **A registry, not a field** (review F10): `PluginState.registerCustomState(key, { get(), set(value | undefined) })`,
  an additive core change. `getSnapshot` writes every registered key into `snapshot.custom`; `setSnapshot` calls each
  `set` after the tree is applied, with `undefined` (reset to defaults) when the snapshot has no value for it; unknown
  keys pass through unchanged; transition frames do not apply it. So the snapshot manager's entries, their playback
  and `.molx` / `.molj` sessions carry the sequence's state with no other hook. The viewer registers
  `'prism-viewer.sequence'`; a value of another `version` is dropped and said once.

### 7.2 Checkpoints

The history's checkpoints use the plugin's state directly, never its snapshot manager (whose `getStateSnapshot` adds
or replaces entries and whose restore clears; review F1; decision 24):

```ts
interface HistoryPort {                                   // given to prism-ui's history manager by the viewer
  checkpoint(): Checkpoint;                               // synchronous, cheap; measured (§7.4)
  restore(c: Checkpoint): Promise<void>;
}
// the viewer's:
checkpoint = () => ({
  plugin: withoutViewState(plugin.state.getSnapshot({ data: true, behavior: true, componentManager: true,
    canvas3d: true, interactivity: false, structureSelection: false, camera: false, canvas3dContext: false,
    animation: false })),                                 // the sequence's state rides in `custom`, minus what it follows
});
restore = (c) => plugin.state.setSnapshot({ ...c.plugin, structureFocus: undefined })
```

- **What undo restores**: the data tree (structures, components, representations, measurements), behavior settings,
  canvas props, the sequence's style and colors, the hidden contents. **What it never touches** (view state, not
  steps, §7.3): the camera, the selection, the selection mode and the picking level, the focus, which structure the
  sequence follows, hover, the Snapshots list, open panels.
- **Without a canvas** (headless): `camera: false` keeps `getSnapshot` from reading `canvas3d`; a headless test runs an
  undoable operation, then Undo, with no canvas.
- **Across a lazy extension's load**: restoring a checkpoint taken before the load keeps the extension loaded (its
  registrations stay); only its cells go.

### 7.3 One history manager

`managers.history` (prism-ui, README §3.3) holds every undoable step, the person's and the agent's (David, 2026-10-15):

- **A step** is one undoable operation run (`api.run`, operations.md §1) or one group. A group is opened with
  `beginGroup(origin)` and closed with `endGroup()`: an agent's turn (keyed by `origin.turn`), a script's
  `viewer.run(label, …)`, a builder's commit. One checkpoint is taken when a step or group opens.
- **Not steps**: operations with `history: false` (`select.*`, `focus.*`, `camera.*`, `highlight.*`, `view.*`,
  `trajectory.frame`, playing animations, reads; Undo leaves the frame shown);
  a result that changed nothing (`unchanged`: "it is already so"); an operation that changes nothing the history could
  restore (`undoable: false`, an upload: no checkpoint, recorded only as a run, so Undo of the steps before it still
  works).
- **Inside a step**: each top-level run or group gets a history context; work that calls back in with it (an
  operation's `ctx.api`, a group's API, `viewer.run`'s `api`) is part of that step, with no lock and no step of its own.
  Calls without it (another origin's) wait their turn and are steps of their own, so a click during an agent's snippet
  is never rolled back with the snippet.
- **While a turn runs**, the person's own calls still run; the first one closes the turn's group (the turn's later
  calls open a new group under the same turn), so steps never interleave inside one entry. Undo turn undoes the
  turn's groups only while they are the latest entries; otherwise it asks before discarding the later steps.
- **Undo while a run is in flight** (or queued) is disabled, never queued behind it (the button and the command say
  why); Stop aborts the run first.
- **A failed or aborted run** restores the checkpoint its step took, so a loader that failed after its download leaves
  nothing half-built; the result says so (critic, run outcomes).
- **No redo** for now: entries hold only the checkpoint before them. (Listed in operations.md as later.)
- **Changes outside the API** (review F5): inside are runs under `api.run` and `viewer.run`, the history's own
  restores, cells owned by behaviors (allow-listed by transformer: the focus representation, volume streaming,
  custom-property and label cells), and the snapshots' playback (a presentation: Undo after it restores the step's
  checkpoint as usual). Any other change to the data tree makes every checkpoint taken before it unsafe:
  those steps can no longer be undone, and History says so on them.
- **The record**: every run and step is also an entry in the session record ([journal.md](journal.md)).
- Mol\*'s own state history is off (`PluginConfig.State.HistoryCapacity` 0) and Seq\*'s too (its document port with
  `history: false`, as Stories sets it): one undo, the viewer's.

### 7.4 Large scenes

The history measures each checkpoint as it takes it: time, and the size of the snapshot's JSON without data cells'
payloads. When one goes over the size budget, or three in a row over the time budget (one slow checkpoint is a cold
start or a pause, not the scene; a first guess: 50 ms or 5 MB, tuned on a mesoscale session), per-step checkpoints
turn off for that scene, and the status column says so once: "This scene is large: undo is off; Settings
to turn it back on." The agent's turn keeps its one checkpoint unless that one is over a hard limit too, which the turn
card says before the turn runs. Clearing the scene, or loading a smaller one, turns per-step undo back on.

## 8. Open questions

None at the moment.
