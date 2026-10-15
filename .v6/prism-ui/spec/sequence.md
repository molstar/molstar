# The sequence view and the 3D view: one state, shared data

How the Prism viewer keeps Seq\* and the Mol\* plugin in step: what each owns, how data is read once, how what the
sequence shows is decided (by behaviors the viewer and extensions add), how a sequence action reaches the 3D view, and
how both document languages, MolViewSpec and SeqViewSpec, fit. Builds on the first slice
([../../seq/prism-sequence.md](../../seq/prism-sequence.md)) and on Seq\*'s host contract
(`packages/sequences/_spec/viewer.md`, `sources.md`, `libraries.md` § The Mol\* adapter).

Status: **accepted** (2026-10-15).

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
  **through the plugin** (a data cell in the state tree, `Download` with the URL), so the 3D view finds it when it
  loads the same file next, and a session saves it once.
- **Two caches, both read-through** (David, 2026-10-15):
  - **Bytes, in core**: a URL-keyed, ref-counted cache in Mol\*'s asset manager (its "TODO: add URL based ref-counted
    cache?"): a miss downloads, a hit returns the held bytes, an entry lives while something holds it. MVS's downloads,
    the plugin's `Download` and the viewer's loaders all read through it, so a URL is fetched once whoever asks first.
  - **Parses, in the viewer**: an index (`dataByUrl`, a source) of the plugin's data cells by URL and host file with
    their parse cells, kept from the state tree's cell events. A miss parses through the plugin (a data cell and its
    parse in the state tree), so the 3D view and a saved session find it.
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
  first loaded; prism-sequence.md §3); the pLDDT track for predicted models; secondary structure; PAE as a dot plot.
- **Extensions add theirs** (`extension.sequence?: SequenceBehavior[]`, README §3.1): a validation report as a track
  (rcsb, pdbe), membrane orientation, symmetry, interactions around the focus. The composition overrides them by id, as
  any part.
- **They write the base layer by patch** (Seq\*'s `applyPatch`), so untouched nodes keep their identity and only the
  changed ones rerun. The person's own edits (a track hidden, a color chosen) are patches too, from the UI, through the
  API.
- **Ordering**: behaviors patch disjoint parts of the document (each its own tracks, by ref prefix); the follow
  behavior owns the data and the complex. A conflict is a bug the behavior's test catches, not a merge rule.

## 4. Interaction: the link

As built in the first slice ([../../seq/prism-sequence.md](../../seq/prism-sequence.md) §2.3), generalized:

- **One link per shown complex**: `StructureLink.attach(structure)` over the structure the 3D view draws (§2's shared
  parse), `linkViews(seq, view3d, link, options)` for hover and selection both ways, highlight and focus from Mol\*'s
  loci; `pluginView3D` lives in `sources/`.
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
| a track's menu: Show in 3D, Color by this track           | `view.show(target)`, `structure.color(target, theme)`     |
| a structure's name (the strip's one icon)                 | `navigation.reveal('structure-tree', structure)`           |
| a SeqQL result                                            | `selection.set(target)`; `focus(target)`                  |
| hover card actions (later)                                | the same calls                                            |

`HostActions` (molstar-poc's `structureTrackActions`) take these calls instead of touching the plugin. A **target**
is the shared address (README §3.3): Seq\* loci become residue addresses on the link's chains (`coordinatesOf`), so a
sequence action, a command and a tool name the same residues the same way.

## 6. Both documents: MolViewSpec and SeqViewSpec

- **SeqQL and MolQL meet in the link**: `molstar.select` runs a MolQL (PyMOL, VMD, Jmol, MolScript) selection on the
  linked structure and returns residues on Seq\* units; a SeqQL result goes back to 3D through the link's mapping. The
  agent's tools build either language and get one target.
- **Loading an MVS** builds the plugin state; the sequence follows it through the behaviors (§3), so an MVS needs no
  sequence part to get one.
- **Loading an SVS** (`.svsj`, `.svsx`) sets the sequence document; its `molstar:structure` complexes name their sources
  by URI, which the loaders (§2) resolve against the plugin first, so a structure the 3D view has is not fetched again.
  Whether an SVS also loads its structures into 3D is the follow behavior's choice (yes, by default, with the default
  representation).
- **Not both at once in the viewer** (David, 2026-10-15): the viewer loads an MVS or an SVS. A document holding both is
  Stories' (its schema-based document), or a new format later (a "Mol states" zip, `.msx`), not part of the viewer's
  first scope.
- **One sequence viewer** (David, 2026-10-15), switched between Text and Tracks as in the mocks; the keyed model
  (prism-sequence.md §2.1) stays, with `'default'` its only entry.

## 7. Snapshots and undo

- **A viewer snapshot is both engines' state**: the plugin's (`PluginState.Snapshot`: the data and behavior trees,
  camera, canvas, interactivity, focus, selection, components) and the sequence's own choices (style, colors, the
  structure followed). The sequence document itself is rebuilt by the behaviors (§3) from the plugin's state, so it
  need not be saved.
- **In `.molx` and `.molj`** (David, 2026-10-15): the plugin's snapshot gets one optional, owner-keyed field,
  `PluginState.Snapshot.custom?: Record<string, unknown>` (an additive core change), written and read through the
  snapshot manager, so every snapshot entry and a whole session carry it. The viewer keeps the sequence's state under
  `'prism-viewer.sequence'` (`{ style, theme, structure, rail }`); the classic UI and older readers ignore the field, and
  a session without it opens with the sequence's defaults.
- **One history manager** (`managers.history`, David, 2026-10-15) holds every undoable step, the person's and the
  agent's alike: an undoable API call is one entry, an agent's turn is one entry however many calls it made. An entry
  keeps a checkpoint taken before it (the plugin's snapshot through the snapshot manager, and the sequence's choices);
  Undo restores the latest. The history is bounded (a ring); a checkpoint is the state tree's transforms and params,
  not the data, and restoring is a tree diff, so unchanged cells stay. It also records the runs (what ran, its origin,
  its result: contexts.md §5), which History and the eval's run log show.
- Mol\*'s own state history is off (`PluginConfig.State.HistoryCapacity` 0) and Seq\*'s too (its document port with
  `history: false`, as Stories sets it): one undo, the viewer's; `canUndo` is not used.

- **Large scenes** (David, 2026-10-15): the history measures each checkpoint as it takes it (time and size). When one
  goes over the budget (a first guess: 50 ms or 5 MB, tuned by measuring a mesoscale session), per-step checkpoints turn
  off for that scene, and the status column says so once: "This scene is large: undo is off; Settings to turn it back
  on." The agent's turn keeps its single checkpoint (one held at a time, released at the next action), unless that one
  is over a hard limit too, which the turn card then says before the turn runs. Clearing the scene, or loading a
  smaller one, turns per-step undo back on.

## 8. Open questions

None at the moment.
