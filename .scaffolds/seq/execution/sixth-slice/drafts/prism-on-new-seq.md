# Prism on the new Seq\*: the sequence band's plumbing swapped, its behaviour kept

**Draft, 2026-11-23**: not agreed. The sixth slice's "Prism on the new Seq\*" entry ([plan](../plan.md)). Read-only
survey of `packages/prism/viewer` (its `_spec/sequence.md` §1–10, `src/sequence/**`, `src/api/sequence.ts`,
`src/sources/{sequence,sequence-scene,data-by-url,hover-lines,checkpoint}.ts`) against `packages/seq/{core,viewer}` and
their `_spec`s.

## Intent

The author, 2026-11-23: "as part of slice 6, we should replace the old seq. viewer in prism with the new one, so we
should plan that. should be fairly straightforward if we did things right." And, the same day: "One thing for the new
seq in prism: it should use the new "data port" to share the representation with the viewer".

So the swap has a centre: **Prism's data path becomes a Mol\* adapter of the seq data port**
([seq-data-ports.md](../../../design/seq-data-ports.md)): the plugin's structures served to Seq\* as sources from the
plugin's in-memory `Structure`, through `StructureDescription` → `structureNodes`, never a parse of their files. That
design deferred the Mol\* adapter ("no Mol\* adapter", decided 2026-11-19); this is where it gets built. Around it the
rest is plumbing ([the engine, not the game](../../../design/README.md)): what the band shows, how it is chosen, menus,
the tree node, the skin, the history stay Prism's, and `sequence.md`'s behaviour is the parity bar.

## Precedents

- **Today's Prism on the old Seq\*** (`@molstar/seq-core`, `@molstar/seq-viewer`, about 6,000 lines with tests). Its
  Mol\* adapter (`src/sequence/molstar/`, 2,400 lines) is a Seq\* *provider* set: a `molstar:structure` complex
  provider (`complex.ts`), computed providers `molstar:plddt`, `molstar:secondary_structure`, `molstar:pae`,
  `molstar:sifts_coverage`, the `afdb_pae` format and `pdbe`/`afdb` services (`registry.ts`), a `Loaders.parsed` hook
  answering `structure:<key>` with the plugin's parse (`model.ts` constructor), and a link that registers Mol\* chain
  instances as external coordinate systems on a `Translator` (`link.ts`, `link-views.ts`, `seq-units.ts`). Taken: one
  parse, the frame rule, the size gates, ligands whole, the hover-label rule, everything Prism-side. Not taken: computed
  providers (the new language has no `computed` kind: dynamic nodes come from sources), Mol\* `Loci` inside Seq\*'s
  coordinate machinery (the new link speaks Locations on units).
- **The new Seq\*'s data port** (`seq-data-ports.md`, core `_spec` "Seq io", "Data", Protocol; `context/adapters.ts`):
  a host adapter is a `FormatRuntime` put in place of the default by `withAdapter`; `structureNodes` gives keys from
  identifiers; `OpenedSource.changes` + the store's `reattach` keep identity across a host's edit; `values.field` and
  `source.ref`; outward is pure (`unitRows`, `colorRuns`). Built and exercised only by the Storybook's fake host
  (`packages/seq/storybook/src/_stories/support/fake-structure-adapter.ts`, values by field only).
- **The new Seq\*'s structure syncing** (`structure-syncing.md`; core `_spec` "Structures", viewer `_spec` "Structure
  link"): `structureFacts`, `toUnits`, `fromUnits`, `linkStructure(viewer, side)`; the fake host's side translation
  (`.../support/host-adapter.ts`) is the reference a real host follows.
- **Mol\*'s layering** (`@molstar/io` → `@molstar/model`): convert at the edge onto one shape. Taken, as the port did.

## Decision

### 1. The Mol\* data adapter (the centre)

- **One format runtime, `molstar:structure`, in place of nothing** (no core runtime reads it; `withAdapter` also takes
  `mmcif`/`bcif` so a loaded spec naming a file the plugin holds is served from the plugin's parse, not read again).
  `open(entry)`: by `entry.ref` (`structure:<key>`, Prism's host key: a structure with no file, a predicted or edited
  model) or by `entry.uri` normalized as the asset manager keys it, looked up in `dataByUrl` (unchanged,
  `sources/data-by-url.ts`) → the plugin's `Structure` of that trajectory → `describeStructure(structure)` →
  `structureNodes(description, { source })`. Why: the port's rule, one path; Seq\*'s nodes from the model the 3D view
  draws, so a unit's residues are its atoms' residues by construction.
- **`describeStructure`** (new, Prism-side, `@molstar/model` only): entities, components, identifiers, `chains` (one row
  a `label_asym_id` of the **deposited model**, no operator: today's rule, `chain-filter.ts`, `complex.ts`; an
  assembly's copies are their chain once), `residues` (`auth_seq_id`, `ins_code`, `observed`), `links` from
  `struct_conn`, compounds, heavy `atoms` for non-polymers, and named per-residue `values`: `b_factor` always, and any
  channel an extension registers (below). Reads `Model` hierarchy and properties, never a CIF frame; SIFTS is the one
  exception (§1 gaps). Water out (`structureNodes`' default).
- **Tables by region**: `values` by `field` (B-factor; `plddt` from model-archive) and `ranges` of the secondary
  structure (`kind` = `helix`/`strand`/`turn`, from `ComputedSecondaryStructure` on the drawn structure, else the
  model's), so the built-in `secondary` layer draws them. Why ranges from the source, not a computed node: the language
  has none, and a host's secondary structure is the port's own example ("By region").
- **Updates**: the adapter's `changes(fn)` fires when the structure under a key gets other atoms or another hierarchy
  (today a new key: §2 "A frame keeps the document") **or** the 3D view's secondary-structure assignment changes
  (`HostStructure.assignment`, today a param re-set on a computed node). The viewer reopens in place a frame later,
  `reattach` keeps node identity by key, the selection survives by key. **A trajectory frame fires nothing** (same
  atoms: nodes, values and secondary structure stay the frame the document was made for, as today). A removed
  structure disposes its source.
- **Value channels are extensible**: `ViewerExtension.sequence.values?: { name, applies(structure), read(structure) →
  Float32Array per residue }`, read into the description's named values. Why: model-archive's pLDDT track today is a
  `computed` node over `molstar:plddt`; the new form is a `values` node `{ source, field: 'plddt' }`, and the viewer
  itself knows nothing of pLDDT (README decision 27).
- **Nothing fetched**: the adapter has no fallback for structures (a `uri` the plugin does not hold is `unavailable`,
  its message "load it in the viewer first", today's `svs.source.load`); other formats a loaded spec names read through
  a `ReaderAssetsPort` over `dataByUrl.bytes` (and an archive's entries), refusing the rest (`AssetRefused` →
  `seq.source.unreadable`). Why: §2's rule, "the sequence fetches nothing".

**Where it lives: the import boundary decides.** `seqx-core` never imports `@molstar/model` (its allowlist,
`packages.md`; seq-data-ports "Layering"); the adapter must. It does not need `@molstar/plugin`: what it reads is a
`Structure` and its custom properties. So two layers:
- `describeStructure`, the Mol\* format runtime and the link's translation (`Loci` ↔ `StructureTarget`): imports
  `@molstar/model` and `@molstar/seqx-core/{data,io/structure,structure}` only; **no plugin, no Prism**, held by a
  boundary test (as `src/_test/boundary.test.ts` holds today's). Today `src/sequence/molstar/`, exported as
  `@molstar/prism-viewer/sequence/molstar`; movable unchanged to a bridge package (`@molstar/seqx-molstar`, the port's
  name for it) when a second consumer appears (Question 1).
- The plugin side (which structure, `dataByUrl`, `changes` from plugin state events, `assignment`): Prism's
  (`sources/sequence.ts`, `sequence/model.ts`).

**Given back as the host's.** The port's outward functions replace today's link code where they can:
- **Hover and selection**: `linkStructure(viewer, side)` (`@molstar/seqx-viewer/link-structure`) does the echoes,
  ownership (EC-141), coalescing, gates and clearing that `link-views.ts` does today. Prism's `side` is a thin
  translation over `structureFacts` (residue rows joined once per document against the structure's residue index):
  a run on a polymer unit → `StructureElement.Loci` of that chain's residues `label_seq_id = position + 1` **on every
  operator of the chain** (deposited units, assembly copies: today's rule); a molecule unit → its residue whole; inward
  the reverse, any operator collapsing onto the chain's unit. Writes stay operations with origin `link` (§4).
- **Colours**: "Color by this track" keeps today's named Mol\* themes (`HostActions.colorStructure`); `colorRuns` is
  what Prism's TODO "Coloring a structure by data" (`selector-colors`) would take for any track. Not in this swap
  (a Mol\* theme, not Seq\* plumbing).
- **Annotations**: `unitRows` for a loaded spec's ranges onto the structure (a future "show domains in 3D"); no user
  today. The hover label's track lines come from the card (§10 below), not `unitRows`.

### 2. What `src/sequence/molstar/*` becomes

| File (lines) | Becomes | Mark |
|---|---|---|
| `source.ts` (143: headless parse of bytes) | gone: no parse; the plugin's `Structure` | replaced |
| `complex.ts` (325: complex provider, numberings) | `describeStructure` + `structureNodes` | replaced (rewritten smaller) |
| `chains.ts` (143: chain instances) | the residue index `describeStructure` and the side share | adapt (kept, trimmed) |
| `seq-units.ts` (105), `link.ts` (309: Translator systems), `link-views.ts` (228) | `structureFacts`/`toUnits`/`fromUnits` + `linkStructure`; Prism's side (~150) | replaced |
| `plddt.ts` (62) | a `values` channel (`b_factor`/`plddt` field) | replaced |
| `secondary.ts` (278) | `ranges` by region from the adapter; `secondaryCategory` kept | adapt |
| `registry.ts` (107: formats, services, providers) | `withAdapter(registry, MolstarAdapter)` | replaced |
| `structure-ref.ts` (49), `coverage.ts` (255), `sifts.ts` (102) | SIFTS: gap G2; coverage, "Show in 3D": stated loss (Q5) | gap / loss |
| `pae.ts` (126) | stated loss in the band (Q5); `pae-matrix.ts` (56, no Seq\*) stays for model-archive's card | loss / kept |
| `query.ts` (130: `molql` in SeqQL) | core query entries (`query-form`); the viewer's address model passes none yet | gap G6 (or loss) |
| `index.ts` | the new exports | adapt |

### 3. The mapping: `sequence.md` and `src/sequence/**`

*direct*: exists, a rename. *adapt*: exists, Prism-side code changes. *gap*: Seq\* lacks it (a Seq\* lane, Gn below,
or a stated loss).

| § / file | Today | On the new Seq\* | Mark |
|---|---|---|---|
| §1 two engines | joined by `Loaders`, link, behaviors | the data port (§1 above), `linkStructure`, spec patches | adapt |
| §1 dependencies | `@molstar/seq-core`, `@molstar/seq-viewer` | `@molstar/seqx-core`, `@molstar/seqx-viewer`, **and transitively `@molstar/ui`** (its `dom` frame tick, reactive, menus) | gap G7 |
| §1 adapter export | `sequence/molstar` | the Mol\* adapter (§1), same export | adapt |
| §1 `createSequenceViewers` (`model.ts`) | one `SeqViewerContext` per id | `createSeqViewer(spec, target?)`; headless without a target | direct |
| §1 `PrismSequenceSpec` (`spec.ts`) | Seq\*'s default spec + registries + render props | spread `DefaultSeqViewerSpec()`, trimmed (no `editing`, `keysPanel`, story, `uiHistory`), `withAdapter`, Prism's scheme, config items (`viewerItem`) | adapt |
| §2 runtime chunk (`runtime.ts`) | dynamic import, `rafScheduler` | same import; the scheduler is `@molstar/ui/dom`'s (G7) | adapt |
| §2 one parse | `Loaders.parsed` → `StructureSource` | the adapter by `ref`/`uri` | adapt |
| §2 a frame keeps the document | rebinding `StructureLink` | no `changes`; the side reads the current frame | adapt (simpler) |
| §2 byte cache | core's asset manager | `ReaderAssetsPort` over `dataByUrl` (refusing) | adapt |
| §2 `dataByUrl` (`sources/data-by-url.ts`) | parses by address | unchanged, read by the adapter | direct |
| §2 `load` refuses | `refuse()` in `Loaders` | adapter `unavailable` + refusing reader | adapt |
| §3 follow, Auto (`sources/sequence.ts`) | Prism | Prism | direct |
| §3 built document (`document.ts`) | old `createSeqViewSpec`, `computed` nodes | `createSeqSpec` (`@molstar/seqx-core/lang`): `source` (`molstar:structure`, `ref`), one view, section `main` `of` the complex | adapt |
| §3 chain filter (`chain-filter.ts`) | blocks, `copies` | section `classes` (`[polymer]`), `copies: one|each`, block children (overrides); None = a block of a unit there is not (check `seq.lang` accepts it) | adapt |
| §3 ligands the 3D view draws | `ligand_comps` | `classes` + block overrides naming those units | adapt |
| §3 track providers (`tracks.ts`) | `build` → old `DocNode`s, `computed` | new nodes (`values` `field`, `ranges` from the source); public break | adapt |
| §3 secondary structure track | `molstar:secondary_structure`, `assignment` param, merged copies one row set | host `ranges` by region, `changes` on assignment; merged copies: G1 | gap G1 |
| §3 pLDDT (model-archive) | `computed` `molstar:plddt` | `values` `{ field: 'plddt' }`, an extension value channel | adapt |
| §3 edits by base patch, origin `host:prism-viewer.sequence.edits` | old `applyPatch` | `SpecModel.applyPatch` (`@molstar/core/repr/patch` ops, `requireBase`), a refusal rebuilds | adapt |
| §3 behaviors (`behaviors.ts`) | `patch(shown, Op[])` | same contract, new `Op`; public break | adapt |
| §3 ids `<owner>.sequence.<name>` | refs | node ids (every child has one); same strings where the id charset takes them | adapt |
| §4 link per document (`link.ts`) | `StructureLink.attach` + `linkViews` | `linkStructure(viewer, side)`, made again per document | adapt |
| §4 3D side as operations, one channel per kind | `view3d` | the side's three calls → Prism ops; selection coalesced by the link | adapt |
| §4 state lives in the engines | — | same | direct |
| §4 ligand whole | `StructureLink.toLoci` | `ligand: 'whole'` | direct |
| §4 emptying from the sequence resets the camera | origin `pointer`/`keyboard` | `side.select(null)` carries no cause | gap G3 |
| §4 marks as `marks:<owner>:<id>` channels | highlight channels | `HighlightModel` channels | direct |
| §4 3D hover owned by `3d` | external hover | `HoverModel.setLocation`, owned by the link | direct |
| §4 size gates 2,000 / 20,000 | `linkViews` defaults | `DefaultSyncLimits` (same bounds, `over: units`) | direct |
| §1/§4 `revealHighlight` (`links:*`, `3d`) | optional behavior | only `revealOnSelect`; nothing reveals a side hover or a channel | gap G4 |
| §5 actions as API calls, colouring by theme | `HostActions` | unchanged | direct |
| §6 MVS needs no sequence part | — | same | direct |
| §6 `load.svs` (`.svs`, `.svsj`, `.svsx`) | old codec, `readSvsx` | new JSON/MessagePack codec, `readArchive` (`seqx` archive); old files refused | adapt (break) |
| §6 `matchStructure` | download address, `pdbe`/`afdb` entry, host file | `uri` and `ref` only: the new language has no services | adapt (loss: entry ids, Q5) |
| §6 override (`svsDocument`, `setOverride`) | a copy re-pointed at the host's parse | `SpecModel.apply` (whole root); its `source` served by the adapter by `uri`; Prism's edits patched over | adapt |
| §6 `.svsx` files in memory | `addArchive`, `archiveLoaders` | `archiveReader` behind the reader port | direct |
| §7.1 saved state | version 4 | version 5 (below) | adapt |
| §7.1 edits overlay | `{ [structure]: { tracks, chains } }` | same shape; track ids, unit ids (`label_asym_id`, no operator) unchanged | direct |
| §7.2 checkpoints, Undo restores the sequence's state | custom state | unchanged: Prism's saved state is the truth, the spec a projection (§History) | adapt |
| §7.3 one history manager | Mol\*'s off, prism-ui's on | unchanged; Seq\* gets no history port | adapt |
| §7.4 large scenes | — | unchanged | direct |
| §8 band, fade, tab (`shell-ui.tsx`) | Prism | Prism | direct |
| §8 control (`sequence-control.tsx`) | `mountSection` on Prism's surface | `mountCanvas(viewer, el)` of a one-section view, `canvas.height: 'content'`, `section.maxPx` | adapt |
| §8 menu and commands (`commands.ts`, `api/sequence.ts`) | Prism | Prism, ops' bodies over the new model | adapt |
| §8 Text (`compact` preset + render props) | a preset, view state | section `layout: 'text'`, `align: center`, `max_width`; a spec patch | adapt |
| §8 Tracks (`presets.ts`, `headerWidth: 0`) | host preset | tracks layout, `labelWidth: 0` (no measured column; check a 0 column lays out) | adapt |
| §8 rail (`ruler`, `strip`, `band`, `rail`, `none`) | old navigators | `navigator` nodes (`rail` form, pictures `axis`/`tracks`) and `navigator.positionsRail`; `strip`/`band` looks differ | adapt (parity) |
| §8 height | `layout.contentHeight` | content height, written when it changes | direct |
| §8 track menus | `interaction.on('menu')` → `openContextMenu` | `headers` behavior's `menu` port (Prism's `openContextMenu`), `contextClick: everywhere`; Prism's `header-action` entries for colourings, Hide | adapt |
| §8 hide/show/add/remove tracks | `track.hidden` | `track.hidden` (a one-key `set`) | direct |
| §8 chains | section patch | section patch | direct |
| §8 hidden tracks hint | Prism over the document | `HiddenTracksModel` (`models/hidden-tracks.ts`), drawn by Prism | adapt |
| §8 track gutter (`track-gutter.tsx`) | `SectionLayout`, `EntryKind` | the last frame's layout (`viewer.frame.layout`) rows, or a header slot entry (Q4) | adapt |
| §8 keyboard cursor said | Prism's live region | `mountCanvas`'s announcer; Prism's dropped | direct |
| §8 framing below the band (`sources/framing.ts`) | Prism | Prism | direct |
| §8.1 Sequence node (`sources/sequence-scene.ts`) | over Prism's model | unchanged | direct |
| §9 skin (`skin.ts`) | `Skin` from the ink | `viewer.setTokens(createTokens(…))` from the ink, `theme.source: 'config'` (no `MutationObserver`, no CSS read) | adapt |
| §10 `prism:residues` (`palette.ts`) | registered scheme | a `scheme` entry through `defineExtension` (`prism` namespace) | adapt |
| §10 Zappo, Taylor, Hydrophobicity, none | `KernelColorSchemes` | `zappo`, `taylor`, `hydrophobicity` built-ins | direct |
| §10 legend (`legend.ts`, `legend-card.tsx`) | `schemeLegend` | core `theme` legends; Prism's card | adapt |
| §10 hover label (`sources/hover-lines.ts`, `link.ts`'s `annotateHover`) | `seq.labels.current`, track lines | `hoverCard` behavior with Prism's port: the composed card's rows go to Prism's `LociLabelProvider` when the hover maps, to `HoverLabel` when not; no Seq\* card drawn | adapt |
| Prism rule: one frame tick | `rafScheduler` over prism-ui's tick | Seq\* asks `@molstar/ui/dom`'s tick: two per frame | gap G7 |
| Prism rule: no `MutationObserver` | — | only under `theme.source: 'root'`; Prism sets `config` | direct |

**Count**: 63 rows: 22 direct, 36 adapt, 5 gap (G1, G3, G4, G7 twice); §2's file table adds G2 (SIFTS) and G6
(query) and the stated losses (coverage, PAE in the band, services).

### 4. The gaps, as Seq\* lanes before the swap

- **G1 Host ranges on units, merged copies.** No host source serves `ranges` today (the fake serves values only; the
  mmCIF runtime reads no `struct_conf`). Build: a conformance case for a host source serving ranges on unit systems,
  read by the `secondary` layer, reopened on `changes`; and say what a merged copies block draws of per-unit ranges
  (today: one merged row set). Seq\* lane, with the fake host gaining secondary structure.
- **G2 SIFTS in the description.** `StructureDescription` has no per-entity mapping onto a reference sequence; today's
  `sifts.ts` gives `uniprot:<acc>` numberings and coverage, which a loaded spec mapping P04637 onto 2OCJ needs. Build: an
  optional `segments` table (entity, accession, runs) → a `structure` map per entity (`map.method: structure`, keys
  from identifiers). Seq\* lane (core `io/structure`); else a stated loss for loaded specs.
- **G3 The cause of a selection sent out.** `side.select(target)` needs the origin (`pointer`, `keyboard`, a command,
  an echo-free carry) so Prism resets the camera only on the person's emptying from the sequence (§4). Seq\* lane:
  `select(target, { origin })`, or the selection model's last change origin public.
- **G4 Reveal what the structure lights.** `linkStructure` option `revealOnHover` (and a channel's reveal for marks), at
  most once a frame, never echoed. Seq\* lane (the link owns the ownership rules; a host doing it would duplicate them).
- **G5 `hostAdapter(fallback, serve)`**: the port's helper, built only inside the Storybook's fake. Promote to
  `@molstar/seqx-core/data` with the conformance suite every adapter passes (the port's Consequences): the Mol\*
  adapter is its first real user. Small.
- **G6 Query extensions in the viewer** (core `_spec` "Context": "the viewer's address model does not yet: a
  followup"): `molql` over the linked structure. Not used by Prism's band today; a stated loss until the followup.
- **G7 One frame tick, one DOM runtime.** The new viewer imports `@molstar/ui` (`dom` frames, `DomHost`, menus,
  `--msui-*` tokens); Prism runs on prism-ui's own runtime and tick (`prism/ui/src/reactive/frame-state.ts`). Two
  `requestAnimationFrame`s a frame breaks Prism's rule (`frames.quick.ts`); two runtimes cost size. Not a Seq\* lane:
  a decision (Q2).
- **Region-level invalidation** (the port's "later refinement"): not needed; a structure's reopen is whole and cheap at
  structure sizes. Stays later.

## The saved state: version 5

`prism-viewer.sequence` today: `{ version: 4, style, colors, rail, structure?, override?, edits? }`.

- **Version 5**: `{ version: 5, style: 'tracks' | 'text', colors: string | null, rail: Rail, structure?, override?: {
  spec, url?, structure, archive? }, edits? }`. `spec` is the new language's JSON (`LanguageVersion` recorded in it).
  The Seq\* spec is not saved: it is built from this (the providers, the edits, the override), as today.
- **Migration** reads every version: 1 → 2 → 3 → 4 as today (`runtime.ts`), then 4 → 5: `style` `compact` → `text`;
  `colors` `prism:residues`, zappo, taylor, hydrophobicity kept by name (the new built-ins' names); `rail` names mapped
  (`ruler` → positions rail, `strip`/`band` → nearest navigator, `rail`, `none`); `edits` unchanged (track ids `ruler`,
  `letters`, `ss-track` and providers' ids kept; chain names are `label_asym_id`, the new units' ids without an
  operator); **`override` dropped** with one notice ("the saved sequence document is of an older format; showing the
  structure's own"): an old-language document has `computed` nodes and providers the new language lacks, and no
  converter exists. A stated break (Q3).
- Unknown or absent versions reset to the defaults, as today.

## Public API change (for `.v6/plans/breaking-v6-changes.md`)

- `ViewerExtension.sequence.tracks: SequenceTrackProvider[]`: `build(ctx)` returns new-language nodes;
  `SequenceTrackContext` gains `source` (the structure source's id) beside `complex`; no `computed` node. New
  `ViewerExtension.sequence.values` (value channels, §1). `menu`, `legend`, `describe`, `applies`, `whenChanged`,
  `shown`, `base` unchanged.
- `ViewerExtension.sequence.behaviors: SequenceBehavior[]`: `patch(shown, ops)` takes `@molstar/core/repr/patch` ops
  (the old `Op` from `@molstar/seq-core/format` goes); `SequenceShown` keeps its fields (`complex`, `data`, `section`,
  `after` are node ids).
- `api/sequence.ts`: operations keep their ids and params (`sequence.style`, `.colors`, `.rail`, `.follow`,
  `.track.*`, `.chains`, `.get`) save `rail`'s values (mapped) and `colors`' names; `showDocument` takes a new-language
  spec; `GraphResult` (old viewer) → the patch result; `runtime()` exposes the new readers (`readArchive`, the codec).
- `@molstar/prism-viewer/sequence/molstar`: the old adapter's exports go (`StructureLink`, providers, `readStructure`);
  the new adapter's come. `sequence/molstar/pae-matrix` stays (model-archive's card).
- `load.svs`: old `.svs`/`.svsj`/`.svsx` files are refused with a message; `svs-url` keeps its name (Q3).

## MolViewSpec beside SeqViewSpec

Unchanged rule (§6): the viewer loads one or the other; an MVS builds the plugin state and the sequence follows. The
document Prism shows is the new SeqViewSpec at the core's `LanguageVersion` (0.10.0); the override holds it whole. No
MolViewSpec types enter Seq\* (structure-syncing, decided); Prism's side translates.

## History

One undo, Prism's (`managers.history`, prism-ui). Seq\* is given **no history port**: its `SpecModel` applies edits
at once. Every change of the band's content goes through Prism's API (style, colours, rail, tracks, chains), which
changes Prism's saved state and patches the spec; Undo restores Prism's state and the model rebuilds or patches the
spec. So Seq\*'s own editing paths are left out (`editing` behavior unlisted) and its header menu's items run Prism
operations (the menu port is Prism's). Why: Prism's checkpoint already holds the sequence's state; a second stack
(`hostHistory` is for `@molstar/ui/history`, not prism-ui's) would be two truths. A patch applied outside that path is a
bug (a test: no Seq\* commit without a Prism run).

## Performance and size

- **Must hold**: Prism's `frames.quick.ts` (idle asks no frame; one tick: G7), `waste.test.tsx`, `subscriptions*`,
  `leaks.heap.tsx` (mount and dispose the band: the new viewer's `retire` releases its canvases, observers, frame), the
  scaling tier (a large assembly's selection through `toUnits`' gate), the axe sweep over the band.
- **The adapter**: `describeStructure` once per document key, linear in residues; a frame costs nothing; a hover a few
  µs on 2OCJ's four units (core "Structures", Performance); the link's quick/heap guards (`link.quick.ts`,
  `link.heap.ts`) port to the side translation.
- **Size**: the new viewer's measured bounds (`packages/seq/_measure/sizes.mjs`): default 1,545 kB, minimal 1,090 kB
  minified, `@molstar/core` included. Prism pays its runtime chunk only (a spec without a sequence still bundles none),
  minus what Prism already ships (`@molstar/core`), plus `@molstar/ui`'s part (G7). Measure the chunk before and after
  in the swap lane; Prism's spec trims what the band never draws (navigators beyond the rail, rows, bridges, genome
  chunks are lazy already).

## Tests

- **Port as they are** (Prism behaviour, no Seq\* type): `_test/commands`, `edits`, `track-hint`, `labels`;
  `src/_test/sequence-api`, `sequence-node`, `sequence-tracks`, `sequence-frames` (assertions kept, fixtures rebuilt);
  `sources/_test/sequence`, `data-by-url`, `framing`, `checkpoint`, `restore`; `molstar/_test/chains`,
  `insertion-codes`, `whole-chains`, `whole-residue`, `ligand-link` (the side translation: same cases, new API).
- **Change**: `model.test.ts`, `link.test.ts`, `sequence-control.test.tsx`, `skin.test.ts`, `track-gutter.test.ts`
  (new APIs); `molstar/_test/structure`, `registry`, `secondary`, `link` → the adapter's conformance (G5's suite):
  equal nodes from the plugin's model and from the mmCIF runtime on the same file (1tup, 2ocj, 1c26, 2xwr, AF-P04637,
  `ins1.cif`: the fixtures already here), keys unchanged across a reopen, a frame fires no `changes`, an assignment
  change does; `features/loaders/_test/svs.test.ts` (new files; old ones refused); the version-5 migration from each
  of 1–4.
- **Drop**: `molstar/_test/pae`, `query`, coverage's (stated losses) unless G2/G6 lanes land.
- **Parity**: each story in `src/sequence/_stories/sequences.stories.tsx` beside today's (light, dark, wide, narrow;
  Text and Tracks; merged and apart copies; a ligand; an AlphaFold model with pLDDT), compared by eye and by a
  headless twin asserting the same track list, chain blocks and hover label lines.

## Lane order

1. **Seq\* gaps** (Seq\* packages, their checks): G5 `hostAdapter` and the adapter conformance suite; G1 host ranges
   (the fake gains secondary structure); G3 selection origin; G4 `revealOnHover`; G2 SIFTS segments if Q5 keeps them.
2. **The Mol\* adapter** (Prism, `sequence/molstar/`): `describeStructure`, the format runtime, `changes`, value
   channels, the side translation; tested headless against the mmCIF runtime on the existing fixtures. Old Seq\* still
   draws the band.
3. **The swap** (Prism): model, document, providers, behaviors, control, skin, palette, legend, hover label, menus,
   saved state v5, API; `package.json` drops `@molstar/seq-*`; `_spec/sequence.md`, `README.md`, extension docs,
   `breaking-v6-changes.md`, model-archive's provider, in the same change; G7's choice applied.
4. **Parity**: the stories side by side, the axe sweep, perf, heap and scaling tiers, the size measured.

## Edge cases

- A structure removed while its source opens: the open aborts (`AbortError`), no issue.
- Two structures of one file (two models): two keys, two sources; a `uri` naming the file maps to the first, as today.
- An assembly's copies of a chain lit from 3D: collapsed onto the chain's one unit; lit back on every copy.
- A loaded spec naming a file the plugin does not hold: refused with today's message, the rest of the spec draws.
- A unit id with `-` in a `label_asym_id`: refused by `structureNodes` (`seq.io.structure.id`): a structure whose
  chains carry it loses those chains in the band. Check the PDB's archive for such ids before the swap.
- None chains: a block of a unit there is not, as today; if the new checks refuse it, the band empties by Prism's
  `empty` state instead.

## Rejected

- **Parse the files again in Seq\*'s mmCIF runtime**: two parses, keys from another reading, and no structure without
  a file (predicted, edited). The author's requirement and the port's rule.
- **Keep `computed` nodes as an extension kind**: the language takes no node kinds from extensions (core "Context",
  Extensions); sources are the one way data arrives.
- **A Seq\* history port over prism-ui's stack**: two truths for one band (History above).
- **The adapter in `seqx-core`**: it would import `@molstar/model`, outside the allowlist.
- **A converter from old SeqViewSpec documents**: no user base beyond Prism's own saved sessions; cost of a second
  language's reader (Q3 asks).

## Consequences

- Prism: `src/sequence/**` rewritten over the new packages; `src/sequence/molstar/` becomes the data-port adapter;
  `_spec/sequence.md` §1–2, §3, §4, §6, §7.1, §9, §10 rewritten; `README.md`'s extension rows; `TODO.md`'s "more tracks
  from the structure" becomes value channels.
- Seq\*: G1–G5 in `core` and `viewer`, their `_spec`s; the port design's "no Mol\* adapter" closed by Prism's.
- model-archive: its pLDDT provider rewritten (a value channel + a `values` track).
- **`packages/sequences/storybook` imports `@molstar/prism-viewer/sequence/molstar`** (`_test/molstar/*.test.ts`,
  `file-formats.stories.tsx`): the swap breaks the old Storybook's Mol\* demo before the seventh slice removes it.
  Remove those stories with the swap, or keep the old adapter's exports until then (Q6).

## Questions for the author

1. **Where the Mol\* adapter lives**: in Prism (`sequence/molstar/`, no plugin import, a boundary test), or a bridge
   package `@molstar/seqx-molstar` now. *Recommend*: Prism now, written to move unchanged; a package when the new mol
   viewer or the Seq\* Storybook needs it.
2. **One tick, one runtime** (G7): Prism patches `DomHost.requestAnimationFrame` onto prism-ui's tick (one line in the
   composition root, global), prism-ui's tick delegates to `@molstar/ui/dom`'s, or Seq\*'s viewer takes a frame source
   from its spec. *Recommend*: prism-ui's frame state on `@molstar/ui/dom`'s tick (one tick for the page, no global
   patch); the size of `@molstar/ui` in Prism measured first.
3. **Old SeqViewSpec files and saved overrides**: refused and dropped (a stated break), or a one-way converter for
   what Prism itself wrote. *Recommend*: break; v6 is unreleased.
4. **The track gutter**: keep Prism's gutter reading the layout's rows, or a Seq\* header slot entry (Prism's
   `header-action`) with a narrow label column. *Recommend*: Prism's gutter (the look is Prism's), over
   `viewer.frame.layout`.
5. **SIFTS, coverage, PAE and services in the band's documents**: G2 now (a Seq\* lane), or stated losses until a
   loaded spec needs them. *Recommend*: G2 now (a loaded UniProt-to-structure spec is the band's main SeqViewSpec
   case); coverage, PAE in the band and `pdbe`/`afdb` entry matching as losses.
6. **The old Storybook's Mol\* demo**: removed with the swap, or the old adapter kept exported until the seventh slice.
   *Recommend*: remove those stories with the swap; the seventh slice removes the rest.

## Where "fairly straightforward" holds

It holds for the Prism-owned three quarters (follow, the node, menus, commands, saved edits, the band's shell, the
API's surface: 22 rows direct), and for the link (`linkStructure` was built from Prism's own link). It does not for:
the data path (a new adapter, the port's first real host, 2–3 lanes), host ranges and SIFTS (Seq\* gaps), the history
(Seq\*'s port does not fit prism-ui's stack: decided around, not through), and `@molstar/ui` arriving in a Prism that
has its own runtime (a decision, not code).
