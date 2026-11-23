# Structure syncing: a library for a sequence view beside a structure viewer

## Intent

The author, 2026-11-19: "let's only do a demo with fake 3D viewer to demonstrate how the syncing will work. There will
be a new mol viewer work after we are done with the sequence viewer which will replace prism, so i don't want to waste
work". On the words that cross: "this should be provided by the application layer anyway". And, 2026-11-20: "this is
just a wrapper around events the viewer should expose normally. we should just provide a library that helps with the
coordinate translation and the wiring, and not call it a port."

So syncing is not part of a viewer's composition. A Seq\* viewer publishes its hover, selection, highlights and jumps
and takes commands, as it does for any host. This design is a **library** over that surface:
- **translation**: pure functions between Seq\*'s Locations on units and a structure's residues, through the maps;
- **wiring**: one helper, `linkStructure(viewer, side)`, that subscribes both ways and returns a disposer.

A small fake structure viewer in the Storybook shows every path: hover and selection both ways, focus, "color by this
track", a structure's chains mapped onto sequences. The Mol\* side (a `Model` adapter, a live 3D view) waits for the new
mol viewer.

It builds on [molecules.md](molecules.md) (units with `chain`, `operator`, `observed`, `auth`), [maps.md](maps.md) (the
copy rule, the hull, the `structure` kind) and [track-header-parts.md](track-header-parts.md) (a host's header entry).
Reading units from a host's structure is the seq data ports design's.

## Precedents

- **The current Seq\*, in Prism and in its Mol\* story.** A link between the viewer and a 3D view (`StructureView3D`:
  `highlight`, `select`, `focus`, `onHover`, `onSelect`). Taken:
  - hover and selection both ways;
  - the size gate (EC-167): residue by residue within bounds, else whole chains or nothing;
  - a 3D hover with an owner of its own, so a late 3D leave never clears the pointer's hover, and a 3D hover older than
    the pointer's last own hover (its `at`) is dropped (EC-141);
  - a non-polymer position lighting its whole residue;
  - chain instances keyed by `label_asym_id` and operator (EC-181).

  Not taken:
  - the link speaking Mol\* `Loci`: it imports Mol\* and the plugin;
  - "show in 3D" composing MolViewSpec inside Seq\*: a state is the host's;
  - the host registering external coordinate systems with their maps: units are already the spec's systems;
  - a selection from 3D focusing the sequence view: it reveals here (Decision).
- **Mol\*'s sequence panel.** Hover and click a residue, and the 3D view lights or selects it; a chain picker. Taken: the
  gestures. Not taken: the picker (a section of a complex shows its blocks).
- **Jalview with a structure.** A sequence mapped to chains by SIFTS or by alignment; hover and selection both ways;
  "Colour by Jalview", kept live. Taken: the paths, the colors, SIFTS and an alignment as two maps of one pair (the
  choice below). Not taken: live colors by default (a snapshot, live planned), a chain mapping chosen in a dialog (a map
  is data; which one is taken is a map's `priority` or a section's `via`).
- **ChimeraX's sequence viewer.** Chains associated by alignment, selection both ways, the structure colored by a
  sequence attribute. Taken: association as a map. Not taken: attribute assignment (the seq data ports').
- **RCSB's 1D feature viewer and PDBe's ProtVista with Mol\*.** Hover a feature, its residues light in 3D; click selects
  and focuses. Taken: a feature's extent crossing as residues. Not taken: click focusing (focus is the double-click).
- **MolViewSpec** (`packages/mvs/builder`). Component expressions (`label_asym_id`, `instance_id`, `label_seq_id`,
  `auth_seq_id` with `pdbx_PDB_ins_code`, `label_atom_id`), `color` and `focus` nodes; state only, no hover or
  selection. Taken: by the fake host's adapter, as its viewer's vocabulary, and as the check that the unit facts Seq\*
  gives are enough to write it. Not taken: any MolViewSpec type in Seq\* (decided).
- **This codebase.** Reused as built:
  - `HoverModel`: owners, and projection over linked sections;
  - `SelectionModel` methods;
  - `AddressModel`: `focus`, `reveal`, `hoverAddress` on a channel of its own (EC-182);
  - core's `StaticModel.project`: onto systems no section shows;
  - `findPath`'s rule set (EC-22 to EC-25), extended by the map choice below;
  - a section's `direction: { along: <map> }` and `pieces: <codon map>`: a section already names a map;
  - the animation port: one clock for the sequence view and the camera ([stories.md](stories.md));
  - `hostHistory`: a host's stack over Seq\*'s spec and a structure's state.


## Decision

### Who is who

- **Three parts: Seq\*, a structure viewer, and the host.** Seq\* imports nothing of a structure viewer and declares
  none of its words (no MolViewSpec type, decided). The host's adapter is the only code that knows both. Why: host
  neutrality; the vocabulary is the application's.
- **No port, no behavior, nothing in the composition** (decided, the author, 2026-11-20). A host that never calls the
  library carries none of it; [ports.md](ports.md) and [viewer-composition.md](viewer-composition.md) name nothing for
  structures. Why: the viewer's events and commands are the surface every host already uses; a port would duplicate
  them for one consumer.
- **Two entries** (proposed): `@molstar/seqx-core/structure` (translation: pure, Node, no viewer) and
  `@molstar/seqx-viewer/link-structure` (the wiring: one file over the viewer's public models). Why apart: a headless
  host (the render CLI, a server colouring a structure) needs the translation alone; the wiring is thin and is the one
  part that knows a live viewer.

### Translation (`@molstar/seqx-core/structure`)

Pure, tested in Node, over a `StaticModel` (and its overlay).
- **`structureFacts(model, complex)`**: per unit its address, `chain`, `auth_chain`, `operator`, molecule and its kind
  (polymer, branched, other), length, `observed` ranges, and its `auth` numbering as a numbering entry
  (`label(position)`, `position(label)`, EC-1 to EC-8). Plain data. Why: an adapter needs no model walk, and every host
  translates from the same facts.
- **`structureFacts`' residue rows**: per unit, its residues as rows (position, `auth` label and insertion code,
  observed): the table an adapter joins against its own residue list once, so per message it looks up, never searches.
- **`unitRows(model, store, annotation, complexes)`**: an annotation's rows as targets on units with their fields, for a
  host to write as its own representation ([seq-data-ports.md](seq-data-ports.md)), beside `colorRuns`.
- **`toUnits(model, loc, complexes, options) → { target, cut }`**: a Location carried onto the units of the bound
  complexes along static or overlay paths, under the copy rule and the use's `via` (a map onto a unit a source's
  description gives, SIFTS per chain, resolves as an overlay edge once the description has arrived, pending until
  then; the author, 2026-11-20, ST1's drift: "Yes, as recommended"); unobserved positions cut (no atoms,
  nothing to light); a molecule unit whole (`ligand: whole`); the size gate (EC-167): past `maxRuns` or `maxLength`,
  whole units (`over: units`) or nothing (`skip`). Why the gate here: a 4,000-unit assembly's selection must not become
  4,000 × n runs a frame, whatever the host.
- **`fromUnits(model, target)`**: the reverse: parts on units, on a molecule (every unit), or on a unit's atoms
  (collapsed onto residues unless a section shows the `atoms` form); the rest counted in `seq.structure.unplaced`.
- **`colorRuns(model, store, layer, complexes) → { complete, pending, runs }`**: the layer's compiled color channel (the
  same compile as drawing it, a dark variant included) as `{ color, target }` runs on units, unobserved cut, uncolored
  positions in no run, data not held `pending`, never drawn as empty.
- Round trip: `fromUnits(toUnits(loc))` is the Location's observed part on the units (a property test); each adapter's
  own round trip (the fake's Location ↔ MolViewSpec) is its own test.

### Wiring (`@molstar/seqx-viewer/link-structure`)

```ts
interface StructureSide {
  /** Seq*'s hover, selection and jumps, as targets on units (`null`: cleared). */
  highlight(target: StructureTarget | null): void;
  select(target: StructureTarget | null): void;
  focus(target: StructureTarget, options: { readonly durationMs: number }): void;
  /** The structure viewer's own events, translated by the host; each returns its unsubscribe. */
  onHover?(fn: (target: StructureTarget | null, at?: number) => void): () => void;
  onSelect?(fn: (target: StructureTarget | null) => void): () => void;
  onFocus?(fn: (target: StructureTarget) => void): () => void;
  /** Aborted when the structure viewer goes: the link clears and disposes itself. */
  readonly signal?: AbortSignal;
}
function linkStructure(viewer: SeqViewer, side: StructureSide, options?: LinkStructureOptions): () => void;
```

`StructureTarget` is a Location on units (structured-cloneable, so it crosses `postMessage`). The helper:
- **out**: on `hover.version`, the hover's Location through `toUnits`, sent when its image changes (an unchanged image
  sends nothing); on `selection.value`, the `shown` Location, coalesced to one message a frame (a drag moves it every
  pointer event); on a jump naming a Location, `focus` with the jump's duration (camera and sections on one clock);
- **in**: the side's hover as a hover owned by the link (drawn as a projected hover, no card; the pointer taking the
  hover takes ownership, so a late `null` clears only the link's own, EC-141); its selection through `fromUnits` into
  `selection.set`, then revealed (the zoom kept); its focus as `addresses.focusLocation` (one jump);
- **echoes**: what came in is never sent back (hover by owner; selection by value, `locationsEqual` against the last
  set);
- **clearing**: the disposer, the side's `signal`, or the viewer retiring sends `highlight(null)` and `select(null)`
  once and drops the link's hover;
- **options**: `complexes` (`{ complex, via? }[]`; absent, every complex), each direction (`hover`, `select`, `focus`:
  `both`, `out`, `in`, `none`), `revealOnSelect`, `followJumps`, `limits`, `ligand`. Options, not `ViewerConfig`: they
  belong to one link, and a viewer may have none or two.

Why a helper and not a behavior: what it adds over the viewer's events is echo suppression, throttling and clearing,
all per link; a host with other needs writes its own wiring over the same translation.

### "Color by this track"

A host's header entry ([track-header-parts.md](track-header-parts.md)), never built into Seq\*: its command calls
`colorRuns` and the host builds and loads its viewer's state (the fake's: a MolViewSpec state, one `color` node per
run). A snapshot: again restores the base state; a second track replaces the first; `live` is the host's to rebuild on
the layer's changes. Undo is the host's (`hostHistory`).

### The fake host's adapter (the reference translation)

In the Storybook beside the fake viewer, on `structureFacts` and its residue rows, as a real host's will be:
- a run on a polymer unit is `{ label_asym_id: chain, instance_id: operator, beg_label_seq_id, end_label_seq_id }`;
  `label_seq_id` is the position plus one (a polymer unit's molecule is the entity's whole sequence);
- a branched unit's residue is `auth_seq_id` (with `pdbx_PDB_ins_code`) through `auth`;
- a molecule unit is `{ label_asym_id, instance_id }`; an atom adds `label_atom_id`;
- inward the reverse; `label_entity_id` and class selectors mean every unit of that molecule or class; MolQL is not
  placed.

### Choosing between two maps of one pair of systems

Today two maps joining the same two systems (SIFTS's `structure` map and an `alignment` map from P04637 to one chain's
molecule; a per-unit map beside a map onto the molecule) make every path through them ambiguous (EC-23, EC-24): the
search finds equal-length candidates that compose differently and takes none. Several maps of one entry that join
different pairs (one per chain) already resolve.

Decided, the author, 2026-11-19: "priority + via", both in the 0.5.0 bump; if usage runs short, `priority` in MK and
`via` later.

- **A map's `priority`** (an integer param of every map node, default 0). In `findPath`, after equal compositions are
  merged, distinct shortest candidates are ranked by their maps' priorities sorted weakest first and compared in order (the
  lowest priority among the maps they step through, then the next; unit steps have none) (decided by the author, 2026-11-20: "Compare sorted priorities": a shared unranked map ties and the rival maps decide); the top rank wins when it holds one composition; a tie stays ambiguous. Why:
  - one place decides for every use (follows, hover projection, the structure wiring, queries' `project`), so two
    sections never disagree about where a residue is;
  - the search stays shortest-first (EC-22): a priority never makes a longer path beat a shorter one;
  - the path cache stays keyed by its pair: the param is in the map entry's hash, so a carried search is dropped when a
    priority changes (paths.ts's `deps`);
  - weakest first, not a sum: a path is as trustworthy as its least trusted step, and among equally weak paths the
    next step decides, so a map both share never hides the rivals.
- **A use's `via`**: a section's `via` (a ref to a map, or a list), and a structure binding's (`linkStructure`'s `via`). For that use the search
  leaves out the rivals of the named maps (other maps joining the same two systems); paths that never reach that pair
  are unchanged. It applies to the section's `follow`, its `direction: { along }`, and what is projected onto it (hover,
  selection, the focus hull). Why: a figure comparing SIFTS against an alignment shows both, one section each; the
  precedent is a section already naming a map (`along`, `pieces`). `StaticModel.path(from, to, via?)` keeps one entry
  per pair and `via` set.
- **Not chosen silently**: no rule ranks kinds (`structure` over `alignment`) or specific over general (a per-unit map
  over a map onto the molecule). Why: a wrong map taken quietly is worse than the ambiguity issue, which teaches;
  `seq.resolve.path.ambiguous`'s hint names `priority` and `via`.
- **Links** carry no `via`: a link names addresses and actions; where they land is the view's, so the sections' `via`
  applies.
- Language: `priority` and `via` join the 0.5.0 bump (one bump, unreleased). EC-22 to EC-27 run on every encoding with
  a ranked pair, a tied pair and a `via`.

### The fake structure viewer (Storybook)

- **A Storybook component** (`packages/seq/storybook`, `_stories/support/`) that imports nothing of Seq\*, on
  `@molstar/ui`; its states are MolViewSpec, validated with `@molstar/mvs-builder` (the Storybook may import it,
  [packages.md](packages.md)). Its adapter, beside it, is the host's: the MolViewSpec translation and its `side`.
- **Its structure is its own recorded chain table** per entry (`label_asym_id`, `auth_asym_id`, entity, operators, and
  per residue `label_seq_id`, `auth_seq_id`, `pdbx_PDB_ins_code`, `label_comp_id`, observed), recorded from the mmCIF by
  the examples' record script, independent of Seq\*'s nodes. Why: the mapping crosses two independently held
  representations, as it will with a real viewer.
- **It draws** each chain instance as a row of residue cells (unobserved as a dashed gap, ligands as chips, copies as
  rows), the state's `color` nodes, a hover outline, the selection, and a "camera" frame that `focus` moves over the
  given duration.
- **Its own grammar** (the host's): hover; click selects; Cmd/Ctrl+click toggles; Esc clears; double-click focuses; a
  "select chain" and a "within 5 Å" button for structure-only selections.
- **`latency`** (0 by default) delays its messages: the stories and twins exercise a stale hover.

### Datasets

- **p53 with 2OCJ**: four copies of the DNA-binding domain, author numbering equal to P04637's, unobserved ends, a zinc
  per chain. Recorded in the fourth slice ([examples.md](examples.md)). P04637 maps to the entity's molecule by a
  `structure` map (SIFTS's segments with provenance and status); a story adds an `alignment` map of the same pair, ranked
  lower, and a section with `via` on it.
- **Complex with 1FIN**: two copies of the complex, cyclin A as a construct (P20248 from 172), ATP as a molecule unit;
  its maps `structure` maps too.
- **The `structure` kind is kept** (decided): segments, each with provenance (SIFTS, an alignment, the author's
  numbering), and per-residue status (paired, conflict, unpaired) as [maps.md](maps.md) says; onto a unit or onto the
  entity's molecule. Unobserved is the unit's `observed` (the structure's own fact); a map's unobserved status that
  disagrees is reported once and the unit wins. Why kept: SIFTS's provenance and status are the map's, and a spec
  states them where a reader of the format looks.


## The interaction matrix's cells

With the columns of [ux-foundations.md](ux-foundations.md) section 8, under a `linkStructure` with its defaults; "→
side" is what the side is told. Nothing here exists without the link.

| | select | extend | add / remove | focus | back | zoom in / out | pan | hide | edit |
|---|---|---|---|---|---|---|---|---|---|
| **an element** on a unit or its reference | as the grammar; → side `select` | as the grammar; → side | as the grammar; → side | as the grammar; → side `focus` | as the grammar; → side `focus` when the step names a Location | no: a sequence zoom has no structural equivalent | no | no | no: data |
| **a residue in the structure viewer** | its click: Seq\*'s selection replaced, revealed | its own; arrives whole | its own; arrives whole | its double-click: one jump | its own | the camera's | the camera's | no | no |
| **a range** | drag; → side once a frame | Shift+drag; → side | Cmd/Ctrl+drag, Alt+drag; → side | Enter; → side | Alt+Left; → side | — | — | — | — |
| **a section of a complex** | as [molecules.md](molecules.md); a merged block → every member | as molecules.md | as molecules.md | → side | Alt+Left | as its layout | as its layout | no | `via` by patch |
| **a map** of one pair | — | — | — | — | — | — | — | — | `priority` by patch: every use re-resolves |
| **a track or layer** | — | — | — | — | — | — | — | the snapshot stays | "color by this track": a host's entry |
| **the view** | — | — | — | 0 fits: nothing sent | Alt+Left; → side when the step names a Location | no | no | — | — |

## Configuration

| Item | Kind | Default | Why | Alternatives |
|---|---|---|---|---|
| `complexes` | link option | every complex, no `via` | the spec names no host ref | `{ complex, via? }[]` |
| `hover`, `select`, `focus` | link options | `both` | the precedents link all three both ways | `out`, `in`, `none` |
| `revealOnSelect` | link option | `reveal` | the zoom kept | `focus`, `none` |
| `followJumps` | link option | `true` | back and forward bring the camera | `false` |
| `limits.hover`, `limits.select` | link options | 2,000 runs / 20,000 positions; 20,000 / 200,000; `over: units` | the precedent's bounds | any; `skip` |
| `ligand` | link and `toUnits` option | `whole` | a ligand's position is its residue | `atoms` |
| a map's `priority` | map param (0.5.0) | 0 | equal maps stay ambiguous until ranked | any integer |
| a section's `via` | section param (0.5.0) | none | a per-use choice | a map or a list |
| unobserved positions | fixed | never sent | nothing to light | — |
| the vocabulary | fixed | Locations on units | the application's (decided) | MolViewSpec types, Mol\* `Loci`: the host's |

## In the specs

| Decision | Where |
|---|---|
| Three parts; Seq\* imports nothing of a structure viewer and names none of its words | core › Structures |
| No port, no behavior, nothing in the composition | viewer › Structure link |
| Two entries: `@molstar/seqx-core/structure` (translation) and `@molstar/seqx-viewer/link-structure` (wiring) | core › Structures |
| `structureFacts` and its residue rows | core › Facts |
| `unitRows` | core › Structures › Rows |
| `toUnits`: along static or overlay paths, the use's `via`, unobserved cut, a molecule whole | core › To units |
| The size gate (EC-167) | core › The size gate |
| `fromUnits`, the rest counted | core › From units |
| `colorRuns` | core › Colors |
| The round trip as a property test | core › From units |
| The side (`StructureSide`), targets structured-cloneable | viewer › The side |
| Out: hover by image, selection once a frame, jumps naming a Location | viewer › Out |
| In: the link's own hover, the selection revealed, focus one jump | viewer › In |
| Echoes never go back | viewer › Echoes |
| Clearing | viewer › Clearing |
| The link's options, not `ViewerConfig` | viewer › Structure link |
| A hover or selection pending is held | viewer › Pending |
| "Color by this track": a host's header entry over `colorRuns` | design-only: a host's entry (the Storybook's), never Seq\*'s |
| The fake host's adapter (the reference translation to MolViewSpec) | design-only: the Storybook's, not a package rule |
| A map's `priority`, sorted weakest first | core › Search |
| A use's `via` | core › Resolved view › Via |
| Not chosen silently; the ambiguity's hint | core › Search |
| Links carry no `via` | core › A section's `via` |
| `priority` and `via` in the 0.5.0 language | core › Versions |
| The fake structure viewer (Storybook) | design-only: the Storybook's test fixture |
| Datasets (p53 with 2OCJ, the complex with 1FIN) | design-only: example fixtures |
| The `structure` map kind kept, with provenance and status | core › Language › Maps |
| The link's options and defaults | viewer › Structure link |

## Tensions

- **A library, against a behavior every host gets.** Each host calls `linkStructure` itself; two hosts may wire
  differently. Decided: the library. Standing.
- **Seq\*'s words, against a shared vocabulary.** Each host translates; `structureFacts` and the fake's adapter are the
  reference. Standing.
- **`priority` global, against `via` per use.** Two sections may then light different residues for one hover: the point
  of a comparison. Decided: both. Standing.
- **A tie stays ambiguous**; the issue names the fix. Standing.
- **The fake, against a real 3D view**: it proves messages and maps, not occlusion or picking. The new mol viewer
  re-runs the same twins against its side. Standing.

## Self-check

The claim: **the library covers what a user does with a structure beside a sequence view.**
1. Assume it holds: every path is a viewer event or command carried onto units by pure functions, wired per link.
2. Negated: (a) the viewer publishes everything the wiring needs (negated: it did not; the four additions above); (b) one link
   (negated: two structure viewers: two links, each its own echoes); (c) the host keeps the link (negated: the viewer or
   the structure goes first: the disposer, the signal, the retire); (d) one map a pair (negated: `priority`, `via`);
   (e) synchronous messages (negated: a frame, latency, a stale hover: ownership and the fake's `latency`); (f) small
   (negated: the gate, coalescing); (g) headless (the translation alone, no viewer).
3. Held: (b) to (g) in the Decision; (a) the four additions, built before the wiring.

## Consequences

- **`seqx-core`**: `structure/` (`structureFacts`, `unitRows`, `toUnits`, `fromUnits`, `colorRuns`,
  `seq.structure.*` issues) as the `./structure` entry; `resolve/`: `priority`, `via`,
  `StaticModel.path(from, to, via?)`; `lang/`: the two params (0.5.0); `context/`: `SeqPorts.structure` removed. Tests:
  round trip, the gate monotone, order irrelevant, a priority never lengthens a path, `via` leaves other paths alone;
  EC-22 to EC-27 with a ranked pair, a tie and a `via`; EC-25, EC-48, EC-65, EC-118, EC-1 to EC-8; a quick
  representative (one 2OCJ hover carried within a bound); heavy: `maxUnits`-scale hover. `_spec`: Structures, Paths,
  Language.
- **`seqx-viewer`**: the four additions; `link-structure.ts`. Tests on a recording side: each path both ways,
  echoes never sent back, a stale `null`, one selection message a frame, focus with the jump's duration on a fake clock,
  clearing on each of the three ends. `_spec`: the models' additions, a Structure link section.
- **`seqx-examples`**: 2OCJ recorded (units, its `structure` map, the fake's chain table, one mmCIF, one script); an
  `alignment` map of the same pair.
- **Storybook**: the fake viewer and its adapter (MolViewSpec, validated against `mvs-builder`); stories with headless
  twins: Hover both ways, Selection both ways, Focus and back, Color by this track, Chain maps (SIFTS against an
  alignment by `priority` and `via`); the Design page.
- **Designs**: [maps.md](maps.md) (`priority`, `via`); [ports.md](ports.md) and
  [viewer-composition.md](viewer-composition.md) (no structure port). EC-24 gains "unless ranked or named" (the
  author's edit).

## What the viewer publishes for the wiring

Decided, the author, 2026-11-20: four public additions, built in ST2 before `linkStructure`, each useful to any host:
- **`HoverModel.setLocation(loc, owner)`**: a hover with no hit, drawn as a projected hover (no card), cleared only by
  its owner. Why: without it a link would fall back to a highlight channel's tint, which reads as a bookmark.
- **`AddressModel.revealLocation(loc)`**: reveal by Location, beside `focusLocation`.
- **`view.lastJump: { n, kind, location?, durationMs }`** (`focus`, `reveal`, `back`, `forward`, `fit`), each history
  step keeping the Location its jump named. Why: `view.jumps` counts but says neither where nor how long; the wiring, a
  story's camera and an analytics host read the same record.
- **`viewer.retiredSignal`** (an `AbortSignal`): `retired` is a getter and the dispose hooks are protected; any helper
  over a viewer needs the end.

Already published: `hover.version`, `selection.value`, `highlight.version`; the selection's announcement is the
viewer's own; frames come from `DomHost`.

## Decided by the author

2026-11-19: the vocabulary is the application's (no MolViewSpec types in Seq\*); a structure hover drawn as the hover;
a selection from the structure revealed; "color by this track" a snapshot; the fake viewer on its own chain table; the
`structure` map kind kept; two maps of one pair by "priority + via" (both in 0.5.0; `via` later if usage runs short).
2026-11-20: "this is just a wrapper around events the viewer should expose normally. we should just provide a library
that helps with the coordinate translation and the wiring, and not call it a port": no port, no behavior; translation
and `linkStructure`. The four viewer additions (`HoverModel.setLocation`, `revealLocation`, `view.lastJump`,
`retiredSignal`) accepted as recommended.
