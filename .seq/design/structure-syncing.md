# Structure syncing: a sequence view and a structure viewer, through one port

## Intent

The author, 2026-11-19: "let's only do a demo with fake 3D viewer to demonstrate how the syncing will work. There will
be a new mol viewer work after we are done with the sequence viewer which will replace prism, so i don't want to waste
work". A small Storybook component stands in for a structure viewer (its chains and residues as a schematic) and shows
every syncing path: hover and selection both ways, focus, "color by this track", a structure's chains mapped onto
sequences. The ports, maps and messages are the real ones. The Mol\* side (`Model` adapters, a live 3D view) waits for the
new mol viewer.

And on the words that cross: "this should be provided by the application layer anyway". Seq\* speaks its own locations;
the host's adapter translates them to its viewer's words (MolViewSpec, Mol\*'s loci) and back. The fake viewer is a host.

This design covers the **interaction side** of a structure: what crosses between a Seq\* viewer and a structure viewer,
in which words, how often, and through which map. It builds on:
- [molecules.md](molecules.md): units of a complex, with `chain`, `operator`, `observed` and `auth` numbering;
- [maps.md](maps.md): one map interface, the copy rule, the hull, the `structure` kind;
- [ports.md](ports.md): the structure port, named but empty so far;
- the seq data ports design (fourth slice): Seq\*'s nodes read from a host's representation, and annotations given back
  as the host's;
- the track header parts design (fourth slice): header actions as registry entries.

What it leaves to those designs:
- reading units from a host's structure (seq data ports);
- the header slot that draws the action (track header parts).

What it leaves to the host: every word of its viewer's vocabulary, MolViewSpec included.

## Precedents

- **The current Seq\*, in Prism and in its Mol\* story.** A link between the viewer and a 3D view (`StructureView3D`:
  `highlight`, `select`, `focus`, `onHover`, `onSelect`). Taken:
  - hover and selection both ways;
  - the size gate (EC-167): residue by residue within bounds, else whole chains or nothing;
  - a 3D hover with an owner of its own, so a late 3D leave never clears the pointer's hover (EC-141);
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
  - the hover card's port, on its behavior ([ports.md](ports.md));
  - `hostHistory`: a host's stack over Seq\*'s spec and a structure's state.

### Interactions in the precedents

| In the precedents | Here |
|---|---|
| hover a residue in the sequence, it lights in 3D (every copy when the sequence is the reference) | built: the `structure` behavior, size-gated |
| hover a residue in 3D, the sequence lights (the reference and that chain, not its sibling copies) | built: a hover with the structure as owner, projected under the copy rule (EC-25) |
| click or drag-select in the sequence, 3D selects | built: sent on change, at most once a frame |
| select in 3D (a click, a chain, "residues within 5 Å") and the sequence selects | built: any selection the adapter can place comes in as a Location |
| double-click a residue, and the camera focuses it | built: `address.focus` sends `focus` |
| focus in 3D, and the sequence jumps there | built: one jump, a history entry, animated through the animation port |
| back in the sequence view, the camera goes back | built: a jump that names an address sends it (`structure.followJumps`) |
| color the structure as a track is colored (Jalview, current Seq\*) | built: "color by this track", a host's header action; Seq\* gives the runs, the host builds its state |
| keep the structure's colors in step as the track's change (Jalview) | planned: `structure.colors.follow: live` (the structure lane), default `snapshot` |
| show an item in 3D from its menu (current Seq\*: its card) | planned: a host's `on: 'element'` entry over the item's Location (track header parts) |
| unobserved residues shown as gaps, author numbers on hover | built: the unit's `observed` and `auth` ([molecules.md](molecules.md)); never sent |
| map a sequence to a chain by SIFTS or by an alignment (Jalview, ChimeraX) | built: a `structure` map, an `alignment` map, or both; which one is taken: a map's `priority`, or a section's `via` |
| choose the chain mapping in a dialog (Jalview, ChimeraX) | refused in the viewer: a mapping is data (a map node by patch or the builder); computed by alignment, planned with services |
| superpose, measure, surfaces, the camera's own history | refused: the structure viewer's |
| the structure's B-factors or pLDDT as a track | not this design: the seq data ports' |

## Decision

### Who is who

- **Three parts: Seq\*, a structure viewer, and the host** that composes both. Seq\* imports nothing of a structure
  viewer and declares none of its words. The structure viewer imports nothing of Seq\*. The host's adapter is the only
  code that knows both. Why: host neutrality, and the vocabulary is the application's (decided).
- **The structure port lives on the `structure` behavior** (its `port` param), as the hover card's does. It is not in
  `spec.ports`, and core's `SeqPorts.structure` stub goes. Why: only the behavior talks to a structure viewer; a viewer
  without the behavior carries nothing; ports.md's rule that a UI port a behavior owns lives on it.
- **The port speaks Seq\*'s Locations on units.** Seq\* defines it and the host's adapter implements it:

  ```ts
  /** A place in the bound structures: a Location whose parts are on units of bound complexes (or a unit's atoms). */
  type StructureTarget = Location;

  interface StructurePort {
    hover(target: StructureTarget | null): void;
    select(target: StructureTarget | null): void;
    focus(target: StructureTarget, options: { readonly duration: number }): void;
    /** The structure viewer's own hover, selection and focus, translated by the adapter; each returns its unsubscribe. */
    onHover(fn: (target: StructureTarget | null) => void): () => void;
    onSelect(fn: (target: StructureTarget | null) => void): () => void;
    onFocus(fn: (target: StructureTarget) => void): () => void;
    readonly disposeWithViewer?: boolean;
    dispose?(): void;
  }
  ```

  A target is structured-cloneable (`Int32Array` ranges), so it crosses `postMessage`. Loading a state, "color by"
  included, is not a port method: the host holds its viewer.
- **What the adapter translates with**: `structureFacts(model, complex)` (`seqx-core`, `structure/`, pure): per unit its
  address, `chain`, `auth_chain`, `operator`, molecule and its kind (polymer, branched, other), length, `observed`
  ranges, and its `auth` numbering as a numbering entry (`label(position)`, `position(label)`, EC-1 to EC-8). Plain
  data, published once per resolve. Why: an adapter needs no Seq\* model walk, and every host translates from the same
  facts.

### Binding: which complex is which structure

- **A complex binds to a structure by the host's ref.** The behavior's `structures` param lists `{ complex, ref, via? }`.
  When absent, each complex binds to the ref equal to its id. Why: the spec stays portable; a host loading `2ocj` under
  its own ref says so once. `via` is the map choice for this binding (below).
- **Several structures in one viewer** (1FIN's crystal and the predicted model) are several bindings. The adapter tells
  them apart by the units' complexes. A unit of no bound complex in a target is ignored and counted in one debug issue.

### Translation: a Location onto units, in Seq\*

Pure, in `seqx-core`, `structure/`, tested in Node.

- **Outward**: carry the Location onto the units of the bound complexes along static or overlay paths, under the copy
  rule (`StaticModel.project`), through the binding's map choice. A complex from a source is overlay data, its units
  reached through the overlay's paths.
- **Unobserved positions are never sent**: cut by the unit's `observed`. Why: no atoms, nothing to light, and a message
  stays as long as what it lights.
- **A molecule unit (a ligand) is sent whole** (`structure.ligand: whole`): its one residue. A branched unit's residue
  is a position like any other; the adapter writes it by `auth` (mmCIF gives it no `label_seq_id`).
- **Inward**: the adapter's Location, read as given: parts on units, or on a molecule (every unit of it), or on a
  unit's atoms (collapsed onto residues unless a section shows that unit's `atoms` form, [molecules.md](molecules.md)).
  What the adapter could not place it leaves out and says once through `seq.structure.unplaced` (a count, the adapter's
  reason). A partial target places what it can.
- **The size gate (EC-167)**, outward: a target past `maxRuns` runs or `maxLength` positions becomes its whole units
  (`over: units`, one range each) or nothing (`skip`), per message kind. Why: a 4,000-unit assembly's selection must
  not become 4,000 × n runs every frame.
- **Round trip**: outward then inward through an identity adapter gives the Location's observed part on the units. A
  property test; each adapter's own round trip (the fake's: Location ↔ MolViewSpec) is its own test.

### The fake host's adapter (the reference translation)

In the Storybook beside the fake viewer, on `structureFacts`, written as the new mol viewer's adapter will be:
- a run on a polymer unit is `{ label_asym_id: chain, instance_id: operator, beg_label_seq_id, end_label_seq_id }`;
  `label_seq_id` is the position plus one, because a polymer unit's molecule is the entity's whole sequence
  (`entity_poly_seq.num`);
- a branched unit's residue is `auth_seq_id` (with `pdbx_PDB_ins_code`) through the unit's `auth` numbering;
- a molecule unit is `{ label_asym_id, instance_id }`; an atom adds `label_atom_id`;
- inward the reverse; `label_entity_id` and class selectors mean every unit of that molecule or class; a missing
  `instance_id` every instance of the chain; MolQL is not placed.

Why in the Storybook: the vocabulary is the application's; the fake's adapter proves `structureFacts` is enough and is
the pattern a real host copies.

### The paths

- **Hover, Seq\* to structure.** The hover's Location goes out through the translation each time its image changes; an
  unchanged image sends nothing (§relevant-rerenders). Hover on the reference (P04637) reaches every unit over its
  molecule (EC-48); a merged block (`×4`) stands for every member (EC-65).
- **Hover, structure to Seq\*.**
  - `HoverModel.setLocation(loc, owner)`: a hover with no hit, owned by `structure`, drawn as a projected hover (the
    outline, lighter) on every linked section; no card, not announced.
  - The pointer taking the hover takes ownership; a late `onHover(null)` clears only the structure's own (EC-141).
  - Under the copy rule, chain B lights B's unit and P04637, never A, C or D (EC-25).
  - Never sent back: echoes are cut by owner.
  - A hover, not a highlight channel's tint: hover looks the same everywhere ([ux-foundations.md](ux-foundations.md)
    section 6); the tint is a config value.
- **Selection, Seq\* to structure.** The selection's `shown` Location on change, coalesced to one message a frame;
  clearing sends `select(null)`.
- **Selection, structure to Seq\*.** Replaces the selection as a click does (`selection.set`, origin `structure`), never
  sent back; **reveals** it (`address.reveal`: nothing when shown, a pan when it fits, else a focus), so the zoom is
  kept; announced ("Selected from the structure: chain B, 248 to 273").
- **Focus, Seq\* to structure.** Every jump that names an address sends `focus(target)`: `address.focus`, a reveal that
  moved, a bookmark run, back or forward onto a step holding an address (`structure.followJumps`). The duration is the
  animation port's plan, so camera and sections move on one clock; `view.motion: cut` cuts both. A fit sends nothing.
- **Focus, structure to Seq\*.** `onFocus(target)` runs `address.focus`: one jump, one history entry. Each section shows
  the hull of the image ([maps.md](maps.md)): copies on P04637 fold to one range, a chimera's two ranges to their hull.
- **Viewport, no.** Zoom and pan never move the camera, nor the camera a section. Why: a sequence zoom has no structural
  equivalent, and following the hand both ways would fight.
- **Highlight channels: no, by default.** A structure viewer's highlight is its hover. With `structure.channels: state`,
  "color by" adds each channel's runs.

### "Color by this track"

- **A host's header action** (a registry entry of the track header parts design, never built into Seq\*). A click calls
  a pure core function and hands its runs to the host, which builds and loads its viewer's state.
- **`colorRuns(model, store, layer, units) → { complete, pending, runs }`** (`seqx-core`, `structure/`):
  - it evaluates the layer's compiled color channel over its data, the same compile as drawing it, so the structure gets
    the colors the track draws, a dark variant included (`structure.colors.tokens`);
  - each run is `{ color, target }`: a resolved CSS color and a `StructureTarget` (projected through the binding's map
    choice, unobserved cut, merged, gated as for select with `over: units` refused: a state is built once);
  - positions the layer leaves uncolored are in no run, so the viewer's own color shows (`structure.colors.missing`);
  - data not held yet is `pending` and asked for, never drawn as empty; the action waits for `complete`.
- **The host's state from the runs.** The fake's adapter writes a MolViewSpec state: the base tree unchanged, one
  `color` node per run per bound structure. A Mol\* host writes the same; another viewer writes its own. Seq\* has no
  writer for it (decided: the application's).
- **A snapshot** (`structure.colors.follow: snapshot`): the colors at the click. Again, or `structure.colors.clear`,
  restores the base state; a second track replaces the first. `live` (rebuilt when theme, data or visibility changes)
  is planned.
- **Undo** is the host's: it keeps the structure state in its stack beside the spec (`hostHistory`) if it wants.

### Choosing between two maps of one pair of systems

Today two maps joining the same two systems (SIFTS's `structure` map and an `alignment` map from P04637 to one chain's
molecule; a per-unit map beside a map onto the molecule) make every path through them ambiguous (EC-23, EC-24): the
search finds equal-length candidates that compose differently and takes none. Several maps of one entry that join
different pairs (one per chain) already resolve.

Decided, the author, 2026-11-19: "priority + via", both in the 0.5.0 bump; if usage runs short, `priority` in MK and
`via` later.

- **A map's `priority`** (an integer param of every map node, default 0). In `findPath`, after equal compositions are
  merged, distinct shortest candidates are ranked by their weakest map (the lowest priority among the maps they step
  through; unit steps have none); the top rank wins when it holds one composition; a tie stays ambiguous. Why:
  - one place decides for every use (follows, hover projection, the structure port, queries' `project`), so two
    sections never disagree about where a residue is;
  - the search stays shortest-first (EC-22): a priority never makes a longer path beat a shorter one;
  - the path cache stays keyed by its pair: the param is in the map entry's hash, so a carried search is dropped when a
    priority changes (paths.ts's `deps`);
  - weakest map, not a sum: a path is as trustworthy as its least trusted step.
- **A use's `via`**: a section's `via` (a ref to a map, or a list), and a structure binding's. For that use the search
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
  [packages.md](packages.md)). Its adapter, beside it, is the host's.
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

Rows this design adds or changes, with the columns of [ux-foundations.md](ux-foundations.md) section 8. "→ port" is
what the structure port is sent; "as the grammar" means unchanged in Seq\*.

| | select | extend | add / remove | focus | back | zoom in / out | pan | hide | edit |
|---|---|---|---|---|---|---|---|---|---|
| **an element** on a unit or its reference (a residue, a chip, an item) | as the grammar; → port `select` | as the grammar; → port | as the grammar; → port | as the grammar; → port `focus` | as the grammar; → port `focus` of the step landed on, when it holds an address | no: a sequence zoom has no structural equivalent | no: as zoom | no: hide a track | no: data |
| **a residue in the structure viewer** (its grammar) | click: Seq\*'s selection replaced, revealed, announced | the structure viewer's; arrives as its whole selection | Cmd/Ctrl+click: the same, its whole selection | double-click: `address.focus`, one jump | the structure viewer's own; Seq\*'s Alt+Left is Seq\*'s | the camera's own | the camera's own | no | no |
| **a range** | drag; → port once a frame while dragging | Shift+drag; → port | Cmd/Ctrl+drag, Alt+drag; → port | Enter; → port | Alt+Left; → port | — | — | — | — |
| **the selection** | — | — | — | Enter; → port. `structure.focus-selection` focuses the camera alone (no default key; palette and menu) | Alt+Left | — | — | — | bookmark it: no port effect until it runs |
| **a section of a complex** (copies merged `×N`) | as [molecules.md](molecules.md); a merged block → every member | as molecules.md | as molecules.md | as molecules.md; → port | Alt+Left | as its layout | as its layout | no | `copies` by patch; `via` by patch (the projection onto it changes) |
| **a map** (a `structure` or `alignment` map of one pair) | — | — | — | — | — | — | — | — | `priority` by patch: every use re-resolves; no gesture (data) |
| **a track or layer** | — | — | — | — | — | — | — | hiding a colored track leaves the snapshot (`live` clears it) | "color by this track": a host's entry; again, or `structure.colors.clear`, restores |
| **the view** | — | — | — | 0 fits: no port message | Alt+Left; → port when the step holds an address | no: viewport, as above | no | — | — |
| **a bookmark or link** | as built; → port `select` | — | — | runs its action; → port `focus` | Alt+Left; → port | — | — | — | — |
| **a follower** | as an element | as an element | as an element | as an element | as an element | moves its leader, nothing sent | moves its leader, nothing sent | — | — |
| **a navigator** | no: a pixel is many positions | — | — | no port effect | — | no: its window is the viewport | no | — | — |

## Configuration

On the `structure` behavior (params, `ViewerConfig` items) unless the table says otherwise.

| Item | Kind | Default | Why | Alternatives |
|---|---|---|---|---|
| `port` | behavior param | none (the behavior does nothing) | a host supplies its adapter, as the hover card's | the fake's, the new mol viewer's |
| `structures` | behavior param | each complex bound to the ref equal to its id, no `via` | the spec names no host ref | explicit `{ complex, ref, via? }[]` |
| `structure.hover` | ViewerConfig | `both` | the precedents link hover both ways | `out`, `in`, `none` |
| `structure.select` | ViewerConfig | `both` | as hover | `out`, `in`, `none` |
| `structure.focus` | ViewerConfig | `both` | jumps are what a reader follows | `out`, `in`, `none` |
| `structure.hoverAs` | ViewerConfig | `hover` | hover looks the same everywhere | `channel` (`links.hoverChannel`'s tint) |
| `structure.selectFrom` | ViewerConfig | `replace` | a click in 3D is a click | `keep-unplaced` |
| `structure.revealOnSelect` | ViewerConfig | `reveal` | the reader sees it, the zoom kept | `focus`, `none` |
| `structure.followJumps` | ViewerConfig | `focus` | back and forward bring the camera along | `none` |
| `structure.limits.hover` | ViewerConfig | 2,000 runs, 20,000 positions, `over: units` | the precedent's bounds | any; `over: skip` |
| `structure.limits.select` | ViewerConfig | 20,000 runs, 200,000 positions, `over: units` | as hover, larger: sent once | any |
| `structure.ligand` | ViewerConfig | `whole` | the precedent's; a glycan's residue, not the whole glycan | `atoms` |
| `structure.channels` | ViewerConfig | `none` | a structure viewer's highlight is its hover | `state` (channels' runs in "color by") |
| `structure.colors.follow` | the action's option | `snapshot` | "ad hoc" at the click | `live` (planned) |
| `structure.colors.tokens` | the action's option | the viewer's current tokens | the structure gets the colors the track draws | `light` |
| `structure.colors.missing` | the action's option | no run (the viewer's color) | Seq\* colors only what it has data for | a token (`seq-neutral`) |
| a map's `priority` | map param (language 0.5.0) | 0 | equal maps stay ambiguous until the writer ranks them | any integer |
| a section's `via` | section param (language 0.5.0) | none (the ranked search) | a per-use choice where a figure shows both maps | a map ref or a list |
| unobserved positions in messages | fixed | never sent | no atoms, nothing to light | — |
| the message vocabulary | fixed | Seq\*'s Locations on units | the vocabulary is the application's (decided) | MolViewSpec types, Mol\* `Loci` (refused: the host's) |
| the fake viewer's `latency` | story arg | 0 | exercises the async path when set | ms |

## Tensions

- **Seq\*'s words, against a shared public vocabulary.** Each host writes its translation, and two hosts may translate a
  branched residue or an atom differently. Decided: the application's. Mitigation: `structureFacts` carries every
  identifier the translation needs, and the fake's adapter is the reference. Standing.
- **`priority` global, against `via` per use.** A global rank is one answer everywhere; `via` lets two sections differ,
  so a hover on one lights another position on the other, which is the point of a comparison and a surprise elsewhere.
  Decided: both. Standing.
- **A tie stays ambiguous.** A writer who adds a second map of equal rank loses the first path. Recommended: kept; the
  issue names the fix. Standing.
- **A snapshot, against Jalview's live colors.** A snapshot goes stale when the theme is edited; live reloads a state per
  edit. Snapshot, `live` planned. Standing.
- **Replace, against keeping the unplaced.** A 3D click clears a gene range Seq\* had selected; `keep-unplaced` is
  configurable. Standing.
- **One structure viewer a Seq\* viewer.** Two 3D views need a host adapter that fans out. Standing until a host asks.
- **The fake, against a real 3D view.** No occlusion, picking ambiguity, camera fidelity, or loader delays beyond
  `latency`: it proves the messages and the maps. The new mol viewer's integration re-runs the same twins against its
  adapter. Standing.
- **The binding by ref, against a spec that knows its structure.** The spec names the data, the host the refs.
  Standing.

## Self-check

The claim: **this design covers what a user does with a structure beside a sequence view.**

1. **Assume it holds.** Every syncing path is a Location on units through one port, translated by the host; colors are
   runs the host turns into its state. Every structure position is a unit's position, reached by one chosen map, under
   the copy rule and the unit's numbering.
2. **What it rests on, negated.**
   - (a) *One structure, one viewer.* Negated: a crystal and a model together; two 3D views.
   - (b) *Seq\* already holds the units.* Negated: the user opens another entry in 3D.
   - (c) *Seq\* acts, 3D follows.* Inverse: 3D acts (a 5 Å shell, a chain).
   - (d) *Positions match one to one.* Negated: unobserved residues, insertion codes, altlocs, a branched residue
     without `label_seq_id`, CRO over three positions, a chimera.
   - (e) *Messages are synchronous and in order.* Negated: a frame, latency, a stale hover.
   - (f) *A section shows the unit.* Negated: only P04637, or only the gene.
   - (g) *Colors are made once.* Negated: theme edited, data late, track hidden, a second track; the inverse, uncoloring.
   - (h) *Focus has no history.* Inverse: back.
   - (i) *Small.* Negated: 60 copies, a whole-chain selection.
   - (j) *One map joins a sequence and a chain.* Negated: SIFTS and an alignment; a per-unit map beside one onto the
     molecule; the inverse, a figure that wants both.
   - (k) *Seq\* knows the structure viewer's words.* Negated: it knows none (decided).
3. **The negated world.** (a) bindings by ref, a fan-out adapter; (b) the seq data ports; (c) every path two-way,
   echoes cut by owner and origin; (d) the cut of the unobserved, `auth` with EC-1 to EC-8, atoms collapsing, CRO the
   molecules design's gap 6; (e) structured-cloneable messages, ownership (EC-141), `latency`; (f) projection onto
   unshown systems, the hull; (g) snapshot, clear, replace, pending, `live` planned; (h) `followJumps`; (i) the gate and
   coalescing; (j) `priority` and `via`, ties kept ambiguous; (k) `structureFacts` and the host's adapter, with the fake
   as the reference.
4. **Held both.**
   - Integrated: (a) refs, (c), (d) except CRO's grouping, (e), (f), (g) snapshot, (h), (i), (j), (k).
   - Planned with a lane: (b) seq data ports; `live` colors and "show in 3D" from an element's menu (structure lane,
     track header parts); CRO's grouping (ML2); the map choice (MK).
   - Standing: two 3D views; the fake's fidelity; replace against keep; each host's translation; ties.

## Consequences

The fourth slice builds it after ML2 (units) and MK (the map choice), beside the track header parts lane.

- **`seqx-core`**:
  - `structure/`: `structureFacts`, the translation (outward and inward on Locations, the gate, `seq.structure.*`
    issues), `colorRuns`;
  - `resolve/`: `priority` in the search (rank by weakest map after merging equal compositions), `via` (rivals left
    out), `StaticModel.path(from, to, via?)`, the ambiguous issue's hint; `lang/`: the two params (0.5.0);
  - `context/`: `SeqPorts.structure` removed.
  - Tests: property tests (the round trip through an identity adapter; the gate monotone; order irrelevant; a priority
    never lengthens a path; raising a map's priority over a tie picks it; `via` leaves paths away from the pair
    unchanged); EC-22 to EC-27 on every encoding with a ranked pair, a tie and a `via`; EC-25; EC-48 and EC-65; EC-118;
    EC-1 to EC-8 through `structureFacts`; an unplaced report counted once; a quick representative of the hover target
    (one message for 2OCJ within a bound, none for an unchanged image); heavy: `maxUnits`-scale hover gated to units.
  - `_spec`: Structures (new), Paths (the rank and `via`), Language, Context (ports).
- **`seqx-viewer`**: the `structure` behavior (paths, echoes, coalescing, `followJumps`, bindings with `via`);
  `HoverModel.setLocation(loc, owner)`; the announcement; commands `structure.focus-selection` and
  `structure.colors.clear`; the config items. Tests: each path both ways on a recording port; echoes never sent back; a
  stale `onHover(null)`; one message a frame while dragging; focus with the animation port on a fake clock. `_spec`:
  HoverModel, Behaviors, Ports (Structure), Commands, Configuration.
- **`seqx-examples`**: 2OCJ recorded (its units, its `structure` map, the fake's chain table, from one mmCIF by one
  script); an `alignment` map of the same pair for the choice; views for p53 and the Complex.
- **Storybook**: the fake viewer and its adapter (Location ↔ MolViewSpec, validated against `mvs-builder`, its own
  round-trip tests); stories with headless twins on a fake clock (`latency` set in one): Hover both ways, Selection both
  ways, Focus and back, Color by this track (conservation on 2OCJ, confidence on 1FIN's model), Chain maps (copies,
  author numbering, unobserved gaps, the construct offset, SIFTS against an alignment by `priority` and by `via`); a
  `Seq* (next)/Design/Structure syncing` page.
- **Drawing against the data**: each bound unit's observed residue in the fake viewer has `colorRuns`' color at its
  preimage; `colorRuns` equals the track's drawn colors at a letters zoom; every lit cell in the fake viewer is in the
  image of Seq\*'s hover, and the reverse.
- **Lanes**: MK (the map choice, with `segments` and the `structure` kind); ST1 (strong): `structureFacts`, the
  translation, `colorRuns`; ST2 (strong, after ST1): the behavior and the hover owner; ST3 (after ST2 and the track
  header parts lane): the fake viewer and its adapter, the stories and twins, the Design page.
- **Designs that change when agreed**: [ports.md](ports.md) (the structure row on the behavior, Seq\*'s Locations);
  [maps.md](maps.md) (the `structure` kind kept onto a unit or molecule, unobserved from the unit, `priority` and
  `via`); [viewer-composition.md](viewer-composition.md) (`Ports.structure` removed, "Inside a Mol\* plugin" naming the
  behavior's port and the host's adapter); [molecules.md](molecules.md) Tensions ("a unit and maps.md's structure
  system": the structure map kept). The edge-case catalog's EC-24 gains "unless ranked or named" (the author's edit).

## Decided by the author

2026-11-19; each answer is folded into the sections named.

1. Port on the behavior, core's stub removed: Who is who.
2. **Changed**: no MolViewSpec types in Seq\*; the port speaks Locations on units; the host's adapter translates (the
   fake's in the Storybook): Who is who, the fake host's adapter, "Color by this track".
3. Structure hover as the hover: The paths. 4. Reveal on select: The paths. 5. Snapshot colors. 6. The fake's own chain
   table.
7. **No**: the `structure` kind kept (Datasets); two maps of one pair: "Choosing between two maps".
8. Two maps of one pair: "priority + via", both in 0.5.0 (`via` later if usage runs short).
