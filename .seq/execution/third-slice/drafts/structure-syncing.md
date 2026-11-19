# Structure syncing: a sequence view and a structure viewer, through one port

A fourth-slice draft. Links are written as they will read once the file moves to `design/`.

## Intent

The author, 2026-11-19: "let's only do a demo with fake 3D viewer to demonstrate how the syncing will work. There will
be a new mol viewer work after we are done with the sequence viewer which will replace prism, so i don't want to waste
work". The plan: a small Storybook component stands in for a structure viewer (its chains and residues as a schematic,
taking MolViewSpec-shaped states and sending hover, selection and focus as a structure viewer would) and shows every
syncing path: hover and selection both ways, focus, "color by this track", a structure's chains mapped onto sequences.
The ports, maps and messages are the real ones. The Mol\* side (`Model` adapters, a live 3D view) waits for the new mol
viewer.

This design covers the **interaction side** of a structure: what crosses between a Seq\* viewer and a structure viewer,
in which words, and how often. It builds on what is already decided elsewhere:
- [molecules.md](molecules.md): units of a complex, with `chain`, `operator` and `auth` numbering;
- [maps.md](maps.md): one map interface, the copy rule, the hull;
- [ports.md](ports.md): the structure port, named but empty so far;
- the seq data ports design (fourth slice): Seq\*'s nodes read from a host's representation, and annotations given back
  as host representations;
- the track header parts design (fourth slice): header actions as registry entries.

What it leaves to those designs:
- reading units from a Mol\* `Model` (seq data ports);
- the MolViewSpec writer as a representation entry (seq data ports);
- the header slot that draws the action (track header parts).

## Precedents

- **The current Seq\*, in Prism and in its Mol\* story.** A link between the viewer and a 3D view (`StructureView3D`:
  `highlight`, `select`, `focus`, `onHover`, `onSelect`). Taken:
  - hover and selection both ways;
  - the size gate (EC-167): residue by residue within `maxItems` and `maxLength`, else whole chains or nothing;
  - a 3D hover with an owner of its own, so a late 3D leave never clears the pointer's hover (EC-141);
  - a non-polymer position lighting its whole residue;
  - chain instances keyed by `label_asym_id` and operator (EC-181);
  - "show in 3D" as MolViewSpec composed from the viewer's state at that moment.

  Not taken:
  - the link speaking Mol\* `Loci`: it imports Mol\* and the plugin. Plain MolViewSpec selectors cross here instead.
  - The host registering external coordinate systems with their maps: units are already the spec's systems
    ([molecules.md](molecules.md)).
  - A selection from 3D focusing the sequence view: it reveals here (Decision).
- **Mol\*'s sequence panel.** Hover and click a residue, and the 3D view lights or selects it; a chain picker. Taken: the
  gestures. Not taken: the picker (a section of a complex shows its blocks).
- **Jalview with a structure.** A sequence is mapped to chains by SIFTS or by alignment. Hover and selection work both
  ways. "Colour by Jalview" colors the structure as the alignment is colored, kept live. Taken: the paths and the
  colors. Not taken: live colors by default (a snapshot here, live as a config value), and choosing a chain mapping in a
  dialog (the mapping is data).
- **ChimeraX's sequence viewer.** Chains are associated by alignment, selection works both ways, and the structure is
  colored by a sequence attribute. Taken: association as a map. Not taken: attribute assignment (the structure's data
  coming back to Seq\* is the seq data ports' job).
- **RCSB's 1D feature viewer and PDBe's ProtVista with Mol\*.** Hover a feature and its residues light in 3D. Click
  selects and focuses. Taken: a feature's extent crossing as residues. Not taken: click focusing (focus is the
  double-click in this grammar).
- **MolViewSpec** (`packages/mvs/builder`, `param-types.ts`). Component expressions use `label_asym_id`, `instance_id`,
  `label_seq_id` and `beg_`/`end_label_seq_id`, `auth_seq_id` with `pdbx_PDB_ins_code`, `label_atom_id`, and the
  selectors `polymer` … `water`. There are `color` nodes with selectors, and a `focus` node. MolViewSpec is state only:
  it has no hover or selection. Taken: its expressions are the messages' vocabulary, and its `color` and `focus` nodes
  are what "color by this track" builds. Not taken: MolQL selectors (Seq\* never parses them; they are reported and
  ignored).
- **This codebase.** These pieces are built and reused:
  - `HoverModel`: owners, and projection over linked sections;
  - `SelectionModel` methods;
  - `AddressModel`: `focus`, `reveal`, and `hoverAddress` on a channel of its own (EC-182);
  - core's `projectOntoView` and `StaticModel.project`: onto systems no section shows;
  - the animation port: one clock for the sequence view and the camera ([stories.md](stories.md));
  - the hover card's port, which lives on its behavior (ports.md, HC);
  - `hostHistory`: a host's stack over Seq\*'s spec and a structure's state.

### Interactions in the precedents

| In the precedents | Here |
|---|---|
| hover a residue in the sequence, it lights in 3D (every copy when the sequence is the reference) | built: the `structure` behavior, size-gated |
| hover a residue in 3D, the sequence lights (the reference and that chain, not its sibling copies) | built: a hover with the structure as owner, projected under the copy rule (EC-25) |
| click or drag-select in the sequence, 3D selects | built: sent on change, at most once a frame |
| select in 3D (a click, a chain, "residues within 5 Å") and the sequence selects | built: any selector set comes in as a Location |
| double-click a residue, and the camera focuses it | built: `address.focus` sends `focus` |
| focus in 3D, and the sequence jumps there | built: one jump, a history entry, animated through the animation port |
| back in the sequence view, the camera goes back | built: a jump that names an address sends it (`structure.followJumps`) |
| color the structure as a track is colored (Jalview, current Seq\*) | built: "color by this track", a host's header action building a MolViewSpec state |
| keep the structure's colors in step as the track's change (Jalview) | planned: `structure.colors.follow: live` (the structure lane), default `snapshot` |
| show an item in 3D from its hover card (current Seq\*) | planned: the same builder from an item's extent, a host's card action (the track header parts lane) |
| unobserved residues shown as gaps, author numbers on hover | built: the unit's `observed` and `auth` ([molecules.md](molecules.md)); never sent to 3D |
| choose which chain a sequence maps to (Jalview, ChimeraX) | refused in the viewer: a mapping is data (a map node, by patch or by the builder); computed by alignment, planned with services |
| superpose, measure, surfaces, the camera's own history | refused: the structure viewer's |
| the structure's B-factors or pLDDT as a track | not this design: the seq data ports' (annotations read from the host's representation) |

## Decision

### Who is who

- **Three parts: Seq\*, a structure viewer, and the host** that composes both. Seq\* imports nothing of a structure
  viewer. The structure viewer (the fake one here, the new mol viewer later) imports nothing of Seq\*. They meet
  through one port of plain data. Why: host neutrality. Messages that are only JSON also cross a `postMessage` to a
  viewer in another frame.
- **The structure port lives on the `structure` behavior** (its `port` param), as the hover card's does. It is not in
  `spec.ports`, and core's `SeqPorts.structure` stub goes. Why:
  - only the behavior talks to a structure viewer; the core and a headless render never do;
  - a viewer without the behavior carries nothing;
  - ports.md's own later rule: a UI port that a behavior owns lives on that behavior (gap 1).
- **The port.** Seq\* defines it and the structure viewer's adapter implements it:

  ```ts
  /** Where a structure viewer is acted on: per MolViewSpec structure (its node's `ref`), MolViewSpec component expressions. */
  type StructureTarget = readonly { readonly ref: string; readonly selector: readonly MvsComponentExpression[] }[];

  interface StructurePort {
    hover(target: StructureTarget | null): void;
    select(target: StructureTarget | null): void;
    focus(target: StructureTarget, options: { readonly duration: number }): void;
    /** The structure viewer's own hover, selection and focus; each returns its unsubscribe. */
    onHover(fn: (target: StructureTarget | null) => void): () => void;
    onSelect(fn: (target: StructureTarget | null) => void): () => void;
    onFocus(fn: (target: StructureTarget) => void): () => void;
    readonly disposeWithViewer?: boolean;
    dispose?(): void;
  }
  ```

  `MvsComponentExpression` is the subset of MolViewSpec's `ComponentExpression` fields named below. It is declared in
  Seq\* as a plain type, and Seq\* does not import `@molstar/mvs-builder` (gap 2). Loading a state, "color by" included,
  is not a port method. The host holds its own viewer and loads states into it, so Seq\* never holds a 3D state.

### Binding: which complex is which structure

- **A complex binds to a MolViewSpec structure by ref.** The behavior's `structures` param lists `{ complex, ref }`
  pairs. When the param is absent, each complex in the spec binds to the ref equal to its id. Why: the spec stays
  portable (it names no host ref), and a host loading `2ocj` under its own ref says so once.
- **Several structures in one structure viewer** (1FIN's crystal and the predicted model) are several bindings. A
  target names each one by `ref`. A target naming a ref with no binding is ignored and counted in one debug issue.

### Translation: a Location on units and MolViewSpec selectors

The translation lives in `seqx-core`, `structure/`. It is pure and tested in Node.

- **Outward**, a Location to a target:
  - Carry the Location onto the units of the bound complexes, along static or overlay paths, under the copy rule (the
    core's `StaticModel.project`, onto systems no section shows). A complex from a source is overlay data
    ([molecules.md](molecules.md)), so its units are reached through the overlay's paths.
  - A run on a polymer unit becomes `{ label_asym_id: chain, instance_id: operator, beg_label_seq_id, end_label_seq_id }`.
    `label_seq_id` is the molecule's position plus one, because a polymer unit's molecule is the entity's whole sequence
    (`entity_poly_seq.num`). The relation is exact, with no table.
  - A single residue becomes `label_seq_id`.
  - A molecule unit (a ligand) shown whole becomes `{ label_asym_id, instance_id }`.
  - A branched unit's residue becomes `auth_seq_id` through the unit's `auth` numbering: mmCIF gives a branched
    residue no `label_seq_id`.
  - An atom of the `atoms` form becomes `label_atom_id` on its residue's expression.
  - **Unobserved positions are never sent**: they are cut by the unit's `observed`. Why: they have no atoms, so the
    structure viewer would light nothing, and a message stays as long as what it lights.
- **Inward**, a target to a Location:
  - An expression with `label_seq_id` or its range maps to positions on the unit with the same `label_asym_id` and
    `instance_id`. A missing `instance_id` means every instance of that chain.
  - `auth_seq_id` with `pdbx_PDB_ins_code` is read through the unit's `auth` numbering, with the labels rules EC-1 to
    EC-8.
  - `label_asym_id` alone means the whole unit. `label_entity_id` and the class selectors (`polymer`, `ligand`, …) mean
    every unit of that molecule or class.
  - An atom-level expression collapses onto its residue, unless a section shows that unit's `atoms` form (the group map
    at a path's end, [molecules.md](molecules.md)).
  - MolQL and fields Seq\* cannot place are reported once (`seq.structure.selector`) and skipped. A partial target
    places what it can.
- **The size gate (EC-167).** Outward, a target past `maxItems` expressions or `maxLength` positions becomes its whole
  chains (`over: chains`) or nothing (`skip`), per message kind (Configuration). Why: a 4,000-unit assembly's selection
  must not become 4,000 × n expressions every frame.
- **Round trip.** Inward then outward gives the same target, normalized: sorted, merged, and the unobserved cut.
  Outward then inward gives the Location's observed part on the units. Both are property tests.

### The paths

- **Hover, Seq\* to structure.**
  - The hover's Location goes out through the translation each time its image changes. An unchanged image sends
    nothing (§relevant-rerenders), so per pointer move it costs a comparison.
  - Hover on the reference sequence (P04637) reaches every unit over its molecule (EC-48): all four 2OCJ chains light.
    Hover on a merged block (`×4`) stands for every member (EC-65).
- **Hover, structure to Seq\*.**
  - `HoverModel` gains `setLocation(loc, owner)`: a hover with no hit, owned by `structure`. It is drawn as a projected
    hover is (the outline, lighter) on every linked section. It shows no card and is not announced.
  - The pointer taking the hover takes ownership. A late `onHover(null)` from the structure clears only its own hover
    (EC-141).
  - Under the copy rule, hovering chain B in 3D lights B's unit and P04637, never units A, C or D (EC-25).
  - A hover that came from the structure is never sent back. Echoes are cut by owner.
  - Why a hover and not a highlight channel: hover looks the same everywhere ([ux-foundations.md](ux-foundations.md),
    section 6), and a channel's tint would mean a bookmark. The alternative stays a config value (gap 3).
- **Selection, Seq\* to structure.** The selection's `shown` Location goes out on change, coalesced to one message a
  frame (a drag moves it every pointer event). Clearing it sends `select(null)`.
- **Selection, structure to Seq\*.**
  - The target's Location replaces the selection, as a click does (`selection.set` with a session origin of
    `structure`). It is never sent back. Unplaceable parts are dropped with the issue above.
  - It **reveals** what it selected (`address.reveal`: nothing when it is shown, a pan when it fits, else a focus), so
    the reader sees it without losing the zoom. The value `focus` is in Configuration (gap 4).
  - It is announced: "Selected from the structure: chain B, 248 to 273".
- **Focus, Seq\* to structure.**
  - Every jump that names an address sends `focus(target)`: `address.focus`, a reveal that moved, a bookmark run, and
    back or forward landing on a step that holds an address (`structure.followJumps`).
  - The duration is the plan's from the animation port, so the camera and the sections move together on the host's
    clock. With `view.motion: cut`, both cut.
  - A fit sends nothing: it names no residues, and resetting the camera is the structure viewer's.
- **Focus, structure to Seq\*.** `onFocus(target)` runs `address.focus` on its Location, which is one jump and a history
  entry. Each section shows the hull of the Location's image in display offsets ([maps.md](maps.md), the hull): copies
  on P04637 fold to one range, and a chimera's two ranges become their hull.
- **Viewport, no.** Zoom and pan of a section never move the camera, and the camera never moves a section. Why: a zoom
  on a sequence has no structural equivalent, and following the hand both ways would fight.
- **Highlight channels: no, by default.** MolViewSpec has no transient emphasis, and a structure viewer's highlight is
  its hover. With `structure.channels: state`, "color by" adds each channel as a color node over its residues.

### "Color by this track"

- **A host's header action** (a registry entry from the track header parts design, never built into Seq\*). A click
  calls a pure core function and loads the state it returns into the host's own viewer.
- **`colorRuns(model, store, layer, units) → { complete, pending, runs }`** (`seqx-core`, `structure/`):
  - It evaluates the layer's compiled color channel over its data. This is the same compile as drawing it, so the
    structure gets the colors the track draws, a scheme's dark variant included (`structure.colors.tokens`).
  - It projects each run onto the units and cuts the unobserved positions.
  - Data that is not held yet is reported as `pending` and asked for, never drawn as empty. The action waits for
    `complete` (a region over a structure is small).
- **A MolViewSpec state from the runs.** The `mvs` representation writer comes from the seq data ports design, and this
  design states what that writer must produce for colors:
  - the host's base tree, unchanged;
  - for each bound structure, one `color` node per distinct color, its selector the runs as expressions (merged, gated
    as for select, `over: chains` refused here: a state is built once);
  - positions the layer leaves uncolored get no node, so the structure viewer's own color shows (`structure.colors.missing`).
- **A snapshot** (`structure.colors.follow: snapshot`). The colors are those at the click. Activating the action again,
  or `structure.colors.clear`, loads the base state again. Coloring a second track replaces the first. Why: the plan's
  "the click builds a MolViewSpec state ad hoc". `live` (rebuilt when the layer's theme, data or visibility changes) is
  the planned alternative (gap 5).
- **Undo** is not Seq\*'s: the state is the host's. A host that wants it undoable keeps the structure state in its
  stack beside the spec (`hostHistory`).

### The fake structure viewer (Storybook)

- **A Storybook component** (`packages/seq/storybook`, under `_stories/support/`) that imports nothing of Seq\*. It is
  written on `@molstar/ui` and validates states with `@molstar/mvs-builder`, which the Storybook may import for
  showcases ([packages.md](packages.md)).
- **Its structure is its own recorded chain table** per entry (`label_asym_id`, `auth_asym_id`, entity, operators, and
  per residue `label_seq_id`, `auth_seq_id`, `pdbx_PDB_ins_code`, `label_comp_id`, observed). It is recorded from the
  mmCIF by the examples' record script, independent of Seq\*'s nodes. Why: the mapping must cross two independently
  held representations, as it will with Mol\*'s `Model` (gap 6).
- **It draws** each chain instance as a row of residue cells, with unobserved residues as a dashed gap, ligands as
  chips, and copies as separate rows. It shows the state's `color` nodes, a hover outline, the selection, and a "camera"
  frame that the `focus` node and `focus()` move (animated over the given duration).
- **Its own grammar** (the host's, not Seq\*'s): hover; click selects; Cmd/Ctrl+click toggles; Esc clears;
  double-click focuses; a "select chain" and a "within 5 Å" button stand for structure-only selections. It sends
  exactly what a real viewer would: targets of expressions.
- **`latency` option** (0 by default): delays its messages, so the stories and twins exercise the asynchronous path (a
  stale hover arriving after the pointer moved).

### Datasets

- **p53 with 2OCJ**: four copies of the DNA-binding domain, author numbering equal to P04637's, unobserved ends, and a
  zinc per chain. It is not recorded yet: the fourth slice records it as [examples.md](examples.md) says.
- **Complex with 1FIN**: two copies of the complex, cyclin A as a construct (offset 172 onto P20248), and ATP as a
  molecule unit.
- **The structure map is no new kind** (gap 7). A reference sequence maps to the entity's molecule by any map kind:
  - `offset` for both datasets;
  - `segments` for SIFTS with holes, built with the services lane;
  - `alignment` for a chain with no SIFTS.

  Per-position status comes from what already holds it: unobserved from the unit's `observed`, conflict from the
  letters compared (as codon status is). maps.md's "structure system is a chain instance" reads as a unit, the noun
  molecules.md recommended. Why: no arithmetic of its own, and one less kind to read, write and convert.

## The interaction matrix's cells

These are the rows this design adds or changes, with the columns of [ux-foundations.md](ux-foundations.md) section 8.
"→ 3D" is what the structure port is sent; "as the grammar" means the cell is unchanged in Seq\*.

| | select | extend | add / remove | focus | back | zoom in / out | pan | hide | edit |
|---|---|---|---|---|---|---|---|---|---|
| **an element** on a unit or its reference (a residue, a chip, an item) | as the grammar; → 3D `select` | as the grammar; → 3D | as the grammar; → 3D | as the grammar; → 3D `focus` | as the grammar; → 3D `focus` of the step landed on, when it holds an address | no: a sequence zoom has no structural equivalent | no: as zoom | no: hide a track | no: data |
| **a residue in the structure viewer** (its grammar) | click: Seq\*'s selection replaced, revealed, announced | the structure viewer's; arrives as its whole selection | Cmd/Ctrl+click: the same, its whole selection | double-click: `address.focus`, one jump | the structure viewer's own; Seq\*'s Alt+Left is Seq\*'s | the camera's own; nothing in Seq\* | the camera's own | no | no |
| **a range** | drag; → 3D once a frame while dragging | Shift+drag; → 3D | Cmd/Ctrl+drag, Alt+drag; → 3D | Enter; → 3D | Alt+Left; → 3D | — | — | — | — |
| **the selection** | — | — | — | Enter; → 3D. `structure.focus-selection` focuses the camera alone (no default key; in the command palette and the context menu) | Alt+Left | — | — | — | bookmark it: no 3D effect until it runs |
| **a section of a complex** (copies merged `×N`) | as [molecules.md](molecules.md); a merged block → 3D every member | as molecules.md | as molecules.md | as molecules.md; → 3D | Alt+Left | as its layout | as its layout | no | `copies` by patch: no 3D effect |
| **a track or layer** | — | — | — | — | — | — | — | hiding a colored track leaves the snapshot (`live` clears it) | "color by this track": the header action, a host's entry, Enter on the focused header's action; again, or `structure.colors.clear`, restores |
| **the view** | — | — | — | 0 fits: no 3D message, because a fit names no residues | Alt+Left; → 3D when the step holds an address | no: viewport, as above | no | — | — |
| **a bookmark or link** | as built; → 3D `select` | — | — | runs its action; → 3D `focus` | Alt+Left; → 3D | — | — | — | — |
| **a follower** | as an element | as an element | as an element | as an element | as an element | moves its leader, nothing in 3D | moves its leader, nothing in 3D | — | — |
| **a navigator** | no: a pixel is many positions | — | — | no 3D effect | — | no: its window is the viewport | no | — | — |

## Configuration

All items are on the `structure` behavior (its params, `ViewerConfig` items) unless the table says otherwise.

| Item | Kind | Default | Why | Alternatives |
|---|---|---|---|---|
| `port` | behavior param | none (the behavior does nothing) | a host supplies its viewer's adapter, as the hover card's | the fake viewer's, the new mol viewer's |
| `structures` | behavior param | each complex bound to the ref equal to its id | the spec names no host ref | explicit `{ complex, ref }[]` |
| `structure.hover` | ViewerConfig | `both` | the precedents all link hover both ways | `out`, `in`, `none` |
| `structure.select` | ViewerConfig | `both` | as hover | `out`, `in`, `none` |
| `structure.focus` | ViewerConfig | `both` | jumps are what a reader follows | `out`, `in`, `none` |
| `structure.hoverAs` | ViewerConfig | `hover` | hover looks the same everywhere | `channel` (`links.hoverChannel`'s tint) |
| `structure.selectFrom` | ViewerConfig | `replace` | a click in 3D is a click | `keep-unplaced` (Seq\*'s parts with no structure image kept) |
| `structure.revealOnSelect` | ViewerConfig | `reveal` | the reader sees it, and the zoom is kept | `focus` (the precedent's), `none` |
| `structure.followJumps` | ViewerConfig | `focus` | back and forward bring the camera along | `none` |
| `structure.limits.hover` | ViewerConfig | 2,000 expressions, 20,000 positions, `over: chains` | the precedent's, measured there | any; `over: skip` |
| `structure.limits.select` | ViewerConfig | 20,000, 200,000, `over: chains` | as hover, larger: sent once | any |
| `structure.ligand` | ViewerConfig | `whole` (a ligand's position is its whole residue) | the precedent's; a glycan's residue, not the whole glycan | `atoms` |
| `structure.channels` | ViewerConfig | `none` | MolViewSpec has no transient emphasis | `state` (channels as color nodes in "color by") |
| `structure.colors.follow` | the action's option | `snapshot` | "ad hoc" at the click | `live` (planned) |
| `structure.colors.tokens` | the action's option | the viewer's current tokens (light or dark) | the structure gets the colors the track draws | `light` (a figure on white) |
| `structure.colors.missing` | the action's option | no node (the structure viewer's color) | Seq\* colors only what it has data for | a token (`seq-neutral`) |
| unobserved positions in messages | fixed | never sent | no atoms, nothing to light | — |
| the message vocabulary | fixed | MolViewSpec component expressions | a public format both sides read; the new mol viewer reads MolViewSpec | Mol\* `Loci` (refused: an import of Mol\*) |
| the fake viewer's `latency` | story arg | 0 | exercises the async path when set | ms |

## Tensions

- **MolViewSpec as the message vocabulary, against the new mol viewer's own words.** The new viewer may select with
  something richer (MolQL, its own loci). Recommended: MolViewSpec's expressions, because they are public, already read
  by Mol\*, and plain JSON. A richer form would be a second target kind added later, never instead. Standing.
- **The port on the behavior, against ports.md's table** (gap 1). Recommended: on the behavior, with ports.md's row
  moved and core's stub removed. Open until decided.
- **A snapshot, against Jalview's live colors** (gap 5). A snapshot goes stale when the theme is edited. Live colors
  rebuild a state on every edit (a structure viewer reloading per keystroke of a theme). Recommended: snapshot, with
  `live` planned. Open.
- **Replace, against keeping the unplaced.** A 3D click clears a gene range that Seq\* had selected. That is the
  grammar's click, and `keep-unplaced` is configurable. Standing.
- **One structure viewer a Seq\* viewer.** The port is one. Two 3D views of one sequence need a host adapter that fans
  out (hover to both, selection from either). Standing until a host asks.
- **The fake, against a real 3D view.** It cannot show occlusion, picking ambiguity (a click on overlapping atoms),
  camera fidelity, or a loader's delays beyond `latency`. What it proves is the messages and the maps, not the look of
  3D. Standing: the new mol viewer's integration re-runs the same headless twins against its adapter.
- **The binding by ref, against a spec that knows its structure.** A spec could name its structure's source (entry,
  model, assembly) so that a host loads it. Recommended: the spec names the data (the complex's units and identifiers),
  and the host names refs. Standing.

## Self-check

The claim: **this design covers what a user does with a structure beside a sequence view.**

1. **Assume it holds.** Every syncing path is a message of MolViewSpec selectors through one port: hover, selection and
   focus both ways, and colors as a state the host loads. Every structure position is a unit's position, under the
   copy rule and the unit's numbering.
2. **What it rests on, negated.**
   - (a) *One structure, one viewer.* Negated: a crystal and a model together, two 3D views.
   - (b) *Seq\* already holds the structure's units.* Negated: the user opens another entry in 3D.
   - (c) *Seq\* acts, 3D follows.* Inverse: 3D acts (a 5 Å shell, a chain), Seq\* follows.
   - (d) *Positions match one to one.* Negated: unobserved residues, insertion codes, altlocs, a branched residue
     without `label_seq_id`, CRO over three positions, a chimera on two sequences.
   - (e) *Messages are synchronous and in order.* Negated: a viewer in a frame, latency, a stale hover.
   - (f) *A section shows the unit.* Negated: only P04637 is shown, or only the gene.
   - (g) *Colors are made once.* Negated: the theme is edited, data arrives later, the track is hidden, a second track
     is colored; the inverse, uncoloring.
   - (h) *Focus is one-way and has no history.* Inverse: back.
   - (i) *Small.* Negated: a capsid's 60 copies, a whole-chain selection.
3. **The negated world.**
   - (a) brings refs in targets and bindings per complex; two viewers mean a fan-out adapter.
   - (b) needs the units from the host's `Model`: the seq data ports' job.
   - (c) makes every path two-way, with a target taking any expression and echoes cut by owner and origin.
   - (d) gives the cut of the unobserved, the `auth` numbering with EC-1 to EC-8, branched residues by `auth_seq_id`, and
     atoms collapsing; CRO's three positions to one monomer is the molecules design's gap 6, and inward it lands on all
     three.
   - (e) gives plain-JSON messages, ownership by owner (EC-141), and the fake's `latency`.
   - (f) uses projection onto systems no section shows, outward, and the hull for focus, inward.
   - (g) gives the snapshot with clear, replace on a second track, a pending result, and `live` planned.
   - (h) gives `followJumps`.
   - (i) gives the size gate and once-a-frame selection.
4. **Held both.**
   - Integrated in the Decision: (a) refs, (c), (d) except CRO's grouping, (e), (f), (g) snapshot and clear, (h), (i).
   - Planned with a named lane: (b) with the seq data ports' lane; `live` colors and "show in 3D" from a card with the
     structure lane and track header parts; CRO's grouping with ML2.
   - Standing tensions: two 3D views; the fake's fidelity; replace against keep; MolViewSpec against richer selectors.

## Consequences

The fourth slice builds it. It comes after ML2 (units and their numbering), after the seq data ports' MolViewSpec
writer, and beside the track header parts lane, which draws the action.

- **`seqx-core`**:
  - `structure/`: the translation (outward and inward, the size gate, `seq.structure.*` issues), `colorRuns`, and the
    color part of the MolViewSpec state (through the seq data ports' writer).
  - `context/`: `SeqPorts.structure` removed (gap 1).
  - Tests:
    - property tests: both round trips; adjacent runs merged; the gate monotone in size; a target's order irrelevant;
    - cases: EC-25 (B lights only B and P04637); EC-48 and EC-65 (the reference and a merged block reach every copy);
      EC-118 (unobserved never sent); EC-1 to EC-8 through `auth` with insertion codes; a branched residue; an atom
      collapsing; an unplaceable selector reported once;
    - a quick representative of the hover target: one message built for 2OCJ within a bound, and no message for an
      unchanged image;
    - heavy: `maxUnits`-scale hover on a molecule gated to chains within a bound.
  - `_spec`: a new Structures section and Context (ports).
- **`seqx-viewer`**:
  - the `structure` behavior (paths, echoes, coalescing, `followJumps`);
  - `HoverModel.setLocation(loc, owner)`;
  - the announcement of a selection from the structure;
  - commands `structure.focus-selection` and `structure.colors.clear`;
  - the config items above.
  - Tests: each path both ways on a recording port; echoes never sent back; a stale `onHover(null)` leaves the pointer's
    hover; a drag sends at most one message a frame; focus and the animation port's duration together on a fake clock.
  - `_spec`: HoverModel, Behaviors, Ports (Structure), Commands, Configuration.
- **`seqx-examples`**: 2OCJ recorded (its `seq` units and the fake viewer's chain table, from one mmCIF by one record
  script); views for p53 and the Complex with a bound structure; "color by" fixtures.
- **Storybook**:
  - the fake structure viewer;
  - stories, each with its headless twin (both viewers on a fake clock, the fake's `latency` set in one): Hover both
    ways, Selection both ways, Focus and back, Color by this track (the conservation track on 2OCJ, the confidence
    track on 1FIN's model), Chain maps (copies merged and split, author numbering, unobserved gaps, the construct
    offset);
  - a `Seq* (next)/Design/Structure syncing` page.
- **Drawing against the data**: for each bound unit's observed residue, the fake viewer's color equals `colorRuns` at
  its preimage. `colorRuns` equals the track's drawn colors at a letters zoom over each tile edge. Every lit cell in the
  fake viewer is in the image of Seq\*'s hover, and the reverse.
- **Lanes** (named for the fourth slice's plan):
  - ST1 (strong): core translation and `colorRuns`;
  - ST2 (strong, after ST1): the behavior and the hover owner;
  - ST3 (after ST2 and the track header parts lane): the fake viewer, the stories and twins, the Design page.
- **Designs that change when agreed**:
  - [ports.md](ports.md): the structure row moved to the behavior;
  - [maps.md](maps.md): the `structure` kind dropped, the structure system read as a unit;
  - [viewer-composition.md](viewer-composition.md): `Ports.structure` removed, and "Inside a Mol\* plugin" naming the
    behavior's port.

## Gaps for the author

1. Move the structure port onto the `structure` behavior (as the hover card's port), and remove core's
   `SeqPorts.structure`? Recommended **yes**.
2. Have the messages speak MolViewSpec component expressions, declared in Seq\* as plain types with no
   `@molstar/mvs-*` import (validated against `mvs-builder` in the Storybook's tests)? Recommended **yes**.
3. Should a hover from the structure show in Seq\* as the hover (an outline, no card) rather than as a highlight
   channel's tint? Recommended **yes**.
4. Should a selection made in the structure reveal it in Seq\* (pan, zoom kept) rather than focus it, as the precedent
   did? Recommended **yes**.
5. Should "color by this track" be a snapshot at the click (cleared by clicking again), with live coloring a later
   option? Recommended **yes**.
6. Should the fake viewer read its own chain table recorded from the mmCIF, not Seq\*'s nodes? Recommended **yes**.
7. Drop `structure` as a map kind, so that a reference maps to the entity's molecule by `offset`, `segments` or
   `alignment`, with status from `observed` and the letters? Recommended **yes** (maps.md changes).

## Decided by the author (2026-11-19)

1. Yes ("rest as recommended").
2. **Changed**: "this should be provided by the application layer anyway": Seq\* declares no MolViewSpec expressions,
   not even as plain types; the structure port speaks Seq\*'s own locations (systems, ranges, later units), and the
   host's adapter translates them to MVS and back. The fake viewer is a host: its adapter does the translation.
3. Yes ("the hover cards will usually be taken over by an app anyway and shown based on context").
4. Yes. 5. Yes. 6. Yes.
7. **No**: "we should still support the structure maps in the spec". The `structure` map kind stays. And "can we have
   more than one map per entry?": several maps from one entry (one per chain or unit) already resolve as paths; two
   different maps between the same pair of systems (SIFTS and an alignment) are ambiguous today (EC-22–27): the
   revision designs a way to name which map a section or link goes through, with a recommendation.
