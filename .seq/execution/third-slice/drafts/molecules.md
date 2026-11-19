# Molecules: small molecules, glycans, conjugates, modified residues, complexes

## Intent

The author, 2026-11-19: "i also wonder we should have a better primitive than `(sequence atp :alphabet protein
:letters X :name ATP)` for this? because these will be for example antibody-drug complexes, and in general, these could
be pointers to smiles, carbohydrate polymers (branched entity in mmCIF), etc. and all these need to be mappable between
each other"; "do we support non-standard residues in our "sequence" nodes? if not, should include this"; "chips is a
prolly good way to do these. contacts should be optional and configurable by the app. remember we are building a
library, not a "one off viewer"."

Today (language 0.4.0) a `sequence` is one letter a position from a fixed alphabet: a ligand is a one-letter protein,
MSE is `M` or `X`, 1EHZ's PSU is `U`, a glycan has no form at all, and nothing says where a drug or a glycan is
attached. This design gives Seq\* its primitives for what is not a linear polymer of letters, a position's chemical
component inside one, the complex the layouts design waits for (layouts.md, Blocks: "a section `of` a complex ... when
the complex kind exists"), and the maps and links between them. Seq\* holds **pointers and labels, never chemistry**: a
component code, a SMILES, a GlyTouCan id are strings it matches and shows; parsing, depicting and comparing chemistry is
a host's.

Built in the fourth slice, in its one language bump (0.4.0 → 0.5.0, pre-release: no migration). Meanwhile the third
slice's DG and DC write ligands as one-position sequences drawn as chips (LB), and DT writes 1EHZ's modified nucleotides
as their parents' letters, the components kept in the fixture; each changes once, in the fourth slice.

## Precedents

- **The previous implementation**: two molecule kinds (`polymer`, `ligand`); a ligand indexed by heavy atoms with their
  elements, hydrogens rejected (EC-64); `chem_comp_ids` a dense column a position; `ligand_form: position | atoms` on a
  section; ligand chips in text (EC-99); default blocks polymers then ligands (EC-101); units, complexes and derived copy
  groups (EC-35, EC-65). Taken: the kinds kept apart, units and complexes, derived copies, the forms, the chips. Not
  taken: a ligand as atoms only (a glycan has residues), the dense component column (sparse here: most positions are
  standard), the form on the section (on the block here: one section shows ATP whole beside a glycan's residues).
- **Mol\*'s model** (`@molstar/model`): `entity.type` (`polymer`, `non-polymer`, `branched`, `water`); `Sequence` keeps
  `code`, `compId` and a `label` a position, the label the component id where the code is `X`, a microheterogeneous
  position `(A|B)`; `StructureSequence` by entity; carbohydrates (SNFG) computed from branched entities. Taken: the
  entity classes, a component a position beside its letter, microheterogeneity as alternatives at one position. Not
  taken: `X` for a modified residue with a parent (the parent's letter here, so search and alignment keep working).
- **MolViewSpec**: component selectors `polymer`, `protein`, `nucleic`, `branched`, `ligand`, `ion`, `water`; component
  expressions by `label_asym_id`, `label_seq_id`, `label_comp_id`, `instance_id`. Taken: the class vocabulary, word for
  word, so a host turns a Seq\* selection into an MVS component (the fourth slice's "color by this track" story) with no
  table between them; a unit's `chain` and `operator` as MVS's `label_asym_id` and `instance_id`.
- **The mmCIF dictionary**: `entity_poly_seq` (`num`, `mon_id`, `hetero`), `entity_poly.pdbx_seq_one_letter_code_can`
  (the parents' letters), `pdbx_entity_nonpoly`, `pdbx_entity_branch` (`type`), `pdbx_entity_branch_list` (`num`,
  `comp_id`, `hetero`), `pdbx_entity_branch_link` (`entity_branch_list_num_1/2`, `atom_id_1/2`, leaving atoms,
  `value_order`), `pdbx_branch_scheme`, `struct_conn` (`covale`, `disulf`, `metalc`), `pdbx_struct_mod_residue`
  (`parent_comp_id`), `chem_comp` (`id`, `name`, `type`). The CCD's `chem_comp` adds `one_letter_code` and
  `mon_nstd_parent_comp_id`, **a list**: a component may stand for several parents (GFP's chromophore CRO for T, Y, G).
  Entries' `chem_comp` (io's mmCIF schema) carries no parent; io's schema has no `pdbx_poly_seq_scheme`. Taken: every
  table above as what the mmCIF runtime reads; positions of a branched entity in `num` order.
- **HELM** (Pistoia's notation for complex biomolecules: peptides and oligos with non-natural monomers, ADCs). Monomers
  with a natural analog, `CHEM` polymers given by SMILES, connections between polymers at attachment points
  (`PEPTIDE1,CHEM1,5:R3-1:R1`), `BLOB` for what has no structure. Taken: the parent letter (HELM's natural analog), a
  molecule by a pointer, a connection as a link between two positions with the atoms named. Not taken: R-group
  attachment points (the structure's atom names instead), a monomer library inside the document (a host's).
- **UniProt and ProtVista**: PTM features `Modified residue`, `Glycosylation`, `Lipidation`, `Disulfide bond`,
  `Cross-link`; a disulfide a feature with two ends. Taken: links as two-ended items, sites as marks at positions on
  the polymer. Not taken: the site without the molecule (here the link names both).
- **RCSB's and PDBe's sequence views**: binding-site, modified-residue and glycosylation tracks along the chain
  (RCSB's 1D feature viewer, PDBe's ProtVista); Mol\*'s sequence panel lists polymers, then branched entities, then
  ligands, each residue labeled by its code; glycans as 2D SNFG diagrams on the entry pages (*to verify the current
  forms when the examples are recorded*). Taken: polymers then branched then ligands (EC-101), a residue of a molecule
  labeled by its code, binding sites as a layer the host adds (decided: contacts optional). Not taken: SNFG drawing in
  the core (it needs a symbol table by component: a host's registry entry, later).
- **SnapGene and Benchling**: features on a construct; Benchling's small molecules by SMILES and its monomer library
  with natural analogs for modified oligos. Taken: a construct is a complex without a structure (units and links over
  molecules); SMILES as an identifier.
- **This codebase**: [maps.md](maps.md) (one map interface; structure systems as chain instances),
  [layouts.md](layouts.md) (blocks, chips, copies, contacts optional), [addressing.md](addressing.md) and
  [locations-and-queries.md](locations-and-queries.md) (`(chain :id A)`, forms over data-made nodes),
  [dynamic-nodes.md](dynamic-nodes.md) (a source answering with nodes), [examples.md](examples.md) (1FIN, 4HHB, 1EHZ).

### Interactions in the precedents

| In the precedents | Here |
|---|---|
| Mol\*'s sequence panel: hover or click a ligand or glycan residue, the 3D view lights or selects it | planned: the structure port (the fourth slice's structure lane), reading units' `chain` and `operator` |
| RCSB, PDBe-KB: hover a binding site, the ligand's contacting residues light; pick a ligand from a list | built here as data: contacts a `pairs` layer the builder adds (decided optional), a ligand's contacts a ref; the Complex view |
| UniProt: hover a disulfide, both ends light | built: links are pairs, hovered end to end |
| Mol\*: show or hide water and ions in the sequence panel | built: the section's `classes` by patch |
| the previous implementation: a ligand as one position or its atoms | built: `form` on a block, by patch |
| find the modified residues of a chain | built: `(component :standard false)`, selected by the query command |
| hover a modified residue: its component and name | built: hover text |
| see a glycan's tree, a compound's 2D drawing | refused in the core: chemistry is the host's (it receives the hovered molecule's identifiers); a glycan's tree layer with a host's SNFG table planned later (Tensions) |
| edit a monomer, draw a connection (HELM editors) | refused: data is not edited in the viewer (ux-foundations section 8) |

## Decision

### A position's chemical component (non-standard residues)

- **A `sequence` keeps its letters and gains `components`**: a sparse table (`position`, `id`; optional `length`, 1
  when absent) naming the chemical component at the positions whose component is not the standard one of their letter.
  A position not listed has its letter's standard component (`A` in protein is ALA, in DNA DA, in RNA A: a table in
  `@molstar/core` beside the one-letter codes, base-library.md). Why sparse: a protein of 300 residues has a handful of
  modified ones; a dense column costs every sequence what one MSE needs.
- **The letter is the parent's letter**: MSE is `M`, SEP `S`, TPO `T`, PTR `Y`, PSU `U`, 2MG and 7MG `G`, OMC and 5MC
  `C`, H2U and 5MU `U`, 1MA `A`; a component with no parent is `X` (protein) or `N` (nucleic). Why: drawing, color
  schemes, search, alignment, codon status and translation read letters; the component adds what the position is and
  changes no arithmetic, so 1EHZ aligns to the Rfam seed with no special case.
- **One component over several positions** (`length` > 1): a component standing for several parents (CRO for T, Y, G)
  covers that many positions, the parents' letters, so the sequence is the canonical one (UniProt's, the aligner's).
  A structure's monomer over three positions is the structure map's to group (Tensions).
- **Several components at one position** (microheterogeneity, `hetero`): several rows with one `position`, in the
  source's order; the first is the position's component, the others its alternatives, shown on hover and marked like a
  modified residue. Its letter is the first's parent.
- **Drawn**: a `letters` layer's `modified` param (`mark`, the default, or `none`) puts a mark on a position whose
  component is not standard (the preset's `letters.modifiedMark`: a dot at the cell's top right, in the letter's ink;
  a run of one component length > 1 one bar over its cells). Where a cell has room for three letters (`code: three`,
  or `auto` with room, EC-117) the position shows its component's code (`MSE`) instead of the parent's (`Met`). Why a
  mark by default: a PSU drawn as a bare `U` says something false by omission, and the writer turns it off.
- **Hover** shows the letter, the component code and its name: `U · PSU (pseudouridine)`, alternatives after a bar.
- **Names**: a component's name comes from the spec's `component` node, else the host's `componentNames` registry
  entry (a CCD lookup Mol\*'s integration supplies), else the code alone. The core ships no CCD (its size, and names are
  labels a host already holds).

### Molecules that are not linear polymers

- **A new data kind, `molecule`**: a small molecule, an ion, a glycan, a drug, a linker-payload. Its positions are its
  **residues** (components) in order: `components` (a list of component ids: `[ATP]`, `[NAG NAG BMA MAN MAN]`); its
  **links** between residues (a table: `a`, `b` positions, `a_atom`, `b_atom`, `order`; a glycosidic bond
  `NAG 2 C1 – NAG 1 O4`); its **atoms**, optional (a table: `residue`, non-decreasing, `name`, `element`: heavy atoms
  only, EC-64); its `class` (`branched`, `ligand`, `ion`, `water`, or any word; absent: `ligand`); its `name`; its
  `identifiers` (a table: `namespace`, `value`: `ccd:ATP`, `smiles:…`, `chembl:…`, `glytoucan:…`, `prd:…`,
  `inchikey:…`); inline or from a `source`. Why a kind and not `sequence`: a sequence's arithmetic (alphabets, letter
  schemes, codon maps, aligners) does not apply, and folding both into one kind makes every reader of letters check a
  flag (the previous implementation's `svs.ligand.unsupported` in the aligners and codon maps).
- **A molecule with no pointer is valid** (an unknown ligand, `UNL`; a payload known by name only): a name, one
  residue, no identifiers.
- **Three systems a molecule**: `/data/glycan` (its residues), `/data/glycan/(whole)` (one position, the molecule), and
  `/data/glycan/(atoms)` when it has atoms. A single-component ligand's residues system is one position, as today.
- **Hydrogens are refused** in `atoms` (`seq.lang.molecule.hydrogen`). Why: they are often absent or model-dependent,
  so a molecule's length would change between two entries of one ligand, and per-atom data would not line up.
- **Components as data**: an optional `components` container under `data` holds `component` nodes (id: the CCD code by
  convention; `name`, `parent` (a list), `class`: the CCD's `type` word, `identifiers`). A code needs no node: the node
  adds a name, parents and identifiers. A writer's own monomer (not in the CCD) is a node with a SMILES.

### Pointers match by string; equivalence is the host's

- Two molecules are the same compound when an identifier of one equals an identifier of the other, namespace and value
  exactly (`(molecule :identifier "ccd:ATP")`). Seq\* never canonicalizes a SMILES or compares structures. A host
  registers an `identifierResolver` (registry entry) that adds equivalents (an InChIKey for a SMILES, a CCD code for a
  ChEMBL id); none by default. Why: matching chemistry is chemistry.

### Complexes and units

- **A new data kind, `complex`**: units shown together, in order, with the links between them. Its `unit` children
  each instantiate a sequence or a molecule: `of`, `range` (on the molecule; absent: all), `observed` (ranges of what
  was modeled; absent: all; unobserved gray, EC-118), `chain` (the structure's `label_asym_id`), `auth_chain`,
  `operator` (MVS's `instance_id`: `1_555`, `ASM-X0-1`), `label`. A unit is a system: its molecule's positions over its
  range, never renumbered (EC-13: its range never wraps). A complex adds no arithmetic (EC-35).
- **Entities are molecules**: an mmCIF entity is a `sequence` (polymer) or a `molecule` (non-polymer, branched, water);
  its chains are units. A **construct** (an ADC design, a fusion with its tags) is a complex with no structure behind
  it: the same units and links.
- **Copies are derived** (EC-65): units of one complex with equal `of` and `range` are a copy group; never declared.
  The copy rule holds (EC-25): no path climbs from a unit to its molecule and down to a sibling.
- **Links**: the complex's `links` table: `a`, `b` (unit ids), `a_at`, `b_at` (positions on them), `a_atom`, `b_atom`,
  `kind` (`covalent`, `disulfide`, `glycosidic`, `metal`, or any word), `occupancy` (absent: 1), `label`. Both ends on
  one unit is allowed (an intra-chain disulfide, a staple); a site carrying one of several glycans (glycoforms) or a
  payload on one of many lysines (stochastic conjugation) is several rows with their occupancies. Inside a molecule,
  links are the molecule's (above).
- **A complex's links are a `pairs` entry** (`/data/adc/(links)`): its rows' ends name their units (the pairs kind
  generalized: `a` and `b` may name a complex, rows then carry `a_unit` and `b_unit`). So a `pairs` layer draws them as
  arcs (between blocks of one section, layouts.md's section-scope pairs), a `ranges` layer over them draws each end as
  a mark at its site labeled by the other end (as a reader without `pairs` draws a pair's ends today; a `ranges`
  layer taking `pairs` data on purpose is part of this change), and hover lights both ends.
  Contacts a host computes are a `pairs` node of the same shape (decided optional; never drawn by default).
- **Links and contacts are never maps**: projection does not cross them. Why: a map carries positions by arithmetic;
  a link relates two places without one, and projecting across it would make a glycan's selection select its Asn.
  What a user expects instead (the site lit when the glycan is hovered) is a behavior, below.
- **Annotations on any positions**: `ranges`, `values`, `pairs`, `grid` and `map` take a sequence, a molecule (any of
  its three systems), or a unit; a target that is a molecule reaches every unit over it (EC-48). A per-ligand score
  (AlphaFill's) is a field of the unit's row in a `ranges` over the molecule, or a `values` over `(whole)`.
- **Bounded**: a complex from a source is refused past `SpecBounds.maxUnits` (65,536, a large capsid's assembly) with
  an issue naming an assembly subset; merged copies cost one block (EC-67).

### Maps between the forms: many-to-one

- **A `group` map**: runs of the source each onto one target position, held as offsets (the target's position `r`
  owns source `[off[r], off[r+1])`, contiguous by construction). Atoms to residues (`residue` non-decreasing), residues
  to `(whole)`, a component over `length` positions to the structure's one monomer. Forward, a range goes to the
  positions it touches (as nucleotides to residues, EC-18); back, a position to its run. Implicit for a molecule's own
  forms (from its tables), inline for a writer's (`method: group`, `spans`).
- **One span step a path**: the codon step (span 3) generalizes to a span step (codon or group); at most one a path
  (EC-19), consecutive group steps composing into one (a grouping of groupings is a grouping). Why: paths stay one
  linear map with a scale, as the codon step made them; nothing in this design needs a codon and a group on one path (a
  polymer's atoms are the structure viewer's, not a form here).
- **Follow across a group map is refused** (EC-113: follow goes at scale 1, 3 or 1/3): a block's form is not followed.

### On an axis: the block's `form`

- **`block.form`**: `whole` (one position: the chip of layouts.md, the molecule's label and `×N`), `residues` (a
  position a component, each labeled by its code; a glycan's links drawn when a `pairs` layer over the molecule's links
  is added), or `atoms` (a position a heavy atom, labeled by element symbol, hovered by atom name). Absent: `whole` for
  a molecule, `residues` for a sequence or a polymer unit. Why `whole`: a molecule's residues are a sliver at the fit
  view (the reason chips were decided), and one rule draws ATP, a heme, a glycan and a payload alike. `whole` on a short
  polymer (a peptide ligand) is allowed. Switched by patch (one undo step), as `copies`.
- **A chip's label**: the molecule's `name`, else its one component's code, else its id; truncated at the chip's
  maximum width (LB's rule), whole on hover. The mmCIF runtime writes a branched molecule's `name` as its composition
  in list order (`NAG2 BMA MAN2`): the entity's description is a paragraph.
- **Hit, hover, selection** read a chip as its one position of `(whole)`; through the group maps it stands for every
  residue and atom (the previous implementation's "a ligand shown as one position is the locus of all its atoms").

### A section of a complex

- **`section.of` takes a complex**: a block a unit (layouts.md, Blocks), copies merged by `copies` (EC-65), ordered by
  the section's `classes` (default `[polymer, branched, ligand, ion]`; water and any class not listed left out),
  polymers in the complex's order, then each class in order (EC-101). Why left out: hundreds of waters say nothing on a
  sequence axis; a writer lists `water` to show them.
- **Explicit `block` children** still win: a writer picks units, forms and order.
- **Follow**: a section of units on several systems neither follows nor is followed (`seq.render.follow.blocks`,
  unchanged).

### Behaviors (the host's to keep or drop)

- **`links.hover`** (a `BehaviorProvider` in the default spec): hovering a molecule (a chip, a residue, an atom)
  lights the polymer positions its links reach, as the `linked` emphasis; hovering a linked site lights the molecule.
  A host leaves it out of its spec. Contacts get no such behavior: they are a layer the host adds (decided).
- **The hover carries identifiers**: the hovered element's molecule and its `identifiers` are in the hover record, so a
  host's card draws a structure from a SMILES or a CCD image. No depiction port: one consumer, the host's card, already
  has what it needs.

### Queries

| Form | Prefix | Result |
|---|---|---|
| `(component :id MSE)`, `(component :id [PSU 2MG])` | sequences, molecules, units | Location: the positions with that component (alternatives included) |
| `(component :standard false)` | sequences, units | Location: every non-standard position |
| `(component :class saccharide)` | as above | Location, through the `component` nodes or the host's names entry; a code with neither is no match |
| `(molecule :identifier "ccd:ATP")`, `(molecule :class ligand)` | `data`, a complex | Nodes: molecules (with the host's resolver, its equivalents) |
| `(entity :id 3)` | a complex from a source | Nodes: the molecule its source calls entity 3 |
| `(chain :id A :operator 1_555)`, `(chain :auth B)` | a complex | Nodes: units (the form exists; `:operator` and `:auth` added) |
| `(atom :name O1G)` | a molecule's `(atoms)`, a unit | Location: one atom |
| `(whole)`, `(atoms)` | a molecule, a unit | Nodes: its derived systems |
| `(links)` | a complex | Nodes: its links' pairs entry |

### Reading

- **mmCIF** through a runtime the application registers (`@molstar/seqx-core/runtimes/mmcif`, importing
  `@molstar/io`'s CIF reader and mmCIF schema on the allowlist, packages.md; never in the default registry). A source
  of format `mmcif` answers with nodes ([dynamic-nodes.md](dynamic-nodes.md)): a sequence per polymer entity (letters
  from `pdbx_seq_one_letter_code_can`, components from `entity_poly_seq` where the code is not its letter's standard,
  `hetero` rows as alternatives, `pdbx_struct_mod_residue` corroborating), a molecule per non-polymer
  (`pdbx_entity_nonpoly`) and branched entity (`pdbx_entity_branch_list` in `num` order, `pdbx_entity_branch_link` as
  its links, leaving atoms dropped), a complex of units (`struct_asym`, observed from `atom_site`, `operator` from the
  assembly asked), its links from `struct_conn` (`covale`, `disulf`, `metalc`), names from `chem_comp`. Its options:
  `assembly` (absent: the deposited model), `water` (absent: not read). Units from a source are reached by forms
  (`(chain :id A)`), inline units by id. A molecule's atoms come from the CCD's `chem_comp_atom` when a CCD source is
  given, else the first unit's `atom_site` (a unit's `observed` marks its modeled atoms).
- **The `seq` format** carries the same nodes and tables (the fixtures' form, converted at recording).

### Language, store, resolved model

- **Language (0.5.0)**: `sequence.components`; the kinds `molecule`, `component` (under a `components` container),
  `complex`, `unit`; `map.method: group` with `spans`; `pairs.a`/`b` naming a complex with `a_unit`, `b_unit`; the
  kinds a data reference takes widened to molecules and units (`ranges.target`, `values.target`, `pairs.a`/`b`,
  `grid.x`/`y`, `map.from`/`to`, `block.of`, `section.of`, which takes a complex too); `block.form`;
  `section.classes`; the `letters` layer's `modified`. Checks (`seq.lang.*`): a component position past the end or
  out of order, a component `length` crossing another, a molecule's `atoms.residue` decreasing or past its residues, a
  hydrogen, a link's position past its unit's range, a unit `range` wrapping, `form: atoms` on a molecule with no atoms.
- **Store**: molecules, components, complexes and links are held whole (as pairs and grids, `resolve/held.ts`): an
  entry's tables in typed arrays, keyed by address and hash; a complex from a source is one request, its nodes the
  overlay's.
- **Resolved model**: systems for units and a molecule's three forms; group maps implicit from the tables; copy groups;
  the links' pairs entry; components per system (a sorted position index, the alternatives); the query model's
  `chain`, `node` and the new forms answering.

## The interaction matrix's cells

Rows and columns of ux-foundations section 8; "as an element" means the grammar's cell, unchanged.

| | select | extend | add / remove | focus | back | zoom in / out | pan | hide | edit |
|---|---|---|---|---|---|---|---|---|---|
| **a section of a complex** | — | — | — | — | — | as its layout | as its layout | no: its tracks | `classes`, `copies`, `arrange` by patch (`section.set`, previewed) |
| **a block: a chip** (a molecule `whole`) | click: its one position (every residue and atom, every merged copy) | Shift+click: each block's part between the anchor and it, in reading order (layouts.md) | Cmd/Ctrl+click | double-click: revealed; a sibling section showing its residues or atoms zooms to them | double-click again, Alt+Left | no: a chip is fixed in px (layouts.md) | no: the window skips chips | no: `classes` or the blocks list by patch | `block.form` by patch (`whole`, `residues`, `atoms`), no default key, in the context menu (as `copies`) |
| **an element**: a modified residue | as an element | as an element | as an element | as an element | as an element | — | — | no: `modified: none` on the layer | no: data |
| **an element**: a molecule's residue or atom | as an element; projected onto `(whole)` and the other form | as an element | as an element | as an element | as an element | — | — | no | no: data |
| **an element**: a link (an arc, or a site's mark) | click: the link item; both ends lit | Shift+click | Cmd/Ctrl+click | double-click: its extent, two parts on two systems, each section revealing its own | Alt+Left | — | — | its layer's items by kind (`layer.filter`: hide disulfides) | no: data |
| **a range** over blocks with chips | drag in reading order across a chip (the Complex's text figure) | Shift+drag | Cmd/Ctrl+drag, Alt+drag | Enter | Alt+Left | — | — | — | — |
| **the selection** | `(component :standard false)` and the other forms through the query command | — | — | Enter | Alt+Left | — | — | — | bookmark it (a glycosylation site, a ligand's pocket) |
| **the view** | — | — | — | 0 fits; chips keep their px | — | — | — | — | — |
| **a track or layer** | — | — | — | — | — | — | — | `track.hide` (a sites track, the contacts) | `modified` by patch; a links layer added by the host or builder |
| **a bookmark or link** | — | — | `refs.add` (a ligand's contacts, a site) | runs its action: the molecule and its sites highlighted | Alt+Left | — | — | — | its label |
| **a follower** | no: a section of several systems is never a follower (`seq.render.follow.blocks`) | | | | | | | | |
| **a navigator** | — | — | — | — | — | chips skipped, as by the window (layouts.md) | as above | — | — |

## Configuration

| Item | Whose | Default | Why | Alternatives |
|---|---|---|---|---|
| `block.form` | spec param | `whole` for a molecule, `residues` for a sequence or unit | a molecule reads at the fit view; one rule for every non-polymer | `residues`, `atoms` |
| `section.classes` | spec param | `[polymer, branched, ligand, ion]` | Mol\*'s and MVS's order (EC-101); waters say nothing on an axis | any list; `water` added |
| `copies` | spec param | `merge` (layouts.md) | unchanged | `each` (poses: the builder's for a pose set) |
| `modified` (`letters`) | spec param | `mark` | a modified residue drawn as its parent alone hides it | `none` |
| `letters.modifiedMark` | `RenderConfig` preset | a dot at the cell's top right in the letter's ink; a bar over a long component | reads at 8 px, never covers the letter | `underline` |
| `code` (`letters`) | spec param | `auto` (exists): the component's code where three letters have room | a cell with room says what the residue is | `one`, `three` |
| `class` (molecule) | spec param | `ligand` | the commonest non-polymer | `branched`, `ion`, `water`, any word |
| `componentNames` | registry entry | none (the code alone) | the core ships no CCD | a host's CCD lookup (Mol\*'s) |
| `identifierResolver` | registry entry | none (exact match only) | equivalence is chemistry | a host's (InChIKey, ChEMBL to CCD) |
| `links.hover` | `BehaviorProvider` | in the default spec | a conjugate read with its site | left out by a host |
| contacts | a layer the host or builder adds | none (decided) | "contacts should be optional and configurable by the app" | — |
| the mmCIF runtime | registry entry | not registered | the core reads no structure file by default; its size | registered by the application |
| `assembly`, `water` (mmCIF source) | source params | the deposited model; waters not read | what the entry is; waters as above | an assembly id; `true` |
| `SpecBounds.maxUnits` | `SpecBounds` | 65,536 | a large capsid's assembly | larger, at its memory |
| hydrogens in `atoms` | fixed | refused | lengths and per-atom data line up across entries (EC-64) | — |
| links and contacts as maps | fixed | never | a glycan's selection never selects its site | — |

## Tensions

- **A new kind against a generalized `sequence`.** One kind with an alphabet `components` would keep the language
  smaller; every letter reader would then check it. Recommended: the new kind (Decision). Open only if the author
  prefers one kind.
- **One position, one letter.** CRO stands for three parents; recommended: the canonical letters, the component over
  three positions (`length` 3), the structure's one monomer grouped by the structure map. The alternative (the
  structure's monomers as positions, CRO one `X`) keeps the structure's numbering and loses UniProt's. The structure
  maps' design (the fourth slice) settles the map side; this stays open until it does.
- **A unit and maps.md's structure system.** maps.md calls a chain instance (entry, model, operator, chain) a
  structure system; here it is a unit of a complex read from a structure. Recommended: one noun, the unit; a structure
  map relates a reference sequence (UniProt) to the entity's molecule, and SIFTS's per-chain segments become per-unit
  `observed` and numbering. The structure lane decides.
- **A conjugate far from its site.** Chips follow the polymers, so a glycan on N297 is drawn after the chains. The
  sites (a ranges layer over the links) and `links.hover` carry the attachment; a chip drawn at its site (a pin above
  the residue) would be a placement of its own in every layout. Open: a later design, if a view needs it.
- **A glycan is a tree, an axis is a line.** `residues` lays a glycan in list order, its links as arcs; the tree
  (SNFG) needs a symbol a component, a host's table (Mol\* has one). Planned later: a `tree` drawing of a molecule's
  residues with a host's `symbols` registry entry; not in the fourth slice unless a view asks.
- **Poses and copies.** Ten docking poses of one ligand are a copy group by the rule and merge to `×10` by default,
  hiding that they differ. The builder for a pose set writes `copies: each`; the rule is not bent (copies are derived,
  never declared). Standing.
- **Stochastic conjugation reads as many links.** An average DAR of 3.5 over 40 lysines is 40 rows with occupancies;
  the chip says `×N` of units, not the average. A label field carries it. Standing until an ADC view is used.
- **Atoms without coordinates.** An atoms form with no bonds or positions is a list; it earns its place by per-atom data
  (fit, occupancy, B-factor) and by matching a structure viewer's atoms. Kept; a map pairing two molecules' atoms by
  name (a pose onto its component) waits for a view that needs it.

## Self-check

The claim: "this design covers what a user does with molecules: small molecules, glycans, conjugates, non-standard
residues and complexes".

1. **Assume it holds.** Everything is a system of positions: a molecule has three (whole, residues, atoms) joined by
   group maps; a modified residue is a component on a position whose letter is its parent; a complex instantiates
   molecules as units; links and contacts are pairs; identity is a pointer matched as a string.
2. **What it rests on, negated.** (a) *One position, one component, one letter.* Negated: a component over three
   letters (CRO), two components at one position (microheterogeneity). (b) *A molecule is named by a pointer.* Negated:
   no pointer (UNL, a payload by name), and one compound under several pointers (CCD, SMILES, ChEMBL): "mappable
   between each other". (c) *A link has one site and one partner.* Negated: one site, several possible glycans
   (glycoforms); one payload, many possible lysines; a link between two polymers (a disulfide), within one (a staple),
   a molecule linked twice (a bivalent linker). (d) *A complex is a structure's chains.* Negated: a construct with no
   structure (an ADC design), a set of poses, an assembly of thousands of units. (e) *Its inverse*: a polymer in a
   molecule's form (a peptide ligand, a cyclic peptide), and a molecule's data on the polymer (a ligand's score where
   its contacts are). (f) *A missing view*: the glycan's tree, the compound's drawing, a table of ligands.
3. **The negated world.** (a) gives the component's `length` and alternatives at one position; (b) gives molecules
   valid without identifiers, identifiers as a table, exact matching and the host's resolver; (c) gives links as rows
   with occupancy, ends on any units of the complex (self-links), several links a molecule; (d) gives units over
   molecules with no structure, `copies: each` for poses by the builder, `maxUnits` and merged copies; (e) gives
   `form: whole` on a polymer, cyclic peptides as `cyclic` sequences (the param exists), annotations on molecules and
   units, and shows that a ligand's score belongs to the unit, not its contacts; (f) surfaces the tree and the drawing
   as the host's, and a ligand table as the Complex's list (a ranges over `(whole)` with scores as fields).
4. **Held both.** Integrated in the Decision: the component's `length`, alternatives, the parent letter; `identifiers`
   and the resolver; links with occupancy and self-links as pairs over a complex's units; constructs as complexes;
   annotations on molecules and units; `form: whole` for any block; the hover carrying identifiers; `maxUnits`.
   Planned with a named lane: the structure map's grouping of a multi-letter monomer and the unit as the structure
   system (the fourth slice's structure lane). Standing tensions: CRO's positions until the structure maps settle; a
   conjugate drawn far from its site; the glycan's tree; poses merged as copies by default; the average DAR.

## Consequences

The fourth slice builds it, after the third slice's DG, DC, DT and LB (chips) and with the structure lane beside it.

- **`seqx-core`**: `lang/` (the kinds, params and checks above; the builder's `molecule`, `complex`, `unit`,
  `components`; issues `seq.lang.molecule.*`, `seq.lang.component.*`, `seq.lang.unit.*`, `seq.lang.link.*`;
  `SeqLanguageVersion` 0.5.0); `resolve/` (the held entries, the three systems, group maps, the span step in
  `findPath`, copy groups, the links' pairs entry, the components index); `view/` (a section of a complex, `classes`,
  `form`, chip labels); `render/` (the `letters` mark and the component code; residue labels of a molecule; element
  symbols on atoms); `query/` (the forms above); `data/` (the `seq` format's molecule, component and complex tables);
  `runtimes/mmcif` (new entry, the `@molstar/io` allowlist entries added). Standard components per letter in
  `@molstar/core` beside the one-letter codes (base-library.md; the author's call).
- **`seqx-viewer`**: hover text with component, alternatives and atom name; the hover record's identifiers;
  `links.hover`; `block.form` in the context menu; the matrix cells above.
- **`seqx-examples`**: DT (1EHZ's components, 11 modified residues, marked and found by query); DC (1FIN as a complex
  from its fixture: units, ATP `whole` and `atoms` joined by the group map, AlphaFill's transplants as molecules with
  scores, contacts a `pairs` layer); DG (4HHB's hemes as molecules); Antibody gains an Fc glycan (an IgG1 entry with its
  N297 glycan, *the entry chosen at recording*, CC0) and an ADC construct (trastuzumab with a payload by SMILES on
  lysines, generated, occupancies), with their views and editing views.
- **Storybook**: a `Seq* (next)/Design/Molecules` page; stories for each of the views above, each with its headless
  twin on a fake clock.
- **Tests**: properties (a group map's forward then back covers what it started from; a modified position's letter is
  its parent's; copies derived equal for any unit order); the mmCIF runtime's conformance suite on 1EHZ, 1FIN, 4HHB and
  the glycoprotein (every entity, chain, branch link and `struct_conn` row read, against counts recorded from the
  file); drawing against the data (every modified position in view has its mark; every chip its label; every link with
  both ends in view its arc and its hit box); a quick representative of the frame target with 200 merged copies and a
  glycan chip (EC-67).
- **`_spec` sections that change**: core's Language (kinds, params), Resolve (Systems, Maps: group and the span step,
  Paths), Resolved view (Blocks: a section of a complex, `form`), Render (Layer types: `letters`), Queries (Forms),
  Data (the seq format), a new Runtimes section (mmCIF); viewer's HoverModel, Behaviors, Commands; both Configuration
  tables.

### The fourth slice's lanes

- **ML** (strong) the language and model: the kinds, checks, builder, the seq format's tables, the store's held
  entries, systems, group maps and the span step, copies, links as pairs; the one bump to 0.5.0, shared with the
  slice's other language changes.
- **MC** (strong, after ML) the mmCIF runtime: the allowlist, the tables read, its conformance suite, the fixtures
  re-recorded from mmCIF; sizes measured (the io reader's cost to an application that registers it).
- **MR** (strong, after ML, beside MC) drawing and interaction: a section of a complex, `form` through LB's chips, the
  letters' mark and code, hover, `links.hover`, `block.form`, the queries, the matrix cells, the Design page.
- **MX** (strong, after MC and MR) the examples changed once: DT, DC, DG to the new kinds; the glycan and ADC views;
  stories, headless twins, drawing against the data.

## Gaps for the author

1. A new `molecule` kind beside `sequence`, rather than one kind with an alphabet of components? Recommend **yes**.
2. A modified residue's letter is its parent's (MSE is `M`), the component beside it? Recommend **yes**.
3. A modified residue marked by default in `letters` (`modified: mark`)? Recommend **yes**.
4. A molecule block `whole` (a chip) by default, a glycan included, its residues and atoms by patch? Recommend **yes**.
5. Water left out of a section of a complex by default (`classes` without `water`)? Recommend **yes**.
6. A component over several positions (CRO as T, Y, G): the canonical letters, rather than the structure's one
   monomer? Recommend **yes**, until the structure maps' design says otherwise.
7. The mmCIF runtime as an opt-in entry of `seqx-core` (`runtimes/mmcif`) on `@molstar/io`'s CIF reader, rather than a
   package of its own? Recommend **yes** (one allowlist entry; never in the default registry).
8. The standard component per letter (A ↔ ALA, DA, A) moved into `@molstar/core` beside the one-letter codes, a
   public addition? Recommend **yes**.
9. This slice's DC row "ATP ... joined by a many-to-one map" (examples.md): 0.4.0 has no such map (`offset`, `codon`
   only). DC writes ATP's two forms as two sequences with no map this slice, hover not crossing until the fourth
   slice? Recommend **yes**, and examples.md's row reworded with DC.
