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

**Words.** A *component* is a chemical component code at a position (the CCD's word: MSE, NAG, ATP). A *compound* is
the description of one (its name, parents, identifiers), a node of the optional `compounds` container. `modified` is a
sequence's table of its non-standard positions; `residues` is a molecule's list of components.

## Precedents

- **The previous implementation**: two molecule kinds (`polymer`, `ligand`); a ligand indexed by heavy atoms with their
  elements, hydrogens rejected (EC-64); `chem_comp_ids` a dense column a position; `ligand_form: position | atoms` on a
  section; ligand chips in text (EC-99); default blocks polymers then ligands (EC-101); units, complexes and derived copy
  groups (EC-35, EC-65). Taken: the kinds kept apart, units and complexes, derived copies, the forms, the chips. Not
  taken: a ligand as atoms only (a glycan has residues), the dense component column (sparse here: most positions are
  standard), the form on the section (on the block here: one section shows ATP whole beside a glycan's residues).
- **Mol\*'s model** (`@molstar/model`): `entity.type` (`polymer`, `non-polymer`, `branched`, `water`); `Sequence` keeps
  `code`, `compId` and a `label` a position, the label the component id where the code is `X`, a microheterogeneous
  position `(A|B)` with every component kept; `StructureSequence` by entity; carbohydrates (SNFG) computed from branched
  entities. Taken: the entity classes, a component a position beside its letter, every alternative kept at one
  position. Not taken: `X` for a modified residue with a parent (the parent's letter here, so search and alignment keep
  working).
- **MolViewSpec**: component selectors `polymer`, `protein`, `nucleic`, `branched`, `ligand`, `ion`, `water`; component
  expressions by `label_asym_id`, `label_seq_id`, `label_comp_id`, `instance_id`. Taken: the class vocabulary, word for
  word, so a host turns a Seq\* selection into an MVS component (the fourth slice's "color by this track" story) with no
  table between them; a unit's `chain` and `operator` as MVS's `label_asym_id` and `instance_id`.
- **The mmCIF dictionary**: `entity_poly_seq` (`num`, `mon_id`, `hetero`), `entity_poly.pdbx_seq_one_letter_code_can`
  (the parents' letters), `pdbx_entity_nonpoly`, `pdbx_entity_branch` (`type`), `pdbx_entity_branch_list` (`num`,
  `comp_id`, `hetero`), `pdbx_entity_branch_link` (`entity_branch_list_num_1/2`, `atom_id_1/2`, leaving atoms,
  `value_order`), `pdbx_entity_branch_descriptor` (WURCS, LINUCS), `pdbx_branch_scheme`, `pdbx_nonpoly_scheme`,
  `struct_conn` (`covale`, `disulf`, `metalc`), `pdbx_struct_mod_residue` (`parent_comp_id`), `pdbx_molecule` (a BIRD
  `prd_id` over several `asym_id`s), `atom_site` (`auth_seq_id`, `pdbx_PDB_ins_code`, `label_alt_id`), `chem_comp`
  (`id`, `name`, `type`). The CCD's `chem_comp` adds `one_letter_code` and `mon_nstd_parent_comp_id`, **a list**: a
  component may stand for several parents (GFP's chromophore CRO for T, Y, G). Entries' `chem_comp` (io's mmCIF schema)
  carries no parent; io's schema has no `pdbx_poly_seq_scheme`. Taken: every table above as what the mmCIF runtime
  reads; positions of a branched entity in `num` order.
- **HELM** (Pistoia's notation for complex biomolecules: peptides and oligos with non-natural monomers, ADCs). Monomers
  with a natural analog, `CHEM` polymers given by SMILES, connections between polymers at attachment points
  (`PEPTIDE1,CHEM1,5:R3-1:R1`), `BLOB` for what has no structure. Taken: the parent letter (HELM's natural analog), a
  molecule by a pointer, a connection as a link between two positions with the atoms named. Not taken: R-group
  attachment points (the structure's atom names instead), a monomer library inside the document (a host's).
- **UniProt and ProtVista**: PTM features `Modified residue`, `Glycosylation`, `Lipidation`, `Disulfide bond`,
  `Cross-link`, known by name and evidence, with no CCD code and no glycan; a disulfide a feature with two ends. Taken:
  links as two-ended items, sites as marks at positions, a modification known by name only as an annotation.
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
  `resolve/paths.ts` (a unit's free `up` and `down` steps to its molecule, under the copy rule: units as systems fit
  what is built), [layouts.md](layouts.md) (blocks, chips, copies, contacts optional), [addressing.md](addressing.md)
  and [locations-and-queries.md](locations-and-queries.md) (`(chain :id A)/(at :label 94A :numbering auth)`),
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
| see which residues are modified, at a glance: Mol\*'s sequence panel prints the component code inline where it has no one-letter code (`PSU`, microheterogeneity as `(MSE\|MET)`), the parent letter otherwise; UniProt's and PDBe's feature viewers (Nightingale) and RCSB's 1D viewer give modifications a feature row of their own, a marker per site, typed by kind | built, each a value of `letters.modified` (below): `mark`, `code` (Mol\*'s), `tint`, `none`; the feature row as a `ranges` layer over `modified`-type items or the sequence's `(component :standard false)` sites |
| address a residue by its author number | built: a unit's `auth` numbering |
| see a glycan's tree, a compound's 2D drawing | refused in the core: chemistry is the host's (it receives the hovered molecule's identifiers); a glycan's tree layer with a host's SNFG table planned later (Tensions) |
| edit a monomer, draw a connection (HELM editors) | refused: data is not edited in the viewer (ux-foundations section 8) |

## Decision

### A position's chemical component (non-standard residues)

- **A `sequence` keeps its letters and gains `modified`**: a sparse table (`position`, `component`; optional `length`,
  1 when absent) naming the component at the positions whose component is not the standard one of their letter. A
  position not listed has its letter's standard component (`A` in protein is ALA, in DNA DA, in RNA A: a table in
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
  source's order, every one kept (as Mol\* keeps them); the first gives the letter, all are shown on hover
  (`A · ALA | GLY`), found by the query and marked like a modified residue (gap 8).
- **Drawn** (the author, 2026-11-19: "we need a visual representation for the modified residues in the tracks. Maybe a
  dot or something? How do people do it? Should be configurable either way"): a `letters` layer's `modified` param
  says how a position whose component is not standard is drawn: `mark` (the default: a dot at the cell's top right in
  the letter's ink, `RenderConfig`'s `letters.modifiedMark` for its size and place; one bar over the cells of a
  component of `length` > 1; below the letters' minimum size, where letters become color strips, a tick above the
  strip, so the mark survives zooming out), `code` (the component's code in place of the letter where the cell has
  room, the dot otherwise: Mol\*'s sequence panel), `tint` (the cell outlined in `seq-modified`, a new token, for
  figures where a dot is too small), or `none`. A modification's color is never the residue scheme's (the scheme
  colors the parent letter); a writer who wants modifications by kind adds the feature row, a `ranges` layer over
  the sites, colored by type like any range. Where a cell has room for three letters (`code: three`, or
  `auto` with room, EC-117) the position shows its component's code (`MSE`) instead of the parent's (`Met`).
  **Compact by default** (the author, 2026-11-19: "I would still like the (default) representation to be compact for
  the modified residues"): no value of `modified` changes a cell's width, the letters' spacing or a line's length; the
  default `mark` adds only the corner dot inside the cell, and `code` prints the code only where the cell already has
  room for it (the dot otherwise), never widening a cell to fit it. Why a
  mark by default: a PSU drawn as a bare `U` says something false by omission, and the writer turns it off.
- **Hover** shows the letter, the component code and its name: `U · PSU (pseudouridine)`.
- **Names**: a component's name comes from the spec's `compound` node, else the host's `componentNames` registry entry
  (a CCD lookup Mol\*'s integration supplies), else the code alone. The core ships no CCD (its size, and names are
  labels a host already holds).
- **A modification known by name only** (UniProt's `Modified residue`, `Glycosylation`, `Lipidation`: no CCD code, no
  glycan) stays an annotation: a `ranges` item whose `type` is `modified` (or `glycosylation`, `lipidation`), its name
  its label. `(component :standard false)` finds those items too, over the annotations targeting the prefix's systems
  (`:annotations`, default true), so "where are the modified residues" has one answer whatever the source. Why not a
  `modified` row without a component: the sequence's table says what a position *is* (a structure's chemistry); a
  database's PTM is a claim with evidence about a position, which is what annotations hold (gap 7).

### Molecules that are not linear polymers

- **A new data kind, `molecule`**: a small molecule, an ion, a glycan, a drug, a linker-payload. Its positions are its
  **residues** (`residues`: a list of component codes, a position each: `[ATP]`, `[NAG NAG BMA MAN MAN]`); its
  **atoms**, optional (a table: `residue`, non-decreasing, `name`, `element`: heavy atoms only, EC-64); its `class`
  (`branched`, `ligand`, `ion`, `water`, or any word; absent: `ligand`); its `name`; its `identifiers` (a table:
  `namespace`, `value`); inline or from a `source`; its links between residues a `pairs` child (Links, below). Why a
  kind and not `sequence`: a sequence's arithmetic (alphabets, letter schemes, codon maps, aligners) does not apply, and
  folding both into one kind makes every reader of letters check a flag (the previous implementation's
  `svs.ligand.unsupported` in the aligners and codon maps).
- **A molecule with no pointer is valid** (an unknown ligand, `UNL`; a payload known by name only): a name, one
  residue, no identifiers.
- **Two systems a molecule**: `/data/glycan` (its residues) and `/data/glycan/(atoms)` when it has atoms. The whole
  molecule is `(all)` of its residues system; "whole" is a drawing (the chip), not a system. A single-component
  ligand's residues system is one position, as today.
- **Hydrogens are refused** in `atoms` (`seq.lang.molecule.hydrogen`). Why: they are often absent or model-dependent,
  so a molecule's length would change between two entries of one ligand, and per-atom data would not line up.
- **Alternate locations**: one atom a name; the atoms form takes the first `label_alt_id` an atom name has (`.` when
  the atom has none), the others dropped and counted in one info issue. Why: a position is one atom; altlocs are the
  structure viewer's.
- **Compounds as data**: an optional `compounds` container under `data` holds `compound` nodes (id: the CCD code by
  convention; `name`, `parents` (a list of codes), `class`: the CCD's `type` word, `identifiers`). A code needs no
  node: the node adds a name, parents and identifiers; a writer's own monomer (not in the CCD) is a node with a SMILES.
  Why a container under `data` and not compounds as data nodes beside sequences: a compound is reference data, not a
  system (it has no positions, nothing is placed on it), and its ids are CCD codes that collide with the writer's ids
  of systems (an `A` compound, an `A` chain); one container keeps them in a namespace of their own, as `refs` keeps refs.

### Identifiers match by string; equivalence is the host's

- Two molecules are the same compound when an identifier of one equals an identifier of the other, namespace and value
  exactly (`(molecule :identifier "ccd:ATP")`). Seq\* never canonicalizes a SMILES or compares structures. A host
  registers an `identifierResolver` (registry entry) that adds equivalents (an InChIKey for a SMILES, a CCD code for a
  ChEMBL id); none by default. Why: matching chemistry is chemistry.
- **Namespaces are open words** (gap 9), with a documented vocabulary the readers write: `ccd`, `prd`, `smiles`,
  `inchi`, `inchikey`, `chembl`, `pubchem`, `glytoucan`, `wurcs`, `linucs`. Why open: a host's registry (a company's
  compound ids) is as valid a pointer as a public one.

### Complexes and units

- **A new data kind, `complex`**: units shown together, in order, with the links between them. Its `unit` children
  each instantiate a sequence or a molecule: `of`, `range` (on the molecule; absent: all), `observed` (ranges of what
  was modeled; absent: all; unobserved gray, EC-118), `chain` (the structure's `label_asym_id`), `auth_chain`,
  `operator` (MVS's `instance_id`: `1_555`, `ASM-X0-1`), `label`. A unit is a system: its molecule's positions over its
  range, never renumbered (EC-13: its range never wraps). A complex adds no arithmetic (EC-35). Units are the paths'
  free `up` and `down` steps that `resolve/paths.ts` already takes under the copy rule.
- **A unit's `auth` numbering**: a `numbering` child of the unit (labels a position: the author's number and insertion
  code), read by the mmCIF runtime from `atom_site` (`auth_seq_id`, `pdbx_PDB_ins_code`), `pdbx_branch_scheme` for a
  glycan's units and `pdbx_nonpoly_scheme` for a ligand's, so `(chain :id A)/(at :label 94A :numbering auth)`
  ([addressing.md](addressing.md)) resolves. An unobserved residue has no author number (EC-1 to EC-8 at the edge). If
  the language has no `numbering` kind by the fourth slice, ML1 adds it with the `labels` method only. (Superseded: the
  author, 2026-11-19, "Add it now (0.4.0)": the third slice's NB built the kind with `labels`; ML1 reuses it.)
- **Entities are molecules**: an mmCIF entity is a `sequence` (polymer) or a `molecule` (non-polymer, branched, water);
  its chains are units. A **construct** (an ADC design, a fusion with its tags) is a complex with no structure behind
  it: the same units and links.
- **A unit instantiates one molecule**. A BIRD/PRD molecule spanning several entities (a peptide antibiotic: a polymer
  and its ligands, `pdbx_molecule`) is its units, each with the `prd:` identifier and the PRD instance in a unit field
  (`group`), drawn as separate blocks (a Tension: a chip for the group is later).
- **Copies are derived** (EC-65): units of one complex with equal `of` and `range` are a copy group; never declared.
  The copy rule holds (EC-25): no path climbs from a unit to its molecule and down to a sibling.
- **Annotations on any positions**: `ranges`, `values`, `pairs`, `grid` and `map` take a sequence, a molecule (either
  system), or a unit; a target that is a molecule reaches every unit over it (EC-48). A per-ligand score (AlphaFill's)
  is a field of the unit's row in a `ranges` over the molecule.
- **Bounded**: a complex from a source is overlay data, refused past `ResolveConfig.maxUnits` (65,536, a large
  capsid's assembly; beside `maxNestedDepth`) with an issue naming an assembly subset; merged copies cost one block
  (EC-67).

### Links: one shape, a `pairs` node

- **Links are a `pairs` node, a child of a complex or a molecule**: between units of the complex (a glycan on an Asn, a
  payload on a Cys, a disulfide, a metal), or between residues of the molecule (a glycosidic bond `NAG 2 C1 – NAG 1
  O4`). Its rows: `a_start`, `a_end`, `b_start`, `b_end` (one position each), and under a complex `a_unit`, `b_unit`
  (unit ids); `type` (`covalent`, `disulfide`, `glycosidic`, `metal`, or any word: `kind` is reserved, as pairs
  already uses `type`), `a_atom`, `b_atom`, `order`, `occupancy` (absent: 1), `label`. A child names no `a` or `b`: its
  parent is both. Both ends on one unit is allowed (an intra-chain disulfide, a staple); a site carrying one of several
  glycans (glycoforms) or a payload on one of many lysines (stochastic conjugation) is several rows with their
  occupancies. The mmCIF runtime names the complex's child `links` (`/data/1fin/links`).
- **Ends on per-row units**: `resolve/held.ts` holds pairs on two systems today; it holds a pairs node under a complex
  with a system a row end (two `Int32Array`s of unit indices beside the ends), and the `pairs` layer reads each end's
  system from its row, placing it on the block showing that unit (two blocks of one section: layouts.md's
  section-scope pairs). A `ranges` layer taking `pairs` data draws each end as a mark at its site labeled by the other
  end (as a reader without `pairs` draws a pair's ends today; taking it on purpose is part of this change). Contacts a
  host computes are a `pairs` node of the same shape (decided optional; never drawn by default).
- **Links and contacts are never maps**: projection does not cross them. Why: a map carries positions by arithmetic;
  a link relates two places without one, and projecting across it would make a glycan's selection select its Asn.
  What a user expects instead (the site lit when the glycan is hovered) is a behavior, below.

### Atoms to residues: a group map at a path's ends

- **One many-to-one map**: a molecule's atoms onto its residues, from its `atoms.residue` column (the residue `r` owns
  atoms `[off[r], off[r+1])`, contiguous by construction). Forward, atoms go to the residues they touch (as nucleotides
  to residues, EC-18); back, a residue to its atoms' run.
- **Not on paths**: a path is one `Linear` of runs at a fixed scale, which the codon step keeps because every residue is
  exactly three places; a group map has variable spans, so it never enters `findPath`. It applies only at a path's
  ends: a location on `(atoms)` collapses onto residues before projection, and a projected location on a molecule's
  residues expands onto `(atoms)` after it, when a section shows that form. Paths, their bounds and the codon step are
  unchanged. A writer's own group map (`map.method: group`) is not added: no view needs one.
- **Focus projects, follow does not** (EC-113: follow goes at scale 1, 3 or 1/3): a focus is one jump, a location
  carried through the ends' expansion once; a follow is the continuous image of a window edge every frame, which needs
  a fixed scale a group map lacks.

### On an axis: the block's `form`

- **`block.form`**: `whole` (the block collapsed over its range into one chip, layouts.md: the label and `×N`),
  `residues` (a position a residue or component, a molecule's each labeled by its code; a glycan's links drawn when a
  `pairs` layer over its links is added), or `atoms` (a position a heavy atom, labeled by element symbol, hovered by
  atom name and residue; a molecule with atoms only). Absent: `whole` for a molecule, `residues` for a sequence or a
  polymer unit. Why `whole`: a molecule's residues are a sliver at the fit view (the reason chips were decided), and
  one rule draws ATP, a heme, a glycan and a payload alike. `whole` on a short polymer (a peptide ligand) is the same
  drawing. Switched by patch (one undo step), as `copies`.
- **Hit, hover, selection** on a chip give the block's whole range on its system (a molecule's `(all)`, a unit's
  range); through the ends' expansion they stand for every atom (the previous implementation's "a ligand shown as one
  position is the locus of all its atoms").
- **A chip's label**: the molecule's `name`, else its one component's code, else its id; truncated at the chip's
  maximum width (LB's rule), whole on hover. The mmCIF runtime writes `name` from the entity's description
  (`pdbx_description`, the entry's own words); its `name: composition` option writes a branched molecule's
  composition in list order instead (`NAG2 BMA MAN2`), never a default.

### A section of a complex

- **`section.of` takes a complex**: a block a unit (layouts.md, Blocks), copies merged by `copies` (EC-65), ordered by
  the section's `classes` (default `[polymer, branched, ligand, ion]`; water and any class not listed left out),
  polymers in the complex's order, then each class in order (EC-101). Why left out: hundreds of waters say nothing on a
  sequence axis; a writer lists `water` to show them.
- **Explicit `block` children** still win: a writer picks units, forms and order.
- **Follow**: a section of units on several systems neither follows nor is followed (`seq.render.follow.blocks`,
  unchanged).

### Behaviors (the host's to keep or drop)

- **`links.hover`** (a `BehaviorProvider`): hovering a molecule (a chip, a residue, an atom) lights the polymer
  positions its links reach, as the `linked` emphasis; hovering a linked site lights the molecule. In the default spec
  (gap 10: a gesture default); a host leaves it out. Contacts get no such behavior: they are a layer the host adds
  (decided).
- **The hover carries identifiers**: the hovered element's molecule and its `identifiers` are in the hover record, so a
  host's card draws a structure from a SMILES or a CCD image. No depiction port: one consumer, the host's card, already
  has what it needs.

### Queries

| Form | Prefix | Result |
|---|---|---|
| `(component :id MSE)`, `(component :id [PSU 2MG])` | sequences, molecules, units | Location: the positions with that component (alternatives included) |
| `(component :standard false :annotations true)` | sequences, units | Location: every non-standard position, and the `modified`-type items of their annotations |
| `(component :class saccharide)` | as above | Location, through the `compound` nodes or the host's names entry; a code with neither is no match |
| `(molecule :identifier "ccd:ATP")`, `(molecule :class ligand)` | `data`, a complex | Nodes: molecules (with the host's resolver, its equivalents) |
| `(entity :id 3)` | a complex from a source | Nodes: the molecule its source calls entity 3 |
| `(chain :id A :operator 1_555)`, `(chain :auth B)` | a complex | Nodes: units (the form exists; `:operator` and `:auth` added) |
| `(atom :name O1G)` | a molecule's `(atoms)`, a unit | Location: one atom |
| `(atoms)` | a molecule, a unit | Nodes: its atoms system |

### Reading

- **mmCIF** through a runtime the application registers (`@molstar/seqx-core/runtimes/mmcif`, importing
  `@molstar/io`'s CIF reader and mmCIF schema on the allowlist, packages.md; never in the default registry). A source
  of format `mmcif` answers with nodes ([dynamic-nodes.md](dynamic-nodes.md)):
  - a sequence per polymer entity: letters from `pdbx_seq_one_letter_code_can`; where an entity has none, from
    `entity_poly_seq.mon_id` through the standard one-letter table, `X` (protein) or `N` (nucleic) for an unknown code,
    with one info issue; `modified` from `entity_poly_seq` where the code is not its letter's standard, `hetero` rows
    as alternatives, `pdbx_struct_mod_residue` corroborating;
  - a molecule per non-polymer (`pdbx_entity_nonpoly`) and branched entity (`pdbx_entity_branch_list` in `num` order,
    `pdbx_entity_branch_link` its links child, leaving atoms dropped, `pdbx_entity_branch_descriptor`'s WURCS and
    LINUCS into `identifiers`);
  - a complex of units (`struct_asym`, observed and the `auth` numbering from `atom_site`, `pdbx_branch_scheme`,
    `pdbx_nonpoly_scheme`; `operator` from the assembly asked; `pdbx_molecule`'s PRD per unit), its `links` child from
    `struct_conn` (`covale`, `disulf`, `metalc`), names from `chem_comp`.
  Its options: `assembly` (absent: the deposited model), `water` (absent: not read), `name` (`description`, absent;
  `composition`). Units from a source are reached by forms (`(chain :id A)`), inline units by id. A molecule's atoms
  come from the CCD's `chem_comp_atom` when a CCD source is given, else the first unit's `atom_site` (altlocs as above;
  a unit's `observed` marks its modeled atoms).
- **The `seq` format** carries the same nodes and tables (the fixtures' form, converted at recording).

### Language, store, resolved model

- **Language (0.5.0)**: `sequence.modified`; the kinds `molecule` (`residues`, `atoms`, `class`, `name`,
  `identifiers`), `compounds` and `compound`, `complex`, `unit` (and a unit's `numbering`, if not there); `pairs` as a
  child of a complex or a molecule, with `a_unit`, `b_unit`, `a_atom`, `b_atom`, `order`, `occupancy`; the kinds a data
  reference takes widened to molecules and units (`ranges.target`, `values.target`, `pairs.a`/`b`, `grid.x`/`y`,
  `map.from`/`to`, `block.of`, `section.of`, which takes a complex too); `block.form`; `section.classes`; the `letters`
  layer's `modified`; a `ranges` layer taking `pairs` data. Checks (`seq.lang.*`): a modified position past the end or
  out of order, a `length` crossing the next row, a molecule's `atoms.residue` decreasing or past its residues, a
  hydrogen, a link's position past its unit's range or naming no unit, a unit `range` wrapping, `form: atoms` on a block
  with no atoms.
- **Store**: molecules, compounds, complexes and their pairs are held whole (`resolve/held.ts`, ends on per-row units as
  above): typed arrays keyed by address and hash; a complex from a source is one request, its nodes the overlay's,
  bounded by `ResolveConfig.maxUnits`.
- **Resolved model**: systems for units and a molecule's residues and atoms; the atoms group map applied at a path's
  ends; copy groups; the modified index a system (a sorted position index with its alternatives); the query model's
  `chain`, `node` and the new forms answering.

## The interaction matrix's cells

Rows and columns of ux-foundations section 8; "as an element" means the grammar's cell, unchanged.

| | select | extend | add / remove | focus | back | zoom in / out | pan | hide | edit |
|---|---|---|---|---|---|---|---|---|---|
| **a section of a complex** | — | — | — | — | — | as its layout | as its layout | no: its tracks | `classes`, `copies`, `arrange` by patch (`section.set`, previewed) |
| **a block: a chip** (`form: whole`) | click: the block's whole range (every atom, every merged copy) | Shift+click: each block's part between the anchor and it, in reading order (layouts.md) | Cmd/Ctrl+click | double-click: revealed; a sibling section showing its residues or atoms zooms to them (focus projects through the atoms map, one jump; follow does not, it needs a fixed scale) | double-click again, Alt+Left | no: a chip is fixed in px (layouts.md) | no: the window skips chips | no: `classes` or the blocks list by patch | `block.form` by patch (`whole`, `residues`, `atoms`), no default key, in the context menu (as `copies`) |
| **an element**: a modified residue | as an element | as an element | as an element | as an element | as an element | — | — | no: `modified: none` on the layer | no: data |
| **an element**: a molecule's residue or atom | as an element; atoms collapse onto residues for other sections | as an element | as an element | as an element | as an element | — | — | no | no: data |
| **an element**: a link (an arc, or a site's mark) | click: the link item; both ends lit | Shift+click | Cmd/Ctrl+click | double-click: its extent, two parts on two systems, each section revealing its own | Alt+Left | — | — | its layer's items by type (`layer.filter`: hide disulfides) | no: data |
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
| `modified` (`letters`) | spec param | `mark` | a modified residue drawn as its parent alone hides it | `code`, `tint`, `none` |
| the mark's size and place | `RenderConfig` (`letters.modifiedMark`) | a dot of 3 px at the cell's top right; a tick above a strip | the presets pattern | — |
| `letters.modifiedMark` | `RenderConfig` preset | a dot at the cell's top right in the letter's ink; a bar over a long component | reads at 8 px, never covers the letter | `underline` |
| `code` (`letters`) | spec param | `auto` (exists): the component's code where three letters have room | a cell with room says what the residue is | `one`, `three` |
| `class` (molecule) | spec param | `ligand` | the commonest non-polymer | `branched`, `ion`, `water`, any word |
| `:annotations` of `(component …)` | query param | true | one answer for "where are the modified residues" | false: the sequences' own only |
| `componentNames` | registry entry | none (the code alone) | the core ships no CCD | a host's CCD lookup (Mol\*'s) |
| `identifierResolver` | registry entry | none (exact match only) | equivalence is chemistry | a host's (InChIKey, ChEMBL to CCD) |
| identifier namespaces | open words | the documented vocabulary | a host's registry ids are pointers too | a fixed list |
| `links.hover` | `BehaviorProvider` | off; a host adds it (decided: the author, 2026-11-19, "for the linked highlight: this is a full app level decision anyway how the events from different parts compose"); on in the examples and their stories (the author: "the demo should show it on") | a conjugate read with its site | added by a host |
| contacts | a layer the host or builder adds | none (decided) | "contacts should be optional and configurable by the app" | — |
| the mmCIF runtime | registry entry | not registered | the core reads no structure file by default; its size | registered by the application |
| `assembly`, `water`, `name` (mmCIF source) | source params | the deposited model; waters not read; `description` | what the entry is and says | an assembly id; `true`; `composition` |
| `ResolveConfig.maxUnits` | `ResolveConfig` | 65,536 | a large capsid's assembly; a source's complex is overlay data | larger, at its memory |
| altlocs in `atoms` | fixed | the first `label_alt_id` an atom name has | a position is one atom | — |
| hydrogens in `atoms` | fixed | refused | lengths and per-atom data line up across entries (EC-64) | — |
| links and contacts as maps | fixed | never | a glycan's selection never selects its site | — |
| group maps on paths | fixed | never: at a path's ends only | a path is runs at one scale | — |

## Tensions

- **A new kind against a generalized `sequence`.** One kind with an alphabet of components would keep the language
  smaller; every letter reader would then check it. Recommended: the new kind (gap 1).
- **One position, one letter.** CRO stands for three parents; recommended: the canonical letters, the component over
  three positions (`length` 3), the structure's one monomer grouped by the structure map. The alternative (the
  structure's monomers as positions, CRO one `X`) keeps the structure's numbering and loses UniProt's. The structure
  maps' design (the fourth slice) settles the map side; open until it does (gap 6).
- **A unit and maps.md's structure system.** maps.md calls a chain instance (entry, model, operator, chain) a
  structure system; here it is a unit of a complex read from a structure. Recommended: one noun, the unit; a structure
  map relates a reference sequence (UniProt) to the entity's molecule, and SIFTS's per-chain segments become per-unit
  `observed` and `auth` numbering. The structure lane decides; maps.md changes when this is agreed (Consequences).
- **The mmCIF runtime against the structure lane.** Mol\*'s `Model` already parses entities, `chem_comp`, branched
  entities and microheterogeneity; a Mol\* host would read the same file twice. Options: (a) the Mol\* host derives the
  nodes from its `Model` through the structure port, MC only for headless use (the render CLI, Node); (b) one converter
  from CIF categories to nodes, which MC's runtime feeds from io's reader and the Mol\* host feeds from its model's
  source frame. Recommended (b): one reading of the dictionary, one conformance suite, the host paying no second parse
  (gap 11).
- **A conjugate far from its site.** Chips follow the polymers, so a glycan on N297 is drawn after the chains. The
  sites (a ranges layer over the links) and `links.hover` carry the attachment; a chip drawn at its site (a pin above
  the residue) would be a placement of its own in every layout. Open: a later design, if a view needs it.
- **A glycan is a tree, an axis is a line.** `residues` lays a glycan in list order, its links as arcs; the tree
  (SNFG) needs a symbol a component, a host's table (Mol\* has one). Planned later: a tree drawing of a molecule's
  residues with a host's `symbols` registry entry; not in the fourth slice unless a view asks.
- **A PRD over several units.** A unit instantiates one molecule, so a BIRD molecule is several blocks tied by a field.
  Standing until a view needs one chip for the group.
- **Poses and copies.** Ten docking poses of one ligand are a copy group by the rule and merge to `×10` by default,
  hiding that they differ. The builder for a pose set writes `copies: each`; the rule is not bent. Standing.
- **Stochastic conjugation reads as many links.** An average DAR of 3.5 over 40 lysines is 40 rows with occupancies;
  the chip says `×N` of units, not the average. A label field carries it. Standing until an ADC view is used.
- **Atoms without coordinates.** An atoms form with no bonds or positions is a list; it earns its place by per-atom data
  (fit, occupancy, B-factor) and by matching a structure viewer's atoms. Kept; a map pairing two molecules' atoms by
  name (a pose onto its component) waits for a view that needs it.

## Self-check

The claim: "this design covers what a user does with molecules: small molecules, glycans, conjugates, non-standard
residues and complexes".

1. **Assume it holds.** Everything is a system of positions: a molecule has residues and atoms joined by a group map at
   a path's ends; a modified residue is a component on a position whose letter is its parent; a complex instantiates
   molecules as units; links and contacts are pairs; identity is a pointer matched as a string.
2. **What it rests on, negated.** (a) *One position, one component, one letter.* Negated: a component over three
   letters (CRO), two components at one position. (b) *A molecule is named by a pointer.* Negated: no pointer (UNL),
   and one compound under several pointers. (c) *A link has one site and one partner.* Negated: glycoforms, stochastic
   lysines, a link between two polymers or within one, a molecule linked twice. (d) *A complex is a structure's
   chains.* Negated: a construct with no structure, a set of poses, thousands of units. (e) *Its inverse*: a polymer
   drawn as a molecule (a peptide ligand, a cyclic peptide); a molecule's data on the polymer. (f) *A missing view*:
   the glycan's tree, the compound's drawing, a table of ligands. (g) *The data comes from a structure.* Negated:
   UniProt's PTMs, known by name with no component and no glycan. (h) *A molecule is one entity.* Negated: a BIRD/PRD
   over several entities and units. (i) *A unit's positions are labelled by index.* Negated: the author's numbers and
   insertion codes every structure user types. (j) *One atom per name.* Negated: alternate locations. (k) *Many-to-one
   is a map like the others.* Negated: a variable span has no fixed scale, so it cannot be a path step.
3. **The negated world.** (a) gives `length` and alternatives; (b) molecules valid without identifiers, open
   namespaces, exact matching and the host's resolver; (c) links as `pairs` rows with occupancy, self-links, ends on
   per-row units; (d) units over molecules with no structure, `copies: each` for poses by the builder, `maxUnits`;
   (e) `whole` as a drawing of any block, `cyclic` sequences (the param exists), annotations on molecules and units;
   (f) the tree and drawing as the host's, a ligand table as a ranges over the molecule; (g) name-only modifications
   as `modified`-type annotation items the component query finds; (h) the PRD as an identifier and a unit field;
   (i) the unit's `auth` numbering; (j) the altloc rule; (k) the group map out of `findPath`, at a path's ends, and
   focus projecting while follow does not.
4. **Held both.** Integrated in the Decision: (a) to (e), (g), (i), (j), (k) as above. Planned with a named lane: the
   structure map's grouping of a multi-letter monomer and the unit as the structure system (the structure lane); one
   converter for MC and the Mol\* host (gap 11). Standing tensions: CRO's positions; a conjugate drawn far from its
   site; the glycan's tree; a PRD's units drawn apart; poses merged as copies by default; the average DAR.

## Consequences

The fourth slice builds it, after the third slice's DG, DC, DT and LB (chips) and with the structure lane beside it.

- **`seqx-core`**: `lang/` (the kinds, params and checks above; the builder's `molecule`, `complex`, `unit`,
  `modified`, `compounds`; issues `seq.lang.molecule.*`, `seq.lang.modified.*`, `seq.lang.unit.*`,
  `seq.lang.link.*`; `SeqLanguageVersion` 0.5.0); `resolve/` (held entries with per-row units, unit and atoms systems,
  the atoms group map at a path's ends, copy groups, the modified index, `maxUnits`); `view/` (a section of a complex,
  `classes`, `form`, chip labels); `render/` (the `letters` mark and the component code; residue labels of a molecule;
  element symbols on atoms; the `pairs` layer's ends by row, a `ranges` layer over pairs); `query/` (the forms above);
  `data/` (the `seq` format's molecule, compound and complex tables); `runtimes/mmcif` (new entry, the `@molstar/io`
  allowlist entries added).
- **`@molstar/core`**: the standard component per letter and alphabet, the inverse half of base-library.md's Promote
  row "residue alphabets and schemes" (whose one-letter codes are `getProteinOneLetterCode` and friends in
  `packages/model/src/model/sequence/constants.ts`); the agreement adds a row to base-library.md: "standard components |
  `@molstar/core`'s sequence constants, beside the one-letter codes | the standard component of each letter by
  alphabet (A ↔ ALA, DA, A)".
- **`seqx-viewer`**: hover text with component, alternatives, atom name and residue; the hover record's identifiers;
  `links.hover`; `block.form` in the context menu; the matrix cells above.
- **`seqx-examples`**: DT (1EHZ's 11 modified residues, marked and found by query); DC (1FIN as a complex from its
  fixture: units with `auth` numbering, ATP `whole` and `atoms`, AlphaFill's transplants as molecules with scores,
  contacts a `pairs` layer; `links.hover` on, both ways); DG (4HHB's hemes as molecules, `links.hover` on); Antibody gains an Fc glycan (an IgG1 entry with its N297
  glycan, *the entry chosen at recording*, CC0) and an ADC construct (trastuzumab with a payload by SMILES on lysines,
  generated, occupancies), with their views and editing views.
- **Storybook**: a `Seq* (next)/Design/Molecules` page; stories for each of the views above, each with its headless
  twin on a fake clock.
- **Tests**: properties (atoms collapsed then expanded cover what they started from; a modified position's letter is
  its parent's; copies derived equal for any unit order; a path's result unchanged by a molecule's atoms form); the mmCIF
  runtime's conformance suite on 1EHZ, 1FIN, 4HHB and the glycoprotein (every entity, chain, branch link, `struct_conn`
  row and author number read, against counts recorded from the file; an entity without canonical letters; an altloc);
  drawing against the data (every modified position in view has its mark; every chip its label; every link with both
  ends in view its arc and its hit box); a quick representative of the frame target with 200 merged copies and a glycan
  chip (EC-67).
- **[maps.md](maps.md) when agreed**: a molecule's atoms onto its residues as a many-to-one map applied at a path's ends,
  never a kind on paths; "a structure system is a chain instance" read as a unit of a complex from a structure (the
  structure lane confirming); units as systems on the paths' free steps, as built.
- **`_spec` sections that change**: core's Language (kinds, params), Resolve (Systems; Systems, sources, annotations:
  pairs held with per-row units; Paths: the group map at a path's ends; `maxUnits` in Configuration), Resolved view
  (Blocks: a section of a complex, `form`), Render (Layer types: `letters`, `pairs`, `ranges`), Queries (Forms), Data
  (the seq format), a new Runtimes section (mmCIF); viewer's HoverModel, Behaviors, Commands; both Configuration tables.

### The fourth slice's lanes

- **ML1** (strong, first) the language: kinds, params, checks, builder, the `seq` format's tables, the 0.5.0 bump;
  the bump's owner for the whole slice (other lanes' language changes go through it).
- **ML2** (strong, after ML1) the model: unit and atoms systems, the atoms map at a path's ends, held pairs with per-row
  units, copies, the modified index, the query forms.
- **MC** (strong, after ML1, beside ML2: it needs only the node shapes) the mmCIF converter and runtime: the allowlist,
  the tables read, its conformance suite, the fixtures re-recorded; sizes measured; shared with the structure lane as
  gap 11 decides.
- **MR** (strong, after ML2 and MC) drawing and interaction: a section of a complex, `form` through LB's chips, the
  letters' mark and code, the pairs layer's per-row ends, hover, `links.hover`, `block.form`, the matrix cells, the
  Design page.
- **MX** (strong, after MR) the examples changed once: DT, DC, DG to the new kinds; the glycan and ADC views; stories,
  headless twins, drawing against the data.

## Decided by the author

2026-11-19: "Other recommendations from the design ok" (gaps 1–9 and 12–14 as recommended); gap 10 off; gap 11
superseded by the seq data ports design.

1. A new `molecule` kind beside `sequence`, rather than one kind with an alphabet of components? Recommended **yes**: **decided yes**.
2. A modified residue's letter is its parent's (MSE is `M`), the component beside it? Recommended **yes**: **decided yes**.
3. A modified residue marked by default in `letters` (`modified: mark`)? Recommended **yes**: **decided yes**.
4. A molecule block drawn `whole` (a chip) by default, a glycan included, its residues and atoms by patch? Recommended
   **yes**: **decided yes**.
5. Water left out of a section of a complex by default (`classes` without `water`)? Recommended **yes**: **decided yes**.
6. A component over several positions (CRO as T, Y, G): the canonical letters, rather than the structure's one
   monomer? Recommended **yes**: **decided yes**, until the structure maps' design says otherwise.
7. A modification known by name only (UniProt's PTMs) as a `ranges` item of type `modified` that
   `(component :standard false)` also finds, rather than a `modified` row with a name and no component? Recommended
   **yes**: **decided yes**.
8. Several components at one position: all kept, the first giving the letter, all on hover and in queries? Recommended
   **yes**: **decided yes**.
9. Identifier namespaces open words, with a documented vocabulary, rather than a fixed list? Recommended **yes**: **decided yes**.
10. ~~`links.hover` on in the default spec~~ **Decided off** (the author, 2026-11-19): how events from different parts
    compose is the application's; the library offers the behavior, a host adds it. The examples and their stories
    turn it on, both ways (the author: "For #10 for the two way highlight.. the demo should show it on").
11. ~~One CIF-to-nodes converter~~ **Superseded** (the author, 2026-11-19: "we should also be able to share data
    representation with mol* ... maybe we could have generic representation ports, similar to mol io, we could write
    seq io"): how Seq\* and a host share one data representation (Seq\*'s nodes read from a Mol\* `Model` or
    `Structure` in memory, not parsed again; Seq\*'s annotations given to Mol\* as representations) is its own design,
    seq data ports (the author's name for it: "the same principle not the same functionality"), in the fourth
    slice before MC; MC builds on it.
12. The mmCIF runtime as an opt-in entry of `seqx-core` (`runtimes/mmcif`) on `@molstar/io`'s CIF reader, rather than
    a package of its own? Recommended **yes**: **decided yes**.
13. The standard component per letter in `@molstar/core` beside the one-letter codes (a public addition, a new
    base-library.md row)? Recommended **yes**: **decided yes**.
14. This slice's DC row "ATP ... joined by a many-to-one map" (examples.md): 0.4.0 has no such map (`offset`, `codon`
    only). DC writes ATP's two forms as two sequences with no map this slice, hover not crossing until the fourth
    slice? Recommended **yes**: **decided yes**, and examples.md's row reworded with DC.
15. Structure integration in the fourth slice (the author, 2026-11-19: "let's only do a demo with fake 3D viewer to
    demonstrate how the syncing will work. There will be a new mol viewer work after we are done with the sequence
    viewer which will replace prism"): the syncing is shown on a fake 3D viewer (a storybook stand-in taking
    MolViewSpec-shaped states and sending hover, selection and focus); the Mol\* host, its `Model` adapters (Tensions,
    the mmCIF runtime against the structure lane, option (a)) and a live 3D view wait for the new mol viewer. The
    mmCIF runtime (MC) reads files through `@molstar/io` as designed.
