# The `(elements …)` form: schema selections in Seq\*'s queries

## Intent

The author, 2026-11-20: "we should not really have (atom), we should have (elements ....) and that gets "schema
selection" etc .. similar to "atom groups" combined with "schema item" in molQL .... could still have two generators."
Earlier the same day: forms compose monadically, as MolQL's (`(atom :name C1)` over a glycan's three residues gives
three atoms).

ML2's branch built one form a level (`entity`, `chain`, `component`, `atom`, `atoms`, `molecule`). A structure user
thinks in one record instead: chain, residue, atom, as MolViewSpec's component expressions and MolQL's schema items
write it, and as a structure viewer sends it. This draft gives Seq\* that record as a form, and says which level forms
it replaces.

## Precedents

- **MolQL `struct.type.schemaItem` + `struct.generator.elementSchema`**: a record of mmCIF fields (`label_entity_id`,
  `label_asym_id`, `auth_asym_id`, `instance_id`, `label_seq_id`, `auth_seq_id`, `pdbx_PDB_ins_code`,
  `beg_/end_label_seq_id`, `beg_/end_auth_seq_id`, `beg_/end_auth_ins_code`, `label_comp_id`, `label_atom_id`,
  `type_symbol`, `atom_id`, `atom_index`, `residue_index`); every given field must match; several items a union; one
  set out, found by index, not by testing every atom. Taken: the fields, all-must-match, the union, the index walk. Not
  taken: a list as a table's column (rows zipped); Seq\*'s lists are alternatives (`(component :id [PSU 2MG])`).
- **MolQL `struct.generator.atomGroups`**: `entity-test`, `chain-test`, `residue-test`, `atom-test` (predicates over
  `struct.atomProperty.*`, each run on a level's first atom), `group-by` (default: an atom a set). Taken: tests by level,
  grouping. Not taken: arbitrary predicates (Seq\* has no expression language; its tests are field equalities and the
  chemistry words `component` already reads).
- **MolViewSpec `ComponentExpression`** (`packages/mvs/builder`, `param-types.ts`): MolQL's schema item minus the
  insertion-code range ends and `operator_name`; a component's `selector` is one expression or a list (a union). This
  is what the fake host's adapter writes (structure-syncing draft: a polymer run is `{ label_asym_id, instance_id,
  beg_label_seq_id, end_label_seq_id }`, an atom adds `label_atom_id`). Taken: its field set as the round-trip set.
- **Neither has `label_alt_id`** in the record (MolQL has it as an atom property only).
- **Interactions**: a user types or pastes a selection (the query command: built), follows a link or bookmark holding
  an address (built), and receives a structure viewer's selection (structure syncing inward: planned, ST2). No new
  gestures; the matrix's *selection* row gains this form through the query command.

## Decision

### Fields: MVS's record, Seq\*'s spelling

One table maps each field to its MVS name, so the inward translation is field for field and the outward one (a
Location to expressions, the adapter's) is unchanged.

| Field | MVS / mmCIF | Level | Value |
|---|---|---|---|
| `:entity` | `label_entity_id` | molecule | word |
| `:chain` | `label_asym_id` | unit | word |
| `:auth-chain` | `auth_asym_id` | unit | word |
| `:operator` | `instance_id` | unit | word (absent: every instance) |
| `:seq`, `:beg`, `:end` | `label_seq_id`, `beg_/end_label_seq_id` | residue | integer, 1-based on the unit's molecule; `:beg`/`:end` inclusive, either alone open |
| `:auth-seq`, `:auth-beg`, `:auth-end`, `:ins` | `auth_seq_id`, `beg_/end_auth_seq_id`, `pdbx_PDB_ins_code` | residue | through the unit's `auth` numbering; ends snap inward (EC-1 to EC-8) |
| `:index` | `residue_index` | residue | integer, 0-based on the unit |
| `:component` | `label_comp_id` | residue | word; alternatives at a position match |
| `:atom` | `label_atom_id` | atom | word |
| `:element` | `type_symbol` | atom | word |
| `:atom-index` | `atom_index` | atom | integer, 0-based on the atoms system |

- **A list is alternatives** (`:component [NAG MAN]`, `:atom [N CA C O]`), as `component :id` today. Why: one rule for
  every form's lists; rows of a table are a union of forms (below).
- **`:seq` is the molecule's ordinal**, not the unit's: a unit counts from 0 over its range, and `label_seq_id` is the
  entity's (`entity_poly_seq.num`); they differ only for a unit whose range starts past 0. On a molecule's residues
  (`label_seq_id` is `.` in mmCIF for ligands and glycans) it is the residue's ordinal; the outward adapter writes
  `auth_seq_id` there, as it already does.
- **No `label_alt_id`, `auth_comp_id`, `auth_atom_id`, `atom_id`**: Seq\* keeps one atom a name (molecules.md) and
  holds neither author component and atom names nor model-file atom ids; the record stays what both sides hold.

### The result: the finest field decides

`(elements …)` is a **Location**. With no residue or atom field, the whole residue system of every unit that matches
(`(elements :chain A)`: positions; `(chain :id A)` stays the Nodes form). With residue fields, those positions of the
units' residue systems. With an atom field, positions of the units' **atoms systems** (`<unit>/(atoms)`), and the
residue fields restrict them through the atoms-to-residues map. Why one rule: the record says how deep it goes, as
MVS's does, and a Location on atoms collapses onto residues for every other section already (molecules.md).

### Two generators

1. **`(elements <fields>)`**: the schema item (MolQL's `element-schema` with one item; MVS's expression); a Location.
   The round-trip form: a structure viewer's selection comes in as one `elements` an expression, `(or :of [ … ])` for
   a list.
2. **`(groups <fields> <tests> :by <level>)`**: atom groups. The same fields, plus the chemistry tests the level forms
   held: `:standard` and `:class` (residue: a component's; from `component`), `:identifier` and `:molecule-class`
   (molecule; from `molecule`), `:annotations` (as today). `:by` is `atom`, `residue` (default), `unit` or `molecule`;
   the result is **Elements**, a key a group (owner the unit's or atoms system's address and a 0-based position; a
   complex's address and a unit id or entity id), each key's extent its Location. Why a second generator: clicking,
   labels and bookmarks act on items (locations-and-queries.md), "each modified residue" is a list, not a range; and
   chemistry tests are not schema fields, so `elements` stays exactly the round-trip record.

### Composition

- **Monadic over its prefix**: a data container (each complex and molecule under it), a complex (each unit), a unit,
  a molecule (its residues; unit fields are absent, below), a sequence (residue fields only), or a **Location**
  (refines: the positions it holds, the atoms of the residues it holds). Each element yields its part; one with
  nothing yields nothing; results concatenate (a union).
- `/data/1fin/(chain :id A)/(elements :seq 33 :atom NZ)` = `/data/1fin/(elements :chain A :seq 33 :atom NZ)`.
- `(range …)`, `(at …)` stay: any numbering on any system (a genome, an alignment); `elements` is the structure
  record. `(and …)`, `(or …)`, `(not …)`, `(extent …)`, `(items …)` combine and convert as before.

| ML2's form | Becomes |
|---|---|
| `(atom :name O1G)` | `(elements :atom O1G)` (removed) |
| `(component :id MSE)` | `(elements :component MSE)` |
| `(component :standard false)`, `:class` | `(groups :standard false)`, `(extent :of …)` for a Location (removed, or kept: gap 3) |
| `(entity :id 3)` | `(elements :entity 3)` (Location); the molecule node is `(molecule :entity 3)` (removed) |
| `(chain :id A :operator 1_555)`, `(chain :auth B)` | stays: Nodes, what `(section :of …)` takes |
| `(molecule :identifier … :class …)` | stays: Nodes; gains `:entity` |
| `(atoms)` | stays: the atoms system's address |

### Examples

- **4HHB, a heme**: `/data/4hhb/(elements :auth-chain A :component HEM)`: the heme on α1 (its unit shares
  `auth_asym_id` A with the globin; `:component` keeps the heme). `/data/4hhb/(elements :component HEM :atom FE)`: four
  irons, one a heme unit's atoms system. Its proximal His: `(elements :auth-chain A :auth-seq 87)`.
- **1FIN, CDK2**: `/data/1fin/(elements :chain A :beg 10 :end 20 :atom CA)`: eleven Cα atoms of the glycine loop.
- **A glycan** `[NAG NAG BMA]`: `/data/glycan/(elements :atom C1)`: three atoms, one a residue (the author's case);
  `(elements :component NAG :atom C1)`: two.
- **A modified residue**: `/data/1ehz/(elements :component PSU)`: its pseudouridines; `/data/1ehz/(groups :standard
  false)`: one key a modified residue (11); an MSE's selenium: `(elements :component MSE :element SE)`.

### Address, printing, issues

- **Address**: `path/(elements :chain A :beg 10 :end 20 :atom CA)`, the grammar's own. **Printed** with the fields in
  the table's order, then the tests, then `:by`, so two equal selections print equal (addresses compare and hash).
- **An unknown field**: a form error (`seq.query.form.param`), its hint listing the fields; an MVS name
  (`label_asym_id`) is hinted with its spelling (`:chain`). `:seq` with `:beg`/`:end`, or `:beg` past `:end`: a form
  error.
- **A field the element lacks** (`:chain` on a lone molecule, `:auth-seq` on a unit with no `auth` numbering, `:atom`
  on a unit without atoms, `:entity` on an inline complex): that element yields nothing, with one warning
  `seq.query.elements.field` naming the field and how many elements lacked it. Why: monadic, never the whole query's
  failure (ML2's `atom` failed the form; this changes it).
- **Labels not arrived** (a source's `auth`): staged, as `range` today (`seq.query.numbering.pending`).

## Options

| Option | Cost | Gives |
|---|---|---|
| A. `elements` only; `component`'s tests and `molecule` stay as forms | smallest (`atom` out, `component :id` out) | one record; no grouping; two places for residue tests |
| **B. `elements` + `groups`** (recommended) | `groups`, the Elements keys a level; `component` folded | the round-trip record exact; chemistry tests and grouping in one place |
| C. one form, `:by` turning it into groups | one kind | a form whose result type depends on a param, refused by "a form that mixes types is an issue" |
| Spelling: MVS names verbatim (`:label_asym_id`) | none in translation | long addresses; `auth`/`label` reads as noise to a sequence user |

## Configuration

| Item | Kind | Default | Why | Alternatives |
|---|---|---|---|---|
| field names | fixed (language 0.5.0) | the table above | one table, MVS field for field | MVS's names (gap 1) |
| list semantics | fixed | alternatives | every form's lists | MolQL's columns |
| `groups :by` | form param | `residue` | what a user clicks and labels | `atom`, `unit`, `molecule` |
| `:annotations` of `groups :standard false` | form param | true | as `component` today | false |
| a field the element lacks | fixed | nothing, one warning | monadic | a form error |

## Tensions

- **The word.** The result type *Elements* (items, rows, pairs: keys) and the form `elements` (a Location) collide.
  Options: rename the type (*Keys*, in 0.5.0, before release) or name the form `schema`. Recommended: rename the type.
- **Two ways to a range**: `(range :numbering auth …)` and `(elements :auth-beg … :auth-end …)`. Kept: the first is any
  numbering, the second the record a structure viewer speaks; both resolve through one numbering code.
- **Self-check** ("`elements` covers what a structure user selects"): negated, a user selects by what the record lacks
  (an altloc, a model, a secondary structure, a distance). Integrated: none; altlocs and models are the structure
  viewer's (molecules.md), and `include-surroundings` is a host's contacts as `pairs`. Standing: a selection with those
  comes in as its residues (the adapter's) and is not addressable here.

## Consequences

- **`seqx-core` `query/`**: `elements.ts` (the form; fields over `QueryUnit`, `QueryComponents`, numberings,
  `QueryAtoms`), `groups`; `molecule-forms.ts` loses `atom`, `entity`, `component` (or its `:id`). **`QueryModel`**:
  `QueryUnit.entity` (replacing ML2's key-of-`of` lookup), `QueryAtoms.element` and `residue` (the map is held; the
  form reads it). Tests: each field on the 1FIN and 4HHB fixtures, the four examples, the monadic cases (a glycan's
  three C1s, a unit without atoms), printing round trip, a property test (an MVS expression inward then outward through
  the fake's adapter gives the expression's observed part).
- **`_spec/README.md` "Queries"**: Forms gains `elements`, `groups`; loses `atom`, `entity`; molecules.md's Queries
  table and locations-and-queries.md's first forms follow. ML2's branch carries the change before it merges.
- **Migration**: none; 0.5.0 is unreleased.

## Gaps for the author

1. Seq\*'s short field names (`:chain`, `:auth-seq`), not MVS's verbatim? Recommended: yes.
2. Two generators as `elements` (Location, the record) and `groups` (Elements, tests and `:by`)? Recommended: yes.
3. `component` removed, its `:standard`, `:class` and `:annotations` moving to `groups`? Recommended: yes.
4. Rename the *Elements* result type (to *Keys*) so the form owns the word? Recommended: yes.
5. A Location as `elements`' prefix (refining, the glycan's C1s over a selection)? Recommended: yes.
6. No altloc field? Recommended: yes (neither MVS nor Seq\*'s atoms hold one).
