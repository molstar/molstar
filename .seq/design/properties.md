# Property access

## Intent

The author, 2026-11-20: "we should have a generic filter, and generic properties, but it should be easy to add own
filter, modified/combinator etc. likely won't need it many cases, but in general. you are mixing two things: property
access and language design"; on field names, "support both short and long aliases"; "also make the language extensible
so people could add generators, filters, etc on their own props (basically just adding a few symbols)".

The query language ([locations-and-queries.md](locations-and-queries.md), The language) says how forms combine; this
design says what a predicate reads: properties, typed symbols read per element, as MolQL's `struct.atomProperty.*`
apart from its generators and filters.

## Decision

- **A property is a typed symbol read per element** of one level (`molecule`, `unit`, `residue`, `atom`): `(prop
  chain)` in a predicate, or `:chain A` as an `elements` field (sugar for `(= (prop chain) A)`). Types: `word`,
  `integer`, `number`, `boolean`, `words` (a set; `=` on a set is membership: `component` holds a position's
  alternatives).
- **Built in**, the short name printed and the long one (MolViewSpec's, mmCIF's) read:

| Level | Properties |
|---|---|
| molecule | `entity` (`label_entity_id`), `molecule-class`, `identifier` (words) |
| unit | `chain` (`label_asym_id`), `auth-chain` (`auth_asym_id`), `operator` (`instance_id`) |
| residue | `seq` (`label_seq_id`: 1-based on the unit's molecule), `auth-seq` (`auth_seq_id`), `ins` (`pdbx_PDB_ins_code`), `index` (`residue_index`, 0-based on the unit), `component` (`label_comp_id`, words), `component-class`, `standard`, `modified` (non-standard, or under a `modified`-type annotation item), `b-factor` (`B_iso_or_equiv`), `occupancy` (values channels) |
| atom | `atom` (`label_atom_id`), `element` (`type_symbol`), `atom-index` (`atom_index`) |
| any | `(value <values annotation> :field f)`: a values channel on the element's system |

  Why short printed: an address is typed and read beside a sequence, and one printed form keeps equal selections equal
  as text; the long names are a structure viewer's, so the inward translation is a copy. One property given twice under
  both names is a form error. Not held: `label_alt_id`, `auth_comp_id`, `auth_atom_id`, `atom_id` (Seq\* keeps one atom
  a name and no author names, [molecules.md](molecules.md)). `seq` is the molecule's ordinal, not the unit's: a unit
  counts from 0 over its range, `label_seq_id` is the entity's.
- **How a property is read**: a column per element on one system, `{ length, get(i) }` (typed where it can be), from
  the `QueryModel` hook `property?(system, name): PropertyRead`, `PropertyRead = column | 'pending' | undefined`. Units,
  components, numberings and atoms answer the built-ins (`QueryUnit.entity`; `QueryAtoms.element`, `residue`); a values
  channel reads the store (`values` annotations by region), `pending` staging the system as a numbering not arrived does;
  `b-factor` and `occupancy` are values channels, read where their channel lies (the mmCIF runtime writes them on a
  unit's residues), absent otherwise. A residue property on an atom reads its residue's (through the atoms map); the
  reverse is `seq.query.property.level`, except one case: an atom field on a residues system with no atoms held reads
  its standard component's heavy atoms (a modified residue's parent's), and the result is residues (the author,
  2026-11-20: "Residues, as a structure viewer reads them").
- **An element that lacks a property** (`chain` on a lone molecule, `auth-seq` with no `auth` numbering, `atom` on a
  component with no standard atoms, `entity` on an inline complex): its test is false and it yields nothing, with one warning
  `seq.query.property.absent` naming the property and the count. Why: monadic composition; never the whole query's
  failure.
- **Extensions add properties** as `query-property` registry entries ([extension-openings.md](extension-openings.md)),
  read by every test and filter and usable as `elements` fields:

  ```ts
  | { kind: 'query-property'; name: string;  // `af:plddt`
      level: 'molecule' | 'unit' | 'residue' | 'atom'; type: PropertyType; aliases?: readonly string[];
      read: ((model: QueryModel, system: string) => PropertyRead)
          | { values: string; field: string } }  // a few symbols: a values annotation's field, by kind or path
  ```

  Namespaced `<id>:<name>` (aliases too) and checked at composition as every extension entry. An unknown namespaced
  property is `seq.query.extension.missing`; an unknown bare one `seq.query.form.param`, its hint listing the names (and
  the short name for a long one misspelled).
- **Outward**: no property reaches a structure viewer; the adapter translates the Location selected
  ([structure-syncing.md](structure-syncing.md)).

## Consequences

- `seqx-core` `query/`: the property table, `prop`, `value`, `QueryModel.property`; `context/`: `query-property`.
- Tests: each built-in and alias on the 1FIN and 4HHB fixtures; a lacking property; a pending values channel; an
  extension's property, registered and missing.
- `_spec/README.md` "Queries" gains Properties.
