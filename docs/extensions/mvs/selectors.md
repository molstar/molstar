# MVS selectors

Selectors are used in MVS to define substructures (components) and apply colors, labels, or tooltips to them. MVS nodes that take a `selector` parameter are `component` (creates a component from the parent `structure` node) and `color` (applies coloring to a part of the parent `representation` node).

There are three kinds of selectors:

- **Static selector** is a string that selects a part of the structure based on entity type. The supported static selectors are these:

  `"all", "polymer", "protein", "nucleic", "branched", "ligand", "ion", "water"`

- **Component expression** is an object that selects a set of atoms based on their properties like chain identifier, residue number, or type symbol. The type of a component expression object is:

  ```ts
  {
      label_entity_id?: str,    // Entity identifier
      label_asym_id?: str,      // Chain identifier in label_* numbering
      auth_asym_id?: str,       // Chain identifier in auth_* numbering
      label_seq_id?: int,       // Residue number in label_* numbering
      auth_seq_id?: int,        // Residue number in auth_* numbering
      pdbx_PDB_ins_code?: str,  // PDB insertion code
      beg_label_seq_id?: int,   // Minimum label_seq_id (inclusive), leave blank to start from the beginning of the chain
      end_label_seq_id?: int,   // Maximum label_seq_id (inclusive), leave blank to go to the end of the chain
      beg_auth_seq_id?: int,    // Minimum auth_seq_id (inclusive), leave blank to start from the beginning of the chain
      end_auth_seq_id?: int,    // Maximum auth_seq_id (inclusive), leave blank to go to the end of the chain
      label_atom_id?: str,      // Atom name like 'CA', 'N', 'O', in label_* numbering
      auth_atom_id?: str,       // Atom name like 'CA', 'N', 'O', in auth_* numbering
      type_symbol?: str,        // Element symbol like 'H', 'HE', 'LI', 'BE'
      atom_id?: int,            // Unique atom identifier (_atom_site.id)
      atom_index?: int,         // 0-based index of the atom in the source data
  }
  ```

  A component expression can include any combination of the fields. An expression with multiple fields selects atoms that fulfill all fields at the same time. Examples:

  ```ts
  // Select whole chain A
  selector: { label_asym_id: 'A' }

  // Select residues 100 to 200 (inclusive) in chain B
  selector: { label_asym_id: 'B', beg_label_seq_id: 100, end_label_seq_id: 200 }

  // Select C-alpha atoms in residue 100 (using auth_* numbering) of any chain
  selector: { auth_seq_id: 100, type_symbol: 'C', auth_atom_id: 'CA' }
  ```

- **Union component expression** is an array of simple component expressions. A union component expression is interpreted as set union, i.e. it selects all atoms that fulfill at least one of the expressions in the array. Example:

  ```ts
  // Select chains A, B, and C
  selector: [{ label_asym_id: 'A' }, { label_asym_id: 'B' }, { label_asym_id: 'C' }]

  // Select residues up to 100 (inclusive) in chain A plus all magnesium atoms
  selector: [{ label_asym_id: 'A', end_label_seq_id: 100 }, { type_symbol: 'MG' }]
  ```

An alternative to using selectors is using [MVS annotations](./annotations.md). This means defining the selections in a separate file and referencing them from the MVS file.
