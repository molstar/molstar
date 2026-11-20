## Seq* (next): language 0.5.0: molecules and complexes

- New kinds: `molecule`, `compounds`, `complex` and `unit`, links as pairs between units, `modified` residues on a
  sequence, a block's `form` and a section's `classes`. A unit's positions count from 0, as every system's; its
  declared numberings (`auth`) name the residues.
- A source's `ref` names a host's data beside or instead of its `uri`; `values.field` picks a column.
- `@molstar/core/util/sequence-constants`: the standard component of each letter by alphabet and its inverse
  (`standardComponent`, `standardLetter`), with the amino-acid names table moved there.
