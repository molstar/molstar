## Seq* (next): molecules drawn

- A section `of` a complex shows a block per unit: polymers first, then branched, ligand and ion (water left out);
  copies merge into one `×N` block.
- A block's `form` is `whole` (one chip, a molecule's default), `residues` (by component code) or `atoms` (by element
  symbol); a `block` child naming a derived unit overrides that unit only, and `block.form` switches it in one undo
  step.
- Modified residues are marked (`letters.modifiedMark`: `mark`, `code`, `tint`, `none`), never widening a cell.
- Links are drawn as arcs between their units, and a `ranges` layer over them marks their ends. Hover names a residue's
  component, its alternatives, or an atom, and carries the molecule's identifiers; `links.hover` (opt-in) lights a
  molecule's linked sites, and a site's molecule.
