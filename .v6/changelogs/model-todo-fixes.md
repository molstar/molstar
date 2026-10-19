## Model: fixes

- Highlighting and selecting cross-link restraints works: `CrossLinkRestraint.isLoci` tested the tag of interactions
- `Bond.areLocationsEqual` compares structures, not their labels, and `Bond.remapLoci` places the remapped bonds on
  the new structure (they kept the old one)
- `Model.getInvertedAtomSourceIndex` reports `isIdentity` for file-ordered models (it was always `false`)
- `findBonds` clears the covalent flag of metal bonds from the component dictionary instead of toggling it (the same
  result for dictionary bonds, which are always covalent); `StructureSequence.fromAtomicHierarchy` loses a loop that
  never ran (an entity's sequence stays its first chain's)
- The empty module `@molstar/model/model/sequence/alignment` is removed; alignment is `model/sequence/alignment/*`
