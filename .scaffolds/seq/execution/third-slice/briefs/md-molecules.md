# Lane MD: the molecules design (a design only)

Follow `briefs/design-drafts.md` (Read first, What to write, Rules, Report) for one topic: **molecules**. Draft
`.scaffolds/seq/execution/third-slice/drafts/molecules.md`. No code; no other file edited. Model: strong (opus). Base:
`packages-seq` at `3675cf38c` (read files there; a worktree is not needed, nothing is committed by you).

## The topic, in the author's words (2026-11-19)

- "i also wonder we should have a better primitive than `(sequence atp :alphabet protein :letters X :name ATP)` for
  this? because these will be for example antibody-drug complexes, and in general, these could be pointers to smiles,
  carbohydrate polymers (branched entity in mmCIF), etc. and all these need to be mappable between each other."
- "do we support non-standard residues in our "sequence" nodes? if not, should include this" (we do not: a `sequence`
  is one letter a position from a fixed alphabet, `packages/seq/core/src/lang/check.ts` `Letters`).
- "chips is a prolly good way to do these. contacts should be optional and configurable by the app. remember we are
  building a library, not a "one off viewer"." (decided: a one-position block is a chip in every layout, layouts.md
  Blocks.)

## What it must cover

- **Small molecules**: by a chemical component id (CCD) or a SMILES (a pointer, not a structure the library parses);
  their forms on an axis (one position, the chip; or one position per heavy atom); their names and labels.
- **Branched entities** (mmCIF `pdbx_entity_branch`, `pdbx_entity_branch_link`, `pdbx_entity_branch_list`): a glycan
  as a tree of residues with linkages; how it sits on an axis (a chip, a small tree, a linear list) and is addressed.
- **Conjugates**: a drug on a residue (antibody–drug conjugates), a glycan on an Asn, a covalent ligand: the link
  between a molecule and a polymer's position, and how it is drawn and hovered.
- **Non-standard residues in polymers**: a position's chemical component (MSE, SEP, TPO, PTR; modified nucleotides as
  in 1EHZ: PSU, 2MG, H2U, OMC, OMG, YYG, 5MC, 7MG, 5MU, 1MA) with its parent's letter for drawing, search, alignment;
  how a track marks it, how a query finds it, how it maps to a structure.
- **The complex kind** the layouts design waits for (layouts.md Blocks: "a section `of` a complex ... when the complex
  kind exists"): units, copies, entities, and what a section of it lays out.
- **Maps between them**: molecule ↔ polymer position (a link, a contact), per-atom ↔ one-position systems (a
  many-to-one map), a component ↔ its parent letter; on the maps design (`.scaffolds/seq/design/maps.md`) and addressing
  (`addressing.md`, `locations-and-queries.md`).
- **Reading**: from mmCIF through `@molstar/io` (`chem_comp`, `entity_poly_seq`, the branch tables) and from the
  `seq` format; the base's per-position component codes and its `ligand` node and `ligand_form` read first
  (`packages/sequences/core/src/model/molecule.ts`, `packages/sequences/_spec/format.md`).
- **Host neutrality**: Seq\* holds pointers and labels, never chemistry; a SMILES is drawn by a host if at all.

## Constraints

- The language stays 0.4.0 in this slice; this design is built in the fourth slice (one bump then). Say what changes
  in the language, the store, the resolved model, and what the fourth slice's lanes would be.
- Precedents to read: the base (`packages/sequences`), Mol\*'s model (`packages/model`: `entity`, `chem_comp`,
  branched entities, `sequence` of a model), MolViewSpec (`packages/mvs`), PDBe's and RCSB's sequence views of
  glycans and modified residues, SnapGene/Benchling for construct features, and the mmCIF dictionary's branch tables.
- Every choice a host could make differently is the host's (the author: "a library, not a one-off viewer"); a design
  default never inferred from data.

## Report

As `design-drafts.md`'s Report, plus: the language changes as a list, and the fourth slice's lanes you would cut.
