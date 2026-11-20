/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The kinds of what is not a linear polymer of letters, and the complexes over them (_spec/README.md, Language,
// Molecules and complexes): `molecule`, `compounds` and `compound`, `complex`, `unit`, and the params the other kinds
// take for them (a sequence's `modified`, a link's columns of `pairs`). Spread into `SeqLanguage`'s kinds; what the
// description cannot say is `check-molecules.ts`.

/** The systems a data reference takes: a sequence, a molecule's residues, a unit of a complex. */
export const SystemKinds = ['sequence', 'molecule', 'unit'] as const;

/** An open lowercase word (an identifier's namespace, a molecule's class, a section's class). */
export const OpenWord = '^[a-z][a-z0-9_-]*$';

/** Identifiers of a molecule or a compound: pointers matched as strings, never as chemistry. */
const identifiers = {
  type: 'table',
  columns: {
    namespace: {
      type: 'column',
      of: 'category',
      required: true,
      doc: 'Each identifier’s namespace, an open lowercase word (`ccd`, `prd`, `smiles`, `inchi`, `inchikey`, `chembl`, `pubchem`, `glytoucan`, `wurcs`, `linucs`, a host’s own).',
    },
    value: {
      type: 'list',
      of: { type: 'string', doc: 'A value.' },
      required: true,
      doc: 'Each identifier’s value, matched exactly (`ATP`, a SMILES).',
    },
  },
  doc: 'Its identifiers, a namespace and a value each: two are the same compound when one identifier of each is equal, namespace and value (a host’s resolver adds equivalents).',
} as const;

/** A sequence's non-standard positions (`sequence.modified`). */
export const ModifiedParam = {
  type: 'table',
  columns: {
    position: {
      type: 'column',
      of: 'int32',
      required: true,
      doc: 'The position, in order; several rows at one position are alternatives (microheterogeneity), the first giving the letter.',
    },
    component: {
      type: 'column',
      of: 'category',
      required: true,
      doc: 'The chemical component there (a CCD code: `MSE`, `PSU`).',
    },
    length: {
      type: 'column',
      of: 'int32',
      doc: 'The positions the component covers (1 when absent): one standing for several parents (CRO for T, Y, G), their letters.',
    },
  },
  doc: 'The positions whose component is not their letter’s standard one, sparse: the letter is the parent’s (MSE is `M`), `X` or `N` where there is none. A position not listed has its letter’s standard component.',
} as const;

/** The columns a link has beyond a pair's (a `pairs` child of a complex or a molecule). */
export const LinkColumns = {
  a_unit: {
    type: 'column',
    of: 'category',
    doc: 'Under a complex: the unit (its id) each first end is on.',
  },
  b_unit: {
    type: 'column',
    of: 'category',
    doc: 'Under a complex: the unit (its id) each second end is on.',
  },
  a_atom: { type: 'column', of: 'category', missing: true, doc: 'The atom of each first end (`SG`, `ND2`).' },
  b_atom: { type: 'column', of: 'category', missing: true, doc: 'The atom of each second end.' },
  order: { type: 'column', of: 'category', missing: true, doc: 'Each bond’s order (`sing`, `doub`).' },
  occupancy: {
    type: 'column',
    of: 'float32',
    doc: 'The share of copies each link is in, 0 to 1 (NaN: 1): a glycoform, a stochastic conjugation.',
  },
} as const;

export const MoleculeKinds = {
  molecule: {
    doc: 'A molecule that is not a linear polymer of letters (a small molecule, an ion, a glycan, a drug, a linker-payload): its positions are its residues, a component each; its atoms optional, a system of their own (`<molecule>/(atoms)`).',
    parents: ['data'],
    ids: 'required',
    params: {
      residues: {
        type: 'column',
        of: 'category',
        doc: 'Its residues inline, a component code each, a position each (`[ATP]`, `[NAG NAG BMA MAN MAN]`); instead of `source`.',
      },
      source: { type: 'ref', kinds: ['source'], doc: 'The source that holds it (instead of `residues`).' },
      atoms: {
        type: 'table',
        columns: {
          residue: {
            type: 'column',
            of: 'int32',
            required: true,
            doc: 'The residue each atom is of, non-decreasing: a residue’s atoms are a run.',
          },
          name: { type: 'column', of: 'category', required: true, doc: 'Each atom’s name (`O1G`).' },
          element: { type: 'column', of: 'category', required: true, doc: 'Each atom’s element; heavy atoms only.' },
        },
        doc: 'Its heavy atoms, optional, a position each of its atoms system: never hydrogens (EC-64), one atom a name.',
      },
      class: {
        type: 'string',
        pattern: OpenWord,
        doc: 'What it is, an open word: `branched`, `ligand`, `ion`, `water`, or any (`ligand` when absent).',
      },
      name: {
        type: 'string',
        doc: 'Its name for people (a chip’s label; else its one component’s code, else its id).',
      },
      identifiers,
    },
  },
  compounds: {
    doc: 'The container of the compounds a spec describes; its id is `compounds`: reference data, no system.',
    parents: ['data'],
    ids: 'required',
  },
  compound: {
    doc: 'A chemical component described (its id the code, a CCD code by convention): a name, its parents, identifiers. A code needs no node; a writer’s own monomer is a node with a SMILES.',
    parents: ['compounds'],
    ids: 'required',
    params: {
      name: { type: 'string', doc: 'Its name (`pseudouridine-5’-monophosphate`).' },
      parents: {
        type: 'list',
        of: { type: 'string', doc: 'A parent’s code.' },
        doc: 'The standard components it stands for (CRO: `THR`, `TYR`, `GLY`).',
      },
      class: {
        type: 'string',
        doc: 'What it is, the CCD’s `type` word (`L-PEPTIDE LINKING`, `NON-POLYMER`, `SACCHARIDE`).',
      },
      identifiers,
    },
  },
  complex: {
    doc: 'Units shown together, in order, with the links between them (a `pairs` child): a structure’s chains, or a construct with no structure behind it. It adds no arithmetic.',
    parents: ['data'],
    ids: 'required',
    params: {
      name: { type: 'string', doc: 'Its name for people (its id when absent).' },
    },
  },
  unit: {
    doc: 'A unit of a complex: an instance of a sequence or a molecule, a coordinate system of the molecule’s positions over `range`, its position i the molecule’s `start + i`. Units of one complex with equal `of` and `range` are copies.',
    parents: ['complex'],
    ids: 'required',
    params: {
      of: {
        type: 'ref',
        kinds: ['sequence', 'molecule'],
        required: true,
        doc: 'The sequence or molecule it instantiates.',
      },
      range: {
        type: 'list',
        of: { type: 'integer', min: 0, doc: 'A position.' },
        'min-length': 2,
        'max-length': 2,
        doc: 'The positions of its molecule it has, `[start, end)`, never wrapping (EC-13); all when absent.',
      },
      observed: {
        type: 'table',
        columns: {
          start: { type: 'column', of: 'int32', required: true, doc: 'The first position of each modeled run.' },
          end: { type: 'column', of: 'int32', required: true, doc: 'The position after its last.' },
        },
        doc: 'The runs of its own positions that were modeled, in order; all when absent (the rest drawn gray).',
      },
      chain: { type: 'string', doc: 'The structure’s chain (mmCIF’s `label_asym_id`).' },
      auth_chain: { type: 'string', doc: 'The author’s chain (`auth_asym_id`).' },
      operator: { type: 'string', doc: 'The instance’s operator (MolViewSpec’s `instance_id`: `1_555`, `ASM-X0-1`).' },
      label: { type: 'string', doc: 'Its label for people (its id when absent).' },
      group: { type: 'string', doc: 'The group of units it is one of (a BIRD/PRD instance).' },
    },
  },
} as const;

/** A block's form, a section's classes, a letters layer's `modified`: the params the other kinds take for molecules. */
export const BlockForm = {
  type: 'string',
  enum: ['whole', 'residues', 'atoms'],
  doc: 'How it is drawn: `whole` (one chip over its range, labeled), `residues` (a position a residue or letter) or `atoms` (a position a heavy atom; a molecule with atoms only). `whole` for a molecule or a unit of one, `residues` for a sequence or a unit of one, when absent.',
} as const;

export const SectionClasses = {
  type: 'list',
  of: { type: 'string', pattern: OpenWord, doc: 'A class: `polymer` (a sequence’s units) or a molecule’s class.' },
  doc: 'For a section of a complex: the classes of units shown, in order (polymers in the complex’s order, then each class); a class not listed is left out. `[polymer, branched, ligand, ion]` when absent (water left out).',
} as const;

export const LettersModified = {
  type: 'string',
  enum: ['mark', 'code', 'tint', 'none'],
  doc: 'For `letters`: how a position whose component is not standard is drawn, never widening a cell: `mark` (when absent: a dot at the cell’s top right; a bar over a component of several positions), `code` (its code where the cell has room, the mark otherwise), `tint` (the cell outlined) or `none`.',
} as const;
