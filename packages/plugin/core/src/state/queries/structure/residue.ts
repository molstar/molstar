/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MolScriptBuilder as MS } from '@molstar/model/script/language/builder';
import { StructureSelectionCategory, StructureSelectionQuery } from './query.js';
import { ResidueQuery } from './dynamic.js';

export const nonStandardPolymer = StructureSelectionQuery(
  'Non-standard Residues in Polymers',
  MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
      'entity-test': MS.core.rel.eq([MS.ammp('entityType'), 'polymer']),
      'chain-test': MS.core.rel.eq([MS.ammp('objectPrimitive'), 'atomistic']),
      'residue-test': MS.ammp('isNonStandard'),
    }),
  ]),
  { category: StructureSelectionCategory.Residue },
);

export const ring = StructureSelectionQuery(
  'Rings in Residues',
  MS.struct.modifier.union([MS.struct.generator.rings()]),
  {
    category: StructureSelectionCategory.Residue,
  },
);

export const aromaticRing = StructureSelectionQuery(
  'Aromatic Rings in Residues',
  MS.struct.modifier.union([MS.struct.generator.rings({ 'only-aromatic': true })]),
  { category: StructureSelectionCategory.Residue },
);

const StandardAminoAcids = [
  [['HIS'], 'Histidine'],
  [['ARG'], 'Arginine'],
  [['LYS'], 'Lysine'],
  [['ILE'], 'Isoleucine'],
  [['PHE'], 'Phenylalanine'],
  [['LEU'], 'Leucine'],
  [['TRP'], 'Tryptophan'],
  [['ALA'], 'Alanine'],
  [['MET'], 'Methionine'],
  [['PRO'], 'Proline'],
  [['CYS'], 'Cysteine'],
  [['ASN'], 'Asparagine'],
  [['VAL'], 'Valine'],
  [['GLY'], 'Glycine'],
  [['SER'], 'Serine'],
  [['GLN'], 'Glutamine'],
  [['TYR'], 'Tyrosine'],
  [['ASP'], 'Aspartic Acid'],
  [['GLU'], 'Glutamic Acid'],
  [['THR'], 'Threonine'],
  [['SEC'], 'Selenocysteine'],
  [['PYL'], 'Pyrrolysine'],
  [['UNK'], 'Unknown'],
].sort((a, b) => (a[1] < b[1] ? -1 : a[1] > b[1] ? 1 : 0)) as [string[], string][];

const StandardNucleicBases = [
  [['A', 'DA'], 'Adenosine'],
  [['C', 'DC'], 'Cytidine'],
  [['T', 'DT'], 'Thymidine'],
  [['G', 'DG'], 'Guanosine'],
  [['I', 'DI'], 'Inosine'],
  [['U', 'DU'], 'Uridine'],
  [['N', 'DN'], 'Unknown'],
].sort((a, b) => (a[1] < b[1] ? -1 : a[1] > b[1] ? 1 : 0)) as [string[], string][];

export const AminoAcidSelectionQueries = StandardAminoAcids.map((v) =>
  ResidueQuery(v, StructureSelectionCategory.AminoAcid),
);

export const NucleicBaseSelectionQueries = StandardNucleicBases.map((v) =>
  ResidueQuery(v, StructureSelectionCategory.NucleicBase),
);

export const ResidueSelectionQueries = [
  nonStandardPolymer,
  ring,
  aromaticRing,
  ...AminoAcidSelectionQueries,
  ...NucleicBaseSelectionQueries,
];
