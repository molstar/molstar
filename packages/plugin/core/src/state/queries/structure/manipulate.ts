/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MolScriptBuilder as MS } from '@molstar/query-language/language/builder';
import { BondType } from '@molstar/model/model/structure/model/types';
import { StructureSelectionCategory, StructureSelectionQuery } from './query.js';

export const surroundings = StructureSelectionQuery(
  'Surrounding Residues (5 \u212B) of Selection',
  MS.struct.modifier.union([
    MS.struct.modifier.exceptBy({
      0: MS.struct.modifier.includeSurroundings({
        0: MS.internal.generator.current(),
        radius: 5,
        'as-whole-residues': true,
      }),
      by: MS.internal.generator.current(),
    }),
  ]),
  {
    description: 'Select residues within 5 \u212B of the current selection.',
    category: StructureSelectionCategory.Manipulate,
    referencesCurrent: true,
  },
);

export const surroundingLigands = StructureSelectionQuery(
  'Surrounding Ligands (5 \u212B) of Selection',
  MS.struct.modifier.union([
    MS.struct.modifier.surroundingLigands({
      0: MS.internal.generator.current(),
      radius: 5,
      'include-water': true,
    }),
  ]),
  {
    description: 'Select ligand components within 5 \u212B of the current selection.',
    category: StructureSelectionCategory.Manipulate,
    referencesCurrent: true,
  },
);

export const surroundingAtoms = StructureSelectionQuery(
  'Surrounding Atoms (5 \u212B) of Selection',
  MS.struct.modifier.union([
    MS.struct.modifier.exceptBy({
      0: MS.struct.modifier.includeSurroundings({
        0: MS.internal.generator.current(),
        radius: 5,
        'as-whole-residues': false,
      }),
      by: MS.internal.generator.current(),
    }),
  ]),
  {
    description: 'Select atoms within 5 \u212B of the current selection.',
    category: StructureSelectionCategory.Manipulate,
    referencesCurrent: true,
  },
);

export const complement = StructureSelectionQuery(
  'Inverse / Complement of Selection',
  MS.struct.modifier.union([
    MS.struct.modifier.exceptBy({
      0: MS.struct.generator.all(),
      by: MS.internal.generator.current(),
    }),
  ]),
  {
    description: 'Select everything not in the current selection.',
    category: StructureSelectionCategory.Manipulate,
    referencesCurrent: true,
  },
);

export const covalentlyBonded = StructureSelectionQuery(
  'Residues Covalently Bonded to Selection',
  MS.struct.modifier.union([
    MS.struct.modifier.includeConnected({
      0: MS.internal.generator.current(),
      'layer-count': 1,
      'as-whole-residues': true,
    }),
  ]),
  {
    description: 'Select residues covalently bonded to current selection.',
    category: StructureSelectionCategory.Manipulate,
    referencesCurrent: true,
  },
);

export const covalentlyBondedComponent = StructureSelectionQuery(
  'Covalently Bonded Component',
  MS.struct.modifier.union([
    MS.struct.modifier.includeConnected({
      0: MS.internal.generator.current(),
      'fixed-point': true,
    }),
  ]),
  {
    description: 'Select covalently bonded component based on current selection.',
    category: StructureSelectionCategory.Manipulate,
    referencesCurrent: true,
  },
);

export const covalentlyOrMetallicBonded = StructureSelectionQuery(
  'Residues with Cov. or Metallic Bond to Selection',
  MS.struct.modifier.union([
    MS.struct.modifier.includeConnected({
      0: MS.internal.generator.current(),
      'layer-count': 1,
      'as-whole-residues': true,
      'bond-test': MS.core.flags.hasAny([
        MS.struct.bondProperty.flags(),
        MS.core.type.bitflags([BondType.Flag.Covalent | BondType.Flag.MetallicCoordination]),
      ]),
    }),
  ]),
  {
    description: 'Select residues with covalent or metallic bond to current selection.',
    category: StructureSelectionCategory.Manipulate,
    referencesCurrent: true,
  },
);

export const wholeResidues = StructureSelectionQuery(
  'Whole Residues of Selection',
  MS.struct.modifier.union([
    MS.struct.modifier.wholeResidues({
      0: MS.internal.generator.current(),
    }),
  ]),
  {
    description: 'Expand current selection to whole residues.',
    category: StructureSelectionCategory.Manipulate,
    referencesCurrent: true,
  },
);

export const ManipulateSelectionQueries = [
  surroundings,
  surroundingLigands,
  surroundingAtoms,
  complement,
  covalentlyBonded,
  covalentlyOrMetallicBonded,
  covalentlyBondedComponent,
  wholeResidues,
];
