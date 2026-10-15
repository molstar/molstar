/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MolScriptBuilder as MS } from '@molstar/query-language/language/builder';
import { BondType, CommonProteinCaps, PolymerNames } from '@molstar/model/model/structure/model/types';
import { SetUtils } from '@molstar/core/util/set';
import { StructureSelectionCategory, StructureSelectionQuery } from './query.js';
import { nonPolymerResidueTest, nucleiEntityTest, proteinEntityTest } from './common.js';

export const polymer = StructureSelectionQuery(
  'Polymer',
  MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
      'entity-test': MS.core.logic.and([
        MS.core.rel.eq([MS.ammp('entityType'), 'polymer']),
        MS.core.str.match([
          MS.re('(polypeptide|cyclic-pseudo-peptide|peptide-like|nucleotide|peptide nucleic acid)', 'i'),
          MS.ammp('entitySubtype'),
        ]),
      ]),
    }),
  ]),
  { category: StructureSelectionCategory.Type },
);

export const protein = StructureSelectionQuery(
  'Protein',
  MS.struct.modifier.union([MS.struct.generator.atomGroups({ 'entity-test': proteinEntityTest })]),
  { category: StructureSelectionCategory.Type },
);

export const nucleic = StructureSelectionQuery(
  'Nucleic',
  MS.struct.modifier.union([MS.struct.generator.atomGroups({ 'entity-test': nucleiEntityTest })]),
  { category: StructureSelectionCategory.Type },
);

export const water = StructureSelectionQuery(
  'Water',
  MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
      'entity-test': MS.core.rel.eq([MS.ammp('entityType'), 'water']),
    }),
  ]),
  { category: StructureSelectionCategory.Type },
);

export const ion = StructureSelectionQuery(
  'Ion',
  MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
      'entity-test': MS.core.rel.eq([MS.ammp('entitySubtype'), 'ion']),
    }),
  ]),
  { category: StructureSelectionCategory.Type },
);

export const lipid = StructureSelectionQuery(
  'Lipid',
  MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
      'entity-test': MS.core.rel.eq([MS.ammp('entitySubtype'), 'lipid']),
    }),
  ]),
  { category: StructureSelectionCategory.Type },
);

export const branched = StructureSelectionQuery(
  'Carbohydrate',
  MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
      'entity-test': MS.core.logic.or([
        MS.core.rel.eq([MS.ammp('entityType'), 'branched']),
        MS.core.logic.and([
          MS.core.rel.eq([MS.ammp('entityType'), 'non-polymer']),
          MS.core.str.match([MS.re('oligosaccharide', 'i'), MS.ammp('entitySubtype')]),
        ]),
      ]),
    }),
  ]),
  { category: StructureSelectionCategory.Type },
);

export const branchedPlusConnected = StructureSelectionQuery(
  'Carbohydrate with Connected',
  MS.struct.modifier.union([
    MS.struct.modifier.includeConnected({
      0: branched.expression,
      'layer-count': 1,
      'as-whole-residues': true,
    }),
  ]),
  { category: StructureSelectionCategory.Internal, isHidden: true },
);

export const branchedConnectedOnly = StructureSelectionQuery(
  'Connected to Carbohydrate',
  MS.struct.modifier.union([
    MS.struct.modifier.exceptBy({
      0: branchedPlusConnected.expression,
      by: branched.expression,
    }),
  ]),
  { category: StructureSelectionCategory.Internal, isHidden: true },
);

export const ligand = StructureSelectionQuery(
  'Ligand',
  MS.struct.modifier.union([
    MS.struct.modifier.exceptBy({
      0: MS.struct.modifier.union([
        MS.struct.combinator.merge([
          MS.struct.modifier.union([
            MS.struct.generator.atomGroups({
              'entity-test': MS.core.logic.and([
                MS.core.logic.or([
                  MS.core.rel.eq([MS.ammp('entityType'), 'non-polymer']),
                  MS.core.rel.neq([MS.ammp('entityPrdId'), '']),
                ]),
                MS.core.logic.not([
                  MS.core.str.match([MS.re('(oligosaccharide|lipid|ion)', 'i'), MS.ammp('entitySubtype')]),
                ]),
              ]),
              'chain-test': MS.core.rel.eq([MS.ammp('objectPrimitive'), 'atomistic']),
              'residue-test': MS.core.logic.not([
                MS.core.str.match([MS.re('saccharide', 'i'), MS.ammp('chemCompType')]),
              ]),
            }),
          ]),
          MS.struct.modifier.union([
            MS.struct.generator.atomGroups({
              'entity-test': MS.core.rel.eq([MS.ammp('entityType'), 'polymer']),
              'chain-test': MS.core.rel.eq([MS.ammp('objectPrimitive'), 'atomistic']),
              'residue-test': nonPolymerResidueTest,
            }),
          ]),
        ]),
      ]),
      by: MS.struct.combinator.merge([
        MS.struct.modifier.union([
          MS.struct.generator.atomGroups({
            'entity-test': MS.core.rel.eq([MS.ammp('entityType'), 'polymer']),
            'chain-test': MS.core.rel.eq([MS.ammp('objectPrimitive'), 'atomistic']),
            'residue-test': MS.core.set.has([MS.set(...SetUtils.toArray(PolymerNames)), MS.ammp('label_comp_id')]),
          }),
        ]),
        MS.struct.generator.atomGroups({
          'chain-test': MS.core.rel.eq([MS.ammp('objectPrimitive'), 'atomistic']),
          'residue-test': MS.core.set.has([MS.set(...SetUtils.toArray(CommonProteinCaps)), MS.ammp('label_comp_id')]),
        }),
      ]),
    }),
  ]),
  { category: StructureSelectionCategory.Type },
);

// don't include branched entities as they have their own link representation
export const ligandPlusConnected = StructureSelectionQuery(
  'Ligand with Connected',
  MS.struct.modifier.union([
    MS.struct.modifier.exceptBy({
      0: MS.struct.modifier.union([
        MS.struct.modifier.includeConnected({
          0: ligand.expression,
          'layer-count': 1,
          'as-whole-residues': true,
          'bond-test': MS.core.flags.hasAny([
            MS.struct.bondProperty.flags(),
            MS.core.type.bitflags([BondType.Flag.Covalent | BondType.Flag.MetallicCoordination]),
          ]),
        }),
      ]),
      by: branched.expression,
    }),
  ]),
  { category: StructureSelectionCategory.Internal, isHidden: true },
);

export const ligandConnectedOnly = StructureSelectionQuery(
  'Connected to Ligand',
  MS.struct.modifier.union([
    MS.struct.modifier.exceptBy({
      0: ligandPlusConnected.expression,
      by: ligand.expression,
    }),
  ]),
  { category: StructureSelectionCategory.Internal, isHidden: true },
);

// residues connected to ligands or branched entities
export const connectedOnly = StructureSelectionQuery(
  'Connected to Ligand or Carbohydrate',
  MS.struct.modifier.union([
    MS.struct.combinator.merge([branchedConnectedOnly.expression, ligandConnectedOnly.expression]),
  ]),
  { category: StructureSelectionCategory.Internal, isHidden: true },
);

export const coarse = StructureSelectionQuery(
  'Coarse Elements',
  MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
      'chain-test': MS.core.set.has([MS.set('sphere', 'gaussian'), MS.ammp('objectPrimitive')]),
    }),
  ]),
  { category: StructureSelectionCategory.Type },
);

export const TypeSelectionQueries = [
  polymer,
  protein,
  nucleic,
  water,
  ion,
  lipid,
  branched,
  branchedPlusConnected,
  branchedConnectedOnly,
  ligand,
  ligandPlusConnected,
  ligandConnectedOnly,
  connectedOnly,
  coarse,
];
