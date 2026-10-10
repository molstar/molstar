/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MolScriptBuilder as MS } from '@molstar/query-language/language/builder';
import type { CustomProperty } from '@molstar/model/props/common/custom-property';
import type { Structure } from '@molstar/model/model/structure';
import {
  NucleicBackboneAtoms,
  ProteinBackboneAtoms,
  SecondaryStructureType,
} from '@molstar/model/model/structure/model/types';
import { SecondaryStructureProvider } from '@molstar/model/props/computed/secondary-structure';
import { SetUtils } from '@molstar/core/util/set';
import { StructureSelectionCategory, StructureSelectionQuery } from './query.js';
import { nonPolymerResidueTest, nucleiEntityTest, proteinEntityTest } from './common.js';

export const trace = StructureSelectionQuery(
  'Trace',
  MS.struct.modifier.union([
    MS.struct.combinator.merge([
      MS.struct.modifier.union([
        MS.struct.generator.atomGroups({
          'entity-test': MS.core.rel.eq([MS.ammp('entityType'), 'polymer']),
          'chain-test': MS.core.set.has([MS.set('sphere', 'gaussian'), MS.ammp('objectPrimitive')]),
        }),
      ]),
      MS.struct.modifier.union([
        MS.struct.generator.atomGroups({
          'entity-test': MS.core.rel.eq([MS.ammp('entityType'), 'polymer']),
          'chain-test': MS.core.rel.eq([MS.ammp('objectPrimitive'), 'atomistic']),
          'atom-test': MS.core.set.has([MS.set('CA', 'P'), MS.ammp('label_atom_id')]),
        }),
      ]),
    ]),
  ]),
  { category: StructureSelectionCategory.Structure },
);

// TODO maybe pre-calculate backbone atom properties
export const backbone = StructureSelectionQuery(
  'Backbone',
  MS.struct.modifier.union([
    MS.struct.combinator.merge([
      MS.struct.modifier.union([
        MS.struct.generator.atomGroups({
          'entity-test': proteinEntityTest,
          'chain-test': MS.core.rel.eq([MS.ammp('objectPrimitive'), 'atomistic']),
          'residue-test': MS.core.logic.not([nonPolymerResidueTest]),
          'atom-test': MS.core.set.has([MS.set(...SetUtils.toArray(ProteinBackboneAtoms)), MS.ammp('label_atom_id')]),
        }),
      ]),
      MS.struct.modifier.union([
        MS.struct.generator.atomGroups({
          'entity-test': nucleiEntityTest,
          'chain-test': MS.core.rel.eq([MS.ammp('objectPrimitive'), 'atomistic']),
          'residue-test': MS.core.logic.not([nonPolymerResidueTest]),
          'atom-test': MS.core.set.has([MS.set(...SetUtils.toArray(NucleicBackboneAtoms)), MS.ammp('label_atom_id')]),
        }),
      ]),
    ]),
  ]),
  { category: StructureSelectionCategory.Structure },
);

// TODO maybe pre-calculate sidechain atom property
export const sidechain = StructureSelectionQuery(
  'Sidechain',
  MS.struct.modifier.union([
    MS.struct.combinator.merge([
      MS.struct.modifier.union([
        MS.struct.generator.atomGroups({
          'entity-test': proteinEntityTest,
          'chain-test': MS.core.rel.eq([MS.ammp('objectPrimitive'), 'atomistic']),
          'residue-test': MS.core.logic.not([nonPolymerResidueTest]),
          'atom-test': MS.core.logic.or([
            MS.core.logic.not([
              MS.core.set.has([MS.set(...SetUtils.toArray(ProteinBackboneAtoms)), MS.ammp('label_atom_id')]),
            ]),
          ]),
        }),
      ]),
      MS.struct.modifier.union([
        MS.struct.generator.atomGroups({
          'entity-test': nucleiEntityTest,
          'chain-test': MS.core.rel.eq([MS.ammp('objectPrimitive'), 'atomistic']),
          'residue-test': MS.core.logic.not([nonPolymerResidueTest]),
          'atom-test': MS.core.logic.or([
            MS.core.logic.not([
              MS.core.set.has([MS.set(...SetUtils.toArray(NucleicBackboneAtoms)), MS.ammp('label_atom_id')]),
            ]),
          ]),
        }),
      ]),
    ]),
  ]),
  { category: StructureSelectionCategory.Structure },
);

// TODO maybe pre-calculate sidechain atom property
export const sidechainWithTrace = StructureSelectionQuery(
  'Sidechain with Trace',
  MS.struct.modifier.union([
    MS.struct.combinator.merge([
      MS.struct.modifier.union([
        MS.struct.generator.atomGroups({
          'entity-test': proteinEntityTest,
          'chain-test': MS.core.rel.eq([MS.ammp('objectPrimitive'), 'atomistic']),
          'residue-test': MS.core.logic.not([nonPolymerResidueTest]),
          'atom-test': MS.core.logic.or([
            MS.core.logic.not([
              MS.core.set.has([MS.set(...SetUtils.toArray(ProteinBackboneAtoms)), MS.ammp('label_atom_id')]),
            ]),
            MS.core.rel.eq([MS.ammp('label_atom_id'), 'CA']),
            MS.core.logic.and([
              MS.core.rel.eq([MS.ammp('auth_comp_id'), 'PRO']),
              MS.core.rel.eq([MS.ammp('label_atom_id'), 'N']),
            ]),
          ]),
        }),
      ]),
      MS.struct.modifier.union([
        MS.struct.generator.atomGroups({
          'entity-test': nucleiEntityTest,
          'chain-test': MS.core.rel.eq([MS.ammp('objectPrimitive'), 'atomistic']),
          'residue-test': MS.core.logic.not([nonPolymerResidueTest]),
          'atom-test': MS.core.logic.or([
            MS.core.logic.not([
              MS.core.set.has([MS.set(...SetUtils.toArray(NucleicBackboneAtoms)), MS.ammp('label_atom_id')]),
            ]),
            MS.core.rel.eq([MS.ammp('label_atom_id'), 'P']),
          ]),
        }),
      ]),
    ]),
  ]),
  { category: StructureSelectionCategory.Structure },
);

export const helix = StructureSelectionQuery(
  'Helix',
  MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
      'entity-test': proteinEntityTest,
      'residue-test': MS.core.flags.hasAny([
        MS.ammp('secondaryStructureFlags'),
        MS.core.type.bitflags([SecondaryStructureType.Flag.Helix]),
      ]),
    }),
  ]),
  {
    category: StructureSelectionCategory.Structure,
    ensureCustomProperties: (ctx: CustomProperty.Context, structure: Structure) => {
      return SecondaryStructureProvider.attach(ctx, structure);
    },
  },
);

export const beta = StructureSelectionQuery(
  'Beta Strand/Sheet',
  MS.struct.modifier.union([
    MS.struct.generator.atomGroups({
      'entity-test': proteinEntityTest,
      'residue-test': MS.core.flags.hasAny([
        MS.ammp('secondaryStructureFlags'),
        MS.core.type.bitflags([SecondaryStructureType.Flag.Beta]),
      ]),
    }),
  ]),
  {
    category: StructureSelectionCategory.Structure,
    ensureCustomProperties: (ctx: CustomProperty.Context, structure: Structure) => {
      return SecondaryStructureProvider.attach(ctx, structure);
    },
  },
);

export const StructurePropertySelectionQueries = [trace, backbone, sidechain, sidechainWithTrace, helix, beta];
