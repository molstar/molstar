/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MolScriptBuilder as MS } from '@molstar/query-language/language/builder';
import { BondType } from '@molstar/model/model/structure/model/types';
import { StructureSelectionCategory, StructureSelectionQuery } from './query.js';

export const disulfideBridges = StructureSelectionQuery(
  'Disulfide Bridges',
  MS.struct.modifier.union([
    MS.struct.combinator.merge([
      MS.struct.modifier.union([
        MS.struct.modifier.wholeResidues([
          MS.struct.filter.isConnectedTo({
            0: MS.struct.generator.atomGroups({
              'residue-test': MS.core.set.has([MS.set('CYS'), MS.ammp('auth_comp_id')]),
              'atom-test': MS.core.set.has([MS.set('SG'), MS.ammp('label_atom_id')]),
            }),
            target: MS.struct.generator.atomGroups({
              'residue-test': MS.core.set.has([MS.set('CYS'), MS.ammp('auth_comp_id')]),
              'atom-test': MS.core.set.has([MS.set('SG'), MS.ammp('label_atom_id')]),
            }),
            'bond-test': true,
          }),
        ]),
      ]),
      MS.struct.modifier.union([
        MS.struct.modifier.wholeResidues([
          MS.struct.modifier.union([
            MS.struct.generator.bondedAtomicPairs({
              0: MS.core.flags.hasAny([
                MS.struct.bondProperty.flags(),
                MS.core.type.bitflags([BondType.Flag.Disulfide]),
              ]),
            }),
          ]),
        ]),
      ]),
    ]),
  ]),
  { category: StructureSelectionCategory.Bond },
);

export const nosBridges = StructureSelectionQuery(
  'NOS Bridges',
  MS.struct.modifier.union([
    MS.struct.modifier.wholeResidues([
      MS.struct.filter.isConnectedTo({
        0: MS.struct.generator.atomGroups({
          'residue-test': MS.core.set.has([MS.set('CSO', 'LYS'), MS.ammp('auth_comp_id')]),
          'atom-test': MS.core.set.has([MS.set('OD', 'NZ'), MS.ammp('label_atom_id')]),
        }),
        target: MS.struct.generator.atomGroups({
          'residue-test': MS.core.set.has([MS.set('CSO', 'LYS'), MS.ammp('auth_comp_id')]),
          'atom-test': MS.core.set.has([MS.set('OD', 'NZ'), MS.ammp('label_atom_id')]),
        }),
        'bond-test': true,
      }),
    ]),
  ]),
  { category: StructureSelectionCategory.Bond },
);

export const BondSelectionQueries = [disulfideBridges, nosBridges];
