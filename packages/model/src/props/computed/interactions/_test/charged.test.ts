/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { DefaultAssetManager } from '@molstar/core/util/assets';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure, Unit } from '@molstar/model/model/structure';
import { ValenceModelProvider } from '../../valence-model.js';
import { NegativChargeProvider } from '../charged.js';
import { FeatureType } from '../common.js';
import { FeaturesBuilder } from '../features.js';

/** Acetamide as a ligand: an amide N (named N1: an N named N is taken for backbone) and no acid group */
const acetamide = `data_T
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
HETATM 1 C C1 ACM A 1 . 0 0 0 1 A 1
HETATM 2 C C2 ACM A 1 . 1.52 0 0 1 A 1
HETATM 3 O O1 ACM A 1 . 2.14 1.07 0 1 A 1
HETATM 4 N N1 ACM A 1 . 2.18 -1.17 0 1 A 1
`;

describe('negative charge features', () => {
  it('a ligand N is not a negative charge (perf-audit-2 06/B11)', async () => {
    const parsed = await CIF.parse(acetamide).run();
    if (parsed.isError) throw new Error(parsed.message);
    const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
    const structure = Structure.ofModel(trajectory.representative);
    await ValenceModelProvider.attach(
      { runtime: SyncRuntimeContext, assetManager: new DefaultAssetManager() },
      structure,
    );

    const unit = structure.units[0] as Unit.Atomic;
    const builder = FeaturesBuilder.create();
    NegativChargeProvider.add(structure, unit, builder);
    const features = builder.getFeatures(unit.elements.length);

    let negative = 0;
    for (let i = 0; i < features.count; i++) {
      if (features.types[i] === FeatureType.NegativeCharge) negative++;
    }
    expect(negative).toBe(0);
  });
});
