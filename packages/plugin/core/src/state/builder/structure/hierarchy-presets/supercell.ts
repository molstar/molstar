/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '@molstar/core/math/linear-algebra';
import { Model } from '@molstar/model/model/structure';
import { applyCrystalSymmetry, CrystalSymmetryParams } from './crystal-symmetry.js';
import { TrajectoryHierarchyPresetProvider } from './types.js';

export const SupercellHierarchyPreset = TrajectoryHierarchyPresetProvider({
  id: 'preset-trajectory-supercell',
  alias: 'supercell',
  display: {
    name: 'Super Cell',
    group: 'Preset',
    description: 'Shows the super cell, i.e. the central unit cell and all adjacent unit cells.',
  },
  isApplicable: (o) => {
    return Model.hasCrystalSymmetry(o.data.representative);
  },
  params: CrystalSymmetryParams,
  async apply(trajectory, params, plugin) {
    return await applyCrystalSymmetry(
      { ijkMin: Vec3.create(-1, -1, -1), ijkMax: Vec3.create(1, 1, 1), theme: 'operator-hkl' },
      trajectory,
      params,
      plugin,
    );
  },
});
