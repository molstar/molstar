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

export const UnitcellHierarchyPreset = TrajectoryHierarchyPresetProvider({
  id: 'preset-trajectory-unitcell',
  alias: 'unitcell',
  display: {
    name: 'Unit Cell',
    group: 'Preset',
    description: 'Shows the fully populated unit cell.',
  },
  isApplicable: (o) => {
    return Model.hasCrystalSymmetry(o.data.representative);
  },
  params: CrystalSymmetryParams,
  async apply(trajectory, params, plugin) {
    return await applyCrystalSymmetry(
      { ijkMin: Vec3.create(0, 0, 0), ijkMax: Vec3.create(0, 0, 0) },
      trajectory,
      params,
      plugin,
    );
  },
});
