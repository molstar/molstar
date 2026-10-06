/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { Task } from '@molstar/core/task';
import { parseXyz } from '@molstar/io/reader/xyz/parser';
import { trajectoryFromXyz } from '@molstar/model/formats/structure/xyz';
import { trajectoryProps } from './helpers.js';
import { TrajectoryFormatProvider, directTrajectory, defaultVisuals } from './provider.js';
import { TrajectoryFormatCategory } from './category.js';

export { TrajectoryFromXYZ };
type TrajectoryFromXYZ = typeof TrajectoryFromXYZ;
const TrajectoryFromXYZ = PluginStateTransform.BuiltIn({
  name: 'trajectory-from-xyz',
  display: { name: 'Parse XYZ', description: 'Parse XYZ string and create trajectory.' },
  from: [SO.Data.String],
  to: SO.Molecule.Trajectory,
})({
  apply({ a }) {
    return Task.create('Parse XYZ', async (ctx) => {
      const parsed = await parseXyz(a.data).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      const models = await trajectoryFromXyz(parsed.result).runInContext(ctx);
      const props = trajectoryProps(models);
      return new SO.Molecule.Trajectory(models, props);
    });
  },
});

export const XyzProvider = TrajectoryFormatProvider({
  name: 'xyz',
  label: 'XYZ',
  description: 'XYZ',
  category: TrajectoryFormatCategory,
  stringExtensions: ['xyz'],
  ...directTrajectory(TrajectoryFromXYZ),
  visuals: defaultVisuals,
});
