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
import { parseSdf } from '@molstar/io/reader/sdf/parser';
import { type Model, ArrayTrajectory } from '@molstar/model/model/structure';
import { trajectoryFromSdf } from '@molstar/model/formats/structure/sdf';
import { trajectoryProps } from './helpers.js';
import { TrajectoryFormatProvider, directTrajectory, defaultVisuals } from './provider.js';
import { TrajectoryFormatCategory } from './category.js';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';

export { TrajectoryFromSDF };
type TrajectoryFromSDF = typeof TrajectoryFromSDF;
const TrajectoryFromSDF = PluginStateTransform.BuiltIn({
  name: 'trajectory-from-sdf',
  display: { name: 'Parse SDF', description: 'Parse SDF string and create trajectory.' },
  from: [SO.Data.String],
  to: SO.Molecule.Trajectory,
})({
  apply({ a }) {
    return Task.create('Parse SDF', async (ctx) => {
      const parsed = await parseSdf(a.data).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);

      const models: Model[] = [];

      for (const compound of parsed.result.compounds) {
        const traj = await trajectoryFromSdf(compound).runInContext(ctx);
        for (let i = 0; i < traj.frameCount; i++) {
          models.push(await Task.resolveInContext(traj.getFrameAtIndex(i), ctx));
        }
      }

      const traj = new ArrayTrajectory(models);

      const props = trajectoryProps(traj);
      return new SO.Molecule.Trajectory(traj, props);
    });
  },
});

export const SdfProvider = TrajectoryFormatProvider({
  name: 'sdf',
  label: 'SDF',
  description: 'SDF',
  category: TrajectoryFormatCategory,
  stringExtensions: ['sdf', 'sd'],
  ...directTrajectory(TrajectoryFromSDF),
  visuals: defaultVisuals,
});

/** The Sdf data format. */
export const Sdf: PluginRegistryEntry = {
  formats: [SdfProvider],
};
