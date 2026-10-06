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
import { parseGRO } from '@molstar/io/reader/gro/parser';
import { trajectoryFromGRO } from '@molstar/model/formats/structure/gro';
import { trajectoryProps } from './helpers.js';
import { type TrajectoryFormatProvider, directTrajectory, defaultVisuals } from './provider.js';
import { TrajectoryFormatCategory } from './category.js';

export { TrajectoryFromGRO };
type TrajectoryFromGRO = typeof TrajectoryFromGRO;
const TrajectoryFromGRO = PluginStateTransform.BuiltIn({
  name: 'trajectory-from-gro',
  display: { name: 'Parse GRO', description: 'Parse GRO string and create trajectory.' },
  from: [SO.Data.String],
  to: SO.Molecule.Trajectory,
})({
  apply({ a }) {
    return Task.create('Parse GRO', async (ctx) => {
      const parsed = await parseGRO(a.data).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      const models = await trajectoryFromGRO(parsed.result).runInContext(ctx);
      const props = trajectoryProps(models);
      return new SO.Molecule.Trajectory(models, props);
    });
  },
});

export const GroProvider: TrajectoryFormatProvider = {
  label: 'GRO',
  description: 'GRO',
  category: TrajectoryFormatCategory,
  stringExtensions: ['gro'],
  binaryExtensions: [],
  ...directTrajectory(TrajectoryFromGRO),
  visuals: defaultVisuals,
};
