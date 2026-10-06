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
import { parseMol } from '@molstar/io/reader/mol/parser';
import { trajectoryFromMol } from '@molstar/model/formats/structure/mol';
import { trajectoryProps } from './helpers.js';
import { type TrajectoryFormatProvider, directTrajectory, defaultVisuals } from './provider.js';
import { TrajectoryFormatCategory } from './category.js';

export { TrajectoryFromMOL };
type TrajectoryFromMOL = typeof TrajectoryFromMOL;
const TrajectoryFromMOL = PluginStateTransform.BuiltIn({
  name: 'trajectory-from-mol',
  display: { name: 'Parse MOL', description: 'Parse MOL string and create trajectory.' },
  from: [SO.Data.String],
  to: SO.Molecule.Trajectory,
})({
  apply({ a }) {
    return Task.create('Parse MOL', async (ctx) => {
      const parsed = await parseMol(a.data).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      const models = await trajectoryFromMol(parsed.result).runInContext(ctx);
      const props = trajectoryProps(models);
      return new SO.Molecule.Trajectory(models, props);
    });
  },
});

export const MolProvider: TrajectoryFormatProvider = {
  label: 'MOL',
  description: 'MOL',
  category: TrajectoryFormatCategory,
  stringExtensions: ['mol'],
  ...directTrajectory(TrajectoryFromMOL),
  visuals: defaultVisuals,
};
