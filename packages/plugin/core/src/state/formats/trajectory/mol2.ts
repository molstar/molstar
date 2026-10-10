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
import { parseMol2 } from '@molstar/io/reader/mol2/parser';
import { trajectoryFromMol2 } from '@molstar/model/formats/structure/mol2';
import { trajectoryProps } from './helpers.js';
import { TrajectoryFormatProvider, directTrajectory, defaultVisuals } from './provider.js';
import { TrajectoryFormatCategory } from './category.js';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';

export { TrajectoryFromMOL2 };
type TrajectoryFromMOL2 = typeof TrajectoryFromMOL2;
const TrajectoryFromMOL2 = PluginStateTransform.BuiltIn({
  name: 'trajectory-from-mol2',
  display: { name: 'Parse MOL2', description: 'Parse MOL2 string and create trajectory.' },
  from: [SO.Data.String],
  to: SO.Molecule.Trajectory,
})({
  apply({ a }) {
    return Task.create('Parse MOL2', async (ctx) => {
      const parsed = await parseMol2(a.data, a.label).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      const models = await trajectoryFromMol2(parsed.result).runInContext(ctx);
      const props = trajectoryProps(models);
      return new SO.Molecule.Trajectory(models, props);
    });
  },
});

export const Mol2Provider = TrajectoryFormatProvider({
  name: 'mol2',
  label: 'MOL2',
  description: 'MOL2',
  category: TrajectoryFormatCategory,
  stringExtensions: ['mol2'],
  ...directTrajectory(TrajectoryFromMOL2),
  visuals: defaultVisuals,
});

/** The Mol2 data format. */
export const Mol2: PluginRegistryEntry = {
  formats: [Mol2Provider],
};
