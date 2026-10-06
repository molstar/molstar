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
import { parseLammpsTrajectory } from '@molstar/io/reader/lammps/traj/parser';
import { coordinatesFromLammpsTrajectory } from '@molstar/model/formats/structure/lammps-trajectory';
import { DataFormatProvider } from '@molstar/plugin/state/formats/provider';
import { CoordinatesFormatCategory } from './category.js';

export { CoordinatesFromLammpstraj };
type CoordinatesFromLammpstraj = typeof CoordinatesFromLammpstraj;
const CoordinatesFromLammpstraj = PluginStateTransform.BuiltIn({
  name: 'coordinates-from-lammpstraj',
  display: { name: 'Parse LAMMPSTRAJ', description: 'Parse LAMMPSTRAJ data.' },
  from: [SO.Data.String],
  to: SO.Molecule.Coordinates,
})({
  apply({ a }) {
    return Task.create('Parse LAMMPSTRAJ', async (ctx) => {
      const parsed = await parseLammpsTrajectory(a.data).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      const coordinates = await coordinatesFromLammpsTrajectory(parsed.result).runInContext(ctx);
      return new SO.Molecule.Coordinates(coordinates, { label: a.label, description: 'Coordinates' });
    });
  },
});

export { LammpsTrajectoryProvider };
const LammpsTrajectoryProvider = DataFormatProvider({
  name: 'lammpstrj',
  label: 'LAMMPSTRAJ',
  description: 'LAMMPSTRAJ',
  category: CoordinatesFormatCategory,
  stringExtensions: ['lammpstrj'],
  parse: (plugin, data) => {
    const coordinates = plugin.state.data.build().to(data).apply(CoordinatesFromLammpstraj);

    return coordinates.commit();
  },
});

type LammpsTrajectoryProvider = typeof LammpsTrajectoryProvider;
