/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { UnitStyles } from '@molstar/io/reader/lammps/schema';
import { Task } from '@molstar/core/task';
import { parseLammpsData } from '@molstar/io/reader/lammps/data/parser';
import { trajectoryFromLammpsData } from '@molstar/model/formats/structure/lammps-data';
import { trajectoryProps } from './helpers.js';
import { parseLammpsTrajectory } from '@molstar/io/reader/lammps/traj/parser';
import { trajectoryFromLammpsTrajectory } from '@molstar/model/formats/structure/lammps-trajectory';
import { TrajectoryFormatProvider, directTrajectory, defaultVisuals } from './provider.js';
import { TrajectoryFormatCategory } from './category.js';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';

export { TrajectoryFromLammpsData };
type TrajectoryFromLammpsData = typeof TrajectoryFromLammpsData;
const TrajectoryFromLammpsData = PluginStateTransform.BuiltIn({
  name: 'trajectory-from-lammps-data',
  display: { name: 'Parse Lammps Data', description: 'Parse Lammps Data from string and create trajectory.' },
  from: [SO.Data.String],
  to: SO.Molecule.Trajectory,
  params: {
    unitsStyle: PD.Select('real', PD.arrayToOptions(UnitStyles)),
  },
})({
  apply({ a, params }) {
    return Task.create('Parse Lammps Data', async (ctx) => {
      const parsed = await parseLammpsData(a.data).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      const models = await trajectoryFromLammpsData(parsed.result, params.unitsStyle).runInContext(ctx);
      const props = trajectoryProps(models);
      return new SO.Molecule.Trajectory(models, props);
    });
  },
});

export { TrajectoryFromLammpsTrajData };
type TrajectoryFromLammpsTrajData = typeof TrajectoryFromLammpsTrajData;
const TrajectoryFromLammpsTrajData = PluginStateTransform.BuiltIn({
  name: 'trajectory-from-lammps-traj-data',
  display: { name: 'Parse Lammps traj Data', description: 'Parse Lammps Traj Data string and create trajectory.' },
  from: [SO.Data.String],
  to: SO.Molecule.Trajectory,
  params: {
    unitsStyle: PD.Select('real', PD.arrayToOptions(UnitStyles)),
  },
})({
  apply({ a, params }) {
    return Task.create('Parse Lammps Data', async (ctx) => {
      const parsed = await parseLammpsTrajectory(a.data).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      const models = await trajectoryFromLammpsTrajectory(parsed.result, params.unitsStyle).runInContext(ctx);
      const props = trajectoryProps(models);
      return new SO.Molecule.Trajectory(models, props);
    });
  },
});

export const LammpsDataProvider = TrajectoryFormatProvider({
  name: 'lammps_data',
  label: 'Lammps Data',
  description: 'Lammps Data',
  category: TrajectoryFormatCategory,
  stringExtensions: ['data'],
  ...directTrajectory(TrajectoryFromLammpsData),
  visuals: defaultVisuals,
});

export const LammpsTrajectoryDataProvider = TrajectoryFormatProvider({
  name: 'lammps_traj_data',
  label: 'Lammps Trajectory Data',
  description: 'Lammps Trajectory Data',
  category: TrajectoryFormatCategory,
  stringExtensions: ['lammpstrj'],
  ...directTrajectory(TrajectoryFromLammpsTrajData),
  visuals: defaultVisuals,
});

/** The LammpsData data format. */
export const LammpsData: PluginRegistryEntry = {
  formats: [LammpsDataProvider],
};

/** The LammpsTrajectoryData data format. */
export const LammpsTrajectoryData: PluginRegistryEntry = {
  formats: [LammpsTrajectoryDataProvider],
};
