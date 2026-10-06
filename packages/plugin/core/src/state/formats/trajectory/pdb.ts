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
import { Task } from '@molstar/core/task';
import { parsePDB } from '@molstar/io/reader/pdb/parser';
import { trajectoryFromPDB } from '@molstar/model/formats/structure/pdb';
import { trajectoryProps } from './helpers.js';
import { TrajectoryFormatProvider, directTrajectory, defaultVisuals } from './provider.js';
import { TrajectoryFormatCategory } from './category.js';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';

export { TrajectoryFromPDB };
type TrajectoryFromPDB = typeof TrajectoryFromPDB;
const TrajectoryFromPDB = PluginStateTransform.BuiltIn({
  name: 'trajectory-from-pdb',
  display: { name: 'Parse PDB', description: 'Parse PDB string and create trajectory.' },
  from: [SO.Data.String],
  to: SO.Molecule.Trajectory,
  params: {
    variant: PD.Select('pdb', PD.arrayToOptions(['pdb', 'pdbqt', 'pqr'] as const)),
  },
})({
  apply({ a, params }) {
    return Task.create('Parse PDB', async (ctx) => {
      const parsed = await parsePDB(a.data, a.label, params.variant).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      const models = await trajectoryFromPDB(parsed.result).runInContext(ctx);
      const props = trajectoryProps(models);
      return new SO.Molecule.Trajectory(models, props);
    });
  },
});

export const PdbProvider = TrajectoryFormatProvider({
  name: 'pdb',
  label: 'PDB',
  description: 'PDB',
  category: TrajectoryFormatCategory,
  stringExtensions: ['pdb', 'ent'],
  ...directTrajectory(TrajectoryFromPDB),
  visuals: defaultVisuals,
});

export const PdbqtProvider = TrajectoryFormatProvider({
  name: 'pdbqt',
  label: 'PDBQT',
  description: 'PDBQT',
  category: TrajectoryFormatCategory,
  stringExtensions: ['pdbqt'],
  ...directTrajectory(TrajectoryFromPDB, { variant: 'pdbqt' }),
  visuals: defaultVisuals,
});

export const PqrProvider = TrajectoryFormatProvider({
  name: 'pqr',
  label: 'PQR',
  description: 'PQR',
  category: TrajectoryFormatCategory,
  stringExtensions: ['pqr'],
  ...directTrajectory(TrajectoryFromPDB, { variant: 'pqr' }),
  visuals: defaultVisuals,
});

/** The Pdb data format with its actions. */
export const Pdb: PluginRegistryEntry = {
  formats: [PdbProvider],
  actions: [TrajectoryFromPDB],
};

/** The Pdbqt data format with its actions. */
export const Pdbqt: PluginRegistryEntry = {
  formats: [PdbqtProvider],
  actions: [TrajectoryFromPDB],
};

/** The Pqr data format with its actions. */
export const Pqr: PluginRegistryEntry = {
  formats: [PqrProvider],
  actions: [TrajectoryFromPDB],
};
