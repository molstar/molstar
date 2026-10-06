/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { MmcifProvider } from './mmcif.js';
import { CifCoreProvider } from './cif-core.js';
import { PdbProvider, PdbqtProvider, PqrProvider } from './pdb.js';
import { GroProvider } from './gro.js';
import { XyzProvider } from './xyz.js';
import { LammpsDataProvider, LammpsTrajectoryDataProvider } from './lammps.js';
import { MolProvider } from './mol.js';
import { SdfProvider } from './sdf.js';
import { Mol2Provider } from './mol2.js';

export const BuiltInTrajectoryFormats = [
  ['mmcif', MmcifProvider] as const,
  ['cifCore', CifCoreProvider] as const,
  ['pdb', PdbProvider] as const,
  ['pdbqt', PdbqtProvider] as const,
  ['pqr', PqrProvider] as const,
  ['gro', GroProvider] as const,
  ['xyz', XyzProvider] as const,
  ['lammps_data', LammpsDataProvider] as const,
  ['lammps_traj_data', LammpsTrajectoryDataProvider] as const,
  ['mol', MolProvider] as const,
  ['sdf', SdfProvider] as const,
  ['mol2', Mol2Provider] as const,
] as const;

export type BuiltInTrajectoryFormat = (typeof BuiltInTrajectoryFormats)[number][0];
