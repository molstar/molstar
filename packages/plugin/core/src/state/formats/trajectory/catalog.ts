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
  MmcifProvider,
  CifCoreProvider,
  PdbProvider,
  PdbqtProvider,
  PqrProvider,
  GroProvider,
  XyzProvider,
  LammpsDataProvider,
  LammpsTrajectoryDataProvider,
  MolProvider,
  SdfProvider,
  Mol2Provider,
] as const;

export type BuiltInTrajectoryFormat = (typeof BuiltInTrajectoryFormats)[number]['name'];
