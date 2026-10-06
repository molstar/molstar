/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { DcdProvider } from './dcd.js';
import { XtcProvider } from './xtc.js';
import { TrrProvider } from './trr.js';
import { NctrajProvider } from './nctraj.js';
import { LammpsTrajectoryProvider } from './lammps.js';

export const BuiltInCoordinatesFormats = [
  ['dcd', DcdProvider] as const,
  ['xtc', XtcProvider] as const,
  ['trr', TrrProvider] as const,
  ['nctraj', NctrajProvider] as const,
  ['lammpstrj', LammpsTrajectoryProvider] as const,
] as const;

export type BuiltInCoordinatesFormat = (typeof BuiltInCoordinatesFormats)[number][0];
