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
import { XscProvider } from './xsc.js';

export const BuiltInCoordinatesFormats = [
  DcdProvider,
  XtcProvider,
  TrrProvider,
  NctrajProvider,
  LammpsTrajectoryProvider,
  // not coordinates: the periodic boxes a DCD's frames lack, given to the last trajectory from coordinates
  XscProvider,
] as const;

export type BuiltInCoordinatesFormat = (typeof BuiltInCoordinatesFormats)[number]['name'];
