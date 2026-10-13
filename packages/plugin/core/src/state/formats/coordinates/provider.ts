/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import type { DcdProvider } from './dcd.js';
import type { XtcProvider } from './xtc.js';
import type { TrrProvider } from './trr.js';
import type { LammpsTrajectoryProvider } from './lammps.js';

export type CoordinatesProvider = DcdProvider | XtcProvider | TrrProvider | LammpsTrajectoryProvider;
