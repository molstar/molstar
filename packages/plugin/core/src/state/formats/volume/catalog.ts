/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Aliaksei Chareshneu <chareshneu.tech@gmail.com>
 */

import { Ccp4Provider } from './ccp4.js';
import { Dsn6Provider } from './dsn6.js';
import { CubeProvider } from './cube.js';
import { DxProvider } from './dx.js';
import { DscifProvider } from './density-server.js';
import { SegcifProvider } from './segmentation.js';
import { SfcifProvider } from './structure-factors.js';
import { MtzProvider } from './mtz.js';

export const BuiltInVolumeFormats = [
  ['ccp4', Ccp4Provider] as const,
  ['dsn6', Dsn6Provider] as const,
  ['cube', CubeProvider] as const,
  ['dx', DxProvider] as const,
  ['dscif', DscifProvider] as const,
  ['segcif', SegcifProvider] as const,
  ['sfcif', SfcifProvider] as const,
  ['mtz', MtzProvider] as const,
] as const;

export type BuildInVolumeFormat = (typeof BuiltInVolumeFormats)[number][0];
