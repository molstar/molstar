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
  Ccp4Provider,
  Dsn6Provider,
  CubeProvider,
  DxProvider,
  DscifProvider,
  SegcifProvider,
  SfcifProvider,
  MtzProvider,
] as const;

export type BuiltInVolumeFormat = (typeof BuiltInVolumeFormats)[number]['name'];

/** @deprecated use BuiltInVolumeFormat */
export type BuildInVolumeFormat = BuiltInVolumeFormat;
