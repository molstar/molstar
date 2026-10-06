/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PsfProvider } from './psf.js';
import { PrmtopProvider } from './prmtop.js';
import { TopProvider } from './top.js';

export const BuiltInTopologyFormats = [
  ['psf', PsfProvider] as const,
  ['prmtop', PrmtopProvider] as const,
  ['top', TopProvider] as const,
] as const;

export type BuiltInTopologyFormat = (typeof BuiltInTopologyFormats)[number][0];
