/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PlyProvider } from './ply.js';
import { ObjProvider } from './obj.js';
import { VtpProvider } from './vtp.js';

export const BuiltInShapeFormats = [
  ['ply', PlyProvider] as const,
  ['obj', ObjProvider] as const,
  ['vtp', VtpProvider] as const,
] as const;

export type BuildInShapeFormat = (typeof BuiltInShapeFormats)[number][0];
