/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task } from '@molstar/core/task';
import type { ReadonlyUint8Array } from '@molstar/core/util/type-helpers';
import { ReaderResult as Result } from '../result.js';
import { type SplatFile, SplatRowBytes } from './schema.js';

/** Parses antimatter15's `.splat`: its rows counted, the bytes kept as they are (no copy) */
export function parseSplat(data: ReadonlyUint8Array) {
  return Task.create<Result<SplatFile>>('Parse SPLAT', async () => {
    const count = Math.floor(data.length / SplatRowBytes);
    const rest = data.length - count * SplatRowBytes;
    const warnings = rest ? [`${rest} bytes past the last whole row of ${SplatRowBytes} ignored`] : [];
    const view = new DataView(data.buffer, data.byteOffset, data.byteLength);
    return Result.success<SplatFile>({ count, data, view }, warnings);
  });
}
