/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The recorded byte ranges (_test/fixtures/slices.json and one .bin per file, written by scripts/record.ts): a transport
// that serves what was recorded and fails, naming the range, for anything else.

import { readFileSync } from 'node:fs';
import type { RangeReader } from '../../../src/stream/gmod/index.js';

export interface SliceMeta {
  url: string;
  size: number;
  modified: string | null;
  /** [file offset, bytes asked for, bytes the server sent] in the order stored in the .bin. */
  chunks: [number, number, number][];
}
const dir = new URL('./fixtures/', import.meta.url);
export const meta: Record<string, SliceMeta> = JSON.parse(readFileSync(new URL('slices.json', dir), 'utf8'));

/** The recorded bytes of `key`, as a transport; `served` counts the reads it answered. */
export function sliceReader(key: string): RangeReader & { served: number } {
  const bin = readFileSync(new URL(`${key}.bin`, dir));
  const chunks = meta[key].chunks.map(([offset, asked, length], i, all) => ({
    offset,
    asked,
    length,
    at: all.slice(0, i).reduce((n, c) => n + c[2], 0),
  }));
  const read = (async (offset: number, length: number) => {
    read.served++;
    const c = chunks.find((c) => offset >= c.offset && offset + length <= c.offset + c.asked);
    if (!c) throw new Error(`${key}: bytes ${offset}..${offset + length} are not in the recorded slices`);
    const from = Math.min(offset - c.offset, c.length);
    return new Uint8Array(bin.subarray(c.at + from, c.at + Math.min(offset - c.offset + length, c.length)));
  }) as RangeReader & { served: number };
  read.served = 0;
  return read;
}

/** A whole recorded file (the .tbi). */
export const whole = async (key: string) => sliceReader(key)(0, meta[key].size);
