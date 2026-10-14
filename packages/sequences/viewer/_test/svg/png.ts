/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A PNG read back in Node, independently of the writer: chunk CRCs against `zlib.crc32`, the zlib stream (and its Adler-32)
// by `inflateSync`, the pixels. For the PNG writer's tests (SV1b).

import { crc32, inflateSync } from 'node:zlib';

export interface DecodedPng {
  w: number;
  h: number;
  /** IHDR's bit depth and color type (8 and 6: RGBA, what the writer makes). */
  depth: number;
  colorType: number;
  /** Chunk types in order. */
  chunks: string[];
  /** The scanlines with their filter bytes, inflated. */
  raw: Uint8Array;
  /** The IDAT data: the zlib stream. */
  zlib: Uint8Array;
}
/** Throws on a bad signature, a chunk whose CRC is wrong, a zlib stream that does not inflate (Adler-32 included). */
export function decodePng(png: Uint8Array): DecodedPng {
  if (Buffer.compare(png.subarray(0, 8), Buffer.from([0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a])) !== 0)
    throw new Error('not a PNG');
  const dv = new DataView(png.buffer, png.byteOffset, png.length);
  const chunks: string[] = [];
  const idat: Uint8Array[] = [];
  let [w, h, at, depth, colorType] = [0, 0, 8, 0, 0];
  while (at < png.length) {
    const len = dv.getUint32(at);
    const type = String.fromCharCode(...png.subarray(at + 4, at + 8));
    if (crc32(png.subarray(at + 4, at + 8 + len)) !== dv.getUint32(at + 8 + len)) throw new Error(`bad CRC in ${type}`);
    if (type === 'IHDR') {
      w = dv.getUint32(at + 8);
      h = dv.getUint32(at + 12);
      [depth, colorType] = [png[at + 16], png[at + 17]];
    }
    if (type === 'IDAT') idat.push(png.subarray(at + 8, at + 8 + len));
    chunks.push(type);
    at += 12 + len;
  }
  const zlib = Buffer.concat(idat);
  return { w, h, depth, colorType, chunks, raw: inflateSync(zlib), zlib };
}
/** The pixels back from the scanlines (filter 0 only), as RGBA bytes. */
export function pixelsOf(d: DecodedPng): Uint8Array {
  if (d.depth !== 8 || d.colorType !== 6) throw new Error('not 8-bit RGBA');
  const out = new Uint8Array(d.w * d.h * 4);
  const row = 1 + 4 * d.w;
  for (let y = 0; y < d.h; y++) {
    if (d.raw[y * row] !== 0) throw new Error(`row ${y} is filtered`);
    out.set(d.raw.subarray(y * row + 1, (y + 1) * row), y * 4 * d.w);
  }
  return out;
}
