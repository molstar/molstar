/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The SVG's PNGs: the signature, the header, the rows in stored blocks (filter `none`) with their checksums, and
// base64 as the standard alphabet pads it.

import { adler, crc } from '@molstar/core/util/zip/checksum';
import { describe, expect, it } from 'vitest';
import { base64Of, encodePng } from '../png.js';

const u32 = (b: Uint8Array, at: number) => new DataView(b.buffer, b.byteOffset).getUint32(at);

describe('encodePng', () => {
  it('writes the signature, IHDR, one IDAT of stored rows and IEND, each with its CRC', () => {
    const rgba = Uint8Array.from([255, 0, 0, 255, 0, 255, 0, 128, 0, 0, 255, 0, 1, 2, 3, 4]);
    const png = encodePng(2, 2, rgba);
    expect([...png.subarray(0, 8)]).toEqual([137, 80, 78, 71, 13, 10, 26, 10]);
    expect(u32(png, 8)).toBe(13);
    expect(String.fromCharCode(...png.subarray(12, 16))).toBe('IHDR');
    expect([u32(png, 16), u32(png, 20), png[24], png[25]]).toEqual([2, 2, 8, 6]);
    expect(u32(png, 29)).toBe(crc(png, 12, 17) >>> 0);
    const at = 33;
    const n = u32(png, at);
    expect(String.fromCharCode(...png.subarray(at + 4, at + 8))).toBe('IDAT');
    const z = png.subarray(at + 8, at + 8 + n);
    // zlib header, one final stored block of 2 rows × (1 + 8) bytes, its length and complement
    expect([z[0], z[1], z[2], z[3] | (z[4] << 8), (z[5] | (z[6] << 8)) ^ 0xffff]).toEqual([0x78, 0x01, 1, 18, 18]);
    const raw = z.subarray(7, 7 + 18);
    expect([...raw]).toEqual([0, ...rgba.subarray(0, 8), 0, ...rgba.subarray(8)]);
    expect(u32(z, 7 + 18)).toBe(adler(raw, 0, raw.length) >>> 0);
    expect(u32(png, at + 8 + n)).toBe(crc(png, at + 4, 4 + n) >>> 0);
    expect(String.fromCharCode(...png.subarray(png.length - 8, png.length - 4))).toBe('IEND');
  });

  it('splits rows past a stored block’s 65,535 bytes into several blocks', () => {
    const png = encodePng(256, 256, new Uint8Array(256 * 256 * 4));
    const n = u32(png, 33);
    // the zlib header, the raw rows in blocks of 5 header bytes each, the adler
    expect(n).toBe(2 + 1025 * 256 + 5 * Math.ceil((1025 * 256) / 65535) + 4);
  });
});

describe('base64Of', () => {
  it('encodes with the standard alphabet and padding', () => {
    const b = (s: string) => base64Of(Uint8Array.from(s, (c) => c.charCodeAt(0)));
    expect([b(''), b('f'), b('fo'), b('foo'), b('foob'), b('fooba'), b('foobar')]).toEqual([
      '',
      'Zg==',
      'Zm8=',
      'Zm9v',
      'Zm9vYg==',
      'Zm9vYmE=',
      'Zm9vYmFy',
    ]);
  });
});
