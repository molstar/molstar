/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { gzipSync } from 'node:zlib';
import { parseSplat } from '../splat/parser.js';
import { SplatRow, SplatRowBytes } from '../splat/schema.js';
import { parseSpz } from '../spz/parser.js';
import { SpzMagic } from '../spz/schema.js';

describe('splat', () => {
  it('counts the rows, the bytes kept as they are (no copy)', async () => {
    const data = new Uint8Array(2 * SplatRowBytes + 5);
    const view = new DataView(data.buffer);
    view.setFloat32(SplatRowBytes + SplatRow.position + 4, 2.5, true);
    const parsed = await parseSplat(data).run();
    if (parsed.isError) throw new Error(parsed.message);
    const file = parsed.result;
    expect(file.count).toBe(2);
    expect(parsed.warnings).toEqual(['5 bytes past the last whole row of 32 ignored']);
    expect(file.data).toBe(data);
    expect(file.view.buffer).toBe(data.buffer);
    expect(file.view.getFloat32(SplatRowBytes + SplatRow.position + 4, true)).toBe(2.5);
  });
});

/** a header of `count` splats and `fill`-valued sections of their sizes (version 2 or 3, degree 1: 9 SH bytes) */
function spz(version: number, count: number, degree = 1) {
  const sizes = [9, 1, 3, 3, version >= 3 ? 4 : 3, degree ? 9 : 0];
  const data = new Uint8Array(16 + count * sizes.reduce((a, b) => a + b, 0));
  const view = new DataView(data.buffer);
  view.setUint32(0, SpzMagic, true);
  view.setUint32(4, version, true);
  view.setUint32(8, count, true);
  data[12] = degree;
  data[13] = 12;
  data[14] = 1;
  // each section filled with its index: its bounds checked by where the values change
  let o = 16;
  sizes.forEach((size, s) => {
    data.fill(s + 1, o, o + count * size);
    o += count * size;
  });
  return data;
}

describe('spz', () => {
  it('reads the header and the sections as views of the bytes (no copy), inflated where gzipped', async () => {
    for (const version of [2, 3]) {
      const raw = spz(version, 3);
      for (const input of [raw, gzipSync(raw)]) {
        const parsed = await parseSpz(input).run();
        if (parsed.isError) throw new Error(parsed.message);
        const f = parsed.result;
        expect([f.version, f.count, f.shDegree, f.fractionalBits, f.flags]).toEqual([version, 3, 1, 12, 1]);
        const sections = [f.positions, f.alphas, f.colors, f.scales, f.rotations, f.sh];
        expect(sections.map((s) => s.length)).toEqual([27, 3, 9, 9, version >= 3 ? 12 : 9, 27]);
        sections.forEach((s, i) => expect(Array.from(s).every((v) => v === i + 1)).toBe(true));
        // not gzipped: the input's own memory
        if (input === raw) for (const s of sections) expect(s.buffer).toBe(raw.buffer);
      }
    }
  });

  it('rejects other data: no header, version 4, too short', async () => {
    expect((await parseSpz(new Uint8Array(32)).run()).isError).toBe(true);
    const v4 = spz(3, 1);
    new DataView(v4.buffer).setUint32(4, 4, true);
    const parsed = await parseSpz(v4).run();
    expect(parsed.isError && parsed.message).toMatch(/version 4/);
    const short = spz(3, 2).subarray(0, 40);
    expect((await parseSpz(short).run()).isError).toBe(true);
  });
});
