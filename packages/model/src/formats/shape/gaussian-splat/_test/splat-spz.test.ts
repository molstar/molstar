/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { gzipSync } from 'node:zlib';
import { Quat } from '@molstar/core/math/linear-algebra';
import { parseSplat } from '@molstar/io/reader/splat/parser';
import { parseSpz } from '@molstar/io/reader/spz/parser';
import { SH_C0 } from '../data.js';
import { gaussianSplatsFromSplat } from '../splat.js';
import { gaussianSplatsFromSpz } from '../spz.js';
import { decode, expectRotation } from './decode.js';

async function readSplat(bytes: Uint8Array) {
  const parsed = await parseSplat(bytes).run();
  if (parsed.isError) throw new Error(parsed.message);
  return gaussianSplatsFromSplat(parsed.result);
}

async function readSpz(bytes: Uint8Array, shDegree = 0) {
  const parsed = await parseSpz(bytes).run();
  if (parsed.isError) throw new Error(parsed.message);
  return gaussianSplatsFromSpz(parsed.result, { shDegree });
}

interface TestSplat {
  position: number[];
  /** linear */
  scale: number[];
  rotation: Quat;
  opacity: number;
  /** the base color, 0 to 1 */
  color: number[];
}

const splats: TestSplat[] = [
  {
    position: [1.5, -2.25, 3],
    scale: [0.5, 1, 2],
    rotation: Quat.normalize(Quat(), Quat.create(0.2, -0.5, 0.1, 0.8)),
    opacity: 0.75,
    color: [1, 0.5, 0.25],
  },
  {
    position: [-4, 0.125, -0.5],
    scale: [0.1, 0.2, 0.3],
    rotation: Quat.normalize(Quat(), Quat.create(-0.9, 0.2, 0.1, 0.3)),
    opacity: 0.25,
    color: [0.1, 0.9, 0.6],
  },
];

/** antimatter15's packing (`convert.py`): float32 position and linear scales, color and opacity bytes, `w x y z` bytes */
function splatFile(list: TestSplat[]) {
  const bytes = new Uint8Array(list.length * 32);
  const view = new DataView(bytes.buffer);
  list.forEach((s, i) => {
    const o = i * 32;
    s.position.forEach((v, k) => view.setFloat32(o + k * 4, v, true));
    s.scale.forEach((v, k) => view.setFloat32(o + 12 + k * 4, v, true));
    [...s.color, s.opacity].forEach((v, k) => (bytes[o + 24 + k] = Math.round(v * 255)));
    const [x, y, z, w] = s.rotation;
    [w, x, y, z].forEach((v, k) => (bytes[o + 28 + k] = Math.min(255, Math.round(v * 128 + 128))));
  });
  return bytes;
}

/** Niantic's packing (`load-spz.cc`) of `list` in its frame, version 2 or 3, with `sh(i, k, c)` of degree 1 */
function spzFile(list: TestSplat[], version: 2 | 3, sh?: (i: number, k: number, c: number) => number) {
  const n = list.length;
  const degree = sh ? 1 : 0;
  const dim = sh ? 3 : 0;
  const rotationBytes = version === 3 ? 4 : 3;
  const bytes = new Uint8Array(16 + n * (9 + 1 + 3 + 3 + rotationBytes + dim * 3));
  const view = new DataView(bytes.buffer);
  view.setUint32(0, 0x5053474e, true);
  view.setUint32(4, version, true);
  view.setUint32(8, n, true);
  bytes[12] = degree;
  bytes[13] = 12;
  let o = 16;
  for (const s of list) {
    for (const v of s.position) {
      const f = Math.round(v * 4096) & 0xffffff;
      bytes[o++] = f & 0xff;
      bytes[o++] = (f >> 8) & 0xff;
      bytes[o++] = (f >> 16) & 0xff;
    }
  }
  for (const s of list) bytes[o++] = Math.round(s.opacity * 255);
  for (const s of list) {
    for (const v of s.color) bytes[o++] = Math.round(((v - 0.5) / SH_C0) * 0.15 * 255 + 127.5);
  }
  for (const s of list) for (const v of s.scale) bytes[o++] = Math.round((Math.log(v) + 10) * 16);
  for (const s of list) {
    const q = Array.from(s.rotation);
    if (version === 2) {
      const sign = q[3] < 0 ? -127.5 : 127.5;
      for (let k = 0; k < 3; ++k) bytes[o++] = Math.round(q[k] * sign + 127.5);
    } else {
      let largest = 0;
      for (let k = 1; k < 4; ++k) if (Math.abs(q[k]) > Math.abs(q[largest])) largest = k;
      const negate = q[largest] < 0;
      let comp = largest;
      for (let k = 0; k < 4; ++k) {
        if (k === largest) continue;
        const neg = q[k] < 0 !== negate ? 1 : 0;
        const mag = Math.floor((511 * Math.abs(q[k])) / Math.SQRT1_2 + 0.5);
        comp = (comp << 10) | (neg << 9) | mag;
      }
      view.setUint32(o, comp >>> 0, true);
      o += 4;
    }
  }
  if (sh) {
    for (let i = 0; i < n; ++i) {
      for (let k = 0; k < dim; ++k) for (let c = 0; c < 3; ++c) bytes[o++] = Math.round(sh(i, k, c) * 128 + 128);
    }
  }
  return bytes;
}

/** the splats as stored */
function expectSplats(
  data: Awaited<ReturnType<typeof readSplat>>,
  tolerance: { position: number; logScale: number; color: number },
) {
  expect(data.count).toBe(splats.length);
  decode(data).forEach((d, i) => {
    const s = splats[i];
    expect(Math.abs(d.x - s.position[0])).toBeLessThan(tolerance.position);
    expect(Math.abs(d.y - s.position[1])).toBeLessThan(tolerance.position);
    expect(Math.abs(d.z - s.position[2])).toBeLessThan(tolerance.position);
    const logScales = [d.logScaleX, d.logScaleY, d.logScaleZ];
    for (let k = 0; k < 3; ++k) expect(Math.abs(logScales[k] - Math.log(s.scale[k]))).toBeLessThan(tolerance.logScale);
    expectRotation(d, s.rotation, 2);
    expect(d.opacity).toBeCloseTo(s.opacity, 2);
    const color = [d.r, d.g, d.b];
    for (let c = 0; c < 3; ++c) expect(Math.abs(color[c] - s.color[c])).toBeLessThan(tolerance.color);
  });
}

describe('gaussian splats of .splat and .spz files', () => {
  it('reads .splat rows as stored, y-down', async () => {
    const data = await readSplat(splatFile(splats));
    expect(data.upright).toBe(false);
    expectSplats(data, { position: 1e-6, logScale: 1e-6, color: 0.003 });
  });

  for (const version of [2, 3] as const) {
    it(`reads .spz version ${version}, gzipped or not, as stored: upright`, async () => {
      const raw = spzFile(splats, version);
      const tolerance = { position: 1 / 4096, logScale: 1 / 32, color: 0.02 };
      const data = await readSpz(raw);
      expect(data.upright).toBe(true);
      expectSplats(data, tolerance);
      expectSplats(await readSpz(gzipSync(raw)), tolerance);
    });
  }

  it('reads the view-dependent color of .spz, coefficient-major', async () => {
    const sh = (i: number, k: number, c: number) => Math.round(Math.sin(i * 3 + k * 1.7 + c) * 100) / 128;
    const raw = spzFile(splats, 3, sh);
    expect((await readSpz(raw)).sh).toBeUndefined();
    const data = await readSpz(raw, 3);
    expect(data.sh?.degree).toBe(1);
    for (let i = 0; i < splats.length; ++i) {
      for (let k = 0; k < 3; ++k) {
        for (let c = 0; c < 3; ++c) expect(data.sh!.get(i, k, c)).toBeCloseTo(sh(i, k, c), 6);
      }
    }
  });
});
