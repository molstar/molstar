/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Quat } from '@molstar/core/math/linear-algebra';
import { parsePly } from '@molstar/io/reader/ply/parser';
import { SH_C0 } from '../data.js';
import { decode, expectRotation } from './decode.js';
import { gaussianSplatsFromPly, isCompressedGaussianSplatPly, isGaussianSplatPly } from '../ply.js';

async function ply(text: string) {
  const parsed = await parsePly(text).run();
  if (parsed.isError) throw new Error(parsed.message);
  return parsed.result;
}

function read(file: Awaited<ReturnType<typeof ply>>, shDegree = 0) {
  return gaussianSplatsFromPly(file, { shDegree });
}

const MeshPly = `ply
format ascii 1.0
element vertex 3
property float x
property float y
property float z
element face 1
property list uchar int vertex_index
end_header
0 0 0
1 0 0
0 1 0
3 0 1 2
`;

const StandardHeader = `ply
format ascii 1.0
element vertex VERTICES
property float x
property float y
property float z
property float f_dc_0
property float f_dc_1
property float f_dc_2
property float opacity
property float scale_0
property float scale_1
property float scale_2
property float rot_0
property float rot_1
property float rot_2
property float rot_3
end_header
`;

/** the packing of PlayCanvas's compressed layout, per its decoder: 11-10-11 bits, smallest three, color bytes */
const pack111011 = (x: number, y: number, z: number) =>
  ((Math.round(x * 2047) << 21) | (Math.round(y * 1023) << 11) | Math.round(z * 2047)) >>> 0;
function packRotation(q: Quat) {
  const [x, y, z, w] = Quat.normalize(Quat(), q);
  // the index of the largest (0 w, 1 x, 2 y, 3 z), made positive; the others in the decoder's order
  const values = [w, x, y, z];
  let largest = 0;
  for (let k = 1; k < 4; ++k) if (Math.abs(values[k]) > Math.abs(values[largest])) largest = k;
  const sign = values[largest] < 0 ? -1 : 1;
  const [a, b, c] = [
    [x, y, z],
    [w, y, z],
    [w, x, z],
    [w, x, y],
  ][largest].map((v) => Math.round(((v * sign) / Math.SQRT2) * 1023 + 511.5));
  return ((largest << 30) | (a << 20) | (b << 10) | c) >>> 0;
}
const packColor = (r: number, g: number, b: number, a: number) =>
  ((Math.round(r * 255) << 24) | (Math.round(g * 255) << 16) | (Math.round(b * 255) << 8) | Math.round(a * 255)) >>> 0;

function compressedPly(rotations: Quat[], colorBounds: boolean) {
  const bounds = `-10 -20 -30 10 20 30 -5 -5 -5 1 1 1${colorBounds ? ' 0.2 0.2 0.2 0.6 0.6 0.6' : ''}`;
  const colorProps = colorBounds
    ? 'property float min_r\nproperty float min_g\nproperty float min_b\nproperty float max_r\nproperty float max_g\nproperty float max_b\n'
    : '';
  const rows = rotations.map(
    (q, i) => `${pack111011(0.5, 1, i / 4)} ${packRotation(q)} ${pack111011(1, 0, 0.5)} ${packColor(1, 0.5, 0, 0.75)}`,
  );
  const chunk = [
    'min_x',
    'min_y',
    'min_z',
    'max_x',
    'max_y',
    'max_z',
    'min_scale_x',
    'min_scale_y',
    'min_scale_z',
    'max_scale_x',
    'max_scale_y',
    'max_scale_z',
  ]
    .map((p) => `property float ${p}`)
    .join('\n');
  const packed = ['packed_position', 'packed_rotation', 'packed_scale', 'packed_color']
    .map((p) => `property uint ${p}`)
    .join('\n');
  return `ply\nformat ascii 1.0\nelement chunk 1\n${chunk}\n${colorProps}element vertex ${rotations.length}\n${packed}\nend_header\n${bounds}\n${rows.join('\n')}\n`;
}

/** the standard layout of one splat with `f_rest` of degree 3 (45: 15 a channel, channel-major), `rest(j)` each */
function shPly(rest: (j: number) => number) {
  const header = StandardHeader.replace('VERTICES', '1').replace(
    'end_header',
    `${Array.from({ length: 45 }, (_, j) => `property float f_rest_${j}`).join('\n')}\nend_header`,
  );
  const values = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, ...Array.from({ length: 45 }, (_, j) => rest(j))];
  return `${header}${values.join(' ')}\n`;
}
const rest = (j: number) => Math.round((Math.sin(j * 1.7) * 0.8 + (j % 15) * 0.01) * 1000) / 1000;

describe('gaussian splats of PLY files', () => {
  it('reads the standard layout: logit opacity, log scales, rot_0 as w, the base color by f_dc', async () => {
    const q = Quat.normalize(Quat(), Quat.create(0.2, -0.5, 0.1, 0.8));
    const file = await ply(
      `${StandardHeader.replace('VERTICES', '2')}1 2 3 0 1 -1 0 -1 -2 -3 1 0 0 0\n-4 5 -6 0.5 0 0 2 0.5 0 -0.5 ${q[3]} ${q[0]} ${q[1]} ${q[2]}\n`,
    );
    expect(isGaussianSplatPly(file)).toBe(true);
    expect(isCompressedGaussianSplatPly(file)).toBe(false);
    const data = read(file);
    expect(data.count).toBe(2);
    const [a, b] = decode(data);
    expect([b.x, b.y, b.z]).toEqual([-4, 5, -6]);
    expect([a.logScaleX, a.logScaleY, a.logScaleZ]).toEqual([-1, -2, -3]);
    expectRotation(a, Quat.identity());
    expectRotation(b, q);
    expect(a.opacity).toBeCloseTo(0.5, 6);
    expect(b.opacity).toBeCloseTo(1 / (1 + Math.exp(-2)), 6);
    expect([a.r, a.g, a.b]).toEqual([0.5, 0.5 + SH_C0, 0.5 - SH_C0]);
    expect(b.r).toBeCloseTo(0.5 + SH_C0 * 0.5, 6);
    expect(data.sh).toBeUndefined();
  });

  it('reads the compressed layout: positions and log scales in the chunk bounds, the smallest three of each rotation', async () => {
    // the largest component w, x, y and z in turn
    const rotations = [
      Quat.create(0.1, 0.2, -0.3, 0.9),
      Quat.create(-0.9, 0.2, 0.1, 0.3),
      Quat.create(0.2, 0.8, -0.1, -0.4),
      Quat.create(0.3, -0.1, -0.85, 0.2),
    ];
    const file = await ply(compressedPly(rotations, false));
    expect(isGaussianSplatPly(file)).toBe(true);
    expect(isCompressedGaussianSplatPly(file)).toBe(true);
    decode(read(file)).forEach((s, i) => {
      // x at half its range, y at its max, z by the splat
      expect(s.x).toBeCloseTo(0, 1);
      expect(s.y).toBeCloseTo(20, 4);
      expect(s.z).toBeCloseTo(-30 + (60 * Math.round((i / 4) * 2047)) / 2047, 4);
      expect(s.logScaleX).toBeCloseTo(1, 4);
      expect(s.logScaleY).toBeCloseTo(-5, 4);
      // 11 bits over [-5, 1]: steps of 0.003
      expect(s.logScaleZ).toBeCloseTo(-2, 2);
      expectRotation(s, rotations[i], 2);
      expect(s.opacity).toBeCloseTo(0.75, 2);
      // the color as it is without color bounds
      expect([s.r, s.g, s.b].map((v) => Math.round(v * 255))).toEqual([255, 128, 0]);
    });
  });

  it("lerps the colors in the chunk's color bounds where it has them", async () => {
    const [s] = decode(read(await ply(compressedPly([Quat.identity()], true))));
    // 1 -> 0.6, 0.5 -> 0.4, 0 -> 0.2
    expect([s.r, s.g, s.b].map((v) => Math.round(v * 255))).toEqual([153, 102, 51]);
  });

  it('reads the splats as stored, y-down (turned upright by the shape, not here)', async () => {
    const q = Quat.normalize(Quat(), Quat.create(0.2, -0.5, 0.1, 0.8));
    const file = await ply(
      `${StandardHeader.replace('VERTICES', '1')}1 2 3 0 0 0 0 0 0 0 ${q[3]} ${q[0]} ${q[1]} ${q[2]}\n`,
    );
    const data = read(file);
    expect(data.upright).toBe(false);
    const [s] = decode(data);
    expect([s.x, s.y, s.z]).toEqual([1, 2, 3]);
    expectRotation(s, q);
  });

  it('is told from a mesh by its properties', async () => {
    expect(isGaussianSplatPly(await ply(MeshPly))).toBe(false);
  });

  it('reads the view-dependent color up to the degree asked, channel-major, none by default', async () => {
    const file = await ply(shPly(rest));
    expect(read(file).sh).toBeUndefined();
    for (const degree of [1, 3] as const) {
      const { sh } = read(file, degree);
      expect(sh?.degree).toBe(degree);
      const n = degree === 1 ? 3 : 15;
      for (let c = 0; c < 3; ++c) {
        for (let k = 0; k < n; ++k) expect(sh!.get(0, k, c)).toBeCloseTo(rest(c * 15 + k), 6);
      }
    }
  });

  it("reads the compressed layout's view-dependent color: an sh element of bytes", async () => {
    // degree 1: 9 bytes, channel-major
    const bytes = [0, 64, 128, 191, 255, 32, 100, 200, 10];
    const shHeader = bytes.map((_, j) => `property uchar f_rest_${j}`).join('\n');
    const text = compressedPly([Quat.identity()], false).replace('end_header', `element sh 1\n${shHeader}\nend_header`);
    const { sh } = read(await ply(`${text}${bytes.join(' ')}\n`), 3);
    expect(sh?.degree).toBe(1);
    for (let c = 0; c < 3; ++c) {
      for (let k = 0; k < 3; ++k) expect(sh!.get(0, k, c)).toBeCloseTo((bytes[c * 3 + k] * 8) / 255 - 4, 6);
    }
  });
});
