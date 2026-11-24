/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Task } from '@molstar/core/task';
import { Mat4, Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { parsePly } from '@molstar/io/reader/ply/parser';
import type { GaussianSplat, GaussianSplatData } from '@molstar/model/formats/shape/gaussian-splat/data';
import {
  getSplatOpacity,
  getSplatRotation,
  getSplatScale,
  getSplatShCoefficient,
} from '@molstar/graphics/geo/geometry/gaussian-splats/splat-data';
import { buildGaussianSplats, GaussianSplatFileParams, getSplatTransforms } from '../gaussian-splat.js';
import { shapeFromPly } from '../ply.js';

/** two splats, the second's view-dependent color of degree 1 */
const Splats = [
  { p: [1, 2, 3], s: [-1, 0, 1], q: Quat.normalize(Quat(), Quat.create(0.2, -0.5, 0.1, 0.8)), o: 0.75, c: [1, 0.5, 0] },
  { p: [-4, 5, -6], s: [0.5, -2, 0], q: Quat.identity(), o: 0.25, c: [0.2, 0.4, 1.2] },
];
const sh = (i: number, k: number, c: number) => (i + 1) * (k - 1) * 0.3 + c * 0.1;

const Data: GaussianSplatData = {
  count: Splats.length,
  upright: false,
  getSplat: (i: number, out: GaussianSplat) => {
    const { p, s, q, o, c } = Splats[i];
    Object.assign(out, {
      x: p[0],
      y: p[1],
      z: p[2],
      logScaleX: s[0],
      logScaleY: s[1],
      logScaleZ: s[2],
      qx: q[0],
      qy: q[1],
      qz: q[2],
      qw: q[3],
      opacity: o,
      r: c[0],
      g: c[1],
      b: c[2],
    });
  },
  sh: { degree: 1, get: sh },
};

describe('gaussian splat shapes', () => {
  it('build the geometry of splats read: each decoded into the layouts the shaders read', async () => {
    const { splats, colors } = await Task.create('build', (ctx) => buildGaussianSplats(ctx, Data)).run();
    expect(splats.splatCount).toBe(2);
    const pg = splats.positionGroup.ref.value.array;
    const data = splats.splatData.ref.value.array;
    Splats.forEach(({ p, s, q, o }, i) => {
      expect(Array.from(pg.subarray(i * 4, i * 4 + 3))).toEqual(p);
      const scale = getSplatScale(Vec3(), data, i);
      for (let k = 0; k < 3; ++k) expect(Math.log(scale[k])).toBeCloseTo(s[k], 3);
      const r = getSplatRotation(Quat(), data, i);
      const sign = Quat.dot(r, q) < 0 ? -1 : 1;
      for (let k = 0; k < 4; ++k) expect(r[k] * sign).toBeCloseTo(q[k], 4);
      expect(getSplatOpacity(data, i)).toBeCloseTo(o, 4);
    });
    // bytes, clamped
    expect(Array.from(colors)).toEqual([255, 128, 0, 51, 102, 255]);
    expect(splats.sh.degree).toBe(1);
    const step = splats.sh.scale / 127.5;
    for (let i = 0; i < 2; ++i) {
      for (let k = 0; k < 3; ++k) {
        for (let c = 0; c < 3; ++c) {
          expect(Math.abs(getSplatShCoefficient(splats.sh, i, k, c) - sh(i, k, c))).toBeLessThanOrEqual(
            step / 2 + 1e-6,
          );
        }
      }
    }
  });

  it('are provided for PLY files of splats, meshes for the others; upright by default', async () => {
    const ply = async (text: string) => {
      const parsed = await parsePly(text).run();
      if (parsed.isError) throw new Error(parsed.message);
      return parsed.result;
    };
    const mesh = await ply(
      'ply\nformat ascii 1.0\nelement vertex 3\nproperty float x\nproperty float y\nproperty float z\nelement face 1\nproperty list uchar int vertex_index\nend_header\n0 0 0\n1 0 0\n0 1 0\n3 0 1 2\n',
    );
    expect((await shapeFromPly(mesh).run()).label).toBe('Mesh');
    const props = ['x', 'y', 'z', 'f_dc_0', 'f_dc_1', 'f_dc_2', 'opacity', 'scale_0', 'scale_1', 'scale_2'];
    const rot = ['rot_0', 'rot_1', 'rot_2', 'rot_3'];
    const header = [...props, ...rot].map((p) => `property float ${p}`).join('\n');
    const splats = await ply(
      `ply\nformat ascii 1.0\nelement vertex 1\n${header}\nend_header\n0 0 0 0 0 0 0 0 0 0 1 0 0 0\n`,
    );
    const provider = await shapeFromPly(splats).run();
    expect(provider.label).toBe('Gaussian Splats');
    expect(provider.geometryUtils.Params).toHaveProperty('kernelSize');
    expect(PD.getDefaultValues(provider.params as typeof GaussianSplatFileParams).upright).toBe(true);
  });

  // the data kept as stored: y-down files (PLY, .splat, SOG) turned where upright, SPZ (stored upright) where not
  it('turn the splats 180 degrees about x by the shape transform where the stored frame is not the one wanted', () => {
    const p = Vec3.create(1, 2, 3);
    const at = (transforms: Mat4[] | undefined) => (transforms ? Vec3.transformMat4(Vec3(), p, transforms[0]) : p);
    expect(getSplatTransforms(undefined, false, false)).toBeUndefined();
    expect(getSplatTransforms(undefined, true, true)).toBeUndefined();
    for (const [stored, upright] of [
      [false, true],
      [true, false],
    ]) {
      const t = getSplatTransforms(undefined, stored, upright);
      expect(t).toHaveLength(1);
      expect(Array.from(at(t)).map((v) => Math.round(v * 1e6) / 1e6)).toEqual([1, -2, -3]);
    }
    // the file's own transforms after the turn
    const move = Mat4.fromTranslation(Mat4(), Vec3.create(10, 0, 0));
    const t = getSplatTransforms([move, Mat4.identity()], false, true)!;
    expect(t).toHaveLength(2);
    expect(Array.from(at(t)).map((v) => Math.round(v * 1e6) / 1e6)).toEqual([11, -2, -3]);
  });
});
