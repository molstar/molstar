/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Mat4 } from '../../math/linear-algebra/3d/mat4.js';
import { Quat } from '../../math/linear-algebra/3d/quat.js';
import { Vec3 } from '../../math/linear-algebra/3d/vec3.js';
import { Clip } from '../clip.js';

/** A point in the frame of object `i` as the shaders computed it before frames: transformed, less the position, rotated back */
function oldLocal(objects: Clip.Objects, i: number, p: Vec3) {
  const c = Vec3.transformMat4(Vec3(), p, Mat4.fromArray(Mat4(), objects.transform, i * 16));
  Vec3.sub(c, c, Vec3.fromArray(Vec3(), objects.position, i * 3));
  const q = Quat.conjugate(Quat(), Quat.fromArray(Quat(), objects.rotation, i * 4));
  return Vec3.transformQuat(c, c, q);
}

describe('clip objects', () => {
  it('give each object its frame and shape: the point as the shaders transformed it, half the scale, the type', () => {
    let seed = 7;
    const random = () => {
      seed = (seed * 16807) % 2147483647;
      return seed / 2147483647 - 0.5;
    };
    const types = ['plane', 'sphere', 'cube', 'cylinder', 'infiniteCone'] as const;
    const props: Clip.Props = {
      variant: 'pixel',
      objects: types.map((type, i) => ({
        type,
        invert: i % 2 === 1,
        position: Vec3.create(random() * 10, random() * 10, random() * 10),
        rotation: { axis: Vec3.create(random(), random(), random()), angle: random() * 360 },
        scale: Vec3.create(1 + random(), 2 + random(), 3 + random()),
        transform: Mat4.fromTranslation(Mat4(), Vec3.create(random(), random(), random())),
      })),
    };
    const { objects } = Clip.getClip(props);
    expect(objects.count).toBe(types.length);
    for (let i = 0; i < objects.count; ++i) {
      const frame = Mat4.fromArray(Mat4(), objects.frame, i * 16);
      for (let k = 0; k < 4; ++k) {
        const p = Vec3.create(random() * 20, random() * 20, random() * 20);
        const expected = oldLocal(objects, i, p);
        const local = Vec3.transformMat4(Vec3(), p, frame);
        for (let c = 0; c < 3; ++c) expect(local[c]).toBeCloseTo(expected[c], 4);
      }
      const shape = objects.shape.slice(i * 4, i * 4 + 4);
      expect(shape.slice(0, 3)).toEqual(objects.scale.slice(i * 3, i * 3 + 3).map((s) => s * 0.5));
      expect(Math.abs(shape[3])).toBe(objects.type[i]);
      expect(shape[3] < 0).toBe(objects.invert[i]);
    }
  });
});
