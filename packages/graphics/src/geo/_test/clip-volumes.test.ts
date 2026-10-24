/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { Sphere3D } from '@molstar/core/math/geometry';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { describe, expect, it } from 'vitest';
import { type ClipObject, clipObjectFrame } from '../clip-volumes.js';

// Where a clip object's drawing goes: its position and rotation, and back through its `transform` (the shader tests a
// point as the transform maps it).

const object = (o: Partial<ClipObject>): ClipObject => ({
  type: 'sphere',
  invert: false,
  position: Vec3.create(1, 2, 3),
  rotation: { axis: Vec3.create(1, 0, 0), angle: 0 },
  scale: Vec3.create(2, 2, 2),
  transform: Mat4.identity(),
  ...o,
});

const at = (frame: Mat4) => Array.from(Vec3.transformMat4(Vec3(), Vec3(), frame));
const around = Sphere3D.create(Vec3(), 10);

describe('clipObjectFrame', () => {
  it('is at the object position without a transform', () => {
    expect(at(clipObjectFrame(object({}), around))).toEqual([1, 2, 3]);
  });

  it('is where the transform takes the object back to (the world point that maps onto its position)', () => {
    // the clip's space shifted by 10 along X: the world point that maps to (1, 2, 3) is (-9, 2, 3)
    const transform = Mat4.fromTranslation(Mat4(), Vec3.create(10, 0, 0));
    at(clipObjectFrame(object({ transform }), around)).forEach((x, i) => expect(x).toBeCloseTo([-9, 2, 3][i], 6));
  });
});
