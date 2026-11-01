/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { Clip } from '../clip.js';

/** A deterministic generator of numbers in [0, 1) */
function random(seed: number) {
  let s = seed >>> 0;
  return () => {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
    return s / 4294967296;
  };
}

const Shapes: { type: keyof typeof Clip.Type; scale: Vec3 }[] = [
  { type: 'plane', scale: Vec3.create(1, 1, 1) },
  { type: 'sphere', scale: Vec3.create(8, 8, 8) },
  { type: 'sphere', scale: Vec3.create(12, 4, 7) },
  { type: 'cube', scale: Vec3.create(6, 9, 4) },
  { type: 'cylinder', scale: Vec3.create(6, 10, 6) },
  { type: 'infiniteCone', scale: Vec3.create(1, 1, 1) },
  { type: 'infiniteCone', scale: Vec3.create(3, 1, 1) },
];

/** The point's clipping by the objects as the shaders test it (`clipTest`) */
function isClipped(objects: Clip.Objects, p: Vec3) {
  for (let i = 0; i < objects.count; ++i) {
    const inside = Clip.getSignedDistance(objects, i, p) <= 0;
    if (inside !== objects.invert[i]) return true;
  }
  return false;
}

describe('Clip.classifySphere', () => {
  const transforms = [
    Mat4.identity(),
    Mat4.fromTranslation(Mat4(), Vec3.create(1, -2, 3)),
    // a scaled transform (image trims): the frame's distances not the points'
    Mat4.fromScaling(Mat4(), Vec3.create(2, 2, 2)),
  ];
  for (const { type, scale } of Shapes) {
    for (const invert of [false, true]) {
      it(`agrees with points sampled in the sphere: ${type} ${scale.join('x')}${invert ? ', inverted' : ''}`, () => {
        const rand = random(7);
        let decided = 0;
        let total = 0;
        for (const transform of transforms) {
          const objects = Clip.getClip({
            variant: 'pixel',
            objects: [
              {
                type,
                invert,
                position: Vec3.create(0.5, -1, 2),
                rotation: { axis: Vec3.normalize(Vec3(), Vec3.create(1, 2, 0.5)), angle: 35 },
                scale,
                transform,
              },
            ],
          } as Clip.Props).objects;
          for (let k = 0; k < 300; ++k) {
            const center = Vec3.create(rand() * 40 - 20, rand() * 40 - 20, rand() * 40 - 20);
            const sphere = Sphere3D.create(center, 0.2 + rand() * 4);
            const c = Clip.classifySphere(objects, sphere);
            total += 1;
            if (c === Clip.SphereClass.Some) continue;
            decided += 1;
            // the center, points within and on the surface
            for (let n = 0; n < 200; ++n) {
              const d = Vec3.normalize(Vec3(), Vec3.create(rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1));
              const r = n % 4 === 0 ? sphere.radius : sphere.radius * Math.cbrt(rand());
              const p = Vec3.scaleAndAdd(Vec3(), center, d, n === 0 ? 0 : r);
              expect(isClipped(objects, p)).toBe(c === Clip.SphereClass.All);
            }
          }
        }
        // decisive for most spheres (a bound, not every sphere `Some`)
        expect(decided / total).toBeGreaterThan(0.6);
      });
    }
  }

  it('ignores the objects of `flag`, and one clipping every point decides', () => {
    const objects = Clip.getClip({
      variant: 'pixel',
      objects: [
        {
          type: 'plane',
          invert: false,
          position: Vec3.create(0, 0, 0),
          rotation: { axis: Vec3.create(1, 0, 0), angle: 0 },
          scale: Vec3.create(1, 1, 1),
          transform: Mat4.identity(),
        },
        {
          type: 'sphere',
          invert: false,
          position: Vec3.create(0, -10, 0),
          rotation: { axis: Vec3.create(1, 0, 0), angle: 0 },
          scale: Vec3.create(4, 4, 4),
          transform: Mat4.identity(),
        },
      ],
    } as Clip.Props).objects;
    // above the plane (clipped by it), far from the sphere
    const above = Sphere3D.create(Vec3.create(0, 5, 0), 1);
    expect(Clip.classifySphere(objects, above)).toBe(Clip.SphereClass.All);
    expect(Clip.classifySphere(objects, above, 1)).toBe(Clip.SphereClass.None);
    // below the plane, inside the sphere
    const inside = Sphere3D.create(Vec3.create(0, -10, 0), 0.5);
    expect(Clip.classifySphere(objects, inside)).toBe(Clip.SphereClass.All);
    expect(Clip.classifySphere(objects, inside, 2)).toBe(Clip.SphereClass.None);
    expect(Clip.classifySphere(objects, Sphere3D.create(Vec3.create(0, -10, 0), 3))).toBe(Clip.SphereClass.Some);
  });
});
