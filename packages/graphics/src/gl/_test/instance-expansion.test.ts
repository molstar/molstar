/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Frustum3D } from '@molstar/core/math/geometry/primitives/frustum3d';
import { Plane3D } from '@molstar/core/math/geometry/primitives/plane3d';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import {
  createExpansionCopies,
  createInstanceExpansion,
  getCopyBuckets,
  getCopyOrderImage,
  getCopyRange,
  getCopyReach,
  isSharedCopyOrder,
} from '../instance-expansion.js';

/** a frustum holding everything within 1e6 of the origin */
function openFrustum() {
  const plane = (x: number, y: number, z: number) => Plane3D.create(Vec3.create(x, y, z), 1e6);
  return Frustum3D.create(
    plane(-1, 0, 0),
    plane(1, 0, 0),
    plane(0, 1, 0),
    plane(0, -1, 0),
    plane(0, 0, -1),
    plane(0, 0, 1),
  );
}

/** copies of a template's sphere (radius 5 about the origin) placed at `zs` along z, turned by `turn` about y */
function createCopies(zs: number[], turn = 0) {
  const transforms = new Float32Array(zs.length * 16);
  const m = Mat4();
  zs.forEach((z, i) => {
    Mat4.fromRotation(m, turn * i, Vec3.create(0, 1, 0));
    Mat4.setTranslation(m, Vec3.create(0, 0, z));
    Mat4.toArray(m, transforms, i * 16);
  });
  return createExpansionCopies(transforms, zs.length, Sphere3D.create(Vec3(), 5), 1, 0);
}

describe('instance expansion', () => {
  // the camera's plane z = 0, looking along z: a copy's distance its z
  const cameraPlane = Plane3D.create(Vec3.create(0, 0, 1), 0);
  const lod = { depthScale: 1, pixelSize: 0, version: 1 };

  it('orders a level of copies by their centers distance', () => {
    const copies = createCopies([50, 10, 30, 20, 40]);
    const buckets = getCopyBuckets(copies, cameraPlane, openFrustum(), lod);
    expect([...buckets.order]).toEqual([1, 3, 2, 4, 0]);
    expect([...buckets.distances]).toEqual([10, 20, 30, 40, 50]);
    // the copies within a reach of 6 of the distances 25 to 35: 20, 30 and 40
    expect(getCopyRange(buckets.distances, 0, 5, 25, 35, 6, [0, 0])).toEqual([1, 4]);
    // none
    expect(getCopyRange(buckets.distances, 0, 5, 100, 200, 6, [0, 0])).toEqual([5, 5]);
    expect(getCopyRange(buckets.distances, 0, 5, -Infinity, Infinity, 0, [0, 0])).toEqual([0, 5]);
  });

  it('takes the copies an object is within distances of, a superset', () => {
    const zs = Array.from({ length: 64 }, (_, i) => (i * 37) % 64).map((k) => k * 3);
    const copies = createCopies(zs, 0.4);
    const sphere = Sphere3D.create(Vec3.create(3, 4, 0), 2);
    const expansion = createInstanceExpansion(copies, sphere, 0);
    // the object's center 5 from the template's, its radius 2
    const reach = getCopyReach(expansion);
    expect(reach).toBeCloseTo(7, 4);
    const buckets = getCopyBuckets(copies, cameraPlane, openFrustum(), { ...lod, version: 2 });
    const { spheres } = expansion;
    for (const [min, max] of [
      [0, 40],
      [41, 90],
      [90, 200],
    ]) {
      const [from, to] = getCopyRange(buckets.distances, 0, copies.count, min, max, reach, [0, 0]);
      const taken = new Set(buckets.order.subarray(from, to));
      for (let i = 0; i < copies.count; ++i) {
        const d = spheres[i * 4 + 2];
        const r = spheres[i * 4 + 3];
        if (d + r >= min && d - r <= max) expect(taken.has(i)).toBe(true);
      }
    }
  });

  it('keeps a template order a view, rewriting the least recent', () => {
    const copies = createCopies([50, 10, 30, 20, 40]);
    const images = [1, 2, 3, 4].map((version) => {
      getCopyBuckets(copies, cameraPlane, openFrustum(), { ...lod, version });
      return getCopyOrderImage(copies);
    });
    expect(new Set(images).size).toBe(4);
    for (const image of images) expect(isSharedCopyOrder(image)).toBe(true);
    expect([...images[0].array.subarray(0, 5)]).toEqual([1, 3, 2, 4, 0]);
    const first = images[0].version;
    // an earlier view's: the same image, not written again
    getCopyBuckets(copies, cameraPlane, openFrustum(), { ...lod, version: 2 });
    expect(getCopyOrderImage(copies)).toBe(images[1]);
    expect(images[0].version).toBe(first);
    // a fifth view's: the first's image rewritten (an object drawing by it culled again)
    getCopyBuckets(copies, cameraPlane, openFrustum(), { ...lod, version: 5 });
    expect(getCopyOrderImage(copies)).toBe(images[0]);
    expect(images[0].version).toBe(first! + 1);
    // an object's own order not a template's
    expect(isSharedCopyOrder({ array: new Float32Array(4), width: 1, height: 1 })).toBe(false);
  });

  it("culls a copy by its objects' reach, not its template's sphere alone (targets past the particles' centers)", () => {
    const copies = createCopies([10, 30]);
    // the view from z = 40 on: the copy at 30's template sphere (radius 5) short of it, the one at 10's too
    const plane = (x: number, y: number, z: number, c = 1e6) => Plane3D.create(Vec3.create(x, y, z), c);
    const beyond40 = Frustum3D.create(
      plane(-1, 0, 0),
      plane(1, 0, 0),
      plane(0, 1, 0),
      plane(0, -1, 0),
      plane(0, 0, -1),
      plane(0, 0, 1, -40),
    );
    const seen = (version: number) => {
      const b = getCopyBuckets(copies, cameraPlane, beyond40, { ...lod, version });
      return [...b.order.subarray(0, b.offsets[1])];
    };
    expect(seen(2)).toEqual([]);
    // an object whose instances reach 12 past the template's center (a target's geometry): the copy at 30 in view, the
    // buckets of the same view culled again
    createInstanceExpansion(copies, Sphere3D.create(Vec3.create(0, 0, 4), 8), 0);
    expect(copies.reach).toBe(12);
    expect(seen(2)).toEqual([1]);
  });

  it('keeps a view of moved copies by its transform', () => {
    const copies = createCopies([50, 10, 30, 20, 40]);
    const view = { ...lod, version: 1 };
    getCopyBuckets(copies, cameraPlane, openFrustum(), view);
    const unmoved = getCopyOrderImage(copies);
    // moved away from the camera, the order reversed
    const flip = Mat4.fromRotation(Mat4(), Math.PI, Vec3.create(0, 1, 0));
    Mat4.setTranslation(flip, Vec3.create(0, 0, 100));
    const buckets = getCopyBuckets(copies, cameraPlane, openFrustum(), view, flip);
    expect([...buckets.order]).toEqual([0, 4, 2, 3, 1]);
    expect([...buckets.distances]).toEqual([50, 60, 70, 80, 90]);
    const moved = getCopyOrderImage(copies);
    expect(moved).not.toBe(unmoved);
    // the same transform (another matrix) the same buckets, not sorted again
    const same = Mat4.copy(Mat4(), flip);
    buckets.distances[0] = -1;
    expect(getCopyBuckets(copies, cameraPlane, openFrustum(), view, same).distances[0]).toBe(-1);
    expect(getCopyOrderImage(copies)).toBe(moved);
    // unmoved again: its image of the view kept
    getCopyBuckets(copies, cameraPlane, openFrustum(), view);
    expect([...buckets.order]).toEqual([1, 3, 2, 4, 0]);
    expect(getCopyOrderImage(copies)).toBe(unmoved);
  });
});
