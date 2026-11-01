/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { PositionData } from '../common.js';
import { OrderedSet } from '@molstar/core/data/int/ordered-set';
import type { NumberArray, PickRequired } from '@molstar/core/util/type-helpers';
import type { Box3D } from './box3d.js';
import type { Axes3D } from './axes3d.js';
import { Vec3 } from '../../linear-algebra/3d/vec3.js';
import { Mat4 } from '../../linear-algebra/3d/mat4.js';
import { EPSILON } from '../../linear-algebra/3d/common.js';

interface Sphere3D {
  center: Vec3;
  radius: number;
  /** points on or inside the bounded set that bound it in many directions: its extreme points or its own points */
  extrema?: Vec3[];
  /**
   * The radius of a ball around each of `extrema` (0 where absent): a set of few points with radii keeps its points as
   * extrema, and their largest radius, rather than extreme points on their balls (never more extrema than points)
   */
  extremaRadius?: number;
}

function Sphere3D() {
  return Sphere3D.zero();
}

namespace Sphere3D {
  export function hasExtrema(sphere: Sphere3D): sphere is PickRequired<Sphere3D, 'extrema'> {
    return sphere.extrema !== undefined;
  }

  export function create(center: Vec3, radius: number): Sphere3D {
    return { center, radius };
  }
  export function zero(): Sphere3D {
    return { center: Vec3(), radius: 0 };
  }

  export function clone(a: Sphere3D): Sphere3D {
    const out = create(Vec3.clone(a.center), a.radius);
    if (hasExtrema(a))
      setExtrema(
        out,
        a.extrema.map((e) => Vec3.clone(e)),
        a.extremaRadius,
      );
    return out;
  }

  /** The radius of a ball around each extremum (`extremaRadius`), 0 where absent */
  export function getExtremaRadius(sphere: Sphere3D) {
    return sphere.extremaRadius ?? 0;
  }

  export function set(out: Sphere3D, center: Vec3, radius: number) {
    Vec3.copy(out.center, center);
    out.radius = radius;
    return out;
  }

  export function setZero(out: Sphere3D) {
    Vec3.set(out.center, 0, 0, 0);
    out.radius = 0;
    out.extrema = void 0;
    out.extremaRadius = void 0;
    return out;
  }

  export function copy(out: Sphere3D, a: Sphere3D) {
    Vec3.copy(out.center, a.center);
    out.radius = a.radius;
    if (hasExtrema(a))
      setExtrema(
        out,
        a.extrema.map((e) => Vec3.clone(e)),
        a.extremaRadius,
      );
    else {
      out.extrema = void 0;
      out.extremaRadius = void 0;
    }
    return out;
  }

  /**
   * Note that `extrema` must not be reused elsewhere; `extremaRadius` the radius of a ball around each (`Sphere3D`), 0
   * where not given
   */
  export function setExtrema(out: Sphere3D, extrema: Vec3[], extremaRadius = 0): Sphere3D {
    if (out.extrema === extrema) {
      // set in place (`expand` of a sphere into itself)
    } else if (out.extrema !== undefined) {
      out.extrema.length = 0;
      for (let i = 0, il = extrema.length; i < il; ++i) out.extrema.push(extrema[i]);
    } else {
      out.extrema = extrema;
    }
    // not left from an earlier set: a reused sphere would cover more than its points
    out.extremaRadius = extremaRadius > 0 ? extremaRadius : void 0;
    return out;
  }

  export function computeBounding(data: PositionData): Sphere3D {
    const { x, y, z, indices } = data;
    let cx = 0,
      cy = 0,
      cz = 0;
    let radiusSq = 0;

    const size = OrderedSet.size(indices);
    for (let t = 0; t < size; t++) {
      const i = OrderedSet.getAt(indices, t);
      cx += x[i];
      cy += y[i];
      cz += z[i];
    }

    if (size > 0) {
      cx /= size;
      cy /= size;
      cz /= size;
    }

    for (let t = 0; t < size; t++) {
      const i = OrderedSet.getAt(indices, t);
      const dx = x[i] - cx,
        dy = y[i] - cy,
        dz = z[i] - cz;
      const d = dx * dx + dy * dy + dz * dz;
      if (d > radiusSq) radiusSq = d;
    }

    return { center: Vec3.create(cx, cy, cz), radius: Math.sqrt(radiusSq) };
  }

  /** Transform sphere with a Mat4 */
  export function transform(out: Sphere3D, sphere: Sphere3D, m: Mat4): Sphere3D {
    Vec3.transformMat4(out.center, sphere.center, m);
    const scale = Mat4.getMaxScaleOnAxis(m);
    out.radius = sphere.radius * scale;
    if (hasExtrema(sphere)) {
      setExtrema(
        out,
        sphere.extrema.map((e) => Vec3.transformMat4(Vec3(), e, m)),
        getExtremaRadius(sphere) * scale,
      );
    }
    return out;
  }

  /** Translate sphere by Vec3 */
  export function translate(out: Sphere3D, sphere: Sphere3D, v: Vec3) {
    Vec3.add(out.center, sphere.center, v);
    if (hasExtrema(sphere)) {
      setExtrema(
        out,
        sphere.extrema.map((e) => Vec3.add(Vec3(), e, v)),
        sphere.extremaRadius,
      );
    }
    return out;
  }

  /** Scale sphere by a number */
  export function scale(out: Sphere3D, sphere: Sphere3D, s: number) {
    Vec3.scale(out.center, sphere.center, s);
    out.radius = sphere.radius * s;
    if (hasExtrema(sphere)) {
      setExtrema(
        out,
        sphere.extrema.map((e) => Vec3.scale(Vec3(), e, s)),
        getExtremaRadius(sphere) * Math.abs(s),
      );
    }
    return out;
  }

  /** Scale sphere by a number but without extrema */
  export function scaleNX(out: Sphere3D, sphere: Sphere3D, s: number) {
    Vec3.scale(out.center, sphere.center, s);
    out.radius = sphere.radius * s;
    return out;
  }

  export function toArray<T extends NumberArray>(s: Sphere3D, out: T, offset: number) {
    Vec3.toArray(s.center, out, offset);
    out[offset + 3] = s.radius;
    return out;
  }

  export function fromArray(out: Sphere3D, array: NumberArray, offset: number) {
    Vec3.fromArray(out.center, array, offset);
    out.radius = array[offset + 3];
    return out;
  }

  export function fromBox3D(out: Sphere3D, box: Box3D) {
    Vec3.scale(out.center, Vec3.add(out.center, box.max, box.min), 0.5);
    out.radius = Vec3.distance(out.center, box.max);

    Sphere3D.setExtrema(out, [
      Vec3.create(box.min[0], box.min[1], box.min[2]),
      Vec3.create(box.max[0], box.max[1], box.max[2]),
      Vec3.create(box.max[0], box.min[1], box.min[2]),
      Vec3.create(box.min[0], box.max[1], box.max[2]),
      Vec3.create(box.min[0], box.min[1], box.max[2]),
      Vec3.create(box.max[0], box.min[1], box.max[2]),
      Vec3.create(box.max[0], box.max[1], box.min[2]),
      Vec3.create(box.min[0], box.max[1], box.min[2]),
    ]);

    return out;
  }

  export function fromAxes3D(out: Sphere3D, axes: Axes3D) {
    Vec3.copy(out.center, axes.origin);
    out.radius = Math.max(Vec3.magnitude(axes.dirA), Vec3.magnitude(axes.dirB), Vec3.magnitude(axes.dirC));
    return out;
  }

  const tmpCenter = Vec3();
  /** Get a tight sphere around a transformed box */
  export function fromDimensionsAndTransform(out: Sphere3D, dimensions: Vec3, transform: Mat4) {
    const [x, y, z] = dimensions;

    const cpA = Vec3.create(0, 0, 0);
    Vec3.transformMat4(cpA, cpA, transform);
    const cpB = Vec3.create(x, y, z);
    Vec3.transformMat4(cpB, cpB, transform);
    const cpC = Vec3.create(x, 0, 0);
    Vec3.transformMat4(cpC, cpC, transform);
    const cpD = Vec3.create(0, y, z);
    Vec3.transformMat4(cpD, cpD, transform);

    const cpE = Vec3.create(0, 0, z);
    Vec3.transformMat4(cpE, cpE, transform);
    const cpF = Vec3.create(x, 0, z);
    Vec3.transformMat4(cpF, cpF, transform);
    const cpG = Vec3.create(x, y, 0);
    Vec3.transformMat4(cpG, cpG, transform);
    const cpH = Vec3.create(0, y, 0);
    Vec3.transformMat4(cpH, cpH, transform);

    Vec3.add(tmpCenter, cpA, cpB);
    Vec3.scale(tmpCenter, tmpCenter, 0.5);
    const d = Math.max(Vec3.distance(cpA, cpB), Vec3.distance(cpC, cpD));
    Sphere3D.set(out, tmpCenter, d / 2);
    Sphere3D.setExtrema(out, [cpA, cpB, cpC, cpD, cpE, cpF, cpG, cpH]);

    return out;
  }

  const tmpAddVec3 = Vec3();
  export function addVec3(out: Sphere3D, s: Sphere3D, v: Vec3) {
    const d = Vec3.distance(s.center, v);
    if (d < s.radius) return Sphere3D.copy(out, s);
    Vec3.sub(tmpAddVec3, s.center, v);
    Vec3.sub(tmpAddVec3, s.center, tmpAddVec3);
    Vec3.setMagnitude(tmpAddVec3, tmpAddVec3, s.radius);
    Vec3.scale(out.center, Vec3.add(tmpAddVec3, tmpAddVec3, v), 0.5);
    out.radius = Vec3.distance(out.center, v);
    return out;
  }

  /** Expand sphere radius by another sphere */
  export function expandBySphere(out: Sphere3D, sphere: Sphere3D, by: Sphere3D) {
    Vec3.copy(out.center, sphere.center);
    out.radius = Math.max(sphere.radius, Vec3.distance(sphere.center, by.center) + by.radius);
    if (hasExtrema(sphere) && hasExtrema(by)) {
      // the larger ball around both's extrema: covers each
      setExtrema(
        out,
        [...sphere.extrema.map((e) => Vec3.clone(e)), ...by.extrema.map((e) => Vec3.clone(e))],
        Math.max(getExtremaRadius(sphere), getExtremaRadius(by)),
      );
    } else {
      out.extrema = void 0;
      out.extremaRadius = void 0;
    }
    return out;
  }

  /**
   * Expand sphere radius by `delta` (not negative): the extrema kept, each one's ball grown by `delta` (`extremaRadius`),
   * which covers the expanded set exactly
   */
  export function expand(out: Sphere3D, sphere: Sphere3D, delta: number): Sphere3D {
    Vec3.copy(out.center, sphere.center);
    out.radius = sphere.radius + delta;
    if (sphere.radius < 1e-12 || (sphere.extrema?.length ?? 0) <= 1) {
      out.extrema = void 0;
      out.extremaRadius = void 0;
      return out;
    }
    if (hasExtrema(sphere)) {
      const extrema = out === sphere ? sphere.extrema : sphere.extrema.map((e) => Vec3.clone(e));
      setExtrema(out, extrema, Math.max(0, getExtremaRadius(sphere) + delta));
    }
    return out;
  }

  /**
   * Returns whether or not the spheres have exactly the same center and radius (when compared with ===)
   */
  export function exactEquals(a: Sphere3D, b: Sphere3D) {
    return a.radius === b.radius && Vec3.exactEquals(a.center, b.center);
  }

  /**
   * Returns whether or not the spheres have approximately the same center and radius.
   */
  export function equals(a: Sphere3D, b: Sphere3D) {
    const ar = a.radius;
    const br = b.radius;
    return Math.abs(ar - br) <= EPSILON * Math.max(1.0, Math.abs(ar), Math.abs(br)) && Vec3.equals(a.center, b.center);
  }

  /**
   * Check if `a` includes `b`, use `extrema` of `b` when available
   */
  export function includes(a: Sphere3D, b: Sphere3D) {
    if (hasExtrema(b)) {
      const r = getExtremaRadius(b);
      for (const e of b.extrema) {
        if (Vec3.distance(a.center, e) + r > a.radius) return false;
      }
      return true;
    } else {
      return Vec3.distance(a.center, b.center) + b.radius <= a.radius;
    }
  }

  /** Check if `a` and `b` are overlapping */
  export function overlaps(a: Sphere3D, b: Sphere3D) {
    return Vec3.distance(a.center, b.center) <= a.radius + b.radius;
  }

  /** Get the signed distance of `a` and `b` */
  export function distance(a: Sphere3D, b: Sphere3D) {
    return Vec3.distance(a.center, b.center) - a.radius + b.radius;
  }

  /** Get the distance of v from sphere. If negative, v is inside sphere */
  export function distanceToVec(sphere: Sphere3D, v: Vec3): number {
    const { center, radius } = sphere;
    return Vec3.distance(v, center) - radius;
  }
}

export { Sphere3D };
