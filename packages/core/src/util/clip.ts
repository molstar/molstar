/**
 * Copyright (c) 2021-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { EPSILON } from '@molstar/core/math/linear-algebra/3d/common';
import { Mat3 } from '@molstar/core/math/linear-algebra/3d/mat3';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { Quat } from '@molstar/core/math/linear-algebra/3d/quat';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { Vec4 } from '@molstar/core/math/linear-algebra/3d/vec4';
import { Plane3D } from '@molstar/core/math/geometry/primitives/plane3d';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { degToRad } from '@molstar/core/math/misc';
import { ParamDefinition as PD } from './param-definition.js';
import { stringToWords } from './string.js';

export interface Clip {
  variant: Clip.Variant;
  objects: Clip.Objects;
}

export function Clip() {}

export namespace Clip {
  /** Clip object types */
  export const Type = {
    none: 0, // to switch clipping off
    plane: 1,
    sphere: 2,
    cube: 3,
    cylinder: 4,
    infiniteCone: 5,
  };

  export type Variant = 'instance' | 'pixel';

  /** additional objects are ignored */
  export const MaxObjectCount = 16;

  export type Objects = {
    count: number;
    type: number[];
    invert: boolean[];
    position: number[];
    rotation: number[];
    scale: number[];
    /** Transform point by this before testing */
    transform: number[];
    /**
     * Per object the transform of a point into the object's frame, `count` 4x4 matrices: `transform`, less `position`,
     * rotated back by `rotation` (`getFrame`); where the shaders test the shape. The first column's `w` (a point is
     * read with `.xyz`) is the frame's largest scale (`classifySphereObject`)
     */
    frame: number[];
    /**
     * Per object the shape in its frame, `count` 4-vectors: half the scale, the type in `w`, negative where inverted
     * (`getShape`)
     */
    shape: number[];
  };

  export type Rotation = { axis: Vec3; angle: number };

  export const Params = {
    variant: PD.Select('pixel', PD.arrayToOptions<Variant>(['instance', 'pixel'])),
    objects: PD.ObjectList(
      {
        type: PD.Select(
          'plane',
          PD.objectToOptions(Type, (t) => stringToWords(t)),
        ),
        invert: PD.Boolean(false),
        position: PD.Vec3(Vec3()),
        rotation: PD.Group<Rotation>(
          {
            axis: PD.Vec3(Vec3.create(1, 0, 0)),
            angle: PD.Numeric(0, { min: -180, max: 180, step: 1 }, { description: 'Angle in Degrees' }),
          },
          { isExpanded: true },
        ),
        scale: PD.Vec3(Vec3.create(1, 1, 1)),
        transform: PD.Mat4(Mat4.identity()),
        /**
         * Who made it, and which of theirs it is (a clip manager's `clip-<uuid>`): the same object on several
         * representations is one cut by it. Not used for clipping; empty: untagged.
         */
        tag: PD.Optional(PD.Text('', { isHidden: true })),
      },
      (o) => stringToWords(o.type),
    ),
  };
  export type Params = typeof Params;
  export type Props = PD.Values<Params>;

  function createClipObjects(count: number) {
    return {
      count: 0,
      type: new Array(count).fill(1),
      invert: new Array(count).fill(false),
      position: new Array(count * 3).fill(0),
      rotation: new Array(count * 4).fill(0),
      scale: new Array(count * 3).fill(1),
      transform: new Array(count * 16).fill(0),
      frame: new Array(count * 16).fill(0),
      shape: new Array(count * 4).fill(0),
    };
  }

  const mF = Mat4();
  const mT = Mat4();
  const m3A = Mat3();
  const m3B = Mat3();
  const v3Eigen = Vec3();
  /**
   * How much the linear part of `m` lengthens a vector at most, its spectral norm (the square root of the largest
   * eigenvalue of its transpose times it): the longest column is less under a shear or a scale after a rotation. A little
   * more, as a bound: what it decides is read as 32-bit floats
   */
  function getLargestScale(m: Mat4) {
    Mat3.fromMat4(m3A, m);
    Mat3.mul(m3B, Mat3.transpose(Mat3(), m3A), m3A);
    Mat3.symmetricEigenvalues(v3Eigen, m3B);
    return Math.sqrt(Math.max(v3Eigen[0], v3Eigen[1], v3Eigen[2], 0)) * (1 + 1e-6);
  }

  /** The frame of the object at `index` of `objects` (see `Objects.frame`) into `out` at `index` */
  function setFrame(objects: Objects, index: number) {
    const { position, rotation, transform, frame } = objects;
    Quat.fromArray(qA, rotation, index * 4);
    Quat.conjugate(qA, qA);
    Mat4.fromQuat(mF, qA);
    Mat4.fromTranslation(mT, Vec3.set(vB, -position[index * 3], -position[index * 3 + 1], -position[index * 3 + 2]));
    Mat4.mul(mF, mF, mT);
    Mat4.mul(mF, mF, Mat4.fromArray(mT, transform, index * 16));
    Mat4.toArray(mF, frame, index * 16);
    frame[index * 16 + 3] = getLargestScale(mF);
  }
  /** The shape of the object at `index` of `objects` (see `Objects.shape`) */
  function setShape(objects: Objects, index: number) {
    const { type, invert, scale, shape } = objects;
    const o = index * 4;
    shape[o] = scale[index * 3] * 0.5;
    shape[o + 1] = scale[index * 3 + 1] * 0.5;
    shape[o + 2] = scale[index * 3 + 2] * 0.5;
    shape[o + 3] = invert[index] ? -type[index] : type[index];
  }

  const qA = Quat();
  const qB = Quat();
  const vA = Vec3();
  const vB = Vec3();
  const mA = Mat4();
  const mB = Mat4();

  export function getClip(props: Props, clip?: Clip): Clip {
    const count = Math.min(props.objects.length, MaxObjectCount);
    const objects = clip?.objects || createClipObjects(count);
    const { type, invert, position, rotation, scale, transform } = objects;
    for (let i = 0; i < count; ++i) {
      const p = props.objects[i];
      type[i] = Type[p.type];
      invert[i] = p.invert;
      Vec3.toArray(p.position, position, i * 3);
      Vec3.normalize(vA, p.rotation.axis);
      Quat.toArray(Quat.setAxisAngle(qA, vA, degToRad(p.rotation.angle)), rotation, i * 4);
      Vec3.toArray(p.scale, scale, i * 3);
      Mat4.toArray(p.transform, transform, i * 16);
      setFrame(objects, i);
      setShape(objects, i);
    }
    objects.count = count;
    return { variant: props.variant, objects };
  }

  const v4A = Vec4();

  /**
   * Plane of the clip object at `index`, in the space of the points
   * given to the clip test (i.e. with `transform` undone) and with
   * the clipped half-space on its positive side
   */
  export function getPlane(out: Plane3D, objects: Objects, index: number): Plane3D {
    const { invert, position, rotation, transform } = objects;
    Vec3.transformQuat(vA, Vec3.unitY, Quat.fromArray(qA, rotation, index * 4));
    Vec4.set(v4A, vA[0], vA[1], vA[2], -Vec3.dot(vA, Vec3.fromArray(vB, position, index * 3)));
    if (invert[index]) Vec4.scale(v4A, v4A, -1);
    Mat4.fromArray(mA, transform, index * 16);
    Vec4.transformMat4(v4A, v4A, Mat4.transpose(mA, mA));
    return Plane3D.setUnnormalized(out, v4A[0], v4A[1], v4A[2], v4A[3]);
  }

  /** the smallest half axis an ellipsoid or an elliptic cylinder is measured with (`getSignedDistance`) */
  const MinHalfAxis = 1e-6;

  /**
   * The signed distance of `p` to the shape of the clip object at `index` in its frame (`frame`, `shape`), as the
   * shaders' `localSD`: at most 0 inside (where an object not inverted clips), not a true distance for all shapes (see
   * `classifySphere`)
   */
  export function getSignedDistance(objects: Objects, index: number, p: Vec3): number {
    const { frame, shape } = objects;
    const f = index * 16;
    const x = frame[f] * p[0] + frame[f + 4] * p[1] + frame[f + 8] * p[2] + frame[f + 12];
    const y = frame[f + 1] * p[0] + frame[f + 5] * p[1] + frame[f + 9] * p[2] + frame[f + 13];
    const z = frame[f + 2] * p[0] + frame[f + 6] * p[1] + frame[f + 10] * p[2] + frame[f + 14];
    const o = index * 4;
    const sx = shape[o],
      sy = shape[o + 1],
      sz = shape[o + 2];
    const type = Math.round(Math.abs(shape[o + 3]));
    if (type === Type.plane) return -y;
    if (type === Type.sphere) {
      // half axes of 0 (a point, a disc) at a tiny size: finite, not 0 / 0 (as the shaders' `localSD`)
      const ax = Math.max(sx, MinHalfAxis),
        ay = Math.max(sy, MinHalfAxis),
        az = Math.max(sz, MinHalfAxis);
      return (Math.hypot(x / ax, y / ay, z / az) - 1) * Math.min(ax, ay, az);
    }
    if (type === Type.cube) {
      const dx = Math.abs(x) - sx,
        dy = Math.abs(y) - sy,
        dz = Math.abs(z) - sz;
      return Math.min(Math.max(dx, dy, dz), 0) + Math.hypot(Math.max(dx, 0), Math.max(dy, 0), Math.max(dz, 0));
    }
    if (type === Type.cylinder) {
      // elliptic (radii `sx`, `sz`): across its axis scaled by the smaller (at most the true distance, exact round)
      const ax = Math.max(sx, MinHalfAxis),
        az = Math.max(sz, MinHalfAxis);
      const dx = (Math.hypot(x / ax, z / az) - 1) * Math.min(ax, az),
        dy = Math.abs(y) - sy;
      return Math.min(Math.max(dx, dy), 0) + Math.hypot(Math.max(dx, 0), Math.max(dy, 0));
    }
    if (type === Type.infiniteCone) return sx * Math.hypot(x, y) + sy * z;
    return 0.1;
  }

  /** How a sphere lies against clip objects (`classifySphere`): no point clipped, some, every point */
  export const SphereClass = { None: 0, Some: 1, All: 2 } as const;
  export type SphereClass = (typeof SphereClass)[keyof typeof SphereClass];

  const vC = Vec3();
  /**
   * How `sphere` lies against the clip object at `index`: from the signed distance of its center, which divided by its
   * largest slope (1; the infinite cone's `|shape.xy|`; the ellipsoid's and the elliptic cylinder's distance already
   * scaled by their smallest axis)
   * bounds the true one, the radius scaled by the frame's largest scale: `All` within the clipped side by more, `None`
   * outside it by more, `Some` else (exact where it decides, conservative elsewhere); as the sphere raster's
   * `clipSphereClass`
   */
  export function classifySphereObject(objects: Objects, index: number, sphere: Sphere3D): SphereClass {
    const { frame, shape } = objects;
    const scale = frame[index * 16 + 3];
    const o = index * 4;
    const type = Math.round(Math.abs(shape[o + 3]));
    if (type === Type.none) return SphereClass.None;
    const slope = type === Type.infiniteCone ? Math.hypot(shape[o], shape[o + 1]) : 1;
    const r = sphere.radius * scale * slope;
    const sd = getSignedDistance(objects, index, Vec3.copy(vC, sphere.center));
    // within the clipped side at least this far (it is `sd <= 0`, inverted `sd > 0`)
    const depth = shape[o + 3] < 0 ? sd : -sd;
    if (depth > r) return SphereClass.All;
    if (depth >= -r) return SphereClass.Some;
    return SphereClass.None;
  }

  /**
   * How `sphere` lies against all clip objects (`classifySphereObject`; a point is clipped by any): `All` where one clips
   * every point, `None` where none clips any, `Some` else; `flag` the objects to ignore, a bit each (a clipping layer's)
   */
  export function classifySphere(objects: Objects, sphere: Sphere3D, flag = 0): SphereClass {
    let some = false;
    for (let i = 0, il = objects.count; i < il; ++i) {
      if (flag & (1 << i)) continue;
      const c = classifySphereObject(objects, i, sphere);
      if (c === SphereClass.All) return SphereClass.All;
      if (c === SphereClass.Some) some = true;
    }
    return some ? SphereClass.Some : SphereClass.None;
  }

  /**
   * Test if the surface of the clip object at `index` can intersect `sphere` (`classifySphereObject`'s `Some`): exact
   * where the distance is (planes, cubes, round cylinders and spheres), conservative elsewhere
   */
  export function canIntersectSphere(objects: Objects, index: number, sphere: Sphere3D): boolean {
    return classifySphereObject(objects, index, sphere) === SphereClass.Some;
  }

  export function areEqual(cA: Clip, cB: Clip) {
    if (cA.variant !== cB.variant) return false;
    if (cA.objects.count !== cB.objects.count) return false;

    const oA = cA.objects,
      oB = cB.objects;
    for (let i = 0, il = oA.count; i < il; ++i) {
      if (oA.invert[i] !== oB.invert[i]) return false;
      if (oA.type[i] !== oB.type[i]) return false;

      Vec3.fromArray(vA, oA.position, i * 3);
      Vec3.fromArray(vB, oB.position, i * 3);
      if (!Vec3.equals(vA, vB)) return false;

      Vec3.fromArray(vA, oA.scale, i * 3);
      Vec3.fromArray(vB, oB.scale, i * 3);
      if (!Vec3.equals(vA, vB)) return false;

      Quat.fromArray(qA, oA.rotation, i * 4);
      Quat.fromArray(qB, oB.rotation, i * 4);
      if (!Quat.equals(qA, qB)) return false;

      Mat4.fromArray(mA, oA.transform, i * 16);
      Mat4.fromArray(mB, oB.transform, i * 16);
      if (!Mat4.areEqual(mA, mB, EPSILON)) return false;
    }
    return true;
  }
}
