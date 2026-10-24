/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { addCylinder, addSimpleCylinder } from '../../geo/geometry/mesh/builder/cylinder.js';
import { addSphere } from '../../geo/geometry/mesh/builder/sphere.js';
import { MeshBuilder } from '../../geo/geometry/mesh/mesh-builder.js';
import { Torus } from '../../geo/primitive/torus.js';

// The parts handles are drawn from, into a `MeshBuilder` state at its current group (one group a handle). Built in the
// item's own frame; the controller places the frame (`facingCamera` for parts that should always face the viewer).

const Up = Vec3.create(0, 1, 0);

/** A knob (a sphere) of radius `r` at `at` */
export function knob(state: MeshBuilder.State, at: Vec3, r: number) {
  addSphere(state, at, r, 2);
}

/** A tube from `a` to `b` of radius `r` */
export function tube(state: MeshBuilder.State, a: Vec3, b: Vec3, r: number) {
  addCylinder(state, a, b, 1, { radiusTop: r, radiusBottom: r, radialSegments: 12 });
}

/** An arrow from `from` to `to`: a shaft of radius `r`, a head at its last fifth */
export function arrow(state: MeshBuilder.State, from: Vec3, to: Vec3, r: number) {
  const shaftEnd = Vec3.lerp(Vec3(), from, to, 0.8);
  addCylinder(state, from, shaftEnd, 1, { radiusTop: r, radiusBottom: r, radialSegments: 12 });
  // as MolViewSpec's arrows: `addSimpleCylinder` puts `radiusBottom` at its start (`addCylinder` does not)
  addSimpleCylinder(state, shaftEnd, to, { radiusTop: 0, radiusBottom: r * 3, radialSegments: 16, bottomCap: true });
}

/** A ring of radius `r` and tube `width` around the frame's up (+Y), at height `y` (or placed by `t`) */
export function ring(state: MeshBuilder.State, r: number, width: number, y = 0, t?: Mat4) {
  const m = t ?? Mat4.mul(Mat4(), Mat4.fromTranslation(Mat4(), Vec3.create(0, y, 0)), Mat4.rotX90 as Mat4);
  MeshBuilder.addPrimitive(state, m, Torus({ radius: r, tube: width, radialSegments: 8, tubularSegments: 64 }));
}

/** A frame at `transform`'s place turned so its up (+Y) points from there to the camera (`eye`) */
export function facingCamera(transform: Mat4, eye: Vec3): Mat4 {
  const at = Vec3.fromArray(Vec3(), transform, 12);
  const toEye = Vec3.sub(Vec3(), eye, at);
  if (Vec3.magnitude(toEye) < 1e-9) return transform;
  Vec3.normalize(toEye, toEye);
  const axis = Vec3.cross(Vec3(), Up, toEye);
  const turn =
    Vec3.magnitude(axis) < 1e-6
      ? Vec3.dot(Up, toEye) > 0
        ? Mat4.identity()
        : Mat4.fromRotation(Mat4(), Math.PI, Vec3.create(1, 0, 0))
      : Mat4.fromRotation(
          Mat4(),
          Math.acos(Math.max(-1, Math.min(1, Vec3.dot(Up, toEye)))),
          Vec3.normalize(axis, axis),
        );
  return Mat4.mul(Mat4(), Mat4.fromTranslation(Mat4(), at), turn);
}
