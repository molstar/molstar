/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { addCylinder, addSimpleCylinder } from '../../geo/geometry/mesh/builder/cylinder.js';
import { addSphere } from '../../geo/geometry/mesh/builder/sphere.js';
import type { LinesBuilder } from '../../geo/geometry/lines/lines-builder.js';
import { MeshBuilder } from '../../geo/geometry/mesh/mesh-builder.js';
import { Torus } from '../../geo/primitive/torus.js';

// The parts handles are drawn from, into a `MeshBuilder` state at its current group (one group a handle). Built in the
// item's own frame; the controller places the frame (what faces the viewer is points: turned to it by nature).

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
  const shaftEnd = arrowShaftEnd(from, to);
  addCylinder(state, from, shaftEnd, 1, { radiusTop: r, radiusBottom: r, radialSegments: 12 });
  arrowHead(state, from, to, r);
}

/** Where an arrow's shaft ends and its head starts */
export function arrowShaftEnd(from: Vec3, to: Vec3, out = Vec3()) {
  return Vec3.lerp(out, from, to, 0.8);
}

/** An arrow's head alone (its shaft drawn otherwise, e.g. as a line), as `arrow` draws it for a shaft of radius `r` */
export function arrowHead(state: MeshBuilder.State, from: Vec3, to: Vec3, r: number) {
  // as MolViewSpec's arrows: `addSimpleCylinder` puts `radiusBottom` at its start (`addCylinder` does not)
  const shaftEnd = arrowShaftEnd(from, to);
  addSimpleCylinder(state, shaftEnd, to, { radiusTop: 0, radiusBottom: r * 3, radialSegments: 16, bottomCap: true });
}

/** A ring of radius `r` around the frame's up (+Y) as a closed polyline of `segments` lines (placed by `t`) */
export function ringLines(lines: LinesBuilder, r: number, group: number, t?: Mat4, segments = 64) {
  ellipseLines(lines, r, r, group, t, segments);
}

/** An ellipse of half axes `rx` (along X) and `rz` (along Z) in the frame's XZ plane, as `ringLines` (placed by `t`) */
export function ellipseLines(lines: LinesBuilder, rx: number, rz: number, group: number, t?: Mat4, segments = 64) {
  const p = Vec3();
  const q = Vec3();
  for (let k = 0; k < segments; k++) {
    const a0 = (2 * Math.PI * k) / segments;
    const a1 = (2 * Math.PI * (k + 1)) / segments;
    Vec3.set(p, rx * Math.cos(a0), 0, rz * Math.sin(a0));
    Vec3.set(q, rx * Math.cos(a1), 0, rz * Math.sin(a1));
    if (t) {
      Vec3.transformMat4(p, p, t);
      Vec3.transformMat4(q, q, t);
    }
    lines.add(p[0], p[1], p[2], q[0], q[1], q[2], group);
  }
}

/** A ring of radius `r` and tube `width` around the frame's up (+Y), at height `y` (or placed by `t`) */
export function ring(state: MeshBuilder.State, r: number, width: number, y = 0, t?: Mat4) {
  const m = t ?? Mat4.mul(Mat4(), Mat4.fromTranslation(Mat4(), Vec3.create(0, y, 0)), Mat4.rotX90 as Mat4);
  MeshBuilder.addPrimitive(state, m, Torus({ radius: r, tube: width, radialSegments: 8, tubularSegments: 64 }));
}

/** A flat disk of radius `r` across the frame's up (+Y), at height `y` */
export function disk(state: MeshBuilder.State, r: number, y: number) {
  const a = Vec3.create(0, y, 0);
  const b = Vec3.create(0, y + r * 0.01, 0);
  addCylinder(state, a, b, 1, { radiusTop: r, radiusBottom: r, radialSegments: 48, topCap: true, bottomCap: true });
}
