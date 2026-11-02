/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import type { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { Vec4 } from '@molstar/core/math/linear-algebra/3d/vec4';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v4set = Vec4.set;

/**
 * Bounding rectangle of a clipped, perspective-projected 3D Sphere.
 * Michael Mara, Morgan McGuire. 2013
 *
 * Specialization by Arseny Kapoulkine, MIT License Copyright (c) 2018
 * https://github.com/zeux/niagara
 */
function perspectiveProjectSphere(out: Vec4, p: Vec3, r: number, projection: Mat4) {
  const prx = p[0] * r;
  const pry = p[1] * r;
  const prz = p[2] * r;
  const pzr2 = p[2] * p[2] - r * r;

  const vx = Math.sqrt(p[0] * p[0] + pzr2);
  const minx = ((vx * p[0] - prz) / (vx * p[2] + prx)) * projection[0];
  const maxx = ((vx * p[0] + prz) / (vx * p[2] - prx)) * projection[0];

  const vy = Math.sqrt(p[1] * p[1] + pzr2);
  const miny = ((vy * p[1] - prz) / (vy * p[2] + pry)) * projection[5];
  const maxy = ((vy * p[1] + prz) / (vy * p[2] - pry)) * projection[5];

  return v4set(out, maxx * -0.5 + 0.5, miny * 0.5 + 0.5, minx * -0.5 + 0.5, maxy * 0.5 + 0.5);
}

function orthographicProjectSphere(out: Vec4, p: Vec3, r: number, projection: Mat4) {
  const sx = projection[0];
  const sy = projection[5];

  const minx = (p[0] + r) * sx;
  const maxx = (p[0] - r) * sx;
  const miny = (p[1] + r) * sy;
  const maxy = (p[1] - r) * sy;

  return v4set(out, maxx * 0.5 + 0.5, miny * -0.5 + 0.5, minx * 0.5 + 0.5, maxy * -0.5 + 0.5);
}

/**
 * The screen rectangle of the sphere at `p` (view space) of radius `r`, as (min x, min y, max x, max y) in [0, 1] of the
 * viewport, y from the top; a perspective projection needs the sphere in front of the eye
 */
export function projectSphere(out: Vec4, p: Vec3, r: number, projection: Mat4) {
  return projection[11] === -1
    ? perspectiveProjectSphere(out, p, r, projection)
    : orthographicProjectSphere(out, p, r, projection);
}
