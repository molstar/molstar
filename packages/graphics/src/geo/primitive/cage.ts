/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Mat4, Vec3 } from '@molstar/core/math/linear-algebra';

export interface Cage {
  /** Float32Array, so that loops over it keep one array kind (`Vec3.fromFloat32Array`) */
  readonly vertices: Float32Array;
  readonly edges: ArrayLike<number>;
}

/** Keeps a Float32Array of `vertices` as it is; copies other arrays into one */
export function createCage(vertices: ArrayLike<number>, edges: ArrayLike<number>): Cage {
  return { vertices: vertices instanceof Float32Array ? vertices : new Float32Array(vertices), edges };
}

export function cloneCage(cage: Cage): Cage {
  return {
    vertices: new Float32Array(cage.vertices),
    edges: new Uint32Array(cage.edges),
  };
}

const tmpV = Vec3.zero();

/** Transform primitive in-place */
export function transformCage(cage: Cage, t: Mat4) {
  const { vertices } = cage;
  for (let i = 0, il = vertices.length; i < il; i += 3) {
    // position
    Vec3.transformMat4(tmpV, Vec3.fromFloat32Array(tmpV, vertices, i), t);
    Vec3.toFloat32Array(tmpV, vertices, i);
  }
  return cage;
}
