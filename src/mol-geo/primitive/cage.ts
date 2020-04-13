/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { NumberArray } from '../../mol-util/type-helpers';

export interface Cage {
    readonly vertices: ArrayLike<number>
    readonly edges: ArrayLike<number>
}

export function createCage(vertices: ArrayLike<number>, edges: ArrayLike<number>): Cage {
    return { vertices, edges };
}

export function cloneCage(cage: Cage): Cage {
    return {
        vertices: new Float32Array(cage.vertices),
        edges: new Uint32Array(cage.edges)
    };
}

const tmpV = Vec3.zero();

/** Transform primitive in-place */
export function transformCage(cage: Cage, t: Mat4) {
    const { vertices } = cage;
    for (let i = 0, il = vertices.length; i < il; i += 3) {
        // position
        Vec3.transformMat4(tmpV, Vec3.fromArray(tmpV, vertices, i), t);
        Vec3.toArray(tmpV, vertices as NumberArray, i);
    }
    return cage;
}