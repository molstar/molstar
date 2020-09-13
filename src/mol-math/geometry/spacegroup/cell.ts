/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '../../linear-algebra';

export { Cell };

interface Cell {
    readonly size: Vec3
    readonly anglesInRadians: Vec3
}

function Cell() {
    return Cell.empty();
}

namespace Cell {
    export function create(size: Vec3, anglesInRadians: Vec3): Cell {
        return { size, anglesInRadians };
    }

    export function empty(): Cell {
        return create(Vec3(), Vec3());
    }

    export function fromBasis(x: Vec3, y: Vec3, z: Vec3) {
        const a = Vec3.magnitude(x);
        const b = Vec3.magnitude(y);
        const c = Vec3.magnitude(z);

        const alpha = Math.acos(Vec3.dot(y, z) / (b * c));
        const beta = Math.acos(Vec3.dot(x, z) / (a * c));
        const gamma = Math.acos(Vec3.dot(x, y) / (a * b));

        if (a <= 0 || b <= 0 || c <= 0 || alpha >= Math.PI || beta >= Math.PI || gamma >= Math.PI) {
            return empty();
        } else {
            return create(Vec3.create(a, b, c), Vec3.create(alpha, beta, gamma));
        }
    }
}