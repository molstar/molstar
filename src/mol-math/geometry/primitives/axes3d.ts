/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3, Mat4, Mat3 } from '../../linear-algebra';

interface Axes3D { origin: Vec3, dirA: Vec3, dirB: Vec3, dirC: Vec3 }

function Axes3D() {
    return Axes3D.empty();
}

namespace Axes3D {
    export function create(origin: Vec3, dirA: Vec3, dirB: Vec3, dirC: Vec3): Axes3D { return { origin, dirA, dirB, dirC }; }
    export function empty(): Axes3D { return { origin: Vec3(), dirA: Vec3(), dirB: Vec3(), dirC: Vec3() }; }

    export function copy(out: Axes3D, a: Axes3D): Axes3D {
        Vec3.copy(out.origin, a.origin);
        Vec3.copy(out.dirA, a.dirA);
        Vec3.copy(out.dirB, a.dirB);
        Vec3.copy(out.dirC, a.dirC);
        return out;
    }

    export function clone(a: Axes3D): Axes3D {
        return copy(empty(), a);
    }

    /** Get size of each direction */
    export function size(size: Vec3, axes: Axes3D): Vec3 {
        return Vec3.set(size, Vec3.magnitude(axes.dirA) * 2, Vec3.magnitude(axes.dirB) * 2, Vec3.magnitude(axes.dirC) * 2);
    }

    const tmpSizeV = Vec3();
    /** Get volume of the oriented box wrapping the axes */
    export function volume(axes: Axes3D): number {
        size(tmpSizeV, axes);
        return tmpSizeV[0] * tmpSizeV[1] * tmpSizeV[2];
    }

    export function normalize(out: Axes3D, a: Axes3D) {
        Vec3.copy(out.origin, a.origin);
        Vec3.normalize(out.dirA, a.dirA);
        Vec3.normalize(out.dirB, a.dirB);
        Vec3.normalize(out.dirC, a.dirC);
        return out;
    }

    const tmpTransformMat3 = Mat3.zero();
    /** Transform axes with a Mat4 */
    export function transform(out: Axes3D, a: Axes3D, m: Mat4): Axes3D {
        Vec3.transformMat4(out.origin, a.origin, m);
        const n = Mat3.directionTransform(tmpTransformMat3, m);
        Vec3.transformMat3(out.dirA, a.dirA, n);
        Vec3.transformMat3(out.dirB, a.dirB, n);
        Vec3.transformMat3(out.dirC, a.dirC, n);
        return out;
    }
}

export { Axes3D };