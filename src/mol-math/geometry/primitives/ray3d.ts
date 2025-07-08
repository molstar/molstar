/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4 } from '../../linear-algebra/3d/mat4';
import { Vec3 } from '../../linear-algebra/3d/vec3';

interface Ray3D { origin: Vec3, direction: Vec3 }

function Ray3D() {
    return Ray3D.create(Vec3.create(0, 0, 0), Vec3.create(1, 0, 0));
}

namespace Ray3D {
    export function create(origin: Vec3, direction: Vec3): Ray3D { return { origin, direction }; }

    export function copy(out: Ray3D, r: Ray3D): Ray3D {
        Vec3.copy(out.origin, r.origin);
        Vec3.copy(out.direction, r.direction);
        return out;
    }

    export function clone(r: Ray3D): Ray3D {
        return copy(Ray3D(), r);
    }

    export function targetTo(out: Ray3D, ray: Ray3D, target: Vec3): Ray3D {
        Vec3.copy(out.origin, ray.origin);
        Vec3.normalize(out.direction, Vec3.sub(out.direction, target, ray.origin));
        return out;
    }

    /** Transform ray with a Mat4 */
    export function transform(out: Ray3D, ray: Ray3D, m: Mat4): Ray3D {
        Vec3.transformMat4(out.origin, ray.origin, m);
        Vec3.transformDirection(out.direction, ray.direction, m);
        return out;
    }
}

export { Ray3D };
