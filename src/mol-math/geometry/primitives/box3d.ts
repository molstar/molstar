/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PositionData } from '../common';
import { OrderedSet } from '../../../mol-data/int';
import { Sphere3D } from './sphere3d';
import { Vec3 } from '../../linear-algebra/3d/vec3';
import { Mat4 } from '../../linear-algebra/3d/mat4';

interface Box3D { min: Vec3, max: Vec3 }

function Box3D() {
    return Box3D.zero();
}

namespace Box3D {
    export function create(min: Vec3, max: Vec3): Box3D { return { min, max }; }
    export function zero(): Box3D { return { min: Vec3(), max: Vec3() }; }

    export function copy(out: Box3D, a: Box3D): Box3D {
        Vec3.copy(out.min, a.min);
        Vec3.copy(out.max, a.max);
        return out;
    }

    export function set(out: Box3D, min: Vec3, max: Vec3): Box3D {
        Vec3.copy(out.min, min);
        Vec3.copy(out.max, max);
        return out;
    }

    export function clone(a: Box3D): Box3D {
        return copy(zero(), a);
    }

    const tmpV = Vec3();

    /** Get box from sphere, uses extrema if available */
    export function fromSphere3D(out: Box3D, sphere: Sphere3D): Box3D {
        if (Sphere3D.hasExtrema(sphere) && sphere.extrema.length >= 14) { // 14 extrema with coarse boundary helper
            return fromVec3Array(out, sphere.extrema);
        }
        Vec3.set(tmpV, sphere.radius, sphere.radius, sphere.radius);
        Vec3.sub(out.min, sphere.center, tmpV);
        Vec3.add(out.max, sphere.center, tmpV);
        return out;
    }

    export function addVec3Array(out: Box3D, array: Vec3[]): Box3D {
        for (let i = 0, il = array.length; i < il; i++) {
            add(out, array[i]);
        }
        return out;
    }

    export function fromVec3Array(out: Box3D, array: Vec3[]): Box3D {
        setEmpty(out);
        addVec3Array(out, array);
        return out;
    }

    export function addSphere3D(out: Box3D, sphere: Sphere3D): Box3D {
        if (Sphere3D.hasExtrema(sphere) && sphere.extrema.length >= 14) { // 14 extrema with coarse boundary helper
            return addVec3Array(out, sphere.extrema);
        }
        add(out, Vec3.subScalar(tmpV, sphere.center, sphere.radius));
        add(out, Vec3.addScalar(tmpV, sphere.center, sphere.radius));
        return out;
    }

    export function addBox3D(out: Box3D, box: Box3D): Box3D {
        add(out, box.min);
        add(out, box.max);
        return out;
    }

    export function intersectsSphere3D(box: Box3D, sphere: Sphere3D) {
        // Find the point on the AABB closest to the sphere center.
        Vec3.clamp(tmpV, sphere.center, box.min, box.max);
        // If that point is inside the sphere, the AABB and sphere intersect.
        return Vec3.squaredDistance(tmpV, sphere.center) <= (sphere.radius * sphere.radius);
    }

    export function computeBounding(data: PositionData): Box3D {
        const min = Vec3.create(Number.MAX_VALUE, Number.MAX_VALUE, Number.MAX_VALUE);
        const max = Vec3.create(-Number.MAX_VALUE, -Number.MAX_VALUE, -Number.MAX_VALUE);

        const { x, y, z, indices } = data;
        for (let t = 0, _t = OrderedSet.size(indices); t < _t; t++) {
            const i = OrderedSet.getAt(indices, t);
            min[0] = Math.min(x[i], min[0]);
            min[1] = Math.min(y[i], min[1]);
            min[2] = Math.min(z[i], min[2]);
            max[0] = Math.max(x[i], max[0]);
            max[1] = Math.max(y[i], max[1]);
            max[2] = Math.max(z[i], max[2]);
        }
        return { min, max };
    }

    /** Get size/extent of the box */
    export function size(size: Vec3, box: Box3D): Vec3 {
        return Vec3.sub(size, box.max, box.min);
    }

    const tmpSizeV = Vec3();
    /** Get volume of the box */
    export function volume(box: Box3D): number {
        size(tmpSizeV, box);
        return tmpSizeV[0] * tmpSizeV[1] * tmpSizeV[2];
    }

    /** Sets min to Number.MAX_VALUE and max to -Number.MAX_VALUE */
    export function setEmpty(box: Box3D): Box3D {
        Vec3.set(box.min, Number.MAX_VALUE, Number.MAX_VALUE, Number.MAX_VALUE);
        Vec3.set(box.max, -Number.MAX_VALUE, -Number.MAX_VALUE, -Number.MAX_VALUE);
        return box;
    }

    /** Add point to box */
    export function add(box: Box3D, point: Vec3): Box3D {
        Vec3.min(box.min, box.min, point);
        Vec3.max(box.max, box.max, point);
        return box;
    }

    /** Expand box by delta */
    export function expand(out: Box3D, box: Box3D, delta: Vec3): Box3D {
        Vec3.sub(out.min, box.min, delta);
        Vec3.add(out.max, box.max, delta);
        return out;
    }

    export function expandUniformly(out: Box3D, box: Box3D, delta: number): Box3D {
        Vec3.subScalar(out.min, box.min, delta);
        Vec3.addScalar(out.max, box.max, delta);
        return out;
    }

    export function scale(out: Box3D, box: Box3D, scale: number) {
        Vec3.scale(out.min, box.min, scale);
        Vec3.scale(out.max, box.max, scale);
        return out;
    }

    const tmpTransformV = Vec3();
    /** Transform box with a Mat4 */
    export function transform(out: Box3D, box: Box3D, m: Mat4): Box3D {
        const [minX, minY, minZ] = box.min;
        const [maxX, maxY, maxZ] = box.max;
        setEmpty(out);
        add(out, Vec3.transformMat4(tmpTransformV, Vec3.set(tmpTransformV, minX, minY, minZ), m));
        add(out, Vec3.transformMat4(tmpTransformV, Vec3.set(tmpTransformV, minX, minY, maxZ), m));
        add(out, Vec3.transformMat4(tmpTransformV, Vec3.set(tmpTransformV, minX, maxY, minZ), m));
        add(out, Vec3.transformMat4(tmpTransformV, Vec3.set(tmpTransformV, minX, maxY, maxZ), m));
        add(out, Vec3.transformMat4(tmpTransformV, Vec3.set(tmpTransformV, maxX, minY, minZ), m));
        add(out, Vec3.transformMat4(tmpTransformV, Vec3.set(tmpTransformV, maxX, minY, maxZ), m));
        add(out, Vec3.transformMat4(tmpTransformV, Vec3.set(tmpTransformV, maxX, maxY, minZ), m));
        add(out, Vec3.transformMat4(tmpTransformV, Vec3.set(tmpTransformV, maxX, maxY, maxZ), m));
        return out;
    }

    export function containsVec3(box: Box3D, v: Vec3) {
        return !(
            v[0] < box.min[0] || v[0] > box.max[0] ||
            v[1] < box.min[1] || v[1] > box.max[1] ||
            v[2] < box.min[2] || v[2] > box.max[2]
        );
    }

    export function overlaps(a: Box3D, b: Box3D) {
        return !(
            a.max[0] < b.min[0] || a.min[0] > b.max[0] ||
            a.max[1] < b.min[1] || a.min[1] > b.max[1] ||
            a.max[2] < b.min[2] || a.min[2] > b.max[2]
        );
    }

    export function containsSphere3D(box: Box3D, s: Sphere3D) {
        const c = s.center;
        const r = s.radius;
        return (
            c[0] - r < box.min[0] || c[0] + r > box.max[0] ||
            c[1] - r < box.min[1] || c[1] + r > box.max[1] ||
            c[2] - r < box.min[2] || c[2] + r > box.max[2]
        ) ? false : true;
    }

    export function center(out: Vec3, box: Box3D): Vec3 {
        return Vec3.center(out, box.max, box.min);
    }

    export function exactEquals(a: Box3D, b: Box3D) {
        return Vec3.exactEquals(a.min, b.min) && Vec3.exactEquals(a.max, b.max);
    }

    export function equals(a: Box3D, b: Box3D) {
        return Vec3.equals(a.min, b.min) && Vec3.equals(a.max, b.max);
    }
}

export { Box3D };
