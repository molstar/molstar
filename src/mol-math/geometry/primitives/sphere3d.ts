/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3, Mat4, EPSILON } from '../../linear-algebra';
import { PositionData } from '../common';
import { OrderedSet } from '../../../mol-data/int';
import { NumberArray, PickRequired } from '../../../mol-util/type-helpers';
import { Box3D } from './box3d';
import { Axes3D } from './axes3d';

interface Sphere3D {
    center: Vec3,
    radius: number,
    extrema?: Vec3[]
}

function Sphere3D() {
    return Sphere3D.zero();
}

namespace Sphere3D {
    export function hasExtrema(sphere: Sphere3D): sphere is PickRequired<Sphere3D, 'extrema'> {
        return sphere.extrema !== undefined;
    }

    export function create(center: Vec3, radius: number): Sphere3D { return { center, radius }; }
    export function zero(): Sphere3D { return { center: Vec3(), radius: 0 }; }

    export function clone(a: Sphere3D): Sphere3D {
        const out = create(Vec3.clone(a.center), a.radius);
        if (hasExtrema(a)) out.extrema = a.extrema;
        return out;
    }

    export function set(out: Sphere3D, center: Vec3, radius: number) {
        Vec3.copy(out.center, center);
        out.radius = radius;
        return out;
    }

    export function copy(out: Sphere3D, a: Sphere3D) {
        Vec3.copy(out.center, a.center);
        out.radius = a.radius;
        if (hasExtrema(a)) setExtrema(out, a.extrema);
        return out;
    }

    export function setExtrema(out: Sphere3D, extrema: Vec3[]): Sphere3D {
        if (out.extrema !== undefined) {
            out.extrema.length = 0;
            out.extrema.push(...extrema);
        } else {
            out.extrema = [...extrema];
        }
        return out;
    }

    export function computeBounding(data: PositionData): Sphere3D {
        const { x, y, z, indices } = data;
        let cx = 0, cy = 0, cz = 0;
        let radiusSq = 0;

        const size = OrderedSet.size(indices);
        for (let t = 0; t < size; t++) {
            const i = OrderedSet.getAt(indices, t);
            cx += x[i];
            cy += y[i];
            cz += z[i];
        }

        if (size > 0) {
            cx /= size;
            cy /= size;
            cz /= size;
        }

        for (let t = 0; t < size; t++) {
            const i = OrderedSet.getAt(indices, t);
            const dx = x[i] - cx, dy = y[i] - cy, dz = z[i] - cz;
            const d = dx * dx + dy * dy + dz * dz;
            if (d > radiusSq) radiusSq = d;
        }

        return { center: Vec3.create(cx, cy, cz), radius: Math.sqrt(radiusSq) };
    }

    /** Transform sphere with a Mat4 */
    export function transform(out: Sphere3D, sphere: Sphere3D, m: Mat4): Sphere3D {
        Vec3.transformMat4(out.center, sphere.center, m);
        out.radius = sphere.radius * Mat4.getMaxScaleOnAxis(m);
        return out;
    }

    export function toArray(s: Sphere3D, out: NumberArray, offset: number) {
        Vec3.toArray(s.center, out, offset);
        out[offset + 3] = s.radius;
    }

    export function fromArray(out: Sphere3D, array: NumberArray, offset: number) {
        Vec3.fromArray(out.center, array, offset);
        out.radius = array[offset + 3];
        return out;
    }

    export function fromBox3D(out: Sphere3D, box: Box3D) {
        Vec3.scale(out.center, Vec3.add(out.center, box.max, box.min), 0.5);
        out.radius = Vec3.distance(out.center, box.max);
        return out;
    }

    export function fromAxes3D(out: Sphere3D, axes: Axes3D) {
        Vec3.copy(out.center, axes.origin);
        out.radius = Math.max(Vec3.magnitude(axes.dirA), Vec3.magnitude(axes.dirB), Vec3.magnitude(axes.dirC));
        return out;
    }

    const tmpAddVec3 = Vec3();
    export function addVec3(out: Sphere3D, s: Sphere3D, v: Vec3) {
        const d = Vec3.distance(s.center, v);
        if (d < s.radius) return Sphere3D.copy(out, s);
        Vec3.sub(tmpAddVec3, s.center, v);
        Vec3.sub(tmpAddVec3, s.center, tmpAddVec3);
        Vec3.setMagnitude(tmpAddVec3, tmpAddVec3, s.radius);
        Vec3.scale(out.center, Vec3.add(tmpAddVec3, tmpAddVec3, v), 0.5);
        out.radius = Vec3.distance(out.center, v);
        return out;
    }

    /** Expand sphere radius by another sphere */
    export function expandBySphere(out: Sphere3D, sphere: Sphere3D, by: Sphere3D) {
        Vec3.copy(out.center, sphere.center);
        out.radius = Math.max(sphere.radius, Vec3.distance(sphere.center, by.center) + by.radius);
        if (hasExtrema(sphere) && hasExtrema(by)) {
            setExtrema(out, [
                ...sphere.extrema.map(e => Vec3.clone(e)),
                ...by.extrema.map(e => Vec3.clone(e))
            ]);
        }
        return out;
    }

    const tmpDir = Vec3();
    /** Expand sphere radius by delta */
    export function expand(out: Sphere3D, sphere: Sphere3D, delta: number): Sphere3D {
        Vec3.copy(out.center, sphere.center);
        out.radius = sphere.radius + delta;
        if (hasExtrema(sphere)) {
            setExtrema(out, sphere.extrema.map(e => {
                Vec3.sub(tmpDir, e, sphere.center);
                const dist = Vec3.distance(sphere.center, e);
                Vec3.normalize(tmpDir, tmpDir);
                return Vec3.scaleAndAdd(Vec3(), sphere.center, tmpDir, dist + delta);
            }));
        }
        return out;
    }

    /**
     * Returns whether or not the spheres have exactly the same center and radius (when compared with ===)
     */
    export function exactEquals(a: Sphere3D, b: Sphere3D) {
        return a.radius === b.radius && Vec3.exactEquals(a.center, b.center);
    }

    /**
     * Returns whether or not the spheres have approximately the same center and radius.
     */
    export function equals(a: Sphere3D, b: Sphere3D) {
        const ar = a.radius;
        const br = b.radius;
        return (Math.abs(ar - br) <= EPSILON * Math.max(1.0, Math.abs(ar), Math.abs(br)) &&
                Vec3.equals(a.center, b.center));
    }

    /**
     * Check if `a` includes `b`, use `extrema` of `b` when available
     */
    export function includes(a: Sphere3D, b: Sphere3D) {
        if (hasExtrema(b)) {
            for (const e of b.extrema) {
                if (Vec3.distance(a.center, e) > a.radius) return false;
            }
            return true;
        } else {
            return Vec3.distance(a.center, b.center) + b.radius <= a.radius;
        }
    }

    /** Check if `a` and `b` are overlapping */
    export function overlaps(a: Sphere3D, b: Sphere3D) {
        return Vec3.distance(a.center, b.center) <= a.radius + b.radius;
    }

    /** Get the signed distance of `a` and `b` */
    export function distance(a: Sphere3D, b: Sphere3D) {
        return Vec3.distance(a.center, b.center) - a.radius + b.radius;
    }
}

export { Sphere3D };