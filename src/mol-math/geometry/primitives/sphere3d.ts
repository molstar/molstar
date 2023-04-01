/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
import { PrincipalAxes } from '../../linear-algebra/matrix/principal-axes';

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
        if (hasExtrema(a)) out.extrema = a.extrema.map(e => Vec3.clone(e));
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
        if (hasExtrema(a)) setExtrema(out, a.extrema.map(e => Vec3.clone(e)));
        return out;
    }

    /** Note that `extrema` must not be reused elsewhere */
    export function setExtrema(out: Sphere3D, extrema: Vec3[]): Sphere3D {
        if (out.extrema !== undefined) {
            out.extrema.length = 0;
            out.extrema.push(...extrema);
        } else {
            out.extrema = extrema;
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
        if (hasExtrema(sphere)) {
            setExtrema(out, sphere.extrema.map(e => Vec3.transformMat4(Vec3(), e, m)));
        }
        return out;
    }

    /** Translate sphere by Vec3 */
    export function translate(out: Sphere3D, sphere: Sphere3D, v: Vec3) {
        Vec3.add(out.center, sphere.center, v);
        if (hasExtrema(sphere)) {
            setExtrema(out, sphere.extrema.map(e => Vec3.add(Vec3(), e, v)));
        }
        return out;
    }

    export function toArray<T extends NumberArray>(s: Sphere3D, out: T, offset: number) {
        Vec3.toArray(s.center, out, offset);
        out[offset + 3] = s.radius;
        return out;
    }

    export function fromArray(out: Sphere3D, array: NumberArray, offset: number) {
        Vec3.fromArray(out.center, array, offset);
        out.radius = array[offset + 3];
        return out;
    }

    export function fromBox3D(out: Sphere3D, box: Box3D) {
        Vec3.scale(out.center, Vec3.add(out.center, box.max, box.min), 0.5);
        out.radius = Vec3.distance(out.center, box.max);

        Sphere3D.setExtrema(out, [
            Vec3.create(box.min[0], box.min[1], box.min[2]),
            Vec3.create(box.max[0], box.max[1], box.max[2]),
            Vec3.create(box.max[0], box.min[1], box.min[2]),
            Vec3.create(box.min[0], box.max[1], box.max[2]),
            Vec3.create(box.min[0], box.min[1], box.max[2]),
            Vec3.create(box.max[0], box.min[1], box.max[2]),
            Vec3.create(box.max[0], box.max[1], box.min[2]),
            Vec3.create(box.min[0], box.max[1], box.min[2]),
        ]);

        return out;
    }

    export function fromAxes3D(out: Sphere3D, axes: Axes3D) {
        Vec3.copy(out.center, axes.origin);
        out.radius = Math.max(Vec3.magnitude(axes.dirA), Vec3.magnitude(axes.dirB), Vec3.magnitude(axes.dirC));
        return out;
    }

    const tmpCenter = Vec3();
    /** Get a tight sphere around a transformed box */
    export function fromDimensionsAndTransform(out: Sphere3D, dimensions: Vec3, transform: Mat4) {
        const [x, y, z] = dimensions;

        const cpA = Vec3.create(0, 0, 0); Vec3.transformMat4(cpA, cpA, transform);
        const cpB = Vec3.create(x, y, z); Vec3.transformMat4(cpB, cpB, transform);
        const cpC = Vec3.create(x, 0, 0); Vec3.transformMat4(cpC, cpC, transform);
        const cpD = Vec3.create(0, y, z); Vec3.transformMat4(cpD, cpD, transform);

        const cpE = Vec3.create(0, 0, z); Vec3.transformMat4(cpE, cpE, transform);
        const cpF = Vec3.create(x, 0, z); Vec3.transformMat4(cpF, cpF, transform);
        const cpG = Vec3.create(x, y, 0); Vec3.transformMat4(cpG, cpG, transform);
        const cpH = Vec3.create(0, y, 0); Vec3.transformMat4(cpH, cpH, transform);

        Vec3.add(tmpCenter, cpA, cpB);
        Vec3.scale(tmpCenter, tmpCenter, 0.5);
        const d = Math.max(Vec3.distance(cpA, cpB), Vec3.distance(cpC, cpD));
        Sphere3D.set(out, tmpCenter, d / 2);
        Sphere3D.setExtrema(out, [cpA, cpB, cpC, cpD, cpE, cpF, cpG, cpH]);

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
        if (sphere.radius < 1e-12 || (sphere.extrema?.length ?? 0) <= 1) {
            out.extrema = void 0;
            return out;
        }
        if (hasExtrema(sphere)) {
            const positions = new Float32Array(sphere.extrema.length * 3);
            for (let i = 0; i < sphere.extrema.length; i++) {
                Vec3.toArray(sphere.extrema[i], positions, i * 3);
            }

            const axes = PrincipalAxes.calculateMomentsAxes(positions);
            Axes3D.scale(axes, Axes3D.normalize(axes, axes), delta);

            setExtrema(out, sphere.extrema.map(e => {
                Vec3.normalize(tmpDir, Vec3.sub(tmpDir, e, sphere.center));
                const o = Vec3.clone(e);

                const sA = Vec3.dot(tmpDir, axes.dirA) < 0 ? -1 : 1;
                Vec3.scaleAndAdd(o, o, axes.dirA, sA);

                const sB = Vec3.dot(tmpDir, axes.dirB) < 0 ? -1 : 1;
                Vec3.scaleAndAdd(o, o, axes.dirB, sB);

                const sC = Vec3.dot(tmpDir, axes.dirC) < 0 ? -1 : 1;
                Vec3.scaleAndAdd(o, o, axes.dirC, sC);

                if (Vec3.distance(out.center, o) > out.radius) {
                    if (sphere.extrema.length >= 14) { // 14 extrema with coarse boundary helper
                        Vec3.normalize(tmpDir, Vec3.sub(tmpDir, o, sphere.center));
                    }
                    Vec3.scaleAndAdd(o, out.center, tmpDir, out.radius);
                }

                return o;
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

    /** Get the distance of v from sphere. If negative, v is inside sphere */
    export function distanceToVec(sphere: Sphere3D, v: Vec3): number {
        const { center, radius } = sphere;
        return Vec3.distance(v, center) - radius;
    }
}

export { Sphere3D };
