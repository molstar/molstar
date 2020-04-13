/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3, Mat4 } from '../../linear-algebra';
import { PositionData } from '../common';
import { OrderedSet } from '../../../mol-data/int';
import { Sphere3D } from './sphere3d';

interface Box3D { min: Vec3, max: Vec3 }

function Box3D() {
    return Box3D.empty();
}

namespace Box3D {
    export function create(min: Vec3, max: Vec3): Box3D { return { min, max }; }
    export function empty(): Box3D { return { min: Vec3(), max: Vec3() }; }

    export function copy(out: Box3D, a: Box3D): Box3D {
        Vec3.copy(out.min, a.min);
        Vec3.copy(out.max, a.max);
        return out;
    }

    export function clone(a: Box3D): Box3D {
        return copy(empty(), a);
    }

    export function fromSphere3D(out: Box3D, sphere: Sphere3D): Box3D {
        const r = Vec3.create(sphere.radius, sphere.radius, sphere.radius);
        Vec3.sub(out.min, sphere.center, r);
        Vec3.add(out.max, sphere.center, r);
        return out;
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

    export function scale(out: Box3D, box: Box3D, scale: number) {
        Vec3.scale(out.min, box.min, scale);
        Vec3.scale(out.max, box.max, scale);
        return out;
    }

    const tmpTransformV = Vec3();
    /** Transform box with a Mat4 */
    export function transform(out: Box3D, box: Box3D, m: Mat4): Box3D {
        const [ minX, minY, minZ ] = box.min;
        const [ maxX, maxY, maxZ ] = box.max;
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
}

export { Box3D };