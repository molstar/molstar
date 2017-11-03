/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Model } from '../model'
import Operator from './operator'
import { Vec3, Mat4 } from 'mol-base/math/linear-algebra-3d'

interface Unit extends Readonly<{
    // Structure-level unique identifier of the unit.
    id: number,

    // Provides access to the underlying data.
    model: Model,

    // Determines the operation applied to this unit.
    // The transform and and inverse are baked into the "getPosition" function
    operator: Operator,

    // Reference some commonly accessed things for faster access.
    residueIndex: ArrayLike<number>,
    chainIndex: ArrayLike<number>,
    hierarchy: Model['hierarchy'],
    conformation: Model['conformation']
}> {
    // returns the untransformed position. Used for spatial queries.
    getInvariantPosition(atom: number, slot: Vec3): Vec3,

    // gets the transformed position of the specified atom.
    getPosition(atom: number, slot: Vec3): Vec3,

    // optimized x/y/z coordinate projections for your query needs.
    x(atom: number): number,
    y(atom: number): number,
    z(atom: number): number
}

namespace Unit {
    export function create(model: Model, operator: Operator): Unit {
        const h = model.hierarchy;
        const { __x: x, __y: y, __z: z } = model.conformation;
        const getInvariantPosition = makeGetInvariantPosition(x, y, z);
        const getPosition = operator.isIdentity ? getInvariantPosition : makeGetPosition(operator.matrix, x, y, z);
        return {
            id: nextUnitId(),
            model,
            operator,
            residueIndex: h.residueSegments.segmentMap,
            chainIndex: h.chainSegments.segmentMap,
            hierarchy: model.hierarchy,
            conformation: model.conformation,
            getInvariantPosition,
            getPosition,
            x: operator.isIdentity ? makeGetCoord(x) : makeX(operator.matrix, x, y, z),
            y: operator.isIdentity ? makeGetCoord(y) : makeY(operator.matrix, x, y, z),
            z: operator.isIdentity ? makeGetCoord(z) : makeZ(operator.matrix, x, y, z)
        };
    }
}

export default Unit;

let _id = 0;
function nextUnitId() {
    const ret = _id;
    _id = (_id + 1) % 0x3fffffff;
    return ret;
}

function makeGetInvariantPosition(x: ArrayLike<number>, y: ArrayLike<number>, z: ArrayLike<number>) {
    return (a: number, r: Vec3): Vec3 => {
        r[0] = x[a];
        r[1] = y[a];
        r[2] = z[a];
        return r;
    }
}

function makeGetCoord(xs: ArrayLike<number>) {
    return (i: number) => xs[i];
}

function isWConst(m: Mat4) {
    return m[3] === 0 && m[7] === 0 && m[11] === 0;
}

function makeX(m: Mat4, xs: ArrayLike<number>, ys: ArrayLike<number>, zs: ArrayLike<number>) {
    const xx = m[0], yy = m[4], zz = m[8], ww = m[12];

    if (isWConst(m)) {
        const w = m[15] || 1.0;
        return (i: number) => (xx * xs[i] + yy * ys[i] + zz * zs[i] + ww) / w;
    }

    return (i: number) => {
        const x = xs[i], y = ys[i], z = zs[i], w = (m[3] * x + m[7] * y + m[11] * z + m[15]) || 1.0;
        return (xx * x + yy * y + zz * z + ww) / w;
    }
}

function makeY(m: Mat4, xs: ArrayLike<number>, ys: ArrayLike<number>, zs: ArrayLike<number>) {
    const xx = m[1], yy = m[5], zz = m[9], ww = m[13];

    if (isWConst(m)) {
        const w = m[15] || 1.0;
        return (i: number) => (xx * xs[i] + yy * ys[i] + zz * zs[i] + ww) / w;
    }

    return (i: number) => {
        const x = xs[i], y = ys[i], z = zs[i], w = (m[3] * x + m[7] * y + m[11] * z + m[15]) || 1.0;
        return (xx * x + yy * y + zz * z + ww) / w;
    }
}

function makeZ(m: Mat4, xs: ArrayLike<number>, ys: ArrayLike<number>, zs: ArrayLike<number>) {
    const xx = m[2], yy = m[6], zz = m[10], ww = m[14];

    if (isWConst(m)) {
        const w = m[15] || 1.0;
        return (i: number) => (xx * xs[i] + yy * ys[i] + zz * zs[i] + ww) / w;
    }

    return (i: number) => {
        const x = xs[i], y = ys[i], z = zs[i], w = (m[3] * x + m[7] * y + m[11] * z + m[15]) || 1.0;
        return (xx * x + yy * y + zz * z + ww) / w;
    }
}

function makeGetPosition(m: Mat4, x: ArrayLike<number>, y: ArrayLike<number>, z: ArrayLike<number>) {
    return (a: number, r: Vec3): Vec3 => {
        r[0] = x[a];
        r[1] = y[a];
        r[2] = z[a];
        Vec3.transformMat4(r, r, m);
        return r;
    }
}