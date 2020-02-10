/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Mat4 from './mat4';
import Vec3 from './vec3';
import { EVD } from '../matrix/evd';
import { CentroidHelper } from '../../../mol-math/geometry/centroid-helper';
import Matrix from '../matrix/matrix';
import { Sphere3D } from '../../geometry/primitives/sphere3d';

export { MinimizeRmsd };
namespace MinimizeRmsd {
    export interface Result {
        bTransform: Mat4,
        rmsd: number
    }

    export interface Positions { x: ArrayLike<number>, y: ArrayLike<number>, z: ArrayLike<number> }
    export namespace Positions {
        export function empty(n: number) {
            return { x: new Float64Array(n), y: new Float64Array(n), z: new Float64Array(n) };
        }
    }

    export interface Input {
        a: Positions,
        b: Positions,
        centerA?: Vec3,
        centerB?: Vec3
    }

    export function compute(data: Input, result?: MinimizeRmsd.Result) {
        if (typeof result === 'undefined') result = { bTransform: Mat4.zero(), rmsd: 0.0 };
        findMinimalRmsdTransformImpl(new RmsdTransformState(data, result));
        return result;
    }
}

class RmsdTransformState {
    a: MinimizeRmsd.Positions;
    b: MinimizeRmsd.Positions;

    centerA: Vec3;
    centerB: Vec3;

    evdCache: EVD.Cache = EVD.createCache(4);

    translateB = Mat4.identity();
    rotateB = Mat4.identity();
    tempMatrix = Mat4.identity();

    result: MinimizeRmsd.Result;

    constructor(data: MinimizeRmsd.Input, into: MinimizeRmsd.Result) {
        this.a = data.a;
        this.b = data.b;

        if (data.centerA) this.centerA = data.centerA;
        else this.centerA = data.centerA = CentroidHelper.fromArrays(data.a, Sphere3D()).center;

        if (data.centerB) this.centerB = data.centerB;
        else this.centerB = data.centerB = CentroidHelper.fromArrays(data.b, Sphere3D()).center;

        this.result = into;
    }
}

function computeN(state: RmsdTransformState) {
    const N = state.evdCache.matrix;

    Matrix.makeZero(N);

    const xsA = state.a.x, ysA = state.a.y, zsA = state.a.z;
    const xsB = state.b.x, ysB = state.b.y, zsB = state.b.z;
    const cA = state.centerA;
    const cB = state.centerB;

    let sizeSq = 0.0;

    for (let i = 0, _l = state.a.x.length; i < _l; i++) {
        const aX = xsA[i] - cA[0], aY = ysA[i] - cA[1], aZ = zsA[i] - cA[2];
        const bX = xsB[i] - cB[0], bY = ysB[i] - cB[1], bZ = zsB[i] - cB[2];

        sizeSq += aX * aX + aY * aY + aZ * aZ + bX * bX + bY * bY + bZ * bZ;

        Matrix.add(N, 0, 0, aX * bX + aY * bY + aZ * bZ);
        Matrix.add(N, 0, 1, -(aZ * bY) + aY * bZ);
        Matrix.add(N, 0, 2, aZ * bX - aX * bZ);
        Matrix.add(N, 0, 3, -(aY * bX) + aX * bY);
        Matrix.add(N, 1, 0, -(aZ * bY) + aY * bZ);
        Matrix.add(N, 1, 1, aX * bX - aY * bY - aZ * bZ);
        Matrix.add(N, 1, 2, aY * bX + aX * bY);
        Matrix.add(N, 1, 3, aZ * bX + aX * bZ);
        Matrix.add(N, 2, 0, aZ * bX - aX * bZ);
        Matrix.add(N, 2, 1, aY * bX + aX * bY);
        Matrix.add(N, 2, 2, -(aX * bX) + aY * bY - aZ * bZ);
        Matrix.add(N, 2, 3, aZ * bY + aY * bZ);
        Matrix.add(N, 3, 0, -(aY * bX) + aX * bY);
        Matrix.add(N, 3, 1, aZ * bX + aX * bZ);
        Matrix.add(N, 3, 2, aZ * bY + aY * bZ);
        Matrix.add(N, 3, 3, -(aX * bX) - aY * bY + aZ * bZ);

        // conjugate instead of transpose.
        // var l = new Quaternion(-a.X, -a.Y, -a.Z, 0).RightMultiplicationToMatrix();
        // l.Transpose();
        // var r = new Quaternion(b.X, b.Y, b.Z, 0).LeftMultiplicationToMatrix();
        // N += l * r;
    }

    return sizeSq;
}

function makeTransformMatrix(state: RmsdTransformState) {
    const ev = state.evdCache.matrix;

    const qX = Matrix.get(ev, 1, 3);
    const qY = Matrix.get(ev, 2, 3);
    const qZ = Matrix.get(ev, 3, 3);
    const qW = Matrix.get(ev, 0, 3);

    const n1 = 2 * qY * qY;
    const n2 = 2 * qZ * qZ;
    const n3 = 2 * qX * qX;
    const n4 = 2 * qX * qY;
    const n5 = 2 * qW * qZ;
    const n6 = 2 * qX * qZ;
    const n7 = 2 * qW * qY;
    const n8 = 2 * qY * qZ;
    const n9 = 2 * qW * qX;

    let m = state.translateB;
    // translation to center
    Mat4.setValue(m, 0, 3, -state.centerB[0]);
    Mat4.setValue(m, 1, 3, -state.centerB[1]);
    Mat4.setValue(m, 2, 3, -state.centerB[2]);

    m = state.rotateB;
    // rotation
    Mat4.setValue(m, 0, 0, 1 - n1 - n2);
    Mat4.setValue(m, 0, 1, n4 + n5);
    Mat4.setValue(m, 0, 2, n6 - n7);
    Mat4.setValue(m, 1, 0, n4 - n5);
    Mat4.setValue(m, 1, 1, 1 - n3 - n2);
    Mat4.setValue(m, 1, 2, n8 + n9);
    Mat4.setValue(m, 2, 0, n6 + n7);
    Mat4.setValue(m, 2, 1, n8 - n9);
    Mat4.setValue(m, 2, 2, 1 - n3 - n1);
    Mat4.setValue(m, 3, 3, 1);

    Mat4.mul(state.tempMatrix, state.rotateB, state.translateB);

    m = state.translateB;
    // translation to center
    Mat4.setValue(m, 0, 3, state.centerA[0]);
    Mat4.setValue(m, 1, 3, state.centerA[1]);
    Mat4.setValue(m, 2, 3, state.centerA[2]);

    Mat4.mul(state.result.bTransform, state.translateB, state.tempMatrix);
}

function findMinimalRmsdTransformImpl(state: RmsdTransformState): void {
    const sizeSq = computeN(state);

    EVD.compute(state.evdCache);
    let rmsd = sizeSq - 2.0 * state.evdCache.eigenValues[3];
    rmsd = rmsd < 0.0 ? 0.0 : Math.sqrt(rmsd / state.a.x.length);
    makeTransformMatrix(state);
    state.result.rmsd = rmsd;
}