/**
 * Copyright (c) 2018-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { Matrix } from './matrix';
import { Vec3 } from '../3d/vec3';
import { svd } from './svd';
import { NumberArray } from '../../../mol-util/type-helpers';
import { Axes3D } from '../../geometry/primitives/axes3d';
import { EPSILON } from '../3d/common';

export { PrincipalAxes };

interface PrincipalAxes {
    momentsAxes: Axes3D
    boxAxes: Axes3D
}

namespace PrincipalAxes {
    export function ofPositions(positions: NumberArray): PrincipalAxes {
        const momentsAxes = calculateMomentsAxes(positions);
        const boxAxes = calculateBoxAxes(positions, momentsAxes);
        return { momentsAxes, boxAxes };
    }

    export function calculateMomentsAxes(positions: NumberArray): Axes3D {
        if (positions.length === 3) {
            return Axes3D.create(Vec3.fromArray(Vec3(), positions, 0), Vec3.create(1, 0, 0), Vec3.create(0, 1, 0), Vec3.create(0, 0, 1));
        }

        const points = Matrix.fromArray(positions, 3, positions.length / 3);
        const n = points.rows;
        const n3 = n / 3;
        const A = Matrix.create(3, 3);
        const W = Matrix.create(1, 3);
        const U = Matrix.create(3, 3);
        const V = Matrix.create(3, 3);

        // calculate
        const mean = Matrix.meanRows(points);
        const pointsM = Matrix.subRows(Matrix.clone(points), mean);
        const pointsT = Matrix.transpose(Matrix.create(n, 3), pointsM);
        Matrix.multiplyABt(A, pointsT, pointsT);
        svd(A, W, U, V);

        // origin
        const origin = Vec3.create(mean[0], mean[1], mean[2]);

        // directions
        const dirA = Vec3.create(U.data[0], U.data[3], U.data[6]);
        const dirB = Vec3.create(U.data[1], U.data[4], U.data[7]);
        const dirC = Vec3.create(U.data[2], U.data[5], U.data[8]);
        Vec3.scale(dirA, dirA, Math.sqrt(W.data[0] / n3));
        Vec3.scale(dirB, dirB, Math.sqrt(W.data[1] / n3));
        Vec3.scale(dirC, dirC, Math.sqrt(W.data[2] / n3));

        return Axes3D.create(origin, dirA, dirB, dirC);
    }

    export function calculateNormalizedAxes(momentsAxes: Axes3D): Axes3D {
        const a = Axes3D.clone(momentsAxes);
        if (Vec3.magnitude(a.dirC) < EPSILON) {
            Vec3.cross(a.dirC, a.dirA, a.dirB);
        }
        return Axes3D.normalize(a, a);
    }

    const tmpBoxVec = Vec3();
    /**
     * Get the scale/length for each dimension for a box around the axes
     * to enclose the given positions
     */
    export function calculateBoxAxes(positions: NumberArray, momentsAxes: Axes3D): Axes3D {
        if (positions.length === 3) {
            return Axes3D.clone(momentsAxes);
        }

        let d1a = -Infinity;
        let d1b = -Infinity;
        let d2a = -Infinity;
        let d2b = -Infinity;
        let d3a = -Infinity;
        let d3b = -Infinity;

        const p = Vec3();
        const t = Vec3();

        const center = momentsAxes.origin;
        const a = calculateNormalizedAxes(momentsAxes);

        for (let i = 0, il = positions.length; i < il; i += 3) {
            Vec3.projectPointOnVector(p, Vec3.fromArray(p, positions, i), a.dirA, center);
            const dp1 = Vec3.dot(a.dirA, Vec3.normalize(t, Vec3.sub(t, p, center)));
            const dt1 = Vec3.distance(p, center);
            if (dp1 > 0) {
                if (dt1 > d1a) d1a = dt1;
            } else {
                if (dt1 > d1b) d1b = dt1;
            }

            Vec3.projectPointOnVector(p, Vec3.fromArray(p, positions, i), a.dirB, center);
            const dp2 = Vec3.dot(a.dirB, Vec3.normalize(t, Vec3.sub(t, p, center)));
            const dt2 = Vec3.distance(p, center);
            if (dp2 > 0) {
                if (dt2 > d2a) d2a = dt2;
            } else {
                if (dt2 > d2b) d2b = dt2;
            }

            Vec3.projectPointOnVector(p, Vec3.fromArray(p, positions, i), a.dirC, center);
            const dp3 = Vec3.dot(a.dirC, Vec3.normalize(t, Vec3.sub(t, p, center)));
            const dt3 = Vec3.distance(p, center);
            if (dp3 > 0) {
                if (dt3 > d3a) d3a = dt3;
            } else {
                if (dt3 > d3b) d3b = dt3;
            }
        }

        const dirA = Vec3.setMagnitude(Vec3(), a.dirA, (d1a + d1b) / 2);
        const dirB = Vec3.setMagnitude(Vec3(), a.dirB, (d2a + d2b) / 2);
        const dirC = Vec3.setMagnitude(Vec3(), a.dirC, (d3a + d3b) / 2);

        const okDirA = Vec3.isFinite(dirA);
        const okDirB = Vec3.isFinite(dirB);
        const okDirC = Vec3.isFinite(dirC);

        const origin = Vec3();
        const addCornerHelper = function (d1: number, d2: number, d3: number) {
            Vec3.copy(tmpBoxVec, center);

            if (okDirA) Vec3.scaleAndAdd(tmpBoxVec, tmpBoxVec, a.dirA, d1);
            if (okDirB) Vec3.scaleAndAdd(tmpBoxVec, tmpBoxVec, a.dirB, d2);
            if (okDirC) Vec3.scaleAndAdd(tmpBoxVec, tmpBoxVec, a.dirC, d3);
            Vec3.add(origin, origin, tmpBoxVec);
        };
        addCornerHelper(d1a, d2a, d3a);
        addCornerHelper(d1a, d2a, -d3b);
        addCornerHelper(d1a, -d2b, -d3b);
        addCornerHelper(d1a, -d2b, d3a);
        addCornerHelper(-d1b, -d2b, -d3b);
        addCornerHelper(-d1b, -d2b, d3a);
        addCornerHelper(-d1b, d2a, d3a);
        addCornerHelper(-d1b, d2a, -d3b);
        Vec3.scale(origin, origin, 1 / 8);

        return Axes3D.create(origin, dirA, dirB, dirC);
    }
}
