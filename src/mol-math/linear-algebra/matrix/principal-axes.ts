/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Matrix from './matrix';
import { Vec3, Mat4 } from '../3d';
import { svd } from './svd';
import { NumberArray } from '../../../mol-util/type-helpers';

export { PrincipalAxes }

interface PrincipalAxes {
    begA: Vec3
    endA: Vec3
    begB: Vec3
    endB: Vec3
    begC: Vec3
    endC: Vec3

    center: Vec3

    vecA: Vec3
    vecB: Vec3
    vecC: Vec3

    normVecA: Vec3
    normVecB: Vec3
    normVecC: Vec3
}

namespace PrincipalAxes {
    /**
     * @param points 3xN matrix
     */
    export function ofPoints(points: Matrix<3, number>): PrincipalAxes {
        const n = points.rows
        const n3 = n / 3
        const pointsT = Matrix.create(n, 3)
        const A = Matrix.create(3, 3)
        const W = Matrix.create(1, 3)
        const U = Matrix.create(3, 3)
        const V = Matrix.create(3, 3)

        // calculate
        const mean = Matrix.meanRows(points)
        Matrix.subRows(points, mean)
        Matrix.transpose(pointsT, points)
        Matrix.multiplyABt(A, pointsT, pointsT)
        svd(A, W, U, V)

        // center
        const center = Vec3.create(mean[0], mean[1], mean[2])

        // normalized
        const normVecA = Vec3.create(U.data[0], U.data[3], U.data[6])
        const normVecB = Vec3.create(U.data[1], U.data[4], U.data[7])
        const normVecC = Vec3.create(U.data[2], U.data[5], U.data[8])

        // scaled
        const vecA = Vec3.scale(Vec3.zero(), normVecA, Math.sqrt(W.data[0] / n3))
        const vecB = Vec3.scale(Vec3.zero(), normVecB, Math.sqrt(W.data[1] / n3))
        const vecC = Vec3.scale(Vec3.zero(), normVecC, Math.sqrt(W.data[2] / n3))

        // points
        const begA = Vec3.sub(Vec3.clone(center), center, vecA)
        const endA = Vec3.add(Vec3.clone(center), center, vecA)
        const begB = Vec3.sub(Vec3.clone(center), center, vecB)
        const endB = Vec3.add(Vec3.clone(center), center, vecB)
        const begC = Vec3.sub(Vec3.clone(center), center, vecC)
        const endC = Vec3.add(Vec3.clone(center), center, vecC)

        return {
            begA, endA, begB, endB, begC, endC,
            center,
            vecA, vecB, vecC,
            normVecA, normVecB, normVecC
        }
    }

    /**
     * Set basis matrix for given axes
     */
    export function setBasisMatrix(out: Mat4, principalAxes: PrincipalAxes) {
        Mat4.setAxes(out, principalAxes.normVecB, principalAxes.normVecA, principalAxes.normVecC)
        if (Mat4.determinant(out) < 0) Mat4.scaleUniformly(out, out, -1)
        return out
    }

    /**
     * Get the scale/length for each dimension for a box around the axes
     * to enclose the given positions
     */
    export function getProjectedScaleForPositions(positions: NumberArray, principalAxes: PrincipalAxes) {
        let d1a = -Infinity
        let d1b = -Infinity
        let d2a = -Infinity
        let d2b = -Infinity
        let d3a = -Infinity
        let d3b = -Infinity

        const p = Vec3()
        const t = Vec3()

        const { center, normVecA, normVecB, normVecC } = principalAxes

        for (let i = 0, il = positions.length; i < il; i += 3) {
            Vec3.projectPointOnVector(p, Vec3.fromArray(p, positions, i), normVecA, center)
            const dp1 = Vec3.dot(normVecA, Vec3.normalize(t, Vec3.sub(t, p, center)))
            const dt1 = Vec3.distance(p, center)
            if (dp1 > 0) {
                if (dt1 > d1a) d1a = dt1
            } else {
                if (dt1 > d1b) d1b = dt1
            }

            Vec3.projectPointOnVector(p, Vec3.fromArray(p, positions, i), normVecB, center)
            const dp2 = Vec3.dot(normVecB, Vec3.normalize(t, Vec3.sub(t, p, center)))
            const dt2 = Vec3.distance(p, center)
            if (dp2 > 0) {
                if (dt2 > d2a) d2a = dt2
            } else {
                if (dt2 > d2b) d2b = dt2
            }

            Vec3.projectPointOnVector(p, Vec3.fromArray(p, positions, i), normVecC, center)
            const dp3 = Vec3.dot(normVecC, Vec3.normalize(t, Vec3.sub(t, p, center)))
            const dt3 = Vec3.distance(p, center)
            if (dp3 > 0) {
                if (dt3 > d3a) d3a = dt3
            } else {
                if (dt3 > d3b) d3b = dt3
            }
        }

        return {
            d1a: d1a,
            d2a: d2a,
            d3a: d3a,
            d1b: -d1b,
            d2b: -d2b,
            d3b: -d3b
        }
    }
}