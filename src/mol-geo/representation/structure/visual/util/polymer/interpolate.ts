/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from 'mol-math/linear-algebra';

const tmpNormal = Vec3.zero()
const tangentVec = Vec3.zero()
const normalVec = Vec3.zero()
const binormalVec = Vec3.zero()
const prevNormal = Vec3.zero()
const firstControlPoint = Vec3.zero()
const lastControlPoint = Vec3.zero()
const firstTangentVec = Vec3.zero()
const lastTangentVec = Vec3.zero()
const firstNormalVec = Vec3.zero()
const lastNormalVec = Vec3.zero()

/**
 * Populate normalVectors by interpolating from firstDirection to lastDirection with
 * resulting vector perpendicular to tangentVectors and binormalVectors
 */
export function interpolateNormals(controlPoints: Helpers.NumberArray, tangentVectors: Helpers.NumberArray, normalVectors: Helpers.NumberArray, binormalVectors: Helpers.NumberArray, firstDirection: Vec3, lastDirection: Vec3) {
    const n = controlPoints.length / 3

    Vec3.fromArray(firstControlPoint, controlPoints, 0)
    Vec3.fromArray(lastControlPoint, controlPoints, (n - 1) * 3)
    Vec3.fromArray(firstTangentVec, tangentVectors, 0)
    Vec3.fromArray(lastTangentVec, tangentVectors,  (n - 1) * 3)

    Vec3.normalize(tmpNormal, Vec3.sub(tmpNormal, firstControlPoint, firstDirection))
    Vec3.orthogonalize(firstNormalVec, firstTangentVec, tmpNormal)

    Vec3.normalize(tmpNormal, Vec3.sub(tmpNormal, lastControlPoint, lastDirection))
    Vec3.orthogonalize(lastNormalVec, lastTangentVec, tmpNormal)

    if (Vec3.dot(firstNormalVec, lastNormalVec) < 0) {
        Vec3.scale(lastNormalVec, lastNormalVec, -1)
    }

    Vec3.copy(prevNormal, firstNormalVec)

    for (let i = 0; i < n; ++i) {
        const t = i === 0 ? 0 : 1 / (n - i)
        Vec3.normalize(tmpNormal, Vec3.slerp(tmpNormal, prevNormal, lastNormalVec, t))

        Vec3.fromArray(tangentVec, tangentVectors, i * 3)

        Vec3.orthogonalize(normalVec, tangentVec, tmpNormal)
        Vec3.toArray(normalVec, normalVectors, i * 3)

        Vec3.copy(prevNormal, normalVec)

        Vec3.normalize(binormalVec, Vec3.cross(binormalVec, tangentVec, normalVec))
        Vec3.toArray(binormalVec, binormalVectors, i * 3)
    }
}