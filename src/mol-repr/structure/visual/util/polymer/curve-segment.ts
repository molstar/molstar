/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from 'mol-math/linear-algebra';
import { NumberArray } from 'mol-util/type-helpers';
import { lerp } from 'mol-math/interpolate';

export interface CurveSegmentState {
    curvePoints: NumberArray,
    tangentVectors: NumberArray,
    normalVectors: NumberArray,
    binormalVectors: NumberArray,
    widthValues: NumberArray,
    heightValues: NumberArray,
    linearSegments: number
}

export interface CurveSegmentControls {
    p0: Vec3, p1: Vec3, p2: Vec3, p3: Vec3, p4: Vec3,
    d12: Vec3, d23: Vec3
}

export function createCurveSegmentState(linearSegments: number): CurveSegmentState {
    const n = linearSegments + 1
    const pn = n * 3
    return {
        curvePoints: new Float32Array(pn),
        tangentVectors: new Float32Array(pn),
        normalVectors: new Float32Array(pn),
        binormalVectors: new Float32Array(pn),
        widthValues: new Float32Array(n),
        heightValues: new Float32Array(n),
        linearSegments
    }
}

export function interpolateCurveSegment(state: CurveSegmentState, controls: CurveSegmentControls, tension: number, shift: number) {
    interpolatePointsAndTangents(state, controls, tension, shift)
    interpolateNormals(state, controls)
}

const tanA = Vec3.zero()
const tanB = Vec3.zero()
const curvePoint = Vec3.zero()

export function interpolatePointsAndTangents(state: CurveSegmentState, controls: CurveSegmentControls, tension: number, shift: number) {
    const { curvePoints, tangentVectors, linearSegments } = state
    const { p0, p1, p2, p3, p4 } = controls

    const shift1 = 1 - shift

    for (let j = 0; j <= linearSegments; ++j) {
        const t = j * 1.0 / linearSegments;
        if (t < shift1) {
            Vec3.spline(curvePoint, p0, p1, p2, p3, t + shift, tension)
            Vec3.spline(tanA, p0, p1, p2, p3, t + shift + 0.01, tension)
            Vec3.spline(tanB, p0, p1, p2, p3, t + shift - 0.01, tension)
        } else {
            Vec3.spline(curvePoint, p1, p2, p3, p4, t - shift1, tension)
            Vec3.spline(tanA, p1, p2, p3, p4, t - shift1 + 0.01, tension)
            Vec3.spline(tanB, p1, p2, p3, p4, t - shift1 - 0.01, tension)
        }
        Vec3.toArray(curvePoint, curvePoints, j * 3)
        Vec3.normalize(tangentVec, Vec3.sub(tangentVec, tanA, tanB))
        Vec3.toArray(tangentVec, tangentVectors, j * 3)
    }
}

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
export function interpolateNormals(state: CurveSegmentState, controls: CurveSegmentControls) {
    const { curvePoints, tangentVectors, normalVectors, binormalVectors } = state
    const { d12: firstDirection, d23: lastDirection } = controls

    const n = curvePoints.length / 3

    Vec3.fromArray(firstControlPoint, curvePoints, 0)
    Vec3.fromArray(lastControlPoint, curvePoints, (n - 1) * 3)
    Vec3.fromArray(firstTangentVec, tangentVectors, 0)
    Vec3.fromArray(lastTangentVec, tangentVectors,  (n - 1) * 3)

    Vec3.orthogonalize(firstNormalVec, firstTangentVec, Vec3.sub(tmpNormal, firstControlPoint, firstDirection))
    Vec3.orthogonalize(lastNormalVec, lastTangentVec, Vec3.sub(tmpNormal, lastControlPoint, lastDirection))

    if (Vec3.dot(firstNormalVec, lastNormalVec) < 0) {
        Vec3.scale(lastNormalVec, lastNormalVec, -1)
    }

    Vec3.copy(prevNormal, firstNormalVec)

    for (let i = 0; i < n; ++i) {
        const t = i === 0 ? 0 : 1 / (n - i)

        Vec3.fromArray(tangentVec, tangentVectors, i * 3)

        Vec3.orthogonalize(normalVec, tangentVec, Vec3.slerp(tmpNormal, prevNormal, lastNormalVec, t))
        Vec3.toArray(normalVec, normalVectors, i * 3)

        Vec3.copy(prevNormal, normalVec)

        Vec3.normalize(binormalVec, Vec3.cross(binormalVec, tangentVec, normalVec))
        Vec3.toArray(binormalVec, binormalVectors, i * 3)
    }
}

export function interpolateSizes(state: CurveSegmentState, w0: number, w1: number, w2: number, h0: number, h1: number, h2: number, shift: number) {
    const { widthValues, heightValues, linearSegments } = state

    const shift1 = 1 - shift

    for (let i = 0; i <= linearSegments; ++i) {
        const t = i * 1.0 / linearSegments;
        if (t < shift1) {
            widthValues[i] = lerp(w0, w1, t + shift)
            heightValues[i] = lerp(h0, h1, t + shift)
        } else {
            widthValues[i] = lerp(w1, w2, t - shift1)
            heightValues[i] = lerp(h1, h2, t - shift1)
        }
    }
}