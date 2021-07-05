/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '../../../../../mol-math/linear-algebra';
import { NumberArray } from '../../../../../mol-util/type-helpers';
import { lerp, smoothstep } from '../../../../../mol-math/interpolate';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3fromArray = Vec3.fromArray;
const v3toArray = Vec3.toArray;
const v3normalize = Vec3.normalize;
const v3sub = Vec3.sub;
const v3spline = Vec3.spline;
const v3slerp = Vec3.slerp;
const v3copy = Vec3.copy;
const v3cross = Vec3.cross;
const v3orthogonalize = Vec3.orthogonalize;
const v3matchDirection = Vec3.matchDirection;
const v3scale = Vec3.scale;
const v3add = Vec3.add;

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
    secStrucFirst: boolean, secStrucLast: boolean
    p0: Vec3, p1: Vec3, p2: Vec3, p3: Vec3, p4: Vec3,
    d12: Vec3, d23: Vec3
}

export function createCurveSegmentState(linearSegments: number): CurveSegmentState {
    const n = linearSegments + 1;
    const pn = n * 3;
    return {
        curvePoints: new Float32Array(pn),
        tangentVectors: new Float32Array(pn),
        normalVectors: new Float32Array(pn),
        binormalVectors: new Float32Array(pn),
        widthValues: new Float32Array(n),
        heightValues: new Float32Array(n),
        linearSegments
    };
}

export function interpolateCurveSegment(state: CurveSegmentState, controls: CurveSegmentControls, tension: number, shift: number) {
    interpolatePointsAndTangents(state, controls, tension, shift);
    interpolateNormals(state, controls);
}

const tanA = Vec3();
const tanB = Vec3();
const curvePoint = Vec3();

export function interpolatePointsAndTangents(state: CurveSegmentState, controls: CurveSegmentControls, tension: number, shift: number) {
    const { curvePoints, tangentVectors, linearSegments } = state;
    const { p0, p1, p2, p3, p4, secStrucFirst, secStrucLast } = controls;

    const shift1 = 1 - shift;

    const tensionBeg = secStrucFirst ? 0.5 : tension;
    const tensionEnd = secStrucLast ? 0.5 : tension;

    for (let j = 0; j <= linearSegments; ++j) {
        const t = j * 1.0 / linearSegments;
        if (t < shift1) {
            const te = lerp(tensionBeg, tension, t);
            v3spline(curvePoint, p0, p1, p2, p3, t + shift, te);
            v3spline(tanA, p0, p1, p2, p3, t + shift + 0.01, tensionBeg);
            v3spline(tanB, p0, p1, p2, p3, t + shift - 0.01, tensionBeg);
        } else {
            const te = lerp(tension, tensionEnd, t);
            v3spline(curvePoint, p1, p2, p3, p4, t - shift1, te);
            v3spline(tanA, p1, p2, p3, p4, t - shift1 + 0.01, te);
            v3spline(tanB, p1, p2, p3, p4, t - shift1 - 0.01, te);
        }
        v3toArray(curvePoint, curvePoints, j * 3);
        v3normalize(tangentVec, v3sub(tangentVec, tanA, tanB));
        v3toArray(tangentVec, tangentVectors, j * 3);
    }
}

const tmpNormal = Vec3();
const tangentVec = Vec3();
const normalVec = Vec3();
const binormalVec = Vec3();
const prevNormal = Vec3();
const nextNormal = Vec3();
const firstTangentVec = Vec3();
const lastTangentVec = Vec3();
const firstNormalVec = Vec3();
const lastNormalVec = Vec3();

/**
 * Populate normalVectors by interpolating from firstDirection to lastDirection with
 * resulting vector perpendicular to tangentVectors and binormalVectors
 */
export function interpolateNormals(state: CurveSegmentState, controls: CurveSegmentControls) {
    const { curvePoints, tangentVectors, normalVectors, binormalVectors } = state;
    const { d12: firstDirection, d23: lastDirection } = controls;

    const n = curvePoints.length / 3;

    v3fromArray(firstTangentVec, tangentVectors, 0);
    v3fromArray(lastTangentVec, tangentVectors, (n - 1) * 3);

    v3orthogonalize(firstNormalVec, firstTangentVec, firstDirection);
    v3orthogonalize(lastNormalVec, lastTangentVec, lastDirection);
    v3matchDirection(lastNormalVec, lastNormalVec, firstNormalVec);

    v3copy(prevNormal, firstNormalVec);

    const n1 = n - 1;
    for (let i = 0; i < n; ++i) {
        const j = smoothstep(0, n1, i) * n1;
        const t = i === 0 ? 0 : 1 / (n - j);

        v3fromArray(tangentVec, tangentVectors, i * 3);

        v3orthogonalize(normalVec, tangentVec, v3slerp(tmpNormal, prevNormal, lastNormalVec, t));
        v3toArray(normalVec, normalVectors, i * 3);

        v3copy(prevNormal, normalVec);

        v3normalize(binormalVec, v3cross(binormalVec, tangentVec, normalVec));
        v3toArray(binormalVec, binormalVectors, i * 3);
    }

    for (let i = 1; i < n1; ++i) {
        v3fromArray(prevNormal, normalVectors, (i - 1) * 3);
        v3fromArray(normalVec, normalVectors, i * 3);
        v3fromArray(nextNormal, normalVectors, (i + 1) * 3);

        v3scale(normalVec, v3add(normalVec, prevNormal, v3add(normalVec, nextNormal, normalVec)), 1 / 3);
        v3toArray(normalVec, normalVectors, i * 3);

        v3fromArray(tangentVec, tangentVectors, i * 3);
        v3normalize(binormalVec, v3cross(binormalVec, tangentVec, normalVec));
        v3toArray(binormalVec, binormalVectors, i * 3);
    }
}

export function interpolateSizes(state: CurveSegmentState, w0: number, w1: number, w2: number, h0: number, h1: number, h2: number, shift: number) {
    const { widthValues, heightValues, linearSegments } = state;

    const shift1 = 1 - shift;

    for (let i = 0; i <= linearSegments; ++i) {
        const t = i * 1.0 / linearSegments;
        if (t < shift1) {
            widthValues[i] = lerp(w0, w1, t + shift);
            heightValues[i] = lerp(h0, h1, t + shift);
        } else {
            widthValues[i] = lerp(w1, w2, t - shift1);
            heightValues[i] = lerp(h1, h2, t - shift1);
        }
    }
}