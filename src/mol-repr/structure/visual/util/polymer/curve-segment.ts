/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { Vec3 } from '../../../../../mol-math/linear-algebra';
import { NumberArray } from '../../../../../mol-util/type-helpers';
import { lerp } from '../../../../../mol-math/interpolate';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3fromArray = Vec3.fromArray;
const v3toArray = Vec3.toArray;
const v3normalize = Vec3.normalize;
const v3sub = Vec3.sub;
const v3spline = Vec3.spline;
const v3copy = Vec3.copy;
const v3cross = Vec3.cross;
const v3orthogonalize = Vec3.orthogonalize;
const v3scale = Vec3.scale;
const v3scaleAndAdd = Vec3.scaleAndAdd;
const v3dot = Vec3.dot;
const v3isZero = Vec3.isZero;

export interface CurveSegmentState {
    curvePoints: NumberArray,
    tangentVectors: NumberArray,
    normalVectors: NumberArray,
    binormalVectors: NumberArray,
    widthValues: NumberArray,
    heightValues: NumberArray,
    linearSegments: number,
    prevSegmentNormal: Vec3
}

export interface CurveSegmentControls {
    first: boolean,
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
        linearSegments,
        prevSegmentNormal: Vec3()
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
const prevTangentVec = Vec3();
const stepVec = Vec3();
const reflectVec = Vec3();
const firstTangentVec = Vec3();
const lastTangentVec = Vec3();
const firstNormalVec = Vec3();
const lastNormalVec = Vec3();

const HalfPi = Math.PI / 2;

/**
 * Populate normalVectors and binormalVectors with a rotation minimizing frame,
 * propagated by the double reflection method (Wang et al. 2008), seeded from the
 * previous segment's end normal (or from firstDirection at a polymer start). The
 * residual twist to reach lastDirection, taken modulo 180° since the profiles are
 * 2-fold symmetric, is distributed evenly along the segment.
 */
export function interpolateNormals(state: CurveSegmentState, controls: CurveSegmentControls) {
    const { curvePoints, tangentVectors, normalVectors, binormalVectors, prevSegmentNormal } = state;
    const { d12: firstDirection, d23: lastDirection, first } = controls;

    const n = curvePoints.length / 3;
    const n1 = n - 1;

    v3fromArray(firstTangentVec, tangentVectors, 0);
    v3fromArray(lastTangentVec, tangentVectors, n1 * 3);

    if (first || v3isZero(prevSegmentNormal)) {
        v3orthogonalize(firstNormalVec, firstTangentVec, firstDirection);
    } else {
        v3orthogonalize(firstNormalVec, firstTangentVec, prevSegmentNormal);
    }
    v3orthogonalize(lastNormalVec, lastTangentVec, lastDirection);

    v3copy(normalVec, firstNormalVec);
    v3copy(prevTangentVec, firstTangentVec);
    v3toArray(normalVec, normalVectors, 0);

    for (let i = 1; i < n; ++i) {
        v3fromArray(tangentVec, tangentVectors, i * 3);
        v3fromArray(tmpNormal, curvePoints, (i - 1) * 3);
        v3fromArray(stepVec, curvePoints, i * 3);
        v3sub(stepVec, stepVec, tmpNormal);

        const c1 = v3dot(stepVec, stepVec);
        if (c1 > 1e-12) {
            const k1 = -2 / c1;
            v3scaleAndAdd(tmpNormal, normalVec, stepVec, k1 * v3dot(stepVec, normalVec));
            v3scaleAndAdd(reflectVec, prevTangentVec, stepVec, k1 * v3dot(stepVec, prevTangentVec));
            v3sub(reflectVec, tangentVec, reflectVec);
            const c2 = v3dot(reflectVec, reflectVec);
            if (c2 > 1e-12) {
                v3scaleAndAdd(tmpNormal, tmpNormal, reflectVec, (-2 / c2) * v3dot(reflectVec, tmpNormal));
            }
        } else {
            v3copy(tmpNormal, normalVec);
        }
        v3orthogonalize(normalVec, tangentVec, tmpNormal);
        v3toArray(normalVec, normalVectors, i * 3);
        v3copy(prevTangentVec, tangentVec);
    }

    v3cross(binormalVec, lastTangentVec, normalVec);
    let twist = Math.atan2(v3dot(binormalVec, lastNormalVec), v3dot(normalVec, lastNormalVec));
    if (twist > HalfPi) twist -= Math.PI;
    else if (twist < -HalfPi) twist += Math.PI;

    for (let i = 0; i < n; ++i) {
        v3fromArray(tangentVec, tangentVectors, i * 3);
        v3fromArray(normalVec, normalVectors, i * 3);
        if (i > 0 && twist !== 0) {
            const a = twist * (i / n1);
            v3normalize(binormalVec, v3cross(binormalVec, tangentVec, normalVec));
            v3scale(normalVec, normalVec, Math.cos(a));
            v3scaleAndAdd(normalVec, normalVec, binormalVec, Math.sin(a));
            v3toArray(normalVec, normalVectors, i * 3);
        }
        v3normalize(binormalVec, v3cross(binormalVec, tangentVec, normalVec));
        v3toArray(binormalVec, binormalVectors, i * 3);
    }

    v3fromArray(prevSegmentNormal, normalVectors, n1 * 3);
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