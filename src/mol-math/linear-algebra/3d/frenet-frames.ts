/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from './vec3';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3fromArray = Vec3.fromArray;
const v3sub = Vec3.sub;
const v3normalize = Vec3.normalize;
const v3isZero = Vec3.isZero;
const v3set = Vec3.set;
const v3cross = Vec3.cross;
const v3toArray = Vec3.toArray;
const v3copy = Vec3.copy;
const v3magnitude = Vec3.magnitude;
const v3rotateAroundAxis = Vec3.rotateAroundAxis;
const v3dot = Vec3.dot;

/**
 * Compute normal and binormal vectors along a curve using parallel transport
 */
export function computeFrenetFrames(curvePoints: Float32Array, normalVectors: Float32Array, binormalVectors: Float32Array, n: number) {
    const tangent = Vec3();
    const prevTangent = Vec3();
    const normal = Vec3();
    const binormal = Vec3();
    const p0 = Vec3();
    const p1 = Vec3();

    // Compute initial tangent
    v3fromArray(p0, curvePoints, 0);
    v3fromArray(p1, curvePoints, 3);
    v3sub(tangent, p1, p0);
    v3normalize(tangent, tangent);
    if (v3isZero(tangent)) {
        v3set(tangent, 1, 0, 0);
    }

    // Find initial normal (perpendicular to tangent)
    // Use the smallest component of tangent to find a perpendicular vector
    const absX = Math.abs(tangent[0]);
    const absY = Math.abs(tangent[1]);
    const absZ = Math.abs(tangent[2]);

    if (absX <= absY && absX <= absZ) {
        v3set(normal, 1, 0, 0);
    } else if (absY <= absZ) {
        v3set(normal, 0, 1, 0);
    } else {
        v3set(normal, 0, 0, 1);
    }

    // Orthogonalize normal against tangent
    v3cross(binormal, tangent, normal);
    v3normalize(binormal, binormal);
    v3cross(normal, binormal, tangent);
    v3normalize(normal, normal);

    // Store first frame
    v3toArray(normal, normalVectors, 0);
    v3toArray(binormal, binormalVectors, 0);

    // Propagate frames along the curve using parallel transport
    v3copy(prevTangent, tangent);
    for (let i = 1; i < n; ++i) {
        // Compute tangent at this point
        if (i < n - 1) {
            v3fromArray(p0, curvePoints, (i - 1) * 3);
            v3fromArray(p1, curvePoints, (i + 1) * 3);
            v3sub(tangent, p1, p0);
        } else {
            v3fromArray(p0, curvePoints, (i - 1) * 3);
            v3fromArray(p1, curvePoints, i * 3);
            v3sub(tangent, p1, p0);
        }
        v3normalize(tangent, tangent);

        // Parallel transport: rotate the previous frame
        const dot = v3dot(prevTangent, tangent);
        if (dot < 0.9999) {
            const axis = Vec3();
            v3cross(axis, prevTangent, tangent);
            if (v3magnitude(axis) > 0.0001) {
                v3normalize(axis, axis);
                const angle = Math.acos(Math.min(1, Math.max(-1, dot)));

                // Rotate normal and binormal around axis by angle
                v3rotateAroundAxis(normal, normal, axis, angle);
                v3rotateAroundAxis(binormal, binormal, axis, angle);
            }
        }

        // Ensure orthogonality
        v3cross(binormal, tangent, normal);
        v3normalize(binormal, binormal);
        v3cross(normal, binormal, tangent);
        v3normalize(normal, normal);

        v3toArray(normal, normalVectors, i * 3);
        v3toArray(binormal, binormalVectors, i * 3);
        v3copy(prevTangent, tangent);
    }
}
