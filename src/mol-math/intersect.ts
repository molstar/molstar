/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec2, Vec3 } from './linear-algebra';

const tmpT = Vec2();
const tmpOrigin = Vec3();
const tmpDir = Vec3();

export function sphereIntersect(result: Vec3, p0: Vec3, p1: Vec3, center: Vec3, radius: number, invert: boolean) {
    Vec3.copy(tmpOrigin, p0);
    Vec3.sub(tmpDir, p1, p0);

    const maxT = Vec3.magnitude(tmpDir);
    Vec3.normalize(tmpDir, tmpDir);

    const ts = raySphere(tmpT, tmpOrigin, tmpDir, center, radius);
    if (ts === undefined || ts[1] < 0.0 || ts[0] > maxT) {
        return undefined;
    }

    const t = invert ? Math.max(ts[0], 0.0) : Math.min(ts[1], maxT);

    Vec3.scaleAndAdd(result, tmpOrigin, tmpDir, t);
    return result;
};

function solveQuadratic(result: Vec2, a: number, b: number, c: number) {
    const det = b * b - 4.0 * a * c;
    if (det < 0.0) {
        return undefined;
    } else if (det > 0.0) {
        const denom = 1.0 / (2.0 * a);
        const disc = Math.sqrt(det);
        const root0 = (-b + disc) * denom;
        const root1 = (-b - disc) * denom;

        if (root0 < root1) {
            result[0] = root0;
            result[1] = root1;
        } else {
            result[0] = root1;
            result[1] = root0;
        }

        return result;
    }

    const root = -b / (2.0 * a);
    if (root === 0.0) {
        return undefined;
    }

    result[0] = result[1] = root;
    return result;
}

const tmpRaySphereRoots = Vec2();
const tmpDiff = Vec3();

function raySphere(result: Vec2, origin: Vec3, direction: Vec3, center: Vec3, radius: number) {
    const radiusSquared = radius * radius;

    Vec3.sub(tmpDiff, origin, center);

    const a = Vec3.dot(direction, direction);
    const b = 2.0 * Vec3.dot(direction, tmpDiff);
    const c = Vec3.squaredMagnitude(tmpDiff) - radiusSquared;

    const roots = solveQuadratic(tmpRaySphereRoots, a, b, c);
    if (roots === undefined) {
        return undefined;
    }

    Vec2.copy(result, roots);
    return result;
}
