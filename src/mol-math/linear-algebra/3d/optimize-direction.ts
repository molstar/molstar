/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from './vec3';
import { EVD } from '../matrix/evd';
import { Matrix } from '../matrix/matrix';

export interface LeastObstructedDirectionOptions {
    /** Optional centroid/origin. If omitted, centroid is computed from the provided points. */
    origin?: Vec3,

    /** Optional Gaussian falloff distance. If omitted, all points have weight 1. */
    sigma?: number,

    /** Ignore points closer than this to the origin. */
    minDistance?: number,
}

function eachPosition(points: ReadonlyArray<Vec3> | { x: ArrayLike<number>, y: ArrayLike<number>, z: ArrayLike<number> }, callback: (x: number, y: number, z: number) => void) {
    if (Array.isArray(points)) {
        for (const p of points) {
            callback(p[0], p[1], p[2]);
        }
    } else {
        const { x, y, z } = points as { x: ArrayLike<number>, y: ArrayLike<number>, z: ArrayLike<number> };
        const n = Math.min(x.length, y.length, z.length);
        for (let i = 0; i < n; i++) {
            callback(x[i], y[i], z[i]);
        }
    }
}

/**
 * Returns a direction from the selection centroid toward the camera.
 *
 * Input points should usually be nearby, non-selected atom positions.
 *
 * The returned direction is a unit Vec3 or undefined if no valid direction could be computed.
 */
export function leastObstructedDirection(
    points: ReadonlyArray<Vec3> | { x: ArrayLike<number>, y: ArrayLike<number>, z: ArrayLike<number> },
    options: LeastObstructedDirectionOptions = {}
): Vec3 | undefined {
    const origin = options.origin;
    const minDistance = options.minDistance ?? 1e-6;
    const minDistanceSq = minDistance * minDistance;

    const sigma = options.sigma;
    const useWeights = sigma !== void 0 && sigma > 0;
    const twoSigmaSq = useWeights ? 2 * sigma * sigma : 1;

    // Directional second moment:
    // M = sum_i w_i v_i v_i^T
    const evd = EVD.createCache(3);
    const M = evd.matrix;
    Matrix.makeZero(M);

    // Weighted mean direction, used only to choose sign.
    const mean = Vec3.zero();

    let count = 0;
    let weightSum = 0;

    eachPosition(points, (x_, y_, z_) => {
        let x = x_, y = y_, z = z_;
        if (origin) {
            x -= origin[0];
            y -= origin[1];
            z -= origin[2];
        }

        const dSq = x * x + y * y + z * z;
        if (dSq <= minDistanceSq) return;

        const d = Math.sqrt(dSq);
        const invD = 1 / d;

        // Unit obstruction direction v.
        x *= invD;
        y *= invD;
        z *= invD;

        const w = useWeights ? Math.exp(-dSq / twoSigmaSq) : 1;

        // Accumulate symmetric matrix.
        //
        // M = [
        //   xx xy xz
        //   xy yy yz
        //   xz yz zz
        // ]
        Matrix.add(M, 0, 0, w * x * x);
        Matrix.add(M, 0, 1, w * x * y);
        Matrix.add(M, 0, 2, w * x * z);

        Matrix.add(M, 1, 0, w * y * x);
        Matrix.add(M, 1, 1, w * y * y);
        Matrix.add(M, 1, 2, w * y * z);

        Matrix.add(M, 2, 0, w * z * x);
        Matrix.add(M, 2, 1, w * z * y);
        Matrix.add(M, 2, 2, w * z * z);

        mean[0] += w * x;
        mean[1] += w * y;
        mean[2] += w * z;

        count++;
        weightSum += w;
    });

    if (count === 0 || weightSum <= 0) {
        return undefined;
    }

    EVD.compute(evd);

    // EVD sorts eigenvalues ascending, so column 0 is the smallest eigenvector.
    const dir = Vec3.create(
        Matrix.get(M, 0, 0),
        Matrix.get(M, 1, 0),
        Matrix.get(M, 2, 0)
    );

    if (Vec3.magnitude(dir) < 1e-6) {
        return undefined;
    }

    Vec3.normalize(dir, dir);

    // Pick the less-obstructed side of the axis:
    // choose the sign opposite the weighted mean obstruction direction.
    if (Vec3.dot(dir, mean) > 0) {
        Vec3.scale(dir, dir, -1);
    }

    return dir;
}