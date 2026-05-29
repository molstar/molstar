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
 * Estimate a visually open camera direction around a selection.
 *
 * Geometric intuition:
 *
 * The selection centroid is treated as the origin. Each nearby obstruction
 * point is converted into a unit direction on the sphere around the selection:
 *
 *     v_i = normalize(p_i - origin)
 *
 * We then build the directional second-moment matrix:
 *
 *     M = sum_i w_i v_i v_i^T
 *
 * For any candidate view direction `u`, the quadratic form
 *
 *     u^T M u
 *
 * expands to:
 *
 *     sum_i w_i (u · v_i)^2
 *
 * Since `u · v_i = cos(theta_i)`, this value is large when `u` is aligned
 * with many obstruction directions and small when `u` is mostly perpendicular
 * to them. Therefore, the eigenvector of `M` with the smallest eigenvalue is
 * the axis that is least aligned, in a least-squares sense, with the nearby
 * obstruction directions.
 *
 * This gives an unoriented axis: `u` and `-u` have the same score because the
 * dot products are squared. To choose the camera-facing side, we compute the
 * weighted mean obstruction direction:
 *
 *     m = sum_i w_i v_i
 *
 * and return the sign of the axis that points away from this mean direction.
 *
 * In short:
 *
 * - project nearby points onto a sphere around the selection;
 * - find the sparsest angular axis using the smallest eigenvector of their
 *   second-moment matrix;
 * - choose the side of that axis opposite the average obstruction direction.
 *
 * This is a fast, deterministic heuristic. It minimizes average squared
 * angular alignment with nearby points; it is not the exact largest-empty-cone
 * or maximum-clearance solution.
 *
 * The returned vector is a unit direction from the selection centroid toward
 * the camera.
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