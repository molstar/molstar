/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from 'mol-math/linear-algebra/3d';
import { Box3D } from './primitives/box3d';
import { Sphere3D } from './primitives/sphere3d';

/**
 * Usage:
 *
 * 1. .reset(tolerance); tolerance plays part in the "extend" step
 * 2. for each point/sphere call boundaryStep()
 * 3. .finishBoundaryStep
 * 4. for each point/sphere call extendStep
 * 5. use .center/.radius or call getSphere/getBox
 */
export class BoundaryHelper {
    private count = 0;
    private extremes = [Vec3.zero(), Vec3.zero(), Vec3.zero(), Vec3.zero(), Vec3.zero(), Vec3.zero()];
    private u = Vec3.zero();
    private v = Vec3.zero();

    tolerance = 0;
    center: Vec3 = Vec3.zero();
    radius = 0;

    reset(tolerance: number) {
        Vec3.set(this.center, 0, 0, 0);
        for (let i = 0; i < 6; i++) {
            const e = i % 2 === 0 ? Number.MAX_VALUE : -Number.MAX_VALUE;
            Vec3.set(this.extremes[i], e, e, e);
        }
        this.radius = 0;
        this.count = 0;
        this.tolerance = tolerance;
    }

    boundaryStep(p: Vec3, r: number) {
        updateExtremeMin(0, this.extremes[0], p, r);
        updateExtremeMax(0, this.extremes[1], p, r);

        updateExtremeMin(1, this.extremes[2], p, r);
        updateExtremeMax(1, this.extremes[3], p, r);

        updateExtremeMin(2, this.extremes[4], p, r);
        updateExtremeMax(2, this.extremes[5], p, r);
        this.count++;
    }

    finishBoundaryStep() {
        if (this.count === 0) return;

        let maxSpan = 0, mI = 0, mJ = 0;

        for (let i = 0; i < 5; i++) {
            for (let j = i + 1; j < 6; j++) {
                const d = Vec3.squaredDistance(this.extremes[i], this.extremes[j]);
                if (d > maxSpan) {
                    maxSpan = d;
                    mI = i;
                    mJ = j;
                }
            }
        }

        Vec3.add(this.center, this.extremes[mI], this.extremes[mJ]);
        Vec3.scale(this.center, this.center, 0.5);
        this.radius = Vec3.distance(this.center, this.extremes[mI]);
    }

    extendStep(p: Vec3, r: number) {
        const d = Vec3.distance(p, this.center);
        if ((1 + this.tolerance) * this.radius >= r + d) return;

        Vec3.sub(this.u, p, this.center);
        Vec3.normalize(this.u, this.u);

        Vec3.scale(this.v, this.u, -this.radius);
        Vec3.add(this.v, this.v, this.center);
        Vec3.scale(this.u, this.u, r + d);
        Vec3.add(this.u, this.u, this.center);

        Vec3.add(this.center, this.u, this.v);
        Vec3.scale(this.center, this.center, 0.5);
        this.radius = 0.5 * (r + d + this.radius);
    }

    getBox(): Box3D {
        Vec3.copy(this.u, this.extremes[0]);
        Vec3.copy(this.v, this.extremes[0]);

        for (let i = 1; i < 6; i++) {
            Vec3.min(this.u, this.u, this.extremes[i]);
            Vec3.max(this.v, this.v, this.extremes[i]);
        }

        return { min: Vec3.clone(this.u), max: Vec3.clone(this.v) };
    }

    getSphere(): Sphere3D {
        return { center: Vec3.clone(this.center), radius: this.radius };
    }

    constructor() {
        this.reset(0);
    }
}

function updateExtremeMin(d: number, e: Vec3, center: Vec3, r: number) {
    if (center[d] - r < e[d]) {
        Vec3.copy(e, center);
        e[d] -= r;
    }
}

function updateExtremeMax(d: number, e: Vec3, center: Vec3, r: number) {
    if (center[d] + r > e[d]) {
        Vec3.copy(e, center);
        e[d] += r;
    }
}