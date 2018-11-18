/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from 'mol-math/linear-algebra/3d';

export class CentroidHelper {
    private count = 0;

    center: Vec3 = Vec3.zero();
    radiusSq = 0;

    reset() {
        Vec3.set(this.center, 0, 0, 0);
        this.radiusSq = 0;
        this.count = 0;
    }

    includeStep(p: Vec3) {
        Vec3.add(this.center, this.center, p);
        this.count++;
    }

    finishedIncludeStep() {
        if (this.count === 0) return;
        Vec3.scale(this.center, this.center, 1 / this.count);
    }

    radiusStep(p: Vec3) {
        const d = Vec3.squaredDistance(p, this.center);
        if (d > this.radiusSq) this.radiusSq = d;
    }

    constructor() {

    }
}