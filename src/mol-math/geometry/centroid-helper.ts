/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '../../mol-math/linear-algebra/3d';
import { Sphere3D } from './primitives/sphere3d';

export { CentroidHelper };

class CentroidHelper {
    private count = 0;

    center: Vec3 = Vec3();
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

    radiusSphereStep(center: Vec3, radius: number) {
        const _d = Vec3.distance(center, this.center) + radius;
        const d = _d * _d;
        if (d > this.radiusSq) this.radiusSq = d;
    }

    getSphere(sphere?: Sphere3D): Sphere3D {
        if (!sphere) sphere = Sphere3D();
        Vec3.copy(sphere.center, this.center);
        sphere.radius = Math.sqrt(this.radiusSq);
        return sphere;
    }

    getCount() {
        return this.count;
    }

    constructor() { }
}

namespace CentroidHelper {
    const helper = new CentroidHelper();
    const posA = Vec3();
    const posB = Vec3();

    export function fromArrays({ x, y, z }: { x: ArrayLike<number>, y: ArrayLike<number>, z: ArrayLike<number> }, to: Sphere3D) {
        helper.reset();
        const n = x.length;
        for (let i = 0; i < n; i++) {
            Vec3.set(posA, x[i], y[i], z[i]);
            helper.includeStep(posA);
        }
        helper.finishedIncludeStep();
        for (let i = 0; i < n; i++) {
            Vec3.set(posA, x[i], y[i], z[i]);
            helper.radiusStep(posA);
        }
        Vec3.copy(to.center, helper.center);
        to.radius = Math.sqrt(helper.radiusSq);
        return to;
    }

    export function fromProvider(count: number, getter: (i: number, pos: Vec3) => void, to: Sphere3D) {
        helper.reset();
        for (let i = 0; i < count; i++) {
            getter(i, posA);
            helper.includeStep(posA);
        }
        helper.finishedIncludeStep();
        for (let i = 0; i < count; i++) {
            getter(i, posA);
            helper.radiusStep(posA);
        }
        Vec3.copy(to.center, helper.center);
        to.radius = Math.sqrt(helper.radiusSq);
        return to;
    }

    export function fromPairProvider(count: number, getter: (i: number, posA: Vec3, posB: Vec3) => void, to: Sphere3D) {
        helper.reset();
        for (let i = 0; i < count; i++) {
            getter(i, posA, posB);
            helper.includeStep(posA);
            helper.includeStep(posB);
        }
        helper.finishedIncludeStep();
        for (let i = 0; i < count; i++) {
            getter(i, posA, posB);
            helper.radiusStep(posA);
            helper.radiusStep(posB);
        }
        Vec3.copy(to.center, helper.center);
        to.radius = Math.sqrt(helper.radiusSq);
        return to;
    }
}