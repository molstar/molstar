/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '../linear-algebra/3d';
import { CentroidHelper } from './centroid-helper';

// implementing http://www.ep.liu.se/ecp/034/009/ecp083409.pdf

export class EposHelper {
    private minDist: number[] = []
    private maxDist: number[] = []
    private extrema: Vec3[] = []
    private centroidHelper = new CentroidHelper()

    private computeExtrema(i: number, p: Vec3) {
        const d = Vec3.dot(this.dir[i], p)

        if (d < this.minDist[i]) {
            this.minDist[i] = d
            Vec3.copy(this.extrema[i * 2], p)
        }
        if (d > this.maxDist[i]) {
            this.maxDist[i] = d
            Vec3.copy(this.extrema[i * 2 + 1], p)
        }
    }

    includeStep(p: Vec3) {
        for (let i = 0, il = this.dir.length; i < il; ++i) {
            this.computeExtrema(i, p)
        }
    }

    finishedIncludeStep() {
        for (let i = 0; i < this.extrema.length; i++) {
            this.centroidHelper.includeStep(this.extrema[i]);
        }
        this.centroidHelper.finishedIncludeStep();
    }

    radiusStep(p: Vec3) {
        this.centroidHelper.radiusStep(p);
    }

    getSphere() {
        return this.centroidHelper.getSphere();
    }

    reset() {
        for (let i = 0, il = this.dir.length; i < il; ++i) {
            this.minDist[i] = Infinity
            this.maxDist[i] = -Infinity
            this.extrema[i * 2] = Vec3()
            this.extrema[i * 2 + 1] = Vec3()
        }
        this.centroidHelper.reset()
    }

    constructor(private dir: Vec3[]) {
        this.reset()
    }
}

const DirEpos6 = [
    Vec3.create(1, 0, 0), Vec3.create(0, 1, 0), Vec3.create(0, 0, 1)
]
export function Epos6() {
    return new EposHelper(DirEpos6)
}

const DirEpos14 = [
    Vec3.create(1, 0, 0), Vec3.create(0, 1, 0), Vec3.create(0, 0, 1),
    Vec3.create(1/3, 1/3, 1/3), Vec3.create(-1/3, 1/3, 1/3),
    Vec3.create(-1/3, -1/3, 1/3), Vec3.create(1/3, -1/3, 1/3)
]
export function Epos14() {
    return new EposHelper(DirEpos14)
}