/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '../linear-algebra/3d';
import { CentroidHelper } from './centroid-helper';
import { Sphere3D } from '../geometry';

// implementing http://www.ep.liu.se/ecp/034/009/ecp083409.pdf

export class EposHelper {
    private dir: Vec3[]

    private minDist: number[] = []
    private maxDist: number[] = []
    private extrema: Vec3[] = []
    centroidHelper = new CentroidHelper()

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

    private computeSphereExtrema(i: number, center: Vec3, radius: number) {
        const d = Vec3.dot(this.dir[i], center)

        if (d - radius < this.minDist[i]) {
            this.minDist[i] = d - radius
            Vec3.scaleAndSub(this.extrema[i * 2], center, this.dir[i], radius)
        }
        if (d + radius > this.maxDist[i]) {
            this.maxDist[i] = d + radius
            Vec3.scaleAndAdd(this.extrema[i * 2 + 1], center, this.dir[i], radius)
        }
    }

    includeStep(p: Vec3) {
        for (let i = 0, il = this.dir.length; i < il; ++i) {
            this.computeExtrema(i, p)
        }
    }

    includeSphereStep(center: Vec3, radius: number) {
        for (let i = 0, il = this.dir.length; i < il; ++i) {
            this.computeSphereExtrema(i, center, radius)
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

    radiusSphereStep(center: Vec3, radius: number) {
        this.centroidHelper.radiusSphereStep(center, radius);
    }

    getSphere(sphere?: Sphere3D) {
        return this.centroidHelper.getSphere(sphere);
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

    constructor(dir: number[][]) {
        this.dir = dir.map(a => {
            const v = Vec3.create(a[0], a[1], a[2])
            return Vec3.normalize(v, v)
        })
        this.reset()
    }
}

const Type001 = [
    [1, 0, 0], [0, 1, 0], [0, 0, 1]
]

const Type111 = [
    [1, 1, 1], [-1, 1, 1], [-1, -1, 1], [1, -1, 1]
]

const Type011 = [
    [1, 1, 0], [1, -1, 0], [1, 0, 1], [1, 0, -1], [0, 1, 1], [0, 1, -1]
]

const Type012 = [
    [0, 1, 2], [0, 2, 1], [1, 0, 2], [2, 0, 1], [1, 2, 0], [2, 1, 0],
    [0, 1, -2], [0, 2, -1], [1, 0, -2], [2, 0, -1], [1, -2, 0], [2, -1, 0]
]

const Type112 = [
    [1, 1, 2], [2, 1, 1], [1, 2, 1], [1, -1, 2], [1, 1, -2], [1, -1, -2],
    [2, -1, 1], [2, 1, -1], [2, -1, -1], [1, -2, 1], [1, 2, -1], [1, -2, -1]
]

const Type122 = [
    [2, 2, 1], [1, 2, 2], [2, 1, 2], [2, -2, 1], [2, 2, -1], [2, -2, -1],
    [1, -2, 2], [1, 2, -2], [1, -2, -2], [2, -1, 2], [2, 1, -2], [2, -1, -2]
]

const DirEpos6 = [ ...Type001 ]
export function Epos6() {
    return new EposHelper(DirEpos6)
}

const DirEpos14 = [ ...Type001, ...Type111 ]
export function Epos14() {
    return new EposHelper(DirEpos14)
}

const DirEpos26 = [ ...Type001, ...Type111, ...Type011 ]
export function Epos26() {
    return new EposHelper(DirEpos26)
}

const DirEpos98 = [ ...Type001, ...Type111, ...Type011, ...Type012, ...Type112, ...Type122 ]
export function Epos98() {
    return new EposHelper(DirEpos98)
}