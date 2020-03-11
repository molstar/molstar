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

    getHierarchyInput() {
        const sphere = this.centroidHelper.getSphere();
        const normal = Vec3()
        const t = sphere.radius * this.hierarchyThresholdFactor

        let maxDist = -Infinity
        let belowThreshold = false

        for (let i = 0; i < this.extrema.length; i += 2) {
            const halfDist = Vec3.distance(this.extrema[i], this.extrema[i + 1]) / 2
            if (halfDist > maxDist) {
                maxDist = halfDist
                Vec3.normalize(normal, Vec3.sub(normal, this.extrema[i], this.extrema[i + 1]))
            }
            if (halfDist < t) belowThreshold = true
        }

        return belowThreshold ? { sphere, normal } : false
    }

    getSphere(sphere?: Sphere3D) {
        return this.centroidHelper.getSphere(sphere)
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

    constructor(quality: EposQuality, private hierarchyThresholdFactor = 0.66) {
        this.dir = getEposDir(quality)
        this.reset()
    }
}

const tmpV = Vec3()

export class HierarchyHelper {
    private sphere = Sphere3D()
    private normal = Vec3()
    private helperA = new EposHelper(this.quality)
    private helperB = new EposHelper(this.quality)

    private checkSide(p: Vec3) {
        return Vec3.dot(this.normal, Vec3.sub(tmpV, this.sphere.center, p)) > 0
    }

    includeStep(p: Vec3) {
        if (this.checkSide(p)) {
            this.helperA.includeStep(p)
        } else {
            this.helperB.includeStep(p)
        }
    }

    finishedIncludeStep() {
        this.helperA.finishedIncludeStep();
        this.helperB.finishedIncludeStep();
    }

    radiusStep(p: Vec3) {
        if (this.checkSide(p)) {
            this.helperA.radiusStep(p)
        } else {
            this.helperB.radiusStep(p)
        }
    }

    getSphere(): Sphere3D.Hierarchy {
        return {
            center: this.sphere.center,
            radius: this.sphere.radius,
            hierarchy: [this.helperA.getSphere(), this.helperB.getSphere()]
        }
    }

    reset(sphere: Sphere3D, normal: Vec3) {
        Sphere3D.copy(this.sphere, sphere)
        Vec3.copy(this.normal, normal)
        this.helperA.reset()
        this.helperB.reset()
    }

    constructor(private quality: EposQuality) { }
}

type EposQuality = '6' | '14' | '26' | '98'

function getEposDir(quality: EposQuality) {
    let dir: number[][]
    switch (quality) {
        case '6': dir = [ ...Type001 ]; break
        case '14': dir = [ ...Type001, ...Type111 ]; break
        case '26': dir = [ ...Type001, ...Type111, ...Type011 ]; break
        case '98': dir = [ ...Type001, ...Type111, ...Type011, ...Type012, ...Type112, ...Type122 ]; break
    }
    return dir.map(a => {
        const v = Vec3.create(a[0], a[1], a[2])
        return Vec3.normalize(v, v)
    })
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