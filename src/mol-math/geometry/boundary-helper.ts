/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ReadonlyVec3, Vec3 } from '../linear-algebra/3d/vec3';
import { CentroidHelper } from './centroid-helper';
import { Sphere3D } from '../geometry/primitives/sphere3d';
import { Box3D } from './primitives/box3d';
import { memoize1 } from '../../mol-util/memoize';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3dot = Vec3.dot;
const v3copy = Vec3.copy;
const v3scaleAndSub = Vec3.scaleAndSub;
const v3scaleAndAdd = Vec3.scaleAndAdd;

// implementing http://www.ep.liu.se/ecp/034/009/ecp083409.pdf

export class BoundaryHelper {
    private readonly dir: ReadonlyArray<Vec3>;
    private readonly dirLength: number;

    private minDist: number[] = [];
    private maxDist: number[] = [];
    private extrema: Vec3[] = [];
    centroidHelper = new CentroidHelper();

    private computeExtrema(i: number, p: Vec3) {
        const d = v3dot(this.dir[i], p);

        if (d < this.minDist[i]) {
            this.minDist[i] = d;
            v3copy(this.extrema[i * 2], p);
        }
        if (d > this.maxDist[i]) {
            this.maxDist[i] = d;
            v3copy(this.extrema[i * 2 + 1], p);
        }
    }

    private computeSphereExtrema(i: number, center: Vec3, radius: number) {
        const di = this.dir[i];
        const d = v3dot(di, center);

        if (d - radius < this.minDist[i]) {
            this.minDist[i] = d - radius;
            v3scaleAndSub(this.extrema[i * 2], center, di, radius);
        }
        if (d + radius > this.maxDist[i]) {
            this.maxDist[i] = d + radius;
            v3scaleAndAdd(this.extrema[i * 2 + 1], center, di, radius);
        }
    }

    includeSphere(s: Sphere3D) {
        if (Sphere3D.hasExtrema(s) && s.extrema.length > 1) {
            for (const e of s.extrema) {
                this.includePosition(e);
            }
        } else {
            this.includePositionRadius(s.center, s.radius);
        }
    }

    includePosition(p: Vec3) {
        for (let i = 0; i < this.dirLength; ++i) {
            this.computeExtrema(i, p);
        }
    }

    includePositionRadius(center: Vec3, radius: number) {
        for (let i = 0; i < this.dirLength; ++i) {
            this.computeSphereExtrema(i, center, radius);
        }
    }

    finishedIncludeStep() {
        if (this.dirLength < 49) {
            this.refineExtrema();
        }

        for (let i = 0; i < this.extrema.length; i++) {
            this.centroidHelper.includeStep(this.extrema[i]);
        }
        this.centroidHelper.finishedIncludeStep();
    }

    private refineExtrema() {
        const dir98 = getEposDir('98');
        const minDist: number[] = [];
        const maxDist: number[] = [];
        const extrema: Vec3[] = [];

        for (let i = 0, il = dir98.length; i < il; ++i) {
            minDist[i] = Infinity;
            maxDist[i] = -Infinity;
            extrema[i * 2] = Vec3();
            extrema[i * 2 + 1] = Vec3();

            const di = dir98[i];
            for (let j = 0, jl = this.extrema.length; j < jl; ++j) {
                const p = this.extrema[j];
                const d = v3dot(di, p);

                if (d < minDist[i]) {
                    minDist[i] = d;
                    v3copy(extrema[i * 2], p);
                }
                if (d > maxDist[i]) {
                    maxDist[i] = d;
                    v3copy(extrema[i * 2 + 1], p);
                }
            }
        }

        this.extrema = extrema;
    }

    radiusSphere(s: Sphere3D) {
        if (Sphere3D.hasExtrema(s) && s.extrema.length > 1) {
            for (const e of s.extrema) {
                this.radiusPosition(e);
            }
        } else {
            this.radiusPositionRadius(s.center, s.radius);
        }
    }

    radiusPosition(p: Vec3) {
        this.centroidHelper.radiusStep(p);
    }

    radiusPositionRadius(center: Vec3, radius: number) {
        this.centroidHelper.radiusSphereStep(center, radius);
    }

    getSphere(sphere?: Sphere3D) {
        return Sphere3D.setExtrema(this.centroidHelper.getSphere(sphere), this.extrema.slice());
    }

    getBox(box?: Box3D) {
        return Box3D.fromVec3Array(box || Box3D(), this.extrema);
    }

    reset() {
        // discard any extrema added by a previous `refineExtrema` call
        this.extrema.length = this.dirLength * 2;

        for (let i = 0; i < this.dirLength; ++i) {
            this.minDist[i] = Infinity;
            this.maxDist[i] = -Infinity;
            this.extrema[i * 2] = Vec3();
            this.extrema[i * 2 + 1] = Vec3();
        }
        this.centroidHelper.reset();
    }

    constructor(quality: EposQuality) {
        this.dir = getEposDir(quality);
        this.dirLength = this.dir.length;
        this.reset();
    }
}

type EposQuality = '6' | '14' | '26' | '98'

const getEposDir = memoize1((quality: EposQuality): ReadonlyArray<ReadonlyVec3> => {
    let dir: number[][];
    switch (quality) {
        case '6': dir = [...Type001]; break;
        case '14': dir = [...Type001, ...Type111]; break;
        case '26': dir = [...Type001, ...Type111, ...Type011]; break;
        case '98': dir = [...Type001, ...Type111, ...Type011, ...Type012, ...Type112, ...Type122]; break;
    }
    return dir.map(a => {
        const v = Vec3.create(a[0], a[1], a[2]);
        return Vec3.normalize(v, v);
    });
});

const Type001 = [
    [1, 0, 0], [0, 1, 0], [0, 0, 1]
];

const Type111 = [
    [1, 1, 1], [-1, 1, 1], [-1, -1, 1], [1, -1, 1]
];

const Type011 = [
    [1, 1, 0], [1, -1, 0], [1, 0, 1], [1, 0, -1], [0, 1, 1], [0, 1, -1]
];

const Type012 = [
    [0, 1, 2], [0, 2, 1], [1, 0, 2], [2, 0, 1], [1, 2, 0], [2, 1, 0],
    [0, 1, -2], [0, 2, -1], [1, 0, -2], [2, 0, -1], [1, -2, 0], [2, -1, 0]
];

const Type112 = [
    [1, 1, 2], [2, 1, 1], [1, 2, 1], [1, -1, 2], [1, 1, -2], [1, -1, -2],
    [2, -1, 1], [2, 1, -1], [2, -1, -1], [1, -2, 1], [1, 2, -1], [1, -2, -1]
];

const Type122 = [
    [2, 2, 1], [1, 2, 2], [2, 1, 2], [2, -2, 1], [2, 2, -1], [2, -2, -1],
    [1, -2, 2], [1, 2, -2], [1, -2, -2], [2, -1, 2], [2, 1, -2], [2, -1, -2]
];
