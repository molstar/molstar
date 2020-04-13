/*
 * Copyright (c) 2016 - now, David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */

import { Mat4, Vec3 } from '../../../../mol-math/linear-algebra';
import { SpacegroupCell } from '../../../../mol-math/geometry';

/** Information about a region sampled in fractional coordinates */
export interface GridInfo {
    /** Origin in fractional coords. */
    origin: Fractional,
    /** Box dimensions in fractional coords. */
    dimensions: Fractional,
    /** Grid delta in fractional coordinates along each axis (in axis order) */
    delta: Fractional,
    /** Sample count of the grid box */
    sampleCount: number[]
}

/**
 * Grid domain with the supplied info and "kind".
 * The "kind" is used so that the TypeScript compiler
 * can distinguish between different types of grids,
 * e.g. GridDomain<'Data'>, GridDomain<'Query'>, GridDomain<'Block'>, etc.
 */
export interface GridDomain<K> extends GridInfo { kind: K, sampleVolume: number }

export const enum Space { Cartesian, Fractional, Grid }
export interface Coord<S extends Space> { kind: S, '0': number, '1': number, '2': number, [index: number]: number }
export interface Cartesian extends Coord<Space.Cartesian> { }
export interface Fractional extends Coord<Space.Fractional> { }
export interface Grid<K> extends Coord<Space.Grid> { domain: GridDomain<K> }

// CONSTRUCTORS

export function domain<K>(kind: K, info: GridInfo): GridDomain<K> {
    const sc = info.sampleCount;
    return {
        kind,
        delta: info.delta,
        dimensions: info.dimensions,
        origin: info.origin,
        sampleCount: info.sampleCount,
        sampleVolume: sc[0] * sc[1] * sc[2]
    };
}

export function cartesian(x: number, y: number, z: number): Cartesian {
    return { 0: x, 1: y, 2: z, kind: Space.Cartesian };
}

export function fractional(x: number, y: number, z: number): Fractional {
    return { 0: x, 1: y, 2: z, kind: Space.Fractional };
}

export function grid<K>(domain: GridDomain<K>, x: number, y: number, z: number): Grid<K> {
    return { 0: x, 1: y, 2: z, kind: Space.Grid, domain };
}

export function withCoord<C extends (Coord<Space> | Grid<any>)>(a: C, x: number, y: number, z: number): C {
    switch (a.kind) {
        case Space.Cartesian: return cartesian(x, y, z) as C;
        case Space.Fractional: return fractional(x, y, z) as C;
        case Space.Grid: return grid((a as Grid<any>).domain, x, y, z) as C;
    }
}

export function clone<C extends (Coord<Space> | Grid<any>)>(a: C): C {
    return withCoord(a, a[0], a[1], a[2]);
}

// CONVERSIONS

export function cartesianToFractional(a: Cartesian, spacegroup: SpacegroupCell): Fractional {
    const coord = Helpers.transform(a, spacegroup.toFractional);
    return fractional(coord[0], coord[1], coord[2]);
}

export function fractionalToGrid<K>(a: Fractional, domain: GridDomain<K>, snap: 'bottom' | 'top'): Grid<K> {
    const { origin, delta } = domain;
    const coord = grid(domain, 0.1, 0.1, 0.1);
    for (let i = 0; i < 3; i++) {
        coord[i] = Helpers.snap((a[i] - origin[i]) / delta[i], snap);
    }
    return coord;
}

export function gridToFractional<K>(a: Grid<K>): Fractional {
    const { origin, delta } = a.domain;
    const coord = fractional(0.1, 0.1, 0.1);
    for (let i = 0; i < 3; i++) {
        coord[i] = a[i] * delta[i] + origin[i];
    }
    return coord;
}

// MISC

export function clampGridToSamples<K>(a: Grid<K>): Grid<K> {
    const { sampleCount } = a.domain;
    const coord = withCoord(a, 0, 0, 0);
    for (let i = 0; i < 3; i++) {
        if (a[i] < 0) coord[i] = 0;
        else if (a[i] > sampleCount[i]) coord[i] = sampleCount[i];
        else coord[i] = a[i];
    }
    return coord;
}

export function add<S extends Space>(a: Coord<S>, b: Coord<S>): Coord<S> {
    return withCoord(a, a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

export function sub<S extends Space>(a: Coord<S>, b: Coord<S>): Coord<S> {
    return withCoord(a, a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

export function invert<S extends Space>(a: Coord<S>): Coord<S> {
    return withCoord(a, -a[0], -a[1], -a[2]);
}

/** Maps each grid point to a unique integer */
export function linearGridIndex<K>(a: Grid<K>) {
    const samples = a.domain.sampleCount;
    return a[0] + samples[0] * (a[1] + a[2] * samples[1]);
}

export function gridMetrics(dimensions: { [i: number]: number }) {
    return {
        sizeX: dimensions[0],
        sizeXY: dimensions[0] * dimensions[1],
        sizeXYZ: dimensions[0] * dimensions[1] * dimensions[2]
    };
}

export function sampleCounts(dimensions: Fractional, delta: Fractional) {
    return [
        Helpers.snap(dimensions[0] / delta[0], 'top'),
        Helpers.snap(dimensions[1] / delta[1], 'top'),
        Helpers.snap(dimensions[2] / delta[2], 'top')
    ];
}

// to prevent floating point rounding errors
export function round(v: number) {
    return Math.round(10000000 * v) / 10000000;
}

namespace Helpers {
    export function transform(x: { [index: number]: number }, matrix: Mat4) {
        return Vec3.transformMat4(Vec3.zero(), x as Vec3, matrix);
    }

    export function snap(v: number, to: 'bottom' | 'top') {
        switch (to) {
            case 'bottom': return Math.floor(round(v)) | 0;
            case 'top': return Math.ceil(round(v)) | 0;
        }
    }
}