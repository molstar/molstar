/*
 * Copyright (c) 2016 - now, David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */

import * as Coords from './coordinate';
import { SpacegroupCell } from '../../../../mol-math/geometry';

export interface Box<C extends Coords.Coord<Coords.Space>> { a: C, b: C }

export interface Cartesian extends Box<Coords.Cartesian> { }
export interface Fractional extends Box<Coords.Fractional> { }
export interface Grid<K> extends Box<Coords.Grid<K>> { }

// CONVERSIONS

export function cartesianToFractional(box: Cartesian, spacegroup: SpacegroupCell): Fractional {
    const { a: l, b: r } = box;
    const corners = [
        [l[0], l[1], l[2]],
        [r[0], l[1], l[2]],
        [l[0], r[1], l[2]],
        [l[0], l[1], r[2]],
        [r[0], r[1], l[2]],
        [r[0], l[1], r[2]],
        [l[0], r[1], r[2]],
        [r[0], r[1], r[2]],
    ].map(c => Coords.cartesianToFractional(Coords.cartesian(c[0], c[1], c[2]), spacegroup));
    return bounding(corners);
}

export function fractionalToGrid<K>(box: Fractional, domain: Coords.GridDomain<K>): Grid<K> {
    return { a: Coords.fractionalToGrid(box.a, domain, 'bottom'), b: Coords.fractionalToGrid(box.b, domain, 'top') };
}

export function gridToFractional<K>(box: Grid<K>): Fractional {
    return { a: Coords.gridToFractional(box.a), b: Coords.gridToFractional(box.b) };
}

export function fractionalBoxReorderAxes(box: Fractional, axisOrder: number[]) {
    const { a, b } = box;
    return {
        a: Coords.withCoord(a, a[axisOrder[0]], a[axisOrder[1]], a[axisOrder[2]]),
        b: Coords.withCoord(b, b[axisOrder[0]], b[axisOrder[1]], b[axisOrder[2]])
    };
}

export function expandGridBox<K>(box: Grid<K>, by: number) {
    const { a, b } = box;
    return {
        a: Coords.withCoord(a, a[0] - by, a[1] - by, a[2] - by),
        b: Coords.withCoord(b, b[0] + by, b[1] + by, b[2] + by)
    };
}

// MISC

export function shift<C extends Coords.Coord<S>, S extends Coords.Space>(box: Box<C>, offset: C): Box<C> {
    return { a: Coords.add(box.a, offset), b: Coords.add(box.b, offset) } as Box<C>;
}

export function clampGridToSamples<C extends Coords.Grid<K>, K>(box: Box<C>): Box<C> {
    return { a: Coords.clampGridToSamples(box.a), b: Coords.clampGridToSamples(box.b) } as Box<C>;
}

export function fractionalToDomain<K>(box: Fractional, kind: K, delta: Coords.Fractional): Coords.GridDomain<K> {
    const ds = Coords.fractional(box.b[0] - box.a[0], box.b[1] - box.a[1], box.b[2] - box.a[2]);
    return Coords.domain(kind, {
        delta,
        origin: box.a,
        dimensions: ds,
        sampleCount: Coords.sampleCounts(ds, delta)
    });
}

export function fractionalFromBlock(block: Coords.Grid<'Block'>): Fractional {
    const { domain } = block;
    const a = Coords.gridToFractional(block);
    const b = Coords.add(a, domain.delta);
    for (let i = 0; i < 3; i++) {
        b[i] = Math.min(b[i], domain.origin[i] + domain.dimensions[i]);
    }
    return { a, b };
}

export function bounding<C extends Coords.Coord<Coords.Space>>(xs: C[]): Box<C> {
    const a = Coords.clone(xs[0]);
    const b = Coords.clone(xs[0]);

    for (const x of xs) {
        for (let i = 0; i < 3; i++) {
            a[i] = Math.min(a[i], x[i]);
            b[i] = Math.max(b[i], x[i]);
        }
    }
    return { a, b };
}

export function areIntersecting<C extends Coords.Coord<S>, S extends Coords.Space>(box1: Box<C>, box2: Box<C>) {
    for (let i = 0; i < 3; i++) {
        const x = box1.a[i], y = box1.b[i];
        const u = box2.a[i], v = box2.b[i];
        if (x > v || y < u) return false;
    }
    return true;
}

export function intersect<C extends Coords.Coord<S>, S extends Coords.Space>(box1: Box<C>, box2: Box<C>): Box<C> | undefined {
    let a = Coords.clone(box1.a);
    let b = Coords.clone(box1.a);

    for (let i = 0; i < 3; i++) {
        const x = box1.a[i], y = box1.b[i];
        const u = box2.a[i], v = box2.b[i];
        if (x > v || y < u) return void 0;
        a[i] = Math.max(x, u);
        b[i] = Math.min(y, v);
    }
    return { a, b };
}

export function dimensions<C extends Coords.Coord<S>, S extends Coords.Space>(box: Box<C>): number[] {
    return [box.b[0] - box.a[0], box.b[1] - box.a[1], box.b[2] - box.a[2]];
}

export function volume<C extends Coords.Coord<S>, S extends Coords.Space>(box: Box<C>) {
    return (box.b[0] - box.a[0]) * (box.b[1] - box.a[1]) * (box.b[2] - box.a[2]);
}