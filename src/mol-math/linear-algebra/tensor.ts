/**
 * Copyright (c) 2017-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Vec3, Vec4, Mat3 } from './3d';

export interface Tensor { data: Tensor.Data, space: Tensor.Space }

export namespace Tensor {
    export type ArrayCtor = { new (size: number): ArrayLike<number> }

    export interface Data extends Array<number> { '@type': 'tensor' }

    export interface Space {
        readonly rank: number,
        readonly dimensions: ReadonlyArray<number>,
        readonly axisOrderSlowToFast: ReadonlyArray<number>,
        create(array?: ArrayCtor): Tensor.Data,
        get(data: Tensor.Data, ...coords: number[]): number
        set(data: Tensor.Data, ...coordsAndValue: number[]): number
        add(data: Tensor.Data, ...coordsAndValue: number[]): number
        dataOffset(...coords: number[]): number,
        getCoords(dataOffset: number, coords: { [i: number]: number }): number[]
    }

    interface Layout {
        dimensions: number[],
        axisOrderSlowToFast: number[],
        axisOrderFastToSlow: number[],
        accessDimensions: number[],
        // if not specified, use Float64Array
        defaultCtor: ArrayCtor
    }

    function Layout(dimensions: number[], axisOrderSlowToFast: number[], ctor?: ArrayCtor): Layout {
        // need to reverse the axis order for better access.
        const axisOrderFastToSlow: number[] = [];
        for (let i = 0; i < axisOrderSlowToFast.length; i++) axisOrderFastToSlow[i] = axisOrderSlowToFast[axisOrderSlowToFast.length - i - 1];

        const accessDimensions = [1];
        for (let i = 1; i < dimensions.length; i++) accessDimensions[i] = dimensions[axisOrderFastToSlow[i - 1]];
        return { dimensions, axisOrderFastToSlow, axisOrderSlowToFast, accessDimensions, defaultCtor: ctor || Float64Array };
    }

    export function create(space: Space, data: Data): Tensor { return { space, data }; }

    export function Space(dimensions: number[], axisOrderSlowToFast: number[], ctor?: ArrayCtor): Space {
        const layout = Layout(dimensions, axisOrderSlowToFast, ctor);
        const { get, set, add, dataOffset, getCoords } = accessors(layout);
        return { rank: dimensions.length, dimensions, axisOrderSlowToFast, create: creator(layout), get, set, add, dataOffset, getCoords };
    }

    export function Data1(values: ArrayLike<number>): Data { return values as Data; }

    export function Vector(d: number, ctor?: ArrayCtor) { return Space([d], [0], ctor); }
    export function ColumnMajorMatrix(rows: number, cols: number, ctor?: ArrayCtor) { return Space([rows, cols], [1, 0], ctor); }
    export function RowMajorMatrix(rows: number, cols: number, ctor?: ArrayCtor) { return Space([rows, cols], [0, 1], ctor); }

    export function toMat4(out: Mat4, space: Space, data: Tensor.Data): Mat4 {
        if (space.rank !== 2) throw new Error('Invalid tensor rank');
        const d0 = Math.min(4, space.dimensions[0]), d1 = Math.min(4, space.dimensions[1]);
        for (let i = 0; i < d0; i++) {
            for (let j = 0; j < d1; j++) Mat4.setValue(out, i, j, space.get(data, i, j));
        }
        return out;
    }

    export function toMat3(out: Mat3, space: Space, data: Tensor.Data): Mat3 {
        if (space.rank !== 2) throw new Error('Invalid tensor rank');
        const d0 = Math.min(3, space.dimensions[0]), d1 = Math.min(3, space.dimensions[1]);
        for (let i = 0; i < d0; i++) {
            for (let j = 0; j < d1; j++) Mat3.setValue(out, i, j, space.get(data, i, j));
        }
        return out;
    }

    export function toVec3(out: Vec3, space: Space, data: Tensor.Data): Vec3 {
        if (space.rank !== 1) throw new Error('Invalid tensor rank');
        const d0 = Math.min(3, space.dimensions[0]);
        for (let i = 0; i < d0; i++) out[i] = data[i];
        return out;
    }

    export function toVec4(out: Vec4, space: Space, data: Tensor.Data): Vec4 {
        if (space.rank !== 1) throw new Error('Invalid tensor rank');
        const d0 = Math.min(4, space.dimensions[0]);
        for (let i = 0; i < d0; i++) out[i] = data[i];
        return out;
    }

    export function areEqualExact(a: Tensor.Data, b: Tensor.Data) {
        const len = a.length;
        if (len !== b.length) return false;
        for (let i = 0; i < len; i++) if (a[i] !== b[i]) return false;
        return true;
    }

    function accessors(layout: Layout): { get: Space['get'], set: Space['set'], add: Space['add'], dataOffset: Space['dataOffset'], getCoords: Space['getCoords'] } {
        const { dimensions, axisOrderFastToSlow: ao } = layout;
        switch (dimensions.length) {
            case 1: return {
                get: (t, d) => t[d],
                set: (t, d, x) => t[d] = x,
                add: (t, d, x) => t[d] += x,
                dataOffset: (d) => d,
                getCoords: (o, c) => { c[0] = o; return c as number[]; }
            };
            case 2: {
                // column major
                if (ao[0] === 0 && ao[1] === 1) {
                    const rows = dimensions[0];
                    return {
                        get: (t, i, j) => t[j * rows + i],
                        set: (t, i, j, x) => t[j * rows + i] = x,
                        add: (t, i, j, x) => t[j * rows + i] += x,
                        dataOffset: (i, j) => j * rows + i,
                        getCoords: (o, c) => { c[0] = o % rows; c[1] = Math.floor(o / rows) ; return c as number[]; }
                    };
                }
                if (ao[0] === 1 && ao[1] === 0) {
                    const cols = dimensions[1];
                    return {
                        get: (t, i, j) => t[i * cols + j],
                        set: (t, i, j, x) => t[i * cols + j] = x,
                        add: (t, i, j, x) => t[i * cols + j] += x,
                        dataOffset: (i, j) => i * cols + j,
                        getCoords: (o, c) => { c[0] = Math.floor(o / cols); c[1] = o % cols; return c as number[]; }
                    };
                }
                throw new Error('bad axis order');
            }
            case 3: {
                if (ao[0] === 0 && ao[1] === 1 && ao[2] === 2) { // 012 ijk
                    const u = dimensions[0], v = dimensions[1], uv = u * v;
                    return {
                        get: (t, i, j, k) => t[i + j * u + k * uv],
                        set: (t, i, j, k, x ) => t[i + j * u + k * uv] = x,
                        add: (t, i, j, k, x ) => t[i + j * u + k * uv] += x,
                        dataOffset: (i, j, k) => i + j * u + k * uv,
                        getCoords: (o, c) => {
                            const p = Math.floor(o / u);
                            c[0] = o % u;
                            c[1] = p % v;
                            c[2] = Math.floor(p / v);
                            return c as number[];
                        }
                    };
                }
                if (ao[0] === 0 && ao[1] === 2 && ao[2] === 1) { // 021 ikj
                    const u = dimensions[0], v = dimensions[2], uv = u * v;
                    return {
                        get: (t, i, j, k) => t[i + k * u + j * uv],
                        set: (t, i, j, k, x ) => t[i + k * u + j * uv] = x,
                        add: (t, i, j, k, x ) => t[i + k * u + j * uv] += x,
                        dataOffset: (i, j, k) => i + k * u + j * uv,
                        getCoords: (o, c) => {
                            const p = Math.floor(o / u);
                            c[0] = o % u;
                            c[1] = Math.floor(p / v);
                            c[2] = p % v;
                            return c as number[];
                        }
                    };
                }
                if (ao[0] === 1 && ao[1] === 0 && ao[2] === 2) { // 102 jik
                    const u = dimensions[1], v = dimensions[0], uv = u * v;
                    return {
                        get: (t, i, j, k) => t[j + i * u + k * uv],
                        set: (t, i, j, k, x ) => t[j + i * u + k * uv] = x,
                        add: (t, i, j, k, x ) => t[j + i * u + k * uv] += x,
                        dataOffset: (i, j, k) => j + i * u + k * uv,
                        getCoords: (o, c) => {
                            const p = Math.floor(o / u);
                            c[0] = p % v;
                            c[1] = o % u;
                            c[2] = Math.floor(p / v);
                            return c as number[];
                        }
                    };
                }
                if (ao[0] === 1 && ao[1] === 2 && ao[2] === 0) { // 120 jki
                    const u = dimensions[1], v = dimensions[2], uv = u * v;
                    return {
                        get: (t, i, j, k) => t[j + k * u + i * uv],
                        set: (t, i, j, k, x ) => t[j + k * u + i * uv] = x,
                        add: (t, i, j, k, x ) => t[j + k * u + i * uv] += x,
                        dataOffset: (i, j, k) => j + k * u + i * uv,
                        getCoords: (o, c) => {
                            const p = Math.floor(o / u);
                            c[0] = Math.floor(p / v);
                            c[1] = o % u;
                            c[2] = p % v;
                            return c as number[];
                        }
                    };
                }
                if (ao[0] === 2 && ao[1] === 0 && ao[2] === 1) { // 201 kij
                    const u = dimensions[2], v = dimensions[0], uv = u * v;
                    return {
                        get: (t, i, j, k) => t[k + i * u + j * uv],
                        set: (t, i, j, k, x ) => t[k + i * u + j * uv] = x,
                        add: (t, i, j, k, x ) => t[k + i * u + j * uv] += x,
                        dataOffset: (i, j, k) => k + i * u + j * uv,
                        getCoords: (o, c) => {
                            const p = Math.floor(o / u);
                            c[0] = p % v;
                            c[1] = Math.floor(p / v);
                            c[2] = o % u;
                            return c as number[];
                        }
                    };
                }
                if (ao[0] === 2 && ao[1] === 1 && ao[2] === 0) { // 210 kji
                    const u = dimensions[2], v = dimensions[1], uv = u * v;
                    return {
                        get: (t, i, j, k) => t[k + j * u + i * uv],
                        set: (t, i, j, k, x ) => t[k + j * u + i * uv] = x,
                        add: (t, i, j, k, x ) => t[k + j * u + i * uv] += x,
                        dataOffset: (i, j, k) => k + j * u + i * uv,
                        getCoords: (o, c) => {
                            const p = Math.floor(o / u);
                            c[0] = Math.floor(p / v);
                            c[1] = p % v;
                            c[2] = o % u;
                            return c as number[];
                        }
                    };
                }
                throw new Error('bad axis order');
            }
            default: return {
                get: (t, ...c) => t[dataOffset(layout, c)],
                set: (t, ...c) => t[dataOffset(layout, c)] = c[c.length - 1],
                add: (t, ...c) => t[dataOffset(layout, c)] += c[c.length - 1],
                dataOffset: (...c) => dataOffset(layout, c),
                getCoords: (o, c) => getCoords(layout, o, c as number[]),
            };
        }
    }

    function creator(layout: Layout): Space['create'] {
        const { dimensions: ds } = layout;
        let size = 1;
        for (let i = 0, _i = ds.length; i < _i; i++) size *= ds[i];
        return ctor => new (ctor || layout.defaultCtor)(size) as Tensor.Data;
    }

    function dataOffset(layout: Layout, coord: number[]) {
        const { accessDimensions: acc, axisOrderFastToSlow: ao } = layout;
        const d = acc.length - 1;
        let o = acc[d] * coord[ao[d]];
        for (let i = d - 1; i >= 0; i--) {
            o = (o + coord[ao[i]]) * acc[i];
        }
        return o;
    }

    function getCoords(layout: Layout, o: number, coords: number[]) {
        const { dimensions: dim, axisOrderFastToSlow: ao } = layout;
        const d = dim.length;

        let c = o;
        for (let i = 0; i < d; i++) {
            const d = dim[ao[i]];
            coords[ao[i]] = c % d;
            c = Math.floor(c / d);
        }
        coords[ao[d + 1]] = c;

        return coords;
    }

    // Convers "slow to fast" axis order to "fast to slow" and vice versa.
    export function invertAxisOrder(v: number[]) {
        const ret: number[] = [];
        for (let i = 0; i < v.length; i++) {
            ret[i] = v[v.length - i - 1];
        }
        return ret;
    }

    function reorder(xs: number[], indices: number[]) {
        const ret: number[] = [];
        for (let i = 0; i < xs.length; i++) ret[i] = xs[indices[i]];
        return ret;
    }

    export function convertToCanonicalAxisIndicesFastToSlow(order: number[]) {
        const indices = new Int32Array(order.length) as any as number[];
        for (let i = 0; i < order.length; i++) indices[order[i]] = i;
        return (xs: number[]) => reorder(xs, indices);
    }

    export function convertToCanonicalAxisIndicesSlowToFast(order: number[]) {
        const indices = new Int32Array(order.length) as any as number[];
        for (let i = 0; i < order.length; i++) indices[order[order.length - i - 1]] = i;
        return (xs: number[]) => reorder(xs, indices);
    }
}