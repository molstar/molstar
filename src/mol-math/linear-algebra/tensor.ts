/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Mat4, Vec3, Vec4 } from './3d'

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
        return { dimensions, axisOrderFastToSlow, axisOrderSlowToFast, accessDimensions, defaultCtor: ctor || Float64Array }
    }

    export function create(space: Space, data: Data): Tensor { return { space, data }; }

    export function Space(dimensions: number[], axisOrderSlowToFast: number[], ctor?: ArrayCtor): Space {
        const layout = Layout(dimensions, axisOrderSlowToFast, ctor);
        const { get, set } = accessors(layout);
        return { rank: dimensions.length, dimensions, axisOrderSlowToFast, create: creator(layout), get, set };
    }

    export function Data1(values: ArrayLike<number>): Data { return values as Data; }

    export function Vector(d: number, ctor?: ArrayCtor) { return Space([d], [0], ctor); }
    export function ColumnMajorMatrix(rows: number, cols: number, ctor?: ArrayCtor) { return Space([rows, cols], [1, 0], ctor); }
    export function RowMajorMatrix(rows: number, cols: number, ctor?: ArrayCtor) { return Space([rows, cols], [0, 1], ctor); }

    export function toMat4(space: Space, data: Tensor.Data): Mat4 {
        if (space.rank !== 2) throw new Error('Invalid tensor rank');
        const mat = Mat4.zero();
        const d0 = Math.min(4, space.dimensions[0]), d1 = Math.min(4, space.dimensions[1]);
        for (let i = 0; i < d0; i++) {
            for (let j = 0; j < d1; j++) Mat4.setValue(mat, i, j, space.get(data, i, j));
        }
        return mat;
    }

    export function toVec3(space: Space, data: Tensor.Data): Vec3 {
        if (space.rank !== 1) throw new Error('Invalid tensor rank');
        const vec = Vec3.zero();
        const d0 = Math.min(3, space.dimensions[0]);
        for (let i = 0; i < d0; i++) vec[i] = data[i];
        return vec;
    }

    export function toVec4(space: Space, data: Tensor.Data): Vec4 {
        if (space.rank !== 1) throw new Error('Invalid tensor rank');
        const vec = Vec4.zero();
        const d0 = Math.min(4, space.dimensions[0]);
        for (let i = 0; i < d0; i++) vec[i] = data[i];
        return vec;
    }

    export function areEqualExact(a: Tensor.Data, b: Tensor.Data) {
        const len = a.length;
        if (len !== b.length) return false;
        for (let i = 0; i < len; i++) if (a[i] !== b[i]) return false;
        return true;
    }

    function accessors(layout: Layout): { get: Space['get'], set: Space['set'] } {
        const { dimensions, axisOrderFastToSlow: ao } = layout;
        switch (dimensions.length) {
            case 1: return { get: (t, d) => t[d], set: (t, d, x) => t[d] = x };
            case 2: {
                // column major
                if (ao[0] === 0 && ao[1] === 1) {
                    const rows = dimensions[0];
                    return { get: (t, i, j) => t[j * rows + i], set: (t, i, j, x) => t[j * rows + i] = x };
                }
                if (ao[0] === 1 && ao[1] === 0) {
                    const cols = dimensions[1];
                    return { get: (t, i, j) => t[i * cols + j], set: (t, i, j, x) => t[i * cols + j] = x };
                }
                throw new Error('bad axis order')
            }
            case 3: {
                if (ao[0] === 0 && ao[1] === 1 && ao[2] === 2) { // 012 ijk
                    const u = dimensions[0], v = dimensions[1], uv = u * v;
                    return { get: (t, i, j, k) => t[i + j * u + k * uv], set: (t, i, j, k, x ) => t[i + j * u + k * uv] = x };
                }
                if (ao[0] === 0 && ao[1] === 2 && ao[2] === 1) { // 021 ikj
                    const u = dimensions[0], v = dimensions[2], uv = u * v;
                    return { get: (t, i, j, k) => t[i + k * u + j * uv], set: (t, i, j, k, x ) => t[i + k * u + j * uv] = x };
                }
                if (ao[0] === 1 && ao[1] === 0 && ao[2] === 2) { // 102 jik
                    const u = dimensions[1], v = dimensions[0], uv = u * v;
                    return { get: (t, i, j, k) => t[j + i * u + k * uv], set: (t, i, j, k, x ) => t[j + i * u + k * uv] = x };
                }
                if (ao[0] === 1 && ao[1] === 2 && ao[2] === 0) { // 120 jki
                    const u = dimensions[1], v = dimensions[2], uv = u * v;
                    return { get: (t, i, j, k) => t[j + k * u + i * uv], set: (t, i, j, k, x ) => t[j + k * u + i * uv] = x };
                }
                if (ao[0] === 2 && ao[1] === 0 && ao[2] === 1) { // 201 kij
                    const u = dimensions[2], v = dimensions[0], uv = u * v;
                    return { get: (t, i, j, k) => t[k + i * u + j * uv], set: (t, i, j, k, x ) => t[k + i * u + j * uv] = x };
                }
                if (ao[0] === 2 && ao[1] === 1 && ao[2] === 0) { // 210 kji
                    const u = dimensions[2], v = dimensions[1], uv = u * v;
                    return { get: (t, i, j, k) => t[k + j * u + i * uv], set: (t, i, j, k, x ) => t[k + j * u + i * uv] = x };
                }
                throw new Error('bad axis order')
            }
            default: return {
                get: (t, ...c) => t[dataOffset(layout, c)],
                set: (t, ...c) => t[dataOffset(layout, c)] = c[c.length - 1]
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
}