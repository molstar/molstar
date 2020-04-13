/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { NumberArray } from '../../../mol-util/type-helpers';
import { Vec } from '../3d';

interface Matrix<N extends number = number, M extends number = number> {
    data: NumberArray,
    size: number,
    cols: N,
    rows: M
}

namespace Matrix {
    export function create<N extends number, M extends number>(cols: N, rows: M, ctor: { new (size: number): NumberArray } = Float32Array): Matrix<N, M> {
        const size = cols * rows;
        return { data: new ctor(size), size, cols, rows };
    }

    /** Get element assuming data are stored in column-major order */
    export function get(m: Matrix, i: number, j: number) {
        return m.data[m.rows * j + i];
    }

    /** Set element assuming data are stored in column-major order */
    export function set<N extends number, M extends number>(m: Matrix<N, M>, i: number, j: number, value: number) {
        m.data[m.rows * j + i] = value;
        return m;
    }

    /** Add to element assuming data are stored in column-major order */
    export function add<N extends number, M extends number>(m: Matrix<N, M>, i: number, j: number, value: number) {
        m.data[m.rows * j + i] += value;
        return m;
    }

    /** Zero out the matrix */
    export function makeZero<N extends number, M extends number>(m: Matrix<N, M>) {
        m.data.fill(0.0);
        return m;
    }

    export function clone<N extends number, M extends number>(m: Matrix<N, M>): Matrix<N, M> {
        return { data: m.data.slice(), size: m.size, cols: m.cols, rows: m.rows };
    }

    export function fromArray<N extends number, M extends number>(data: NumberArray, cols: N, rows: M): Matrix<N, M> {
        return { data, size: cols * rows, cols, rows };
    }

    export function transpose<N extends number, M extends number>(out: Matrix<M, N>, mat: Matrix<N, M>): Matrix<M, N> {
        if (out.cols !== mat.rows || out.rows !== mat.cols) {
            throw new Error('transpose: matrix dimensions incompatible');
        }
        if (out.data === mat.data) {
            throw new Error('transpose: matrices share memory');
        }
        const nrows = mat.rows, ncols = mat.cols;
        const md = mat.data, mtd = out.data;
        for (let i = 0, mi = 0, mti = 0; i < nrows; mti += 1, mi += ncols, ++i) {
            let ri = mti;
            for (let j = 0; j < ncols; ri += nrows, j++) mtd[ri] = md[mi + j];
        }
        return out;
    }

    /** out = matA * matB' */
    export function multiplyABt<NA extends number, NB extends number, M extends number>(out: Matrix<M, M>, matA: Matrix<NA, M>, matB: Matrix<NB, M>): Matrix<M, M> {
        const ncols = matA.cols, nrows = matA.rows, mrows = matB.rows;
        const ad = matA.data, bd = matB.data, cd = out.data;

        for (let i = 0, matAp = 0, outP = 0; i < nrows; matAp += ncols, i++) {
            for (let pB = 0, j = 0; j < mrows; outP++, j++) {
                let sum = 0.0;
                let pMatA = matAp;
                for (let k = 0; k < ncols; pMatA++, pB++, k++) {
                    sum += ad[pMatA] * bd[pB];
                }
                cd[outP] = sum;
            }
        }
        return out;
    }

    /** Get the mean of rows in `mat` */
    export function meanRows<N extends number, M extends number, V extends Vec<N>>(mat: Matrix<N, M>): V {
        const nrows = mat.rows, ncols = mat.cols;
        const md = mat.data;
        const mean = new Array(ncols) as V;

        for (let j = 0; j < ncols; ++j) mean[ j ] = 0.0;
        for (let i = 0, p = 0; i < nrows; ++i) {
            for (let j = 0; j < ncols; ++j, ++p) mean[ j ] += md[ p ];
        }
        for (let j = 0; j < ncols; ++j) mean[ j ] /= nrows;

        return mean;
    }

    /** Subtract `row` from all rows in `mat` */
    export function subRows<N extends number, M extends number>(mat: Matrix<N, M>, row: NumberArray) {
        const nrows = mat.rows, ncols = mat.cols;
        const md = mat.data;

        for (let i = 0, p = 0; i < nrows; ++i) {
            for (let j = 0; j < ncols; ++j, ++p) md[ p ] -= row[ j ];
        }
        return mat;
    }
}

export default Matrix;