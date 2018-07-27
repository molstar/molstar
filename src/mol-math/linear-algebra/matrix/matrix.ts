/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

interface Matrix { data: Helpers.NumberArray, size: number, cols: number, rows: number }

namespace Matrix {
    export function create(cols: number, rows: number, ctor: { new (size: number): Helpers.NumberArray } = Float32Array): Matrix {
        const size = cols * rows
        return { data: new ctor(size), size, cols, rows }
    }

    export function fromArray(data: Helpers.NumberArray, cols: number, rows: number): Matrix {
        return { data, size: cols * rows, cols, rows }
    }

    export function transpose(out: Matrix, mat: Matrix): Matrix {
        const nrows = mat.rows, ncols = mat.cols
        const md = mat.data, mtd = out.data

        for (let i = 0, mi = 0, mti = 0; i < nrows; mti += 1, mi += ncols, ++i) {
            let ri = mti
            for (let j = 0; j < ncols; ri += nrows, j++) mtd[ri] = md[mi + j]
        }
        return mat
    }

    /** out = matA * matB' */
    export function multiplyABt (out: Matrix, matA: Matrix, matB: Matrix) {
        const ncols = matA.cols, nrows = matA.rows, mrows = matB.rows
        const ad = matA.data, bd = matB.data, cd = out.data

        for (let i = 0, matAp = 0, outP = 0; i < nrows; matAp += ncols, i++) {
            for (let pB = 0, j = 0; j < mrows; outP++, j++) {
                let sum = 0.0
                let pMatA = matAp
                for (let k = 0; k < ncols; pMatA++, pB++, k++) {
                    sum += ad[pMatA] * bd[pB]
                }
                cd[outP] = sum
            }
        }
        return out
    }

    /** Get the mean of rows in `mat` */
    export function meanRows (mat: Matrix) {
        const nrows = mat.rows, ncols = mat.cols
        const md = mat.data
        const mean = new Array(ncols)

        for (let j = 0; j < ncols; ++j) mean[ j ] = 0.0
        for (let i = 0, p = 0; i < nrows; ++i) {
            for (let j = 0; j < ncols; ++j, ++p) mean[ j ] += md[ p ]
        }
        for (let j = 0; j < ncols; ++j) mean[ j ] /= nrows

        return mean
    }

    /** Subtract `row` from all rows in `mat` */
    export function subRows (mat: Matrix, row: Helpers.NumberArray) {
        const nrows = mat.rows, ncols = mat.cols
        const md = mat.data

        for (let i = 0, p = 0; i < nrows; ++i) {
            for (let j = 0; j < ncols; ++j, ++p) md[ p ] -= row[ j ]
        }
        return mat
    }
}

export default Matrix