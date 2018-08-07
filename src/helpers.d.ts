/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

declare module Helpers {
    export type Mutable<T> = {
        -readonly [P in keyof T]: T[P]
    }
    export type TypedArray = Int8Array | Uint8Array | Int16Array | Uint16Array | Int32Array | Uint32Array | Uint8ClampedArray | Float32Array | Float64Array
    export type NumberArray = TypedArray | number[]
    export type UintArray = Uint8Array | Uint16Array | Uint32Array | number[]
    export type ValueOf<T> = T[keyof T]
    export type ArrayCtor<T> = { new(size: number): { [i: number]: T, length: number } }
    /** assignable ArrayLike version */
    export type ArrayLike<T> =  { [i: number]: T, length: number }
}