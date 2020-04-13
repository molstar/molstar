/**
 * Copyright (c) 2017-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { hash2 } from '../util';

/**
 * Represents a pair of two integers as a double,
 * Caution: === does not work, because of NaN, use IntTuple.areEqual for equality
 */
interface IntTuple { '@type': 'int-tuple' }

namespace IntTuple {
    export const Zero: IntTuple = 0 as any;

    const { _int32, _float64, _int32_1, _float64_1 } = (function() {
        const data = new ArrayBuffer(8);
        const data_1 = new ArrayBuffer(8);
        return {
            _int32: new Int32Array(data),
            _float64: new Float64Array(data),
            _int32_1: new Int32Array(data_1),
            _float64_1: new Float64Array(data_1)
        };
    }());

    export function is(x: any): x is IntTuple {
        return typeof x === 'number';
    }

    export function create(fst: number, snd: number): IntTuple {
        _int32[0] = fst;
        _int32[1] = snd;
        return _float64[0] as any;
    }

    export function fst(t: IntTuple): number {
        _float64[0] = t as any;
        return _int32[0];
    }

    export function snd(t: IntTuple): number {
        _float64[0] = t as any;
        return _int32[1];
    }

    /** Normal equality does not work, because NaN === NaN ~> false */
    export function areEqual(a: IntTuple, b: IntTuple) {
        _float64[0] = a as any;
        _float64_1[0] = b as any;
        return _int32[0] === _int32_1[0] && _int32[1] === _int32_1[1];
    }

    export function compare(a: IntTuple, b: IntTuple) {
        _float64[0] = a as any;
        _float64_1[0] = b as any;
        const x = _int32[0] - _int32_1[0];
        if (x !== 0) return x;
        return _int32[1] - _int32_1[1];
    }

    export function compareInArray(xs: ArrayLike<IntTuple>, i: number, j: number) {
        _float64[0] = xs[i] as any;
        _float64_1[0] = xs[j] as any;
        const x = _int32[0] - _int32_1[0];
        if (x !== 0) return x;
        return _int32[1] - _int32_1[1];
    }

    export function hashCode(t: IntTuple) {
        _float64[0] = t as any;
        return hash2(_int32[0], _int32[1]);
    }

    export function toString(t: IntTuple) {
        _float64[0] = t as any;
        return `(${_int32[0]}, ${_int32[1]})`;
    }
}

export default IntTuple;