/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { hash2 } from './hash-functions'

/**
 * Represents a pair of two integers as a double,
 * Caution: === does not work, because of NaN, use IntTuple.areEqual for equality
 */
interface IntTuple { '@type': 'int-tuple' }

namespace IntTuple {
    export interface Unpacked { fst: number, snd: number }

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

    export function is(x: any): x is Unpacked {
        return !!x && typeof x.fst === 'number' && typeof x.snd === 'number';
    }

    export function create(fst: number, snd: number) { return { fst, snd }; }
    export function zero(): Unpacked { return { fst: 0, snd: 0 }; }

    export function pack(fst: number, snd: number): IntTuple {
        _int32[0] = fst;
        _int32[1] = snd;
        return _float64[0] as any;
    }

    export function pack1(t: Unpacked): IntTuple {
        _int32[0] = t.fst;
        _int32[1] = t.snd;
        return _float64[0] as any;
    }

    export function unpack(t: IntTuple, target: Unpacked): Unpacked {
        _float64[0] = t as any;
        target.fst = _int32[0];
        target.snd = _int32[1];
        return target;
    }

    export function unpack1(packed: IntTuple): Unpacked {
        return unpack(packed, zero());
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
}

export default IntTuple