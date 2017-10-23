/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { hash2 } from './hash-functions'

interface IntPair { fst: number, snd: number }

namespace IntPair {
    export interface Packed extends Number { }

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

    export function is(x: any): x is IntPair {
        return !!x && typeof x.fst === 'number' && typeof x.snd === 'number';
    }

    export function create(fst: number, snd: number) { return { fst, snd }; }
    export function zero(): IntPair { return { fst: 0, snd: 0 }; }

    export function pack1(fst: number, snd: number): number {
        _int32[0] = fst;
        _int32[1] = snd;
        return _float64[0];
    }

    export function pack(p: IntPair): number {
        _int32[0] = p.fst;
        _int32[1] = p.snd;
        return _float64[0];
    }

    export function unpack(packed: Packed, target: IntPair): IntPair {
        _float64[0] = packed as number;
        target.fst = _int32[0];
        target.snd = _int32[1];
        return target;
    }

    export function unpack1(packed: Packed): IntPair {
        return unpack(packed, zero());
    }

    export function fst(packed: Packed): number {
        _float64[0] = packed as number;
        return _int32[0];
    }

    export function snd(packed: Packed): number {
        _float64[0] = packed as number;
        return _int32[1];
    }

    export function areEqual(a: Packed, b: Packed) {
        _float64[0] = a as number;
        _float64_1[0] = b as number;
        return _int32[0] === _int32_1[0] && _int32[1] === _int32_1[1];
    }

    export function compare(a: number, b: number) {
        _float64[0] = a;
        _float64_1[0] = b;
        const x = _int32[0] - _int32_1[0];
        if (x !== 0) return x;
        return _int32[1] - _int32_1[1];
    }

    export function compareInArray(xs: ArrayLike<number>, i: number, j: number) {
        _float64[0] = xs[i];
        _float64_1[0] = xs[j];
        const x = _int32[0] - _int32_1[0];
        if (x !== 0) return x;
        return _int32[1] - _int32_1[1];
    }

    export function packedHashCode(packed: Packed) {
        _float64[0] = packed as number;
        return hash2(_int32[0], _int32[1]);
    }
}

export default IntPair