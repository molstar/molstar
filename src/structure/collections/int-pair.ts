/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

interface IntPair { fst: number, snd: number }

namespace IntPair {
    const { _int32, _float64 } = (function() {
        const data = new ArrayBuffer(8);
        return { _int32: new Int32Array(data), _float64: new Float64Array(data) };
    }());

    export function zero(): IntPair { return { fst: 0, snd: 0 }; }

    export function set(fst: number, snd: number): number {
        _int32[0] = fst;
        _int32[1] = snd;
        return _float64[0];
    }

    export function get(packed: number, target: IntPair): IntPair {
        _float64[0] = packed;
        target.fst = _int32[0];
        target.snd = _int32[1];
        return target;
    }
}

export default IntPair