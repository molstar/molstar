/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * from https://github.com/dsehnal/CIFTools.js
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * A helper for building a typed array of token indices.
 */
export interface Tokens {
    indicesLenMinus2: number,
    count: number,
    indices: Int32Array
}

export namespace Tokens {
    function resize(tokens: Tokens) {
        // scale the size using golden ratio, because why not.
        const newBuffer = new Int32Array((1.61 * tokens.indices.length) | 0);
        newBuffer.set(tokens.indices);
        tokens.indices = newBuffer;
        tokens.indicesLenMinus2 = (newBuffer.length - 2) | 0;
    }

    export function add(tokens: Tokens, start: number, end: number) {
        if (tokens.count > tokens.indicesLenMinus2) {
            resize(tokens);
        }
        tokens.indices[tokens.count++] = start;
        tokens.indices[tokens.count++] = end;
    }

    export function addUnchecked(tokens: Tokens, start: number, end: number) {
        tokens.indices[tokens.count++] = start;
        tokens.indices[tokens.count++] = end;
    }

    export function create(size: number): Tokens {
        return {
            indicesLenMinus2: (size - 2) | 0,
            count: 0,
            indices: new Int32Array(size)
        }
    }
}
