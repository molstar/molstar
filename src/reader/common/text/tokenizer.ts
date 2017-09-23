/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * mostly from https://github.com/dsehnal/CIFTools.js
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export interface State<TokenType = any> {
    data: string

    position: number
    length: number

    currentLineNumber: number
    currentTokenStart: number
    currentTokenEnd: number

    currentTokenType: TokenType
}

export function State<TokenType>(data: string, initialTokenType?: TokenType): State<TokenType> {
    return {
        data,
        position: 0,
        length: data.length,
        currentLineNumber: 1,
        currentTokenStart: 0,
        currentTokenEnd: 0,
        currentTokenType: initialTokenType!
    };
}

export function getTokenString(state: State) {
    return state.data.substring(state.currentTokenStart, state.currentTokenEnd);
}

/**
 * Eat everything until a newline occurs.
 */
export function eatLine(state: State) {
    const { data } = state;
    while (state.position < state.length) {
        switch (data.charCodeAt(state.position)) {
            case 10: // \n
                state.currentTokenEnd = state.position;
                ++state.position;
                ++state.currentLineNumber;
                return;
            case 13: // \r
                state.currentTokenEnd = state.position;
                ++state.position;
                ++state.currentLineNumber;
                if (data.charCodeAt(state.position) === 10) {
                    ++state.position;
                }
                return;
            default:
                ++state.position;
                break;
        }
    }
    state.currentTokenEnd = state.position;
}

/** Sets the current token start to the current position */
export function markStart(state: State) {
    state.currentTokenStart = state.position;
}

/** Sets the current token start to current position and moves to the next line. */
export function markLine(state: State) {
    state.currentTokenStart = state.position;
    eatLine(state);
}

/**
 * Eat everything until a whitespace/newline occurs.
 */
export function eatValue(state: State) {
    while (state.position < state.length) {
        switch (state.data.charCodeAt(state.position)) {
            case 9:  // \t
            case 10: // \n
            case 13: // \r
            case 32: // ' '
                state.currentTokenEnd = state.position;
                return;
            default:
                ++state.position;
                break;
        }
    }
    state.currentTokenEnd = state.position;
}

/**
 * Skips all the whitespace - space, tab, newline, CR
 * Handles incrementing line count.
 */
export function skipWhitespace(state: State): number {
    let prev = 10;
    while (state.position < state.length) {
        let c = state.data.charCodeAt(state.position);
        switch (c) {
            case 9: // '\t'
            case 32: // ' '
                prev = c;
                ++state.position;
                break;
            case 10: // \n
                // handle \r\n
                if (prev !== 13) {
                    ++state.currentLineNumber;
                }
                prev = c;
                ++state.position;
                break;
            case 13: // \r
                prev = c;
                ++state.position;
                ++state.currentLineNumber;
                break;
            default:
                return prev;
        }
    }
    return prev;
}

/** Trims spaces and tabs */
export function trim(state: State, start: number, end: number) {
    const { data } = state;
    let s = start, e = end - 1;

    let c = data.charCodeAt(s);
    while ((c === 9 || c === 32) && s <= e) c = data.charCodeAt(++s);
    c = data.charCodeAt(e);
    while ((c === 9 || c === 32) && e >= s) c = data.charCodeAt(--e);

    state.currentTokenStart = s;
    state.currentTokenEnd = e + 1;
    state.position = end;
}

export function trimStr(data: string, start: number, end: number) {
    let s = start, e = end - 1;
    let c = data.charCodeAt(s);
    while ((c === 9 || c === 32) && s <= e) c = data.charCodeAt(++s);
    c = data.charCodeAt(e);
    while ((c === 9 || c === 32) && e >= s) c = data.charCodeAt(--e);
    return data.substring(s, e + 1);
}

export interface Tokens {
    indicesLenMinus2: number,
    count: number,
    indices: Uint32Array
}

export namespace Tokens {
    function resize(tokens: Tokens) {
        // scale the size using golden ratio, because why not.
        const newBuffer = new Uint32Array((1.61 * tokens.indices.length) | 0);
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
            indices: new Uint32Array(size)
        }
    }
}


/**
 * A helper for building a typed array of token indices.
 */
export default Tokens