/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * mostly from https://github.com/dsehnal/CIFTools.js
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { chunkedSubtask, RuntimeContext } from '../../../../mol-task';

export { Tokenizer };

interface Tokenizer {
    data: string,

    position: number,
    length: number,

    lineNumber: number,
    tokenStart: number,
    tokenEnd: number
}

export interface Tokens {
    data: string,
    count: number,
    indices: ArrayLike<number>
}

function Tokenizer(data: string): Tokenizer {
    return {
        data,
        position: 0,
        length: data.length,
        lineNumber: 1,
        tokenStart: 0,
        tokenEnd: 0
    };
}

namespace Tokenizer {
    export function getTokenString(state: Tokenizer) {
        return state.data.substring(state.tokenStart, state.tokenEnd);
    }

    /** Resets the state */
    export function reset (state: Tokenizer) {
        state.position = 0;
        state.lineNumber = 1;
        state.tokenStart = 0;
        state.tokenEnd = 0;
    }

    /**
     * Eat everything until a newline occurs.
     */
    export function eatLine(state: Tokenizer): boolean {
        const { data } = state;
        while (state.position < state.length) {
            switch (data.charCodeAt(state.position)) {
                case 10: // \n
                    state.tokenEnd = state.position;
                    ++state.position;
                    ++state.lineNumber;
                    return true;
                case 13: // \r
                    state.tokenEnd = state.position;
                    ++state.position;
                    ++state.lineNumber;
                    if (data.charCodeAt(state.position) === 10) {
                        ++state.position;
                    }
                    return true;
                default:
                    ++state.position;
                    break;
            }
        }
        state.tokenEnd = state.position;
        return state.tokenStart !== state.tokenEnd;
    }

    /** Sets the current token start to the current position */
    export function markStart(state: Tokenizer) {
        state.tokenStart = state.position;
    }

    /** Sets the current token start to current position and moves to the next line. */
    export function markLine(state: Tokenizer) {
        state.tokenStart = state.position;
        return eatLine(state);
    }

    /** Advance the state by the given number of lines and return line as string. */
    export function readLine(state: Tokenizer): string {
        markLine(state);
        return getTokenString(state);
    }

    function readLinesChunk(state: Tokenizer, count: number, tokens: Tokens) {
        let read = 0;
        for (let i = 0; i < count; i++) {
            if (!markLine(state)) return read;
            TokenBuilder.addUnchecked(tokens, state.tokenStart, state.tokenEnd);
            read++;
        }
        return read;
    }

    /** Advance the state by the given number of lines and return them*/
    export function markLines(state: Tokenizer, count: number): Tokens {
        const lineTokens = TokenBuilder.create(state.data, count * 2);
        readLinesChunk(state, count, lineTokens);
        return lineTokens;
    }

    /** Advance the state by the given number of lines and return them */
    export function readLines(state: Tokenizer, count: number): string[] {
        const ret: string[] = [];
        for (let i = 0; i < count; i++) {
            ret.push(Tokenizer.readLine(state));
        }
        return ret;
    }

    /** Advance the state by the given number of lines and return line starts/ends as tokens. */
    export async function readLinesAsync(state: Tokenizer, count: number, ctx: RuntimeContext, initialLineCount = 100000): Promise<Tokens> {
        const { length } = state;
        const lineTokens = TokenBuilder.create(state.data, count * 2);

        let linesAlreadyRead = 0;
        await chunkedSubtask(ctx, initialLineCount, state, (chunkSize, state) => {
            const linesToRead = Math.min(count - linesAlreadyRead, chunkSize);
            readLinesChunk(state, linesToRead, lineTokens);
            linesAlreadyRead += linesToRead;
            return linesToRead;
        }, (ctx, state) => ctx.update({ message: 'Parsing...', current: state.position, max: length }));

        return lineTokens;
    }

    export function readAllLines(data: string) {
        const state = Tokenizer(data);
        const tokens = TokenBuilder.create(state.data, Math.max(data.length / 80, 2));
        while (markLine(state)) {
            TokenBuilder.add(tokens, state.tokenStart, state.tokenEnd);
        }
        return tokens;
    }

    function readLinesChunkChecked(state: Tokenizer, count: number, tokens: Tokens) {
        let read = 0;
        for (let i = 0; i < count; i++) {
            if (!markLine(state)) return read;
            TokenBuilder.add(tokens, state.tokenStart, state.tokenEnd);
            read++;
        }
        return read;
    }

    export async function readAllLinesAsync(data: string, ctx: RuntimeContext, chunkSize = 100000) {
        const state = Tokenizer(data);
        const tokens = TokenBuilder.create(state.data, Math.max(data.length / 80, 2));

        await chunkedSubtask(ctx, chunkSize, state, (chunkSize, state) => {
            readLinesChunkChecked(state, chunkSize, tokens);
            return state.position < state.length ? chunkSize : 0;
        }, (ctx, state) => ctx.update({ message: 'Parsing...', current: state.position, max: length }));

        return tokens;
    }

    /**
     * Eat everything until a whitespace/newline occurs.
     */
    export function eatValue(state: Tokenizer) {
        while (state.position < state.length) {
            switch (state.data.charCodeAt(state.position)) {
                case 9:  // \t
                case 10: // \n
                case 13: // \r
                case 32: // ' '
                    state.tokenEnd = state.position;
                    return;
                default:
                    ++state.position;
                    break;
            }
        }
        state.tokenEnd = state.position;
    }

    /**
     * Skips all the whitespace - space, tab, newline, CR
     * Handles incrementing line count.
     */
    export function skipWhitespace(state: Tokenizer): number {
        let prev = -1;
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
                        ++state.lineNumber;
                    }
                    prev = c;
                    ++state.position;
                    break;
                case 13: // \r
                    prev = c;
                    ++state.position;
                    ++state.lineNumber;
                    break;
                default:
                    return prev;
            }
        }
        return prev;
    }

    /** Trims spaces and tabs */
    export function trim(state: Tokenizer, start: number, end: number) {
        const { data } = state;
        let s = start, e = end - 1;

        let c = data.charCodeAt(s);
        while ((c === 9 || c === 32) && s <= e) c = data.charCodeAt(++s);
        c = data.charCodeAt(e);
        while ((c === 9 || c === 32) && e >= s) c = data.charCodeAt(--e);

        state.tokenStart = s;
        state.tokenEnd = e + 1;
        state.position = end;
        return state;
    }
}

export function trimStr(data: string, start: number, end: number) {
    let s = start, e = end - 1;
    let c = data.charCodeAt(s);
    while ((c === 9 || c === 32) && s <= e) c = data.charCodeAt(++s);
    c = data.charCodeAt(e);
    while ((c === 9 || c === 32) && e >= s) c = data.charCodeAt(--e);
    return data.substring(s, e + 1);
}

export namespace TokenBuilder {
    interface Builder extends Tokens {
        offset: number,
        indices: Uint32Array,
        indicesLenMinus2: number
    }

    function resize(builder: Builder) {
        // scale the size using golden ratio, because why not.
        const newBuffer = new Uint32Array((1.61 * builder.indices.length) | 0);
        newBuffer.set(builder.indices);
        builder.indices = newBuffer;
        builder.indicesLenMinus2 = (newBuffer.length - 2) | 0;
    }

    export function add(tokens: Tokens, start: number, end: number) {
        const builder = tokens as Builder;
        if (builder.offset > builder.indicesLenMinus2) {
            resize(builder);
        }
        builder.indices[builder.offset++] = start;
        builder.indices[builder.offset++] = end;
        tokens.count++;
    }

    export function addToken(tokens: Tokens, tokenizer: Tokenizer) {
        add(tokens, tokenizer.tokenStart, tokenizer.tokenEnd);
    }

    export function addUnchecked(tokens: Tokens, start: number, end: number) {
        (tokens as Builder).indices[(tokens as Builder).offset++] = start;
        (tokens as Builder).indices[(tokens as Builder).offset++] = end;
        tokens.count++;
    }

    export function create(data: string, size: number): Tokens {
        size = Math.max(10, size);
        return <Builder>{
            data,
            indicesLenMinus2: (size - 2) | 0,
            count: 0,
            offset: 0,
            indices: new Uint32Array(size)
        };
    }
}