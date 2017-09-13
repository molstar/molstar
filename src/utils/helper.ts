/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * from https://github.com/dsehnal/CIFTools.js
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { TokenizerState } from './tokenizer-state'

/**
 * Eat everything until a newline occurs.
 */
export function eatLine(state: TokenizerState) {
    while (state.position < state.length) {
        switch (state.data.charCodeAt(state.position)) {
            case 10: // \n
                state.currentTokenEnd = state.position
                ++state.position
                ++state.currentLineNumber
                return
            case 13: // \r
                state.currentTokenEnd = state.position
                ++state.position
                ++state.currentLineNumber
                if (state.data.charCodeAt(state.position) === 10) {
                    ++state.position
                }
                return
            default:
                ++state.position
        }
    }
    state.currentTokenEnd = state.position;
}

/**
 * Eat everything until a whitespace/newline occurs.
 */
export function eatValue(state: TokenizerState) {
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
export function skipWhitespace(state: TokenizerState): number {
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
