/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * from https://github.com/dsehnal/CIFTools.js
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export interface TokenizerState {
    data: string

    position: number
    length: number

    currentLineNumber: number
    currentTokenStart: number
    currentTokenEnd: number

    currentTokenType?: number
}
