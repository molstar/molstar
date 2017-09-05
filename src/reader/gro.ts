/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { parseInt } from '../utils/number-parser'
import { eatLine, eatValue, skipWhitespace } from '../utils/helper'
import { Tokens } from '../utils/tokens'
import { TokenizerState } from '../utils/tokenizer-state'

import { TextTable } from '../relational/text-table'

import { ParserResult } from '../parser'

/**
 * http://manual.gromacs.org/current/online/gro.html
 */

export interface GroFile {
    data: string;
    blocks: GroBlock[];
}

export interface GroBlock {
    getTable(name: string): TextTable
    addTable(table: TextTable): void
}

export class GroFile implements GroFile {
    data: string;
    blocks: GroBlock[] = [];

    constructor(data: string) {
        this.data = data;
    }
}

export class GroBlock implements GroBlock {
    private tableMap: Map<string, TextTable>;
    private tableList: TextTable[];

    data: string;

    /**
     * Gets a table by its name.
     */
    getTable(name: string) {
        return this.tableMap.get(name);
    }

    /**
     * Adds a table.
     */
    addTable(table: TextTable) {
        this.tableList[this.tableList.length] = table;
        this.tableMap.set(table.name, table);
    }

    constructor(data: string) {
        this.data = data;

        this.tableMap = new Map()
        this.tableList = []
    }
}

export interface GroState extends TokenizerState {
    numberOfAtoms: number
    hasVelocities: boolean
    numberOfDecimalPlaces: number
}

export function createTokenizer(data: string): GroState {
    return {
        data,

        position: 0,
        length: data.length,

        currentLineNumber: 1,
        currentTokenStart: 0,
        currentTokenEnd: 0,

        numberOfAtoms: 0,
        hasVelocities: false,
        numberOfDecimalPlaces: 3
    };
}

/**
 * title string (free format string, optional time in ps after 't=')
 */
function handleTitleString (state: GroState, tokens: Tokens) {
    eatLine(state)
    // console.log('title', state.data.substring(state.currentTokenStart, state.currentTokenEnd))
    let start = state.currentTokenStart
    let end = state.currentTokenEnd
    let valueStart = state.currentTokenStart
    let valueEnd = start

    while (valueEnd < end && !isTime(state.data, valueEnd)) ++valueEnd;

    if (isTime(state.data, valueEnd)) {
        let timeStart = valueEnd + 2

        while (valueEnd > start && isSpaceOrComma(state.data, valueEnd - 1)) --valueEnd;
        Tokens.add(tokens, valueStart, valueEnd)  // title

        while (timeStart < end && state.data.charCodeAt(timeStart) === 32) ++timeStart;
        while (valueEnd > timeStart && state.data.charCodeAt(valueEnd - 1) === 32) --valueEnd;
        Tokens.add(tokens, timeStart, end)  // time
    } else {
        Tokens.add(tokens, valueStart, valueEnd)  // title
        Tokens.add(tokens, valueEnd, valueEnd)  // empty token for time
    }
}

function isSpaceOrComma(data: string, position: number): boolean {
    const c = data.charCodeAt(position);
    return c === 32 || c === 44
}

function isTime(data: string, position: number): boolean {
    // T/t
    const c = data.charCodeAt(position);
    if (c !== 84 && c !== 116) return false;
    // =
    if (data.charCodeAt(position + 1) !== 61) return false;

    return true;
}

// function isDot(state: TokenizerState): boolean {
//     // .
//     if (state.data.charCodeAt(state.currentTokenStart) !== 46) return false;

//     return true;
// }

// function numberOfDecimalPlaces (state: TokenizerState) {
//     // var ndec = firstLines[ 2 ].length - firstLines[ 2 ].lastIndexOf('.') - 1
//     const start = state.currentTokenStart
//     const end = state.currentTokenEnd
//     for (let i = end; start < i; --i) {
//         // .
//         if (state.data.charCodeAt(i) === 46) return end - start - i
//     }
//     throw new Error('Could not determine number of decimal places')
// }

/**
 * number of atoms (free format integer)
 */
function handleNumberOfAtoms (state: GroState, tokens: Tokens) {
    skipWhitespace(state)
    state.currentTokenStart = state.position
    eatValue(state)
    state.numberOfAtoms = parseInt(state.data, state.currentTokenStart, state.currentTokenEnd)
    Tokens.add(tokens, state.currentTokenStart, state.currentTokenEnd)
    eatLine(state)
}

// function checkForVelocities (state: GroState) {

// }

/**
 * This format is fixed, ie. all columns are in a fixed position.
 * Optionally (for now only yet with trjconv) you can write gro files
 * with any number of decimal places, the format will then be n+5
 * positions with n decimal places (n+1 for velocities) in stead
 * of 8 with 3 (with 4 for velocities). Upon reading, the precision
 * will be inferred from the distance between the decimal points
 * (which will be n+5). Columns contain the following information
 * (from left to right):
 *     residue number (5 positions, integer)
 *     residue name (5 characters)
 *     atom name (5 characters)
 *     atom number (5 positions, integer)
 *     position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
 *     velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions with 4 decimal places)
 */
function handleAtoms (state: GroState, block: GroBlock) {
    const name = 'atoms'

    const columns = [ 'residueNumber', 'residueName', 'atomName', 'atomNumber', 'x', 'y', 'z' ]
    if (state.hasVelocities) {
        columns.push('vx', 'vy', 'vz')
    }
    const fieldSizes = [ 5, 5, 5, 5, 8, 8, 8, 8, 8, 8 ]

    const columnCount = columns.length
    const tokens = Tokens.create(state.numberOfAtoms * 2 * columnCount)

    let start: number
    let end: number
    let valueStart: number
    let valueEnd: number = state.position

    for (let i = 0; i < state.numberOfAtoms; ++i) {
        state.currentTokenStart = state.position
        end = state.currentTokenStart
        for (let j = 0; j < columnCount; ++j) {
            start = end
            end = start + fieldSizes[j]

            // trim
            valueStart = start
            valueEnd = end
            while (valueStart < valueEnd && state.data.charCodeAt(valueStart) === 32) ++valueStart;
            while (valueEnd > valueStart && state.data.charCodeAt(valueEnd - 1) === 32) --valueEnd;

            Tokens.addUnchecked(tokens, valueStart, valueEnd)
        }
        state.position = valueEnd
        eatLine(state)
    }

    block.addTable(new TextTable(state.data, name, columns, tokens));
}

/**
 * box vectors (free format, space separated reals), values:
 * v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y),
 * the last 6 values may be omitted (they will be set to zero).
 * Gromacs only supports boxes with v1(y)=v1(z)=v2(z)=0.
 */
function handleBoxVectors (state: GroState, tokens: Tokens) {
    // just read the first three values, ignore any remaining
    for (let i = 0; i < 3; ++i) {
        skipWhitespace(state)
        state.currentTokenStart = state.position
        eatValue(state)
        Tokens.add(tokens, state.currentTokenStart, state.currentTokenEnd)
    }
}

/**
 * Creates an error result.
 */
// function error(line: number, message: string) {
//     return ParserResult.error<GroFile>(message, line);
// }

/**
 * Creates a data result.
 */
function result(data: GroFile) {
    return ParserResult.success(data);
}

function parseInternal(data: string): ParserResult<GroFile> {
    const state = createTokenizer(data)
    const file = new GroFile(data)

    let block = new GroBlock(data)
    file.blocks.push(block)

    const headerColumns = ['title', 'timeInPs', 'numberOfAtoms', 'boxX', 'boxY', 'boxZ']
    const headerTokens = Tokens.create(2 * headerColumns.length)
    let header = new TextTable(state.data, 'header', headerColumns, headerTokens)
    block.addTable(header)

    handleTitleString(state, headerTokens)
    handleNumberOfAtoms(state, headerTokens)
    handleAtoms(state, block)
    handleBoxVectors(state, headerTokens)

    return result(file);
}

export function parse(data: string) {
    return parseInternal(data);
}
