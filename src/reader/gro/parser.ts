/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State as TokenizerState, Tokens, eatLine, skipWhitespace, eatValue, trim } from '../common/text/tokenizer'
import { parseInt } from '../common/text/number-parser'
import { createTokenFields } from '../common/text/token-field'
import * as Data from '../../data/data'
import Result from '../result'

interface StateInfo {
    numberOfAtoms: number
    hasVelocities: boolean
    numberOfDecimalPlaces: number
}

type State = TokenizerState<StateInfo>

function createState(data: string): State {
    return TokenizerState(data, { numberOfAtoms: 0, hasVelocities: false, numberOfDecimalPlaces: 3 });
}

/**
 * title string (free format string, optional time in ps after 't=')
 */
function handleTitleString(state: State, tokens: Tokens) {
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

/**
 * number of atoms (free format integer)
 */
function handleNumberOfAtoms(state: State, tokens: Tokens) {
    skipWhitespace(state)
    state.currentTokenStart = state.position
    eatValue(state)
    state.info.numberOfAtoms = parseInt(state.data, state.currentTokenStart, state.currentTokenEnd)
    Tokens.add(tokens, state.currentTokenStart, state.currentTokenEnd)
    eatLine(state)
}

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
function handleAtoms(state: State) {
    const fieldSizes = [ 5, 5, 5, 5, 8, 8, 8, 8, 8, 8 ];
    const fields = [ 'residueNumber', 'residueName', 'atomName', 'atomNumber', 'x', 'y', 'z' ]
    if (state.info.hasVelocities) {
        fields.push('vx', 'vy', 'vz')
    }

    const fieldCount = fields.length
    const tokens = Tokens.create(state.info.numberOfAtoms * 2 * fieldCount)

    let start: number;
    let end: number;

    for (let i = 0, _i = state.info.numberOfAtoms; i < _i; ++i) {
        state.currentTokenStart = state.position;
        end = state.currentTokenStart;
        for (let j = 0; j < fieldCount; ++j) {
            start = end;
            end = start + fieldSizes[j];

            trim(state, start, end);
            Tokens.addUnchecked(tokens, state.currentTokenStart, state.currentTokenEnd);
        }
        eatLine(state)
    }

    return Data.Category(state.info.numberOfAtoms, createTokenFields(state.data, fields, tokens));
}

/**
 * box vectors (free format, space separated reals), values:
 * v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y),
 * the last 6 values may be omitted (they will be set to zero).
 * Gromacs only supports boxes with v1(y)=v1(z)=v2(z)=0.
 */
function handleBoxVectors(state: State, tokens: Tokens) {
    // just read the first three values, ignore any remaining
    for (let i = 0; i < 3; ++i) {
        skipWhitespace(state);
        state.currentTokenStart = state.position;
        eatValue(state);
        Tokens.add(tokens, state.currentTokenStart, state.currentTokenEnd);
    }
}

function parseInternal(data: string): Result<Data.File> {
    const state = createState(data);

    const headerFields = ['title', 'timeInPs', 'numberOfAtoms', 'boxX', 'boxY', 'boxZ'];
    const headerTokens = Tokens.create(2 * headerFields.length);

    handleTitleString(state, headerTokens);
    handleNumberOfAtoms(state, headerTokens);
    const atoms = handleAtoms(state);
    handleBoxVectors(state, headerTokens);

    const block = Data.Block({
        header: Data.Category(1, createTokenFields(data, headerFields, headerTokens)),
        atoms
    });

    return Result.success(Data.File([block]));
}

export default function parse(data: string) {
    return parseInternal(data);
}