/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Schäfer, Marco <marco.schaefer@uni-tuebingen.de>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Tokens, TokenBuilder, Tokenizer } from '../common/text/tokenizer'
import * as Data from './schema'
import{ ReaderResult } from '../result'
import {Task, RuntimeContext, chunkedSubtask } from 'mol-task'
import { parseInt as fastParseInt, parseFloat as fastParseFloat } from '../common/text/number-parser'

const enum PlyTokenType {
    Value = 0,
    Comment = 1,
    End = 2,
    property = 3,
    element = 4
}

interface State {
    data: string;
    tokenizer: Tokenizer,

    tokenType: PlyTokenType;
    runtimeCtx: RuntimeContext,
    tokens: Tokens[],

    fieldCount: number,

    columnCount: number,
    propertyCount: number,
    vertexCount: number,
    currentVertex: number,
    currentProperty: number,
    currentFace: number,
    currentFaceElement: number,
    faceCount: number,
    endHeader: number,

    initialHead: string[],
    properties: number[],
    vertices: number[],
    colors: number[],
    normals: number[],
    faces: number[],
    propertyNames: string[],
    check: string[],

    commentCharCode: number,
    propertyCharCode: number,
    elementCharCode: number
}

function State(data: string, runtimeCtx: RuntimeContext, opts: PlyOptions): State {
    const tokenizer = Tokenizer(data)
    return {
        data,
        tokenizer,

        tokenType: PlyTokenType.End,
        runtimeCtx,
        tokens: [],

        fieldCount: 0,

        columnCount: 0,
        propertyCount: 0,
        vertexCount: 0,
        currentVertex: 0,
        currentProperty: 0,
        currentFace: 0,
        currentFaceElement: 0,
        faceCount: 0,
        endHeader: 0,

        initialHead: [],
        properties: [],
        vertices: [],
        colors: [],
        normals: [],
        faces: [],
        propertyNames: [],
        check: [],

        commentCharCode: opts.comment.charCodeAt(0),
        propertyCharCode: opts.property.charCodeAt(0),
        elementCharCode: opts.element.charCodeAt(0)
    };
}

/**
 * Eat everything until a delimiter (whitespace) or newline occurs.
 * Ignores whitespace at the end of the value, i.e. trim right.
 * Returns true when a newline occurs after the value.
 */
function eatValue(state: Tokenizer) {
    while (state.position < state.length) {
        const c = state.data.charCodeAt(state.position);
        ++state.position
        switch (c) {
            case 10:  // \n
            case 13:  // \r
                return true;
            case 32: // ' ' Delimeter of ply is space (Unicode 32)
                return true;
            case 9:  // \t
            case 32:  // ' '
                break;
            default:
                ++state.tokenEnd;
                break;
        }
    }
}

function eatLine (state: Tokenizer) {
    while (state.position < state.length) {
        const c = state.data.charCodeAt(state.position);
        ++state.position
        switch (c) {
            case 10:  // \n
            case 13:  // \r
                return true;
            case 9:  // \t
                break;
            default:
                ++state.tokenEnd;
                break;
        }
    }

}

function skipLine(state: Tokenizer) {
    while (state.position < state.length) {
        const c = state.data.charCodeAt(state.position);
        if (c === 10 || c === 13) return  // \n or \r
        ++state.position
    }
}

function getColumns(state: State, numberOfColumns: number) {
    eatLine(state.tokenizer);
    let tmp = Tokenizer.getTokenString(state.tokenizer)
    let split = tmp.split(' ', numberOfColumns);
    return split;
}

/**
 * Move to the next token.
 * Returns true when the current char is a newline, i.e. indicating a full record.
 */
function moveNextInternal(state: State) {
    const tokenizer = state.tokenizer

    if (tokenizer.position >= tokenizer.length) {
        state.tokenType = PlyTokenType.End;
        return true;
    }

    tokenizer.tokenStart = tokenizer.position;
    tokenizer.tokenEnd = tokenizer.position;
    const c = state.data.charCodeAt(tokenizer.position);
    switch (c) {
        case state.commentCharCode:
            state.tokenType = PlyTokenType.Comment;
            skipLine(tokenizer);
            break;
        case state.propertyCharCode: // checks all line beginning with 'p'
            state.check = getColumns(state, 3);
            if (state.check[0] !== 'ply' && state.faceCount === 0) {
                state.propertyNames.push(state.check[1]);
                state.propertyNames.push(state.check[2]);
                state.propertyCount++;
            }
            return;
        case state.elementCharCode: // checks all line beginning with 'e'
            state.check = getColumns(state, 3);
            if (state.check[1] === 'vertex') state.vertexCount= Number(state.check[2]);
            if (state.check[1] === 'face') state.faceCount = Number(state.check[2]);
            if (state.check[0] === 'end_header') state.endHeader = 1;
            return;
        default:                    // for all the other lines
            state.tokenType = PlyTokenType.Value;
            let return_value = eatValue(tokenizer);

            if (state.endHeader === 1) {
                if (state.currentVertex < state.vertexCount) {
                    // TODO the numbers are parsed twice
                    state.properties[state.currentVertex * state.propertyCount + state.currentProperty] = Number(Tokenizer.getTokenString(state.tokenizer));
                    if (state.currentProperty < 3) {
                        state.vertices[state.currentVertex * 3 + state.currentProperty] = fastParseFloat(state.tokenizer.data, state.tokenizer.tokenStart, state.tokenizer.tokenEnd);
                    }
                    if (state.currentProperty >= 3 && state.currentProperty < 6) {
                        state.colors[state.currentVertex * 3 + state.currentProperty - 3] = fastParseInt(state.tokenizer.data, state.tokenizer.tokenStart, state.tokenizer.tokenEnd);
                    }
                    if (state.currentProperty >= 7 && state.currentProperty < 10) {
                        state.normals[state.currentVertex * 3 + state.currentProperty - 7] = fastParseFloat(state.tokenizer.data, state.tokenizer.tokenStart, state.tokenizer.tokenEnd);
                    }
                    state.currentProperty++;
                    if (state.currentProperty === state.propertyCount) {
                        state.currentProperty = 0;
                        state.currentVertex++;
                    }
                    return return_value;
                }
                if (state.currentFace < state.faceCount && state.currentVertex === state.vertexCount) {
                    state.faces[state.currentFace * 4 + state.currentFaceElement] = fastParseInt(state.tokenizer.data, state.tokenizer.tokenStart, state.tokenizer.tokenEnd);
                    state.currentFaceElement++;
                    if (state.currentProperty === 4) {
                        state.currentFaceElement = 0;
                        state.currentFace++;
                    }
                }
            }
            return return_value;
    }
}

/**
 * Moves to the next non-comment token/line.
 * Returns true when the current char is a newline, i.e. indicating a full record.
 */
function moveNext(state: State) {
    let newRecord = moveNextInternal(state);
    while (state.tokenType === PlyTokenType.Comment) { // skip comment lines (marco)
        newRecord = moveNextInternal(state);
    }
    return newRecord
}

function readRecordsChunk(chunkSize: number, state: State) {
    if (state.tokenType === PlyTokenType.End) return 0

    moveNext(state);
    const { tokens, tokenizer } = state;
    let counter = 0;
    while (state.tokenType === PlyTokenType.Value && counter < chunkSize) {
        TokenBuilder.add(tokens[state.fieldCount % state.columnCount], tokenizer.tokenStart, tokenizer.tokenEnd);
        ++state.fieldCount
        moveNext(state);
        ++counter;
    }
    return counter;
}

function readRecordsChunks(state: State) {
    return chunkedSubtask(state.runtimeCtx, 100000, state, readRecordsChunk,
        (ctx, state) => ctx.update({ message: 'Parsing...', current: state.tokenizer.position, max: state.data.length }));
}

function addHeadEntry (state: State) {
    const head = Tokenizer.getTokenString(state.tokenizer)
    state.initialHead.push(head)
    state.tokens.push(TokenBuilder.create(head, state.data.length / 80))
}

function init(state: State) { // only for first two lines to get the format and the coding! (marco)
    let newRecord = moveNext(state)
    while (!newRecord) {  // newRecord is only true when a newline occurs (marco)
        addHeadEntry(state)
        newRecord = moveNext(state);
    }
    addHeadEntry(state)
    newRecord = moveNext(state);
    while (!newRecord) {
        addHeadEntry(state)
        newRecord = moveNext(state);
    }
    addHeadEntry(state)
    if (state.initialHead[0] !== 'ply') {
        console.log('ERROR: this is not a .ply file!')
        throw new Error('this is not a .ply file!');
        return 0;
    }
    if (state.initialHead[2] !== 'ascii') {
        console.log('ERROR: only ASCII-DECODING is supported!');
        throw new Error('only ASCII-DECODING is supported!');
        return 0;
    }
    state.columnCount = state.initialHead.length
    return 1;
}

async function handleRecords(state: State): Promise<Data.PlyData> {
    if (!init(state)) {
        console.log('ERROR: parsing file (PLY) failed!')
        throw new Error('arsing file (PLY) failed!');
    }
    await readRecordsChunks(state)

    return Data.PlyData(state.vertexCount, state.faceCount, state.propertyCount, state.initialHead, state.propertyNames, state.properties, state.vertices, state.colors, state.normals, state.faces)
}

async function parseInternal(data: string, ctx: RuntimeContext, opts: PlyOptions): Promise<ReaderResult<Data.PlyFile>> {
    const state = State(data, ctx, opts);

    ctx.update({ message: 'Parsing...', current: 0, max: data.length });
    const PLYdata = await handleRecords(state)
    const result = Data.PlyFile(PLYdata)

    return ReaderResult.success(result);
}

interface PlyOptions {
    comment: string;
    property: string;
    element: string;
}

export function parse(data: string, opts?: Partial<PlyOptions>) {
    const completeOpts = Object.assign({}, { comment: 'c', property: 'p', element: 'e' }, opts)
    return Task.create<ReaderResult<Data.PlyFile>>('Parse PLY', async ctx => {
        return await parseInternal(data, ctx, completeOpts);
    });
}

export default parse;