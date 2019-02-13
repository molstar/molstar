/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// import { Column } from 'mol-data/db'
import { Tokens, TokenBuilder, Tokenizer } from '../../common/text/tokenizer'
import * as Data from './data-model'
import Field from './field'
import Result from '../../result'
import { Task, RuntimeContext, chunkedSubtask, } from 'mol-task'


const enum PlyTokenType {
    Value = 0,
    Comment = 1,
    End = 2,
    property = 3
}

interface State {
    data: string;
    tokenizer: Tokenizer,

    tokenType: PlyTokenType;
    runtimeCtx: RuntimeContext,
    tokens: Tokens[],

    fieldCount: number,
    recordCount: number,

    columnCount: number,
    initialHead: string[],
    propertyNames: string[],

    commentCharCode: number,
    propertyCharCode: number
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
        recordCount: 0,

        columnCount: 0,
        initialHead: [],
        propertyNames: [],

        commentCharCode: opts.comment.charCodeAt(0),
        propertyCharCode: opts.property.charCodeAt(0)
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
                return;
            case 9:  // \t
            case 32:  // ' '
                break;
            default:
                ++state.tokenEnd;
                break;
        }
    }
}



function skipWhitespace(state: Tokenizer) {
    let prev = -1;
    while (state.position < state.length) {
        const c = state.data.charCodeAt(state.position);
        switch (c) {
            case 9:  // '\t'
            //case 32:  // ' '
                prev = c;
                ++state.position;
                break;
            case 10:  // \n
                // handle \r\n
                if (prev !== 13) {
                    ++state.lineNumber;
                }
                prev = c;
                ++state.position;
                break;
            case 13:  // \r
                prev = c;
                ++state.position;
                ++state.lineNumber;
                break;
            default:
                return;
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

/**
 * Move to the next token.
 * Returns true when the current char is a newline, i.e. indicating a full record.
 */
function moveNextInternal(state: State) {
    const tokenizer = state.tokenizer
    //skipWhitespace(tokenizer);

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
        case state.propertyCharCode:
            state.tokenType = PlyTokenType.property;
            //return eatProperty(tokenizer);
        default:
            state.tokenType = PlyTokenType.Value;
            return eatValue(tokenizer);
    }
}

/**
 * Moves to the next non-comment token/line.
 * Returns true when the current char is a newline, i.e. indicating a full record.
 */
function moveNext(state: State) {
    let newRecord = moveNextInternal(state);
    while (state.tokenType === PlyTokenType.Comment) { // skip comment lines (marco
        newRecord = moveNextInternal(state);
    }
    return newRecord
}

function readRecordsChunk(chunkSize: number, state: State) {
    if (state.tokenType === PlyTokenType.End) return 0

    let newRecord = moveNext(state);
    if (newRecord) ++state.recordCount

    const { tokens, tokenizer } = state;
    let counter = 0;
    while (state.tokenType === PlyTokenType.Value && counter < chunkSize) {
        TokenBuilder.add(tokens[state.fieldCount % state.columnCount], tokenizer.tokenStart, tokenizer.tokenEnd);
        ++state.fieldCount
        newRecord = moveNext(state);
        if (newRecord) ++state.recordCount
        ++counter;
    }
    return counter;
}

function readRecordsChunks(state: State) {
    return chunkedSubtask(state.runtimeCtx, 100000, state, readRecordsChunk,
        (ctx, state) => ctx.update({ message: 'Parsing...', current: state.tokenizer.position, max: state.data.length }));
}

function addColumn (state: State) {
    state.initialHead.push(Tokenizer.getTokenString(state.tokenizer))
    state.tokens.push(TokenBuilder.create(state.tokenizer, state.data.length / 80))
}

function init(state: State) { // only for first line to get the columns! (marco)
    let newRecord = moveNext(state)
    while (!newRecord) {  // newRecord is only true when a newline occurs (marco)
        addColumn(state)
        newRecord = moveNext(state);
    }
    addColumn(state)
    newRecord = moveNext(state);
    while (!newRecord) {
        addColumn(state)
        newRecord = moveNext(state);
    }
    addColumn(state)
    if(state.initialHead[0] !== 'ply'){
        console.log("ERROR: this is not a .ply file!")
        throw new Error("this is not a .ply file!");
        return 0;
    }
    if(state.initialHead[2] !== 'ascii'){
        console.log("ERROR: only ASCII-DECODING is supported!");
        throw new Error("only ASCII-DECODING is supported!");
        return 0;
    }
    state.columnCount = state.initialHead.length
    return 1;
}

async function handleRecords(state: State): Promise<Data.ply_form> {
    if(!init(state)){
        console.log("ERROR: parsing file (PLY) failed!")
        throw new Error("parsing file (PLY) failed!");
    }
    await readRecordsChunks(state)

    const columns: Data.CsvColumns = Object.create(null);
    for (let i = 0; i < state.columnCount; ++i) {
        columns[state.initialHead[i]] = Field(state.tokens[i], state.recordCount);
    }


    return Data.CsvTable(state.recordCount,0,0,0, state.initialHead, columns)
}

async function parseInternal(data: string, ctx: RuntimeContext, opts: PlyOptions): Promise<Result<Data.PlyFile>> {
    const state = State(data, ctx, opts);

    ctx.update({ message: 'Parsing...', current: 0, max: data.length });
    const table = await handleRecords(state)
    const result = Data.CsvFile(table)
    console.log(result);
    return Result.success(result);
}

interface PlyOptions {
    comment: string;
    property: string;
}

export function parse(data: string, opts?: Partial<PlyOptions>) {
    const completeOpts = Object.assign({}, { comment: 'c', property: 'p' }, opts)
    return Task.create<Result<Data.PlyFile>>('Parse PLY', async ctx => {
        return await parseInternal(data, ctx, completeOpts);
    });
}

export default parse;