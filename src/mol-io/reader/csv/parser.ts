/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// import { Column } from 'mol-data/db'
import { Tokens, TokenBuilder, Tokenizer } from '../common/text/tokenizer';
import * as Data from './data-model';
import { Field } from './field';
import { ReaderResult as Result } from '../result';
import { Task, RuntimeContext, chunkedSubtask, } from '../../../mol-task';

const enum CsvTokenType {
    Value = 0,
    Comment = 1,
    End = 2
}

interface State {
    data: string;
    tokenizer: Tokenizer,

    tokenType: CsvTokenType;
    runtimeCtx: RuntimeContext,
    tokens: Tokens[],

    fieldCount: number,
    recordCount: number,

    columnCount: number,
    columnNames: string[],

    quoteCharCode: number,
    commentCharCode: number,
    delimiterCharCode: number,

    noColumnNamesRecord: boolean
}

function State(data: string, runtimeCtx: RuntimeContext, opts: CsvOptions): State {

    const tokenizer = Tokenizer(data);
    return {
        data,
        tokenizer,

        tokenType: CsvTokenType.End,
        runtimeCtx,
        tokens: [],

        fieldCount: 0,
        recordCount: 0,

        columnCount: 0,
        columnNames: [],

        quoteCharCode: opts.quote.charCodeAt(0),
        commentCharCode: opts.comment.charCodeAt(0),
        delimiterCharCode: opts.delimiter.charCodeAt(0),
        noColumnNamesRecord: opts.noColumnNames
    };
}

/**
 * Eat everything until a delimiter or newline occurs.
 * Ignores whitespace at the end of the value, i.e. trim right.
 * Returns true when a newline occurs after the value.
 */
function eatValue(state: Tokenizer, delimiterCharCode: number) {
    while (state.position < state.length) {
        const c = state.data.charCodeAt(state.position);
        ++state.position;
        switch (c) {
            case 10:  // \n
            case 13:  // \r
                return true;
            case delimiterCharCode:
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

/**
 * Eats a quoted value. Can contain a newline.
 * Returns true when a newline occurs after the quoted value.
 *
 * Embedded quotes are represented by a pair of double quotes:
 * - ""xx"" => "xx"
 */
function eatQuoted(state: Tokenizer, quoteCharCode: number, delimiterCharCode: number) {
    ++state.position;
    while (state.position < state.length) {
        const c = state.data.charCodeAt(state.position);
        if (c === quoteCharCode) {
            const next = state.data.charCodeAt(state.position + 1);
            if (next !== quoteCharCode) {
                // get rid of the quotes.
                state.tokenStart++;
                state.tokenEnd = state.position;
                ++state.position;
                return skipEmpty(state, delimiterCharCode);
            }
        }
        ++state.position;
    }
    state.tokenEnd = state.position;
}

/**
 * Skips empty chars.
 * Returns true when the current char is a newline.
 */
function skipEmpty(state: Tokenizer, delimiterCharCode: number) {
    while (state.position < state.length) {
        const c = state.data.charCodeAt(state.position);
        if (c !== 9 && c !== 32 && c !== delimiterCharCode) {  // \t or ' '
            return c === 10 || c === 13;  // \n or \r
        }
        ++state.position;
    }
}

function skipWhitespace(state: Tokenizer) {
    let prev = -1;
    while (state.position < state.length) {
        const c = state.data.charCodeAt(state.position);
        switch (c) {
            case 9:  // '\t'
            case 32:  // ' '
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
        if (c === 10 || c === 13) return;  // \n or \r
        ++state.position;
    }
}

/**
 * Move to the next token.
 * Returns true when the current char is a newline, i.e. indicating a full record.
 */
function moveNextInternal(state: State) {
    const tokenizer = state.tokenizer;
    skipWhitespace(tokenizer);

    if (tokenizer.position >= tokenizer.length) {
        state.tokenType = CsvTokenType.End;
        return false;
    }

    tokenizer.tokenStart = tokenizer.position;
    tokenizer.tokenEnd = tokenizer.position;
    const c = state.data.charCodeAt(tokenizer.position);
    switch (c) {
        case state.commentCharCode:
            state.tokenType = CsvTokenType.Comment;
            skipLine(tokenizer);
            break;
        case state.quoteCharCode:
            state.tokenType = CsvTokenType.Value;
            return eatQuoted(tokenizer, state.quoteCharCode, state.delimiterCharCode);
        default:
            state.tokenType = CsvTokenType.Value;
            return eatValue(tokenizer, state.delimiterCharCode);
    }
}

/**
 * Moves to the next non-comment token/line.
 * Returns true when the current char is a newline, i.e. indicating a full record.
 */
function moveNext(state: State) {
    let newRecord = moveNextInternal(state);
    while (state.tokenType === CsvTokenType.Comment) {
        newRecord = moveNextInternal(state);
    }
    return newRecord;
}

function readRecordsChunk(chunkSize: number, state: State) {
    if (state.tokenType === CsvTokenType.End) return 0;

    let counter = 0;
    let newRecord: boolean | undefined;

    const { tokens, tokenizer } = state;

    while (state.tokenType === CsvTokenType.Value && counter < chunkSize) {
        TokenBuilder.add(tokens[state.fieldCount % state.columnCount], tokenizer.tokenStart, tokenizer.tokenEnd);
        ++state.fieldCount;
        newRecord = moveNext(state);
        if (newRecord) {
            ++state.recordCount;
            ++counter;
        }
    }
    return counter;
}

function readRecordsChunks(state: State) {
    let newRecord = moveNext(state);
    if (newRecord) ++state.recordCount;
    return chunkedSubtask(state.runtimeCtx, 100000, state, readRecordsChunk,
        (ctx, state) => ctx.update({ message: 'Parsing...', current: state.tokenizer.position, max: state.data.length }));
}

function addColumn (state: State) {
    state.columnNames.push(Tokenizer.getTokenString(state.tokenizer));
    state.tokens.push(TokenBuilder.create(state.tokenizer.data, state.data.length / 80));
}

function init(state: State) {
    let newRecord = moveNext(state);
    while (!newRecord) {
        addColumn(state);
        newRecord = moveNext(state);
    }
    addColumn(state);
    state.columnCount = state.columnNames.length;
    if (state.noColumnNamesRecord) {
        state.columnNames.forEach((x, i, arr) => arr[i] = i + '');
        Tokenizer.reset(state.tokenizer);
    }
}

async function handleRecords(state: State): Promise<Data.CsvTable> {
    init(state);
    await readRecordsChunks(state);

    const columns: Data.CsvColumns = Object.create(null);
    for (let i = 0; i < state.columnCount; ++i) {
        columns[state.columnNames[i]] = Field(state.tokens[i]);
    }

    return Data.CsvTable(state.recordCount, state.columnNames, columns);
}

async function parseInternal(data: string, ctx: RuntimeContext, opts: CsvOptions): Promise<Result<Data.CsvFile>> {
    const state = State(data, ctx, opts);

    ctx.update({ message: 'Parsing...', current: 0, max: data.length });
    const table = await handleRecords(state);
    const result = Data.CsvFile(table);
    return Result.success(result);
}

interface CsvOptions {
    quote: string;
    comment: string;
    delimiter: string;
    noColumnNames: boolean;
}

export function parseCsv(data: string, opts?: Partial<CsvOptions>) {
    const completeOpts = Object.assign({}, { quote: '"', comment: '#', delimiter: ',', noColumnNames: false }, opts);
    return Task.create<Result<Data.CsvFile>>('Parse CSV', async ctx => {
        return await parseInternal(data, ctx, completeOpts);
    });
}