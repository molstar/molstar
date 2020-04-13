/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ReaderResult as Result } from '../result';
import { Task, RuntimeContext } from '../../../mol-task';
import { PlyFile, PlyType, PlyElement } from './schema';
import { Tokenizer, TokenBuilder, Tokens } from '../common/text/tokenizer';
import { Column } from '../../../mol-data/db';
import { TokenColumn } from '../common/text/column/token';

// TODO add support for binary ply files
// TODO parse elements asynchronously
// TODO handle lists with appended properties

interface State {
    data: string
    tokenizer: Tokenizer
    runtimeCtx: RuntimeContext

    comments: string[]
    elementSpecs: ElementSpec[]
    elements: PlyElement[]
}

function State(data: string, runtimeCtx: RuntimeContext): State {
    const tokenizer = Tokenizer(data);
    return {
        data,
        tokenizer,
        runtimeCtx,

        comments: [],
        elementSpecs: [],
        elements: []
    };
}

type ColumnProperty = { kind: 'column', type: PlyType, name: string }
type ListProperty = { kind: 'list', countType: PlyType, dataType: PlyType, name: string }
type Property = ColumnProperty | ListProperty

type TableElementSpec = { kind: 'table', name: string, count: number, properties: ColumnProperty[] }
type ListElementSpec = { kind: 'list', name: string, count: number, property: ListProperty }
type ElementSpec = TableElementSpec | ListElementSpec

function markHeader(tokenizer: Tokenizer) {
    const endHeaderIndex = tokenizer.data.indexOf('end_header', tokenizer.position);
    if (endHeaderIndex === -1) throw new Error(`no 'end_header' record found`);
    // TODO set `tokenizer.lineNumber` correctly
    tokenizer.tokenStart = tokenizer.position;
    tokenizer.tokenEnd = endHeaderIndex;
    tokenizer.position = endHeaderIndex;
    Tokenizer.eatLine(tokenizer);
}

function parseHeader(state: State) {
    const { tokenizer, comments, elementSpecs } = state;

    markHeader(tokenizer);
    const headerLines = Tokenizer.getTokenString(tokenizer).split(/\r?\n/);

    if (headerLines[0] !== 'ply') throw new Error(`data not starting with 'ply'`);
    if (headerLines[1] !== 'format ascii 1.0') throw new Error(`format not 'ascii 1.0'`);

    let currentName: string | undefined;
    let currentCount: number | undefined;
    let currentProperties: Property[] | undefined;


    function addCurrentElementSchema() {
        if (currentName !== undefined && currentCount !== undefined && currentProperties !== undefined) {
            let isList = false;
            for (let i = 0, il = currentProperties.length; i < il; ++i) {
                const p = currentProperties[i];
                if (p.kind === 'list') {
                    isList = true;
                    break;
                }
            }
            if (isList && currentProperties.length !== 1) {
                // TODO handle lists with appended properties
                //      currently only the list part will be accessible
            }
            if (isList) {
                elementSpecs.push({
                    kind: 'list',
                    name: currentName,
                    count: currentCount,
                    property: currentProperties[0] as ListProperty
                });
            } else {
                elementSpecs.push({
                    kind: 'table',
                    name: currentName,
                    count: currentCount,
                    properties: currentProperties as ColumnProperty[]
                });
            }
        }
    }

    for (let i = 2, il = headerLines.length; i < il; ++i) {
        const l = headerLines[i];
        const ls = l.split(' ');
        if (l.startsWith('comment')) {
            comments.push(l.substr(8));
        } else if (l.startsWith('element')) {
            addCurrentElementSchema();
            currentProperties = [];
            currentName = ls[1];
            currentCount = parseInt(ls[2]);
        } else if (l.startsWith('property')) {
            if (currentProperties === undefined) throw new Error(`properties outside of element`);
            if (ls[1] === 'list') {
                currentProperties.push({
                    kind: 'list',
                    countType: PlyType(ls[2]),
                    dataType: PlyType(ls[3]),
                    name: ls[4]
                });
            } else {
                currentProperties.push({
                    kind: 'column',
                    type: PlyType(ls[1]),
                    name: ls[2]
                });
            }
        } else if (l.startsWith('end_header')) {
            addCurrentElementSchema();
        } else {
            console.warn('unknown header line');
        }
    }
}

function parseElements(state: State) {
    const { elementSpecs } = state;
    for (let i = 0, il = elementSpecs.length; i < il; ++i) {
        const spec = elementSpecs[i];
        if (spec.kind === 'table') parseTableElement(state, spec);
        else if (spec.kind === 'list') parseListElement(state, spec);
    }
}

function getColumnSchema(type: PlyType): Column.Schema {
    switch (type) {
        case 'char': case 'uchar': case 'int8': case 'uint8':
        case 'short': case 'ushort': case 'int16': case 'uint16':
        case 'int': case 'uint': case 'int32': case 'uint32':
            return Column.Schema.int;
        case 'float': case 'double': case 'float32': case 'float64':
            return Column.Schema.float;
    }
}

function parseTableElement(state: State, spec: TableElementSpec) {
    const { elements, tokenizer } = state;
    const { count, properties } = spec;
    const propertyCount = properties.length;
    const propertyNames: string[] = [];
    const propertyTypes: PlyType[] = [];
    const propertyTokens: Tokens[] = [];
    const propertyColumns = new Map<string, Column<number>>();

    for (let i = 0, il = propertyCount; i < il; ++i) {
        const tokens = TokenBuilder.create(tokenizer.data, count * 2);
        propertyTokens.push(tokens);
    }

    for (let i = 0, il = count; i < il; ++i) {
        for (let j = 0, jl = propertyCount; j < jl; ++j) {
            Tokenizer.skipWhitespace(tokenizer);
            Tokenizer.markStart(tokenizer);
            Tokenizer.eatValue(tokenizer);
            TokenBuilder.addUnchecked(propertyTokens[j], tokenizer.tokenStart, tokenizer.tokenEnd);
        }
    }

    for (let i = 0, il = propertyCount; i < il; ++i) {
        const { type, name } = properties[i];
        const column = TokenColumn(propertyTokens[i], getColumnSchema(type));
        propertyNames.push(name);
        propertyTypes.push(type);
        propertyColumns.set(name, column);
    }

    elements.push({
        kind: 'table',
        rowCount: count,
        propertyNames,
        propertyTypes,
        getProperty: (name: string) => propertyColumns.get(name)
    });
}

function parseListElement(state: State, spec: ListElementSpec) {
    const { elements, tokenizer } = state;
    const { count, property } = spec;

    // initial tokens size assumes triangle index data
    const tokens = TokenBuilder.create(tokenizer.data, count * 2 * 3);

    const offsets = new Uint32Array(count + 1);
    let entryCount = 0;

    for (let i = 0, il = count; i < il; ++i) {
        Tokenizer.skipWhitespace(tokenizer);
        Tokenizer.markStart(tokenizer);
        while (Tokenizer.skipWhitespace(tokenizer) !== 10) {
            ++entryCount;
            Tokenizer.markStart(tokenizer);
            Tokenizer.eatValue(tokenizer);
            TokenBuilder.addToken(tokens, tokenizer);
        }
        offsets[i + 1] = entryCount;
    }

    /** holds row value entries transiently */
    const listValue = {
        entries: [] as number[],
        count: 0
    };

    const column = TokenColumn(tokens, getColumnSchema(property.dataType));

    elements.push({
        kind: 'list',
        rowCount: count,
        name: property.name,
        type: property.dataType,
        value: (row: number) => {
            const offset = offsets[row] + 1;
            const count = column.value(offset - 1);
            for (let i = offset, il = offset + count; i < il; ++i) {
                listValue.entries[i - offset] = column.value(i);
            }
            listValue.count = count;
            return listValue;
        }
    });
}

async function parseInternal(data: string, ctx: RuntimeContext): Promise<Result<PlyFile>> {
    const state = State(data, ctx);
    ctx.update({ message: 'Parsing...', current: 0, max: data.length });
    parseHeader(state);
    parseElements(state);
    const { elements, elementSpecs, comments } = state;
    const elementNames = elementSpecs.map(s => s.name);
    const result = PlyFile(elements, elementNames, comments);
    return Result.success(result);
}

export function parse(data: string) {
    return Task.create<Result<PlyFile>>('Parse PLY', async ctx => {
        return await parseInternal(data, ctx);
    });
}

export default parse;