/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ReaderResult as Result } from '../result';
import { Task, RuntimeContext } from '../../../mol-task';
import { PlyFile, PlyType, PlyElement, PlyTypeByteLength } from './schema';
import { Tokenizer, TokenBuilder, Tokens } from '../common/text/tokenizer';
import { Column } from '../../../mol-data/db';
import { TokenColumn } from '../common/text/column/token';
import { StringLike } from '../../common/string-like';
import { utf8Read, utf8ReadLong } from '../../common/utf8';


// TODO parse elements asynchronously
// TODO handle lists with appended properties

type PlyFormat = 'ascii' | 'binary_little_endian' | 'binary_big_endian'

interface State {
    data: StringLike
    tokenizer: Tokenizer
    runtimeCtx: RuntimeContext
    format: PlyFormat
    formatPreset: boolean
    littleEndian: boolean
    binaryData: Uint8Array | undefined
    binaryOffset: number

    comments: string[]
    elementSpecs: ElementSpec[]
    elements: PlyElement[]
}

function State(data: StringLike, runtimeCtx: RuntimeContext, presetFormat?: PlyFormat): State {
    const tokenizer = Tokenizer(data);
    return {
        data,
        tokenizer,
        runtimeCtx,
        format: presetFormat ?? 'ascii',
        formatPreset: presetFormat !== undefined,
        littleEndian: presetFormat !== 'binary_big_endian',
        binaryData: undefined,
        binaryOffset: 0,

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
    if (!state.formatPreset) {
        const formatLine = headerLines[1];
        if (formatLine === 'format ascii 1.0') {
            state.format = 'ascii';
            state.littleEndian = true;
        } else if (formatLine === 'format binary_little_endian 1.0') {
            state.format = 'binary_little_endian';
            state.littleEndian = true;
        } else if (formatLine === 'format binary_big_endian 1.0') {
            state.format = 'binary_big_endian';
            state.littleEndian = false;
        } else {
            throw new Error(`unsupported PLY format '${formatLine}'`);
        }
    }

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
    if (state.format === 'ascii') {
        for (let i = 0, il = elementSpecs.length; i < il; ++i) {
            const spec = elementSpecs[i];
            if (spec.kind === 'table') parseTableElement(state, spec);
            else if (spec.kind === 'list') parseListElement(state, spec);
        }
    } else {
        for (let i = 0, il = elementSpecs.length; i < il; ++i) {
            const spec = elementSpecs[i];
            if (spec.kind === 'table') parseBinaryTableElement(state, spec);
            else if (spec.kind === 'list') parseBinaryListElement(state, spec);
        }
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

function readBinaryValue(dv: DataView, offset: number, type: PlyType, littleEndian: boolean): number {
    switch (type) {
        case 'char': case 'int8': return dv.getInt8(offset);
        case 'uchar': case 'uint8': return dv.getUint8(offset);
        case 'short': case 'int16': return dv.getInt16(offset, littleEndian);
        case 'ushort': case 'uint16': return dv.getUint16(offset, littleEndian);
        case 'int': case 'int32': return dv.getInt32(offset, littleEndian);
        case 'uint': case 'uint32': return dv.getUint32(offset, littleEndian);
        case 'float': case 'float32': return dv.getFloat32(offset, littleEndian);
        case 'double': case 'float64': return dv.getFloat64(offset, littleEndian);
    }
}

type PlyTypedArray = Int8Array | Uint8Array | Int16Array | Uint16Array | Int32Array | Uint32Array | Float32Array | Float64Array;

function makeTypedArray(type: PlyType, count: number): PlyTypedArray {
    switch (type) {
        case 'char': case 'int8': return new Int8Array(count);
        case 'uchar': case 'uint8': return new Uint8Array(count);
        case 'short': case 'int16': return new Int16Array(count);
        case 'ushort': case 'uint16': return new Uint16Array(count);
        case 'int': case 'int32': return new Int32Array(count);
        case 'uint': case 'uint32': return new Uint32Array(count);
        case 'float': case 'float32': return new Float32Array(count);
        case 'double': case 'float64': return new Float64Array(count);
    }
}

function parseBinaryTableElement(state: State, spec: TableElementSpec) {
    const { elements, binaryData, littleEndian } = state;
    const { count, properties } = spec;
    const propertyCount = properties.length;
    const propertyNames: string[] = [];
    const propertyTypes: PlyType[] = [];
    const propertyColumns = new Map<string, Column<number>>();

    const arrays: PlyTypedArray[] = [];
    for (let i = 0; i < propertyCount; i++) arrays.push(makeTypedArray(properties[i].type, count));

    const dv = new DataView(binaryData!.buffer, binaryData!.byteOffset);
    for (let i = 0; i < count; i++) {
        for (let j = 0; j < propertyCount; j++) {
            const { type } = properties[j];
            arrays[j][i] = readBinaryValue(dv, state.binaryOffset, type, littleEndian);
            state.binaryOffset += PlyTypeByteLength[type];
        }
    }

    for (let i = 0; i < propertyCount; i++) {
        const { type, name } = properties[i];
        const schema = getColumnSchema(type);
        const column = Column.ofArray({ array: arrays[i], schema });
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

function parseBinaryListElement(state: State, spec: ListElementSpec) {
    const { elements, binaryData, littleEndian } = state;
    const { count, property } = spec;
    const { countType, dataType } = property;
    const dv = new DataView(binaryData!.buffer, binaryData!.byteOffset);

    const allEntries: number[] = [];
    const offsets = new Uint32Array(count + 1);
    for (let i = 0; i < count; i++) {
        const listCount = readBinaryValue(dv, state.binaryOffset, countType, littleEndian);
        state.binaryOffset += PlyTypeByteLength[countType];
        offsets[i + 1] = offsets[i] + listCount;
        for (let k = 0; k < listCount; k++) {
            allEntries.push(readBinaryValue(dv, state.binaryOffset, dataType, littleEndian));
            state.binaryOffset += PlyTypeByteLength[dataType];
        }
    }

    /** holds row value entries transiently */
    const listValue = {
        entries: [] as number[],
        count: 0
    };

    elements.push({
        kind: 'list',
        rowCount: count,
        name: property.name,
        type: property.dataType,
        value: (row: number) => {
            const start = offsets[row];
            const end = offsets[row + 1];
            listValue.count = end - start;
            for (let i = start; i < end; i++) {
                listValue.entries[i - start] = allEntries[i];
            }
            return listValue;
        }
    });
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

async function parseInternal(data: StringLike | Uint8Array, ctx: RuntimeContext): Promise<Result<PlyFile>> {
    let stateData: StringLike;
    let binaryData: Uint8Array | undefined;
    let binaryOffset = 0;
    let presetFormat: PlyFormat | undefined;

    if (data instanceof Uint8Array) {
        binaryOffset = findEndHeaderByteOffset(data);
        const headerText = utf8Read(data, 0, binaryOffset);
        const format = detectPlyFormat(headerText);
        presetFormat = format;
        if (format !== 'ascii') {
            stateData = headerText;
            binaryData = data;
        } else {
            // ASCII PLY: decode the full byte array as text
            stateData = utf8ReadLong(data, 0, data.length);
        }
    } else {
        stateData = data;
    }

    const state = State(stateData, ctx, presetFormat);
    if (binaryData !== undefined) {
        state.binaryData = binaryData;
        state.binaryOffset = binaryOffset;
    }

    ctx.update({ message: 'Parsing...', current: 0, max: data.length });
    parseHeader(state);
    parseElements(state);
    const { elements, elementSpecs, comments } = state;
    const elementNames = elementSpecs.map(s => s.name);
    const result = PlyFile(elements, elementNames, comments);
    return Result.success(result);
}

/** Extracts the PLY format from the header text (second line). */
function detectPlyFormat(headerText: string): PlyFormat {
    const newline = headerText.indexOf('\n');
    const formatLine = newline === -1 ? '' : headerText.substring(newline + 1, headerText.indexOf('\n', newline + 1)).replace(/\r$/, '');
    if (formatLine === 'format binary_little_endian 1.0') return 'binary_little_endian';
    if (formatLine === 'format binary_big_endian 1.0') return 'binary_big_endian';
    return 'ascii';
}

/** Returns the byte offset of the first byte after the `end_header\n` line. */
function findEndHeaderByteOffset(data: Uint8Array): number {
    const target = [101, 110, 100, 95, 104, 101, 97, 100, 101, 114]; // 'end_header'
    const tlen = target.length;
    for (let i = 0; i <= data.length - tlen; i++) {
        // Only match at the start of a line
        if (i !== 0 && data[i - 1] !== 10 /* \n */) continue;
        let matched = true;
        for (let j = 0; j < tlen; j++) {
            if (data[i + j] !== target[j]) {
                matched = false;
                break;
            }
        }
        if (!matched) continue;

        let k = i + tlen;
        if (k < data.length && data[k] === 13) k++; // skip \r
        if (k < data.length && data[k] === 10) k++; // skip \n
        return k;
    }
    throw new Error(`no 'end_header' record found`);
}

export function parsePly(data: StringLike | Uint8Array) {
    return Task.create<Result<PlyFile>>('Parse PLY', async ctx => {
        return await parseInternal(data, ctx);
    });
}
