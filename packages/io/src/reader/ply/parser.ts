/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ReaderResult as Result } from '../result.js';
import { Task, type RuntimeContext } from '@molstar/core/task';
import { PlyFile, PlyType, type PlyElement, PlyTypeByteLength } from './schema.js';
import { TextByteTokenizer, TextByteTokenBuilder, type TextByteTokens } from '../common/text/text-byte-tokenizer.js';
import { Column } from '@molstar/core/data/db';
import { TextByteTokenColumn } from '../common/text/column/text-byte-token.js';
import type { StringLike } from '@molstar/core/util/string-like';
import { utf8ReadValue } from '@molstar/core/util/utf8';
import type { ReadonlyUint8Array } from '@molstar/core/util/type-helpers';

// TODO parse elements asynchronously
// TODO handle lists with appended properties

type PlyFormat = 'ascii' | 'binary_little_endian' | 'binary_big_endian';

interface State {
  /** The bytes of the file: an ASCII body is tokenized from them, a binary body read from them. */
  data: ReadonlyUint8Array;
  tokenizer: TextByteTokenizer;
  runtimeCtx: RuntimeContext;
  format: PlyFormat;
  littleEndian: boolean;
  /** The header, `ply` to `end_header` (decoded: it is small). */
  headerText: string;
  binaryOffset: number;

  comments: string[];
  elementSpecs: ElementSpec[];
  elements: PlyElement[];
}

function State(data: ReadonlyUint8Array, runtimeCtx: RuntimeContext): State {
  const tokenizer = TextByteTokenizer(data);
  return {
    data,
    tokenizer,
    runtimeCtx,
    format: 'ascii',
    littleEndian: true,
    headerText: '',
    binaryOffset: 0,

    comments: [],
    elementSpecs: [],
    elements: [],
  };
}

type ColumnProperty = { kind: 'column'; type: PlyType; name: string };
type ListProperty = { kind: 'list'; countType: PlyType; dataType: PlyType; name: string };
type Property = ColumnProperty | ListProperty;

type TableElementSpec = { kind: 'table'; name: string; count: number; properties: ColumnProperty[] };
type ListElementSpec = { kind: 'list'; name: string; count: number; property: ListProperty };
type ElementSpec = TableElementSpec | ListElementSpec;

/**
 * Reads the header, `ply` to the `end_header` line (which starts a line), into `state.headerText`, and moves the
 * tokenizer to the first byte of the body.
 */
function markHeader(state: State) {
  const { data, tokenizer } = state;
  const start = findEndHeader(data);
  let end = start + EndHeader.length;
  // the rest of the `end_header` line
  while (end < data.length && data[end] !== 10 && data[end] !== 13) end++;
  state.headerText = utf8ReadValue(data, 0, end);
  tokenizer.position = end;
  if (data[tokenizer.position] === 13) tokenizer.position++;
  if (data[tokenizer.position] === 10) tokenizer.position++;
  state.binaryOffset = tokenizer.position;
}

function parseHeader(state: State) {
  const { comments, elementSpecs } = state;

  markHeader(state);
  const headerLines = state.headerText.split(/\r?\n/);

  if (headerLines[0].trim() !== 'ply') throw new Error(`data not starting with 'ply'`);
  const formatLine = headerLines[1];
  // by its words, not the exact line: any version, spaces around them (files with a trailing space read before)
  const [keyword, kind] = formatLine?.trim().split(/\s+/) ?? [];
  if (keyword === 'format' && kind === 'ascii') {
    state.format = 'ascii';
    state.littleEndian = true;
  } else if (keyword === 'format' && kind === 'binary_little_endian') {
    state.format = 'binary_little_endian';
    state.littleEndian = true;
  } else if (keyword === 'format' && kind === 'binary_big_endian') {
    state.format = 'binary_big_endian';
    state.littleEndian = false;
  } else {
    throw new Error(`unsupported PLY format '${formatLine}'`);
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
          property: currentProperties[0] as ListProperty,
        });
      } else {
        elementSpecs.push({
          kind: 'table',
          name: currentName,
          count: currentCount,
          properties: currentProperties as ColumnProperty[],
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
          name: ls[4],
        });
      } else {
        currentProperties.push({
          kind: 'column',
          type: PlyType(ls[1]),
          name: ls[2],
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
    case 'char':
    case 'uchar':
    case 'int8':
    case 'uint8':
    case 'short':
    case 'ushort':
    case 'int16':
    case 'uint16':
    case 'int':
    case 'uint':
    case 'int32':
    case 'uint32':
      return Column.Schema.int;
    case 'float':
    case 'double':
    case 'float32':
    case 'float64':
      return Column.Schema.float;
  }
}

function readBinaryValue(dv: DataView, offset: number, type: PlyType, littleEndian: boolean): number {
  switch (type) {
    case 'char':
    case 'int8':
      return dv.getInt8(offset);
    case 'uchar':
    case 'uint8':
      return dv.getUint8(offset);
    case 'short':
    case 'int16':
      return dv.getInt16(offset, littleEndian);
    case 'ushort':
    case 'uint16':
      return dv.getUint16(offset, littleEndian);
    case 'int':
    case 'int32':
      return dv.getInt32(offset, littleEndian);
    case 'uint':
    case 'uint32':
      return dv.getUint32(offset, littleEndian);
    case 'float':
    case 'float32':
      return dv.getFloat32(offset, littleEndian);
    case 'double':
    case 'float64':
      return dv.getFloat64(offset, littleEndian);
  }
}

type PlyTypedArray =
  | Int8Array
  | Uint8Array
  | Int16Array
  | Uint16Array
  | Int32Array
  | Uint32Array
  | Float32Array
  | Float64Array;

function makeTypedArray(type: PlyType, count: number): PlyTypedArray {
  switch (type) {
    case 'char':
    case 'int8':
      return new Int8Array(count);
    case 'uchar':
    case 'uint8':
      return new Uint8Array(count);
    case 'short':
    case 'int16':
      return new Int16Array(count);
    case 'ushort':
    case 'uint16':
      return new Uint16Array(count);
    case 'int':
    case 'int32':
      return new Int32Array(count);
    case 'uint':
    case 'uint32':
      return new Uint32Array(count);
    case 'float':
    case 'float32':
      return new Float32Array(count);
    case 'double':
    case 'float64':
      return new Float64Array(count);
  }
}

function parseBinaryTableElement(state: State, spec: TableElementSpec) {
  const { elements, data: binaryData, littleEndian } = state;
  const { count, properties } = spec;
  const propertyCount = properties.length;
  const propertyNames: string[] = [];
  const propertyTypes: PlyType[] = [];
  const propertyColumns = new Map<string, Column<number>>();

  const arrays: PlyTypedArray[] = [];
  for (let i = 0; i < propertyCount; i++) arrays.push(makeTypedArray(properties[i].type, count));

  const dv = new DataView(binaryData.buffer, binaryData.byteOffset);
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
    getProperty: (name: string) => propertyColumns.get(name),
  });
}

function parseBinaryListElement(state: State, spec: ListElementSpec) {
  const { elements, data: binaryData, littleEndian } = state;
  const { count, property } = spec;
  const { countType, dataType } = property;
  const dv = new DataView(binaryData.buffer, binaryData.byteOffset);

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
    count: 0,
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
    },
  });
}

function parseTableElement(state: State, spec: TableElementSpec) {
  const { elements, tokenizer } = state;
  const { count, properties } = spec;
  const propertyCount = properties.length;
  const propertyNames: string[] = [];
  const propertyTypes: PlyType[] = [];
  const propertyTokens: TextByteTokens[] = [];
  const propertyColumns = new Map<string, Column<number>>();

  for (let i = 0, il = propertyCount; i < il; ++i) {
    const tokens = TextByteTokenBuilder.create(tokenizer.data, count * 2);
    propertyTokens.push(tokens);
  }

  for (let i = 0, il = count; i < il; ++i) {
    for (let j = 0, jl = propertyCount; j < jl; ++j) {
      TextByteTokenizer.skipWhitespace(tokenizer);
      TextByteTokenizer.markStart(tokenizer);
      TextByteTokenizer.eatValue(tokenizer);
      TextByteTokenBuilder.addUnchecked(propertyTokens[j], tokenizer.tokenStart, tokenizer.tokenEnd);
    }
  }

  for (let i = 0, il = propertyCount; i < il; ++i) {
    const { type, name } = properties[i];
    const column = TextByteTokenColumn(propertyTokens[i], getColumnSchema(type));
    propertyNames.push(name);
    propertyTypes.push(type);
    propertyColumns.set(name, column);
  }

  elements.push({
    kind: 'table',
    rowCount: count,
    propertyNames,
    propertyTypes,
    getProperty: (name: string) => propertyColumns.get(name),
  });
}

function parseListElement(state: State, spec: ListElementSpec) {
  const { elements, tokenizer } = state;
  const { count, property } = spec;

  // initial tokens size assumes triangle index data
  const tokens = TextByteTokenBuilder.create(tokenizer.data, count * 2 * 3);

  const offsets = new Uint32Array(count + 1);
  let entryCount = 0;

  for (let i = 0, il = count; i < il; ++i) {
    TextByteTokenizer.skipWhitespace(tokenizer);
    TextByteTokenizer.markStart(tokenizer);
    while (true) {
      // a row ends at its line break (\n, or \r alone) or at the end of the data (a last line without one hung)
      const ws = TextByteTokenizer.skipWhitespace(tokenizer);
      if (ws === 10 || ws === 13 || tokenizer.position >= tokenizer.length) break;
      ++entryCount;
      TextByteTokenizer.markStart(tokenizer);
      TextByteTokenizer.eatValue(tokenizer);
      TextByteTokenBuilder.addToken(tokens, tokenizer);
    }
    offsets[i + 1] = entryCount;
  }

  /** holds row value entries transiently */
  const listValue = {
    entries: [] as number[],
    count: 0,
  };

  const column = TextByteTokenColumn(tokens, getColumnSchema(property.dataType));

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
    },
  });
}

async function parseInternal(data: StringLike | ReadonlyUint8Array, ctx: RuntimeContext): Promise<Result<PlyFile>> {
  const state = State(TextByteTokenizer.bytesOf(data), ctx);
  ctx.update({ message: 'Parsing...', current: 0, max: state.data.length });
  parseHeader(state);
  parseElements(state);
  const { elements, elementSpecs, comments } = state;
  const elementNames = elementSpecs.map((s) => s.name);
  const result = PlyFile(elements, elementNames, comments);
  return Result.success(result);
}

const EndHeader = 'end_header';

/** Returns the byte offset of the `end_header` record, at the start of a line. */
function findEndHeader(data: ReadonlyUint8Array): number {
  const target = [101, 110, 100, 95, 104, 101, 97, 100, 101, 114]; // 'end_header'
  const tlen = target.length;
  for (let i = 0; i <= data.length - tlen; i++) {
    // Only match at the start of a line
    if (i !== 0 && data[i - 1] !== 10 /* \n */ && data[i - 1] !== 13 /* \r */) continue;
    let matched = true;
    for (let j = 0; j < tlen; j++) {
      if (data[i + j] !== target[j]) {
        matched = false;
        break;
      }
    }
    if (matched) return i;
  }
  throw new Error(`no 'end_header' record found`);
}

/** Parses a PLY file, ASCII or binary (by its header), from its bytes or its text (a string is encoded once). */
export function parsePly(data: StringLike | ReadonlyUint8Array) {
  return Task.create<Result<PlyFile>>('Parse PLY', async (ctx) => {
    return await parseInternal(data, ctx);
  });
}
