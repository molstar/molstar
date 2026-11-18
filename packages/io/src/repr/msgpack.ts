/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The MessagePack codec of the tree (core's repr spec "Codecs"; `_spec/repr.md`), on `common/msgpack`'s encoder and
// decoder and their extension hooks, which it does not change. The shape is the JSON codec's (a node is a map
// `{ kind, id?, params?, custom?, children?, args? }`, a list an array, a record a map without `kind`); what
// MessagePack lacks is extensions: NaN, the infinities and -0 (`0x10`), a typed column (`0x11`-`0x18`: its elements'
// bytes, little-endian), a category column (`0x19`) and a bigint (`0x1a`). A column is always typed, whatever the
// language: there is no plain form to decode. The encoder writes a column's bytes where they are (no copy); the decoder
// gives a view into its input where the payload is aligned for the element type, else one copy (`repr.md` §2).

import type { Codec, DecodeOptions, Decoded } from '@molstar/core/repr/codec/index';
import {
  assertKeys,
  assertPrintable,
  CategoryColumn,
  type Column,
  canonicalColumn,
  categoryFromCodes,
  columnConstructor,
  columnType,
  defineIssues,
  IdHint,
  IdPattern,
  type Issue,
  isColumn,
  isNode,
  isReservedKey,
  issueOf,
  KindHint,
  KindPattern,
  MaxDepth,
  type Node,
  type NumericColumnType,
  TreeIssues,
  type TreePath,
  type Value,
  type ValueRecord,
} from '@molstar/core/repr/tree/index';
import { utf8ByteCount, utf8Read, utf8Write } from '@molstar/core/util/utf8';
import { createMsgPackState, parseMsgPackValue } from '../common/msgpack/decode.js';
import { encodeMsgPack, type MsgPackExtData } from '../common/msgpack/encode.js';

/** The extension types. */
export const MsgPackExt = {
  /** NaN, +Infinity, -Infinity, -0: one byte, 0 to 3. */
  number: 0x10,
  int8: 0x11,
  uint8: 0x12,
  int16: 0x13,
  uint16: 0x14,
  int32: 0x15,
  uint32: 0x16,
  float32: 0x17,
  float64: 0x18,
  /** A category: u32 dictionary size, each entry a u32 byte length and its UTF-8, then the int32 codes. */
  category: 0x19,
  /** A bigint: its decimal digits, UTF-8. */
  bigint: 0x1a,
} as const;

/** The codes of reading bytes that are no tree. */
export const MsgPackIssues = defineIssues('repr.codec.msgpack', {
  'repr.codec.msgpack.syntax': {
    severity: 'error',
    description: 'The bytes are not one MessagePack value: truncated, malformed, or followed by more bytes',
  },
});

const ColumnExt: Readonly<Record<NumericColumnType, number>> = {
  int8: MsgPackExt.int8,
  uint8: MsgPackExt.uint8,
  int16: MsgPackExt.int16,
  uint16: MsgPackExt.uint16,
  int32: MsgPackExt.int32,
  uint32: MsgPackExt.uint32,
  float32: MsgPackExt.float32,
  float64: MsgPackExt.float64,
};
const ColumnOfExt = new Map<number, NumericColumnType>(
  (Object.keys(ColumnExt) as NumericColumnType[]).map((t) => [ColumnExt[t], t]),
);
const NodeKeys = new Set(['kind', 'id', 'params', 'custom', 'children', 'args']);
const Specials = [Number.NaN, Number.POSITIVE_INFINITY, Number.NEGATIVE_INFINITY, -0];

// ---------------------------------------------------------------------------------------------------------------
// Encoding

/** An atom or column the encoder cannot write itself, made an extension before it is encoded. */
class Ext {
  constructor(readonly data: MsgPackExtData) {}
}
const SpecialExts = Specials.map((_, i) => new Ext({ type: MsgPackExt.number, data: Uint8Array.of(i) }));

/** The tree as plain maps, arrays and atoms, with `Ext` for what MessagePack lacks; columns are not walked. */
function toPlainNode(node: Node): unknown {
  assertPrintable(node);
  const out: Record<string, unknown> = { kind: node.kind };
  if (node.id !== undefined) out.id = node.id;
  if (node.params !== undefined) out.params = toPlainRecord(node.params);
  if (node.custom !== undefined) out.custom = toPlainRecord(node.custom);
  if (node.children !== undefined) out.children = node.children.map(toPlainNode);
  if (node.args?.length) out.args = node.args.map(toPlainValue);
  return out;
}

function toPlainRecord(record: ValueRecord): Record<string, unknown> {
  assertKeys(record);
  const out: Record<string, unknown> = {};
  // an `undefined` value is absent, as in the other encodings and the hash
  for (const key of Object.keys(record)) if (record[key] !== undefined) out[key] = toPlainValue(record[key]);
  return out;
}

function specialIndex(value: number): number {
  if (Number.isNaN(value)) return 0;
  if (value === Number.POSITIVE_INFINITY) return 1;
  return value === Number.NEGATIVE_INFINITY ? 2 : 3;
}

function toPlainValue(value: Value): unknown {
  switch (typeof value) {
    case 'number':
      return Number.isFinite(value) && !Object.is(value, -0) ? value : SpecialExts[specialIndex(value)];
    case 'bigint': {
      const digits = value.toString();
      const data = new Uint8Array(utf8ByteCount(digits));
      utf8Write(data, 0, digits);
      return new Ext({ type: MsgPackExt.bigint, data });
    }
  }
  if (value === null || typeof value !== 'object') return value;
  if (isColumn(value)) return columnExt(value);
  if (Array.isArray(value)) return value.map(toPlainValue);
  if (isNode(value)) return toPlainNode(value);
  return toPlainRecord(value as ValueRecord);
}

function columnExt(column: Column): Ext {
  canonicalColumn(column);
  if (column instanceof CategoryColumn) return new Ext({ type: MsgPackExt.category, data: categoryBytes(column) });
  const data = new Uint8Array(column.buffer, column.byteOffset, column.byteLength);
  return new Ext({ type: ColumnExt[columnType(column) as NumericColumnType], data });
}

function categoryBytes(column: CategoryColumn): Uint8Array {
  const { dict, codes } = column;
  let size = 4;
  for (const d of dict) size += 4 + utf8ByteCount(d);
  const bytes = new Uint8Array(size + codes.byteLength);
  const view = new DataView(bytes.buffer);
  view.setUint32(0, dict.length, true);
  let at = 4;
  for (const d of dict) {
    const length = utf8ByteCount(d);
    view.setUint32(at, length, true);
    utf8Write(bytes, at + 4, d);
    at += 4 + length;
  }
  bytes.set(new Uint8Array(codes.buffer, codes.byteOffset, codes.byteLength), at);
  return bytes;
}

export function encodeMsgPackTree(node: Node): Uint8Array {
  return encodeMsgPack(toPlainNode(node), (value) => (value instanceof Ext ? value.data : undefined));
}

// ---------------------------------------------------------------------------------------------------------------
// Decoding

/** What an extension decoder gives for a payload that does not read; the walk reports it with its path. */
class Malformed {
  constructor(readonly why: string) {}
}

/** Decodes an extension payload; the columns it makes are noted, so that a `bin` is never taken for a column. */
function extensionDecoder(made: WeakSet<object>) {
  return (type: number, data: Uint8Array): unknown => {
    if (type === MsgPackExt.number) {
      return data.length === 1 && data[0] < Specials.length ? Specials[data[0]] : new Malformed('a special number');
    }
    if (type === MsgPackExt.bigint) {
      const digits = utf8Read(data, 0, data.length);
      return /^-?(0|[1-9][0-9]*)$/.test(digits) ? BigInt(digits) : new Malformed('a bigint');
    }
    const numeric = ColumnOfExt.get(type);
    if (numeric !== undefined) {
      const column = numericColumn(numeric, data);
      if (column) made.add(column);
      return column ?? new Malformed(`a ${numeric} column`);
    }
    if (type === MsgPackExt.category) {
      const column = categoryColumnOf(data);
      if (column) made.add(column);
      return column ?? new Malformed('a category column');
    }
    return new Malformed(`an unknown extension (${type})`);
  };
}

/** A view into the payload when it is aligned for the element type, else one copy. */
function numericColumn(type: NumericColumnType, data: Uint8Array): Column | undefined {
  const Ctor = columnConstructor(type);
  const size = Ctor.BYTES_PER_ELEMENT;
  if (data.length % size !== 0) return undefined;
  if (data.byteOffset % size === 0) return new Ctor(data.buffer as ArrayBuffer, data.byteOffset, data.length / size);
  const copy = new Uint8Array(data.length);
  copy.set(data);
  return new Ctor(copy.buffer, 0, data.length / size);
}

function categoryColumnOf(data: Uint8Array): CategoryColumn | undefined {
  if (data.length < 4) return undefined;
  const view = new DataView(data.buffer, data.byteOffset, data.length);
  const count = view.getUint32(0, true);
  if (count > data.length) return undefined;
  const dict: string[] = new Array(count);
  let at = 4;
  for (let i = 0; i < count; i++) {
    if (at + 4 > data.length) return undefined;
    const length = view.getUint32(at, true);
    at += 4;
    if (length > data.length - at) return undefined;
    dict[i] = utf8Read(data, at, length);
    at += length;
  }
  const rest = data.length - at;
  if (rest % 4 !== 0) return undefined;
  const codes = new Int32Array(rest / 4);
  new Uint8Array(codes.buffer).set(data.subarray(at));
  try {
    return categoryFromCodes(dict, codes);
  } catch {
    return undefined;
  }
}

export function decodeMsgPackTree(bytes: Uint8Array, _options?: DecodeOptions): Decoded {
  const made = new WeakSet<object>();
  let value: unknown;
  try {
    const state = createMsgPackState(bytes);
    state.ext = extensionDecoder(made);
    value = parseMsgPackValue(state);
    if (state.offset !== bytes.length) throw new Error(`${bytes.length - state.offset} byte(s) after the value`);
  } catch (e) {
    // the parser recurses: a document nested past the engine's stack is a depth issue
    if (e instanceof RangeError)
      return { issues: [issueOf(TreeIssues, 'repr.tree.depth_exceeded', `Nested deeper than ${MaxDepth} levels`)] };
    return {
      issues: [
        issueOf(MsgPackIssues, 'repr.codec.msgpack.syntax', 'Not MessagePack', { hint: String((e as Error)?.message) }),
      ],
    };
  }
  const issues: Issue[] = [];
  const node = decodeNode(value, [], issues, 1, made);
  return issues.length ? { issues } : { node, issues };
}

type Code = keyof typeof TreeIssues.codes;
function report(issues: Issue[], code: Code, path: TreePath, message: string, hint?: string) {
  issues.push(issueOf(TreeIssues, code, message, { path, hint }));
}

const isMap = (v: unknown): v is Record<string, unknown> =>
  typeof v === 'object' && v !== null && !Array.isArray(v) && !(v instanceof Malformed) && !ArrayBuffer.isView(v);

function tooDeep(depth: number, path: TreePath, issues: Issue[]) {
  if (depth <= MaxDepth) return false;
  if (!issues.some((i) => i.code === 'repr.tree.depth_exceeded'))
    report(issues, 'repr.tree.depth_exceeded', path, `Nested deeper than ${MaxDepth} levels`);
  return true;
}

function decodeNode(value: unknown, path: TreePath, issues: Issue[], depth: number, made: WeakSet<object>): Node {
  if (tooDeep(depth, path, issues)) return { kind: '' };
  if (!isMap(value) || typeof value.kind !== 'string') {
    report(issues, 'repr.tree.node.shape', path, 'A node must be a map with a string `kind`');
    return { kind: '' };
  }
  const { kind, id, params, custom, children, args } = value;
  for (const key of Object.keys(value)) {
    if (!NodeKeys.has(key))
      report(
        issues,
        'repr.tree.node.shape',
        [...path, key],
        'Not a node key',
        'kind, id, params, custom, children, args',
      );
  }
  if (!KindPattern.test(kind)) report(issues, 'repr.tree.node.kind', [...path, 'kind'], 'Not a kind', KindHint);
  const node: {
    kind: string;
    id?: string;
    params?: ValueRecord;
    custom?: ValueRecord;
    children?: Node[];
    args?: Value[];
  } = { kind };
  if (id !== undefined) {
    if (typeof id !== 'string' || !IdPattern.test(id))
      report(issues, 'repr.tree.node.id', [...path, 'id'], 'Not an id', IdHint);
    else node.id = id;
  }
  for (const [key, v] of [
    ['params', params],
    ['custom', custom],
  ] as const) {
    if (v === undefined) continue;
    if (!isMap(v) || Object.hasOwn(v, 'kind'))
      report(issues, 'repr.tree.node.shape', [...path, key], `\`${key}\` must be a record (a map without \`kind\`)`);
    else node[key] = decodeRecord(v, [...path, key], issues, depth, made);
  }
  if (children !== undefined) {
    if (!Array.isArray(children))
      report(issues, 'repr.tree.node.shape', [...path, 'children'], '`children` must be a list of nodes');
    else node.children = children.map((c, i) => decodeNode(c, [...path, 'children', i], issues, depth + 1, made));
  }
  if (args !== undefined) {
    if (!Array.isArray(args) || node.id !== undefined || node.children !== undefined)
      report(
        issues,
        'repr.tree.node.shape',
        [...path, 'args'],
        '`args` must be a list of values, on a node without id or children',
      );
    else if (args.length) node.args = args.map((a, i) => decodeValue(a, [...path, 'args', i], issues, depth + 1, made));
  }
  return node;
}

function decodeRecord(
  value: Record<string, unknown>,
  path: TreePath,
  issues: Issue[],
  depth: number,
  made: WeakSet<object>,
): ValueRecord {
  const out: Record<string, Value> = {};
  for (const key of Object.keys(value)) {
    if (isReservedKey(key))
      report(issues, 'repr.tree.value.reserved_key', [...path, key], 'A reserved key', 'kind, __proto__, $…');
    else out[key] = decodeValue(value[key], [...path, key], issues, depth, made);
  }
  return out;
}

function decodeValue(value: unknown, path: TreePath, issues: Issue[], depth: number, made: WeakSet<object>): Value {
  switch (typeof value) {
    case 'string':
    case 'number':
    case 'bigint':
    case 'boolean':
      return value;
  }
  if (value === null || value === undefined) return null;
  if (value instanceof Malformed) {
    report(issues, 'repr.tree.value.tag', path, `A malformed extension: not ${value.why}`);
    return null;
  }
  if (isColumn(value) && made.has(value)) return value;
  if ((Array.isArray(value) || isMap(value)) && tooDeep(depth + 1, path, issues)) return null;
  if (Array.isArray(value)) return value.map((v, i) => decodeValue(v, [...path, i], issues, depth + 1, made));
  if (isMap(value)) {
    if (Object.hasOwn(value, 'kind')) return decodeNode(value, path, issues, depth + 1, made);
    return decodeRecord(value, path, issues, depth + 1, made);
  }
  report(issues, 'repr.tree.value.shape', path, 'Not a value');
  return null;
}

/** The MessagePack codec: typed columns of every type, whatever the language. */
export const MsgPackCodec: Codec<Uint8Array> = {
  id: 'msgpack',
  output: 'binary',
  typedColumns: true,
  elides: false,
  encode: encodeMsgPackTree,
  decode: decodeMsgPackTree,
};
