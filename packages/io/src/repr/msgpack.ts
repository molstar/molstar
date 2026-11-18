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
import { utf8ByteCount, utf8Write } from '@molstar/core/util/utf8';
import {
  createMsgPackState,
  type MsgPackState,
  parseMsgPackValue,
  readMsgPackArrayLength,
  readMsgPackMapLength,
  skipMsgPackValue,
} from '../common/msgpack/decode.js';
import { encodeMsgPack, type MsgPackExtData } from '../common/msgpack/encode.js';

const HostLittleEndian = new Uint8Array(new Uint16Array([1]).buffer)[0] === 1;

/** A copy of `bytes` with each `size`-byte element reversed: a column between the wire (little-endian) and a big-endian host. */
export function swapBytes(bytes: Uint8Array, size: number): Uint8Array {
  const out = new Uint8Array(bytes.length);
  for (let i = 0; i + size <= bytes.length; i += size)
    for (let k = 0; k < size; k++) out[i + k] = bytes[i + size - 1 - k];
  return out;
}

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
  'repr.codec.msgpack.duplicate_key': {
    severity: 'error',
    description: 'A map names the same key twice (the MessagePack decoder would keep the last silently)',
  },
  'repr.codec.msgpack.key_type': {
    severity: 'error',
    description: 'A map has a key that is not a string (the decoder would make it one)',
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

/** Throws on a lone surrogate: UTF-8 cannot hold it, and the string would not come back (repr.md §3). */
function assertWellFormed(text: string): string {
  if (!(text as unknown as { isWellFormed(): boolean }).isWellFormed())
    throw new Error('A string holds a lone surrogate, which UTF-8 cannot encode');
  return text;
}

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
  for (const key of Object.keys(record))
    if (record[key] !== undefined) out[assertWellFormed(key)] = toPlainValue(record[key]);
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
  if (typeof value === 'string') return assertWellFormed(value);
  if (value === null || typeof value !== 'object') return value;
  if (isColumn(value)) return columnExt(value);
  if (Array.isArray(value)) return value.map(toPlainValue);
  if (isNode(value)) return toPlainNode(value);
  return toPlainRecord(value as ValueRecord);
}

function columnExt(column: Column): Ext {
  canonicalColumn(column);
  if (column instanceof CategoryColumn) return new Ext({ type: MsgPackExt.category, data: categoryBytes(column) });
  const view = new Uint8Array(column.buffer, column.byteOffset, column.byteLength);
  const data = HostLittleEndian ? view : swapBytes(view, column.BYTES_PER_ELEMENT);
  return new Ext({ type: ColumnExt[columnType(column) as NumericColumnType], data });
}

function categoryBytes(column: CategoryColumn): Uint8Array {
  const { dict, codes } = column;
  let size = 4;
  for (const d of dict) size += 4 + utf8ByteCount(assertWellFormed(d));
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
  const codeBytes = new Uint8Array(codes.buffer, codes.byteOffset, codes.byteLength);
  bytes.set(HostLittleEndian ? codeBytes : swapBytes(codeBytes, 4), at);
  return bytes;
}

export function encodeMsgPackTree(node: Node): Uint8Array {
  return encodeMsgPack(toPlainNode(node), (value) => (value instanceof Ext ? value.data : undefined));
}

// ---------------------------------------------------------------------------------------------------------------
// Decoding

/** Strict UTF-8 that keeps a leading U+FEFF (a value, not a byte order mark): invalid bytes throw, not U+FFFD. */
const strictUtf8 = new TextDecoder('utf-8', { fatal: true, ignoreBOM: true });

/** The text of an extension's bytes; undefined when they are not UTF-8. */
function utf8Strict(data: Uint8Array, offset: number, length: number): string | undefined {
  // short ASCII, the usual entry, without a decoder call
  if (length <= 24) {
    let text = '';
    let i = 0;
    for (; i < length; i++) {
      const c = data[offset + i];
      if (c >= 0x80) break;
      text += String.fromCharCode(c);
    }
    if (i === length) return text;
  }
  try {
    return strictUtf8.decode(data.subarray(offset, offset + length));
  } catch {
    return undefined;
  }
}

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
      const digits = utf8Strict(data, 0, data.length);
      return digits !== undefined && /^-?(0|[1-9][0-9]*)$/.test(digits) ? BigInt(digits) : new Malformed('a bigint');
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

/** Whether a float view holds a NaN: one may carry a payload, which making it canonical would write into the input. */
function hasNaN(view: Float32Array | Float64Array): boolean {
  for (let i = 0; i < view.length; i++) if (Number.isNaN(view[i])) return true;
  return false;
}

/**
 * A view into the payload when it is aligned for the element type, the host is little-endian and a float column holds
 * no NaN (the tree writes the canonical NaN into a column it owns: never into the input); else one copy.
 */
function numericColumn(type: NumericColumnType, data: Uint8Array): Column | undefined {
  const Ctor = columnConstructor(type);
  const size = Ctor.BYTES_PER_ELEMENT;
  if (data.length % size !== 0) return undefined;
  if (HostLittleEndian && data.byteOffset % size === 0) {
    const view = new Ctor(data.buffer as ArrayBuffer, data.byteOffset, data.length / size);
    if (!(view instanceof Float32Array || view instanceof Float64Array) || !hasNaN(view)) return view;
  }
  const copy = new Uint8Array(data.length);
  copy.set(data);
  return new Ctor((HostLittleEndian ? copy : swapBytes(copy, size)).buffer as ArrayBuffer, 0, data.length / size);
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
    const entry = utf8Strict(data, at, length);
    if (entry === undefined) return undefined;
    dict[i] = entry;
    at += length;
  }
  const rest = data.length - at;
  if (rest % 4 !== 0) return undefined;
  const codes = new Int32Array(rest / 4);
  const raw = data.subarray(at);
  new Uint8Array(codes.buffer).set(HostLittleEndian ? raw : swapBytes(raw, 4));
  try {
    return categoryFromCodes(dict, codes);
  } catch {
    return undefined;
  }
}

export function decodeMsgPackTree(bytes: Uint8Array, _options?: DecodeOptions): Decoded {
  const made = new WeakSet<object>();
  let value: unknown;
  const state = createMsgPackState(bytes);
  try {
    state.ext = extensionDecoder(made);
    value = parseMsgPackValue(state);
    if (state.offset !== bytes.length) throw new Error(`${bytes.length - state.offset} byte(s) after the value`);
  } catch (e) {
    // the parser recurses: only a stack overflow is a depth issue; a read past the end is a RangeError too
    if (e instanceof RangeError && /call stack/i.test(e.message))
      return { issues: [issueOf(TreeIssues, 'repr.tree.depth_exceeded', `Nested deeper than ${MaxDepth} levels`)] };
    const offset = Math.min(state.offset, bytes.length);
    return {
      issues: [
        issueOf(MsgPackIssues, 'repr.codec.msgpack.syntax', 'Not MessagePack', {
          offset,
          hint: String((e as Error)?.message),
        }),
      ],
    };
  }
  const keyIssues = checkKeys(bytes);
  if (keyIssues.length) return { issues: keyIssues };
  const issues: Issue[] = [];
  const node = decodeNode(value, [], issues, 1, made);
  return issues.length ? { issues } : { node, issues };
}

/**
 * A second pass over a document that read: every map's keys are strings and given once. The decoder keeps the last of
 * a key and stringifies the others, so it cannot say; extension payloads are skipped as they are. Each issue has the
 * offset of the key and the path of its map.
 */
function checkKeys(bytes: Uint8Array): Issue[] {
  const issues: Issue[] = [];
  const state: MsgPackState = createMsgPackState(bytes);
  const path: Stack = [];
  const scan = () => {
    const type = state.buffer[state.offset];
    if ((type & 0xf0) === 0x90 || type === 0xdc || type === 0xdd) {
      const length = readMsgPackArrayLength(state);
      for (let i = 0; i < length; i++) {
        path.push(i);
        scan();
        path.pop();
      }
    } else if ((type & 0xf0) === 0x80 || type === 0xde || type === 0xdf) {
      const length = readMsgPackMapLength(state);
      const keys: string[] = [];
      const set = length > 8 ? new Set<string>() : undefined;
      for (let i = 0; i < length; i++) {
        const at = state.offset;
        const t = state.buffer[at];
        if (!((t & 0xe0) === 0xa0 || t === 0xd9 || t === 0xda || t === 0xdb)) {
          issues.push(
            issueOf(MsgPackIssues, 'repr.codec.msgpack.key_type', 'A key that is not a string', {
              offset: at,
              path: [...path],
            }),
          );
          skipMsgPackValue(state);
          path.push(String(i));
          scan();
          path.pop();
          continue;
        }
        const key = parseMsgPackValue(state) as string;
        if (set ? set.has(key) : keys.includes(key))
          issues.push(
            issueOf(MsgPackIssues, 'repr.codec.msgpack.duplicate_key', 'A key given twice', {
              offset: at,
              path: [...path],
            }),
          );
        if (set) set.add(key);
        else keys.push(key);
        path.push(key);
        scan();
        path.pop();
      }
    } else skipMsgPackValue(state);
  };
  try {
    scan();
  } catch {
    // the document read once: a stack too shallow for a second walk reports nothing here
  }
  return issues;
}

type Stack = (string | number)[];
type Code = keyof typeof TreeIssues.codes;
function report(issues: Issue[], code: Code, path: TreePath, message: string, hint?: string) {
  issues.push(issueOf(TreeIssues, code, message, { path: [...path], hint }));
}

const isMap = (v: unknown): v is Record<string, unknown> =>
  typeof v === 'object' && v !== null && !Array.isArray(v) && !(v instanceof Malformed) && !ArrayBuffer.isView(v);

function tooDeep(depth: number, path: Stack, issues: Issue[]) {
  if (depth <= MaxDepth) return false;
  if (!issues.some((i) => i.code === 'repr.tree.depth_exceeded'))
    report(issues, 'repr.tree.depth_exceeded', path, `Nested deeper than ${MaxDepth} levels`);
  return true;
}

function decodeNode(value: unknown, path: Stack, issues: Issue[], depth: number, made: WeakSet<object>): Node {
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
    else {
      path.push(key);
      node[key] = decodeRecord(v, path, issues, depth, made);
      path.pop();
    }
  }
  if (children !== undefined) {
    if (!Array.isArray(children))
      report(issues, 'repr.tree.node.shape', [...path, 'children'], '`children` must be a list of nodes');
    else {
      path.push('children', 0);
      const out: Node[] = new Array(children.length);
      for (let i = 0; i < children.length; i++) {
        path[path.length - 1] = i;
        out[i] = decodeNode(children[i], path, issues, depth + 1, made);
      }
      path.length -= 2;
      node.children = out;
    }
  }
  if (args !== undefined) {
    if (!Array.isArray(args) || node.id !== undefined || node.children !== undefined)
      report(
        issues,
        'repr.tree.node.shape',
        [...path, 'args'],
        '`args` must be a list of values, on a node without id or children',
      );
    else if (args.length) {
      path.push('args', 0);
      const out: Value[] = new Array(args.length);
      for (let i = 0; i < args.length; i++) {
        path[path.length - 1] = i;
        out[i] = decodeValue(args[i], path, issues, depth + 1, made);
      }
      path.length -= 2;
      node.args = out;
    }
  }
  return node;
}

function decodeRecord(
  value: Record<string, unknown>,
  path: Stack,
  issues: Issue[],
  depth: number,
  made: WeakSet<object>,
): ValueRecord {
  const out: Record<string, Value> = {};
  for (const key of Object.keys(value)) {
    if (isReservedKey(key))
      report(issues, 'repr.tree.value.reserved_key', [...path, key], 'A reserved key', 'kind, __proto__, $…');
    else {
      path.push(key);
      out[key] = decodeValue(value[key], path, issues, depth, made);
      path.pop();
    }
  }
  return out;
}

function decodeValue(value: unknown, path: Stack, issues: Issue[], depth: number, made: WeakSet<object>): Value {
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
  if (Array.isArray(value)) {
    path.push(0);
    const out: Value[] = new Array(value.length);
    for (let i = 0; i < value.length; i++) {
      path[path.length - 1] = i;
      out[i] = decodeValue(value[i], path, issues, depth + 1, made);
    }
    path.pop();
    return out;
  }
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
