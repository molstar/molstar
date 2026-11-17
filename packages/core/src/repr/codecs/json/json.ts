/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The JSON encoding of the tree (_spec/README.md "Codecs"; format.md): the default interchange, readable by any tool.
// A node is `{ kind, id?, params?, custom?, children?, args? }`, a list an array, a record an object without `kind`
// (the keys `kind`, `__proto__` and `$…` are reserved); the atoms JSON lacks are tagged objects,
// `{ "$number": "NaN" }` (also `Infinity`, `-Infinity`, `-0`) and `{ "$bigint": "123" }`. A typed column in a declared
// column param (the encoder is given the language) is a plain array: a float column's NaN is `null`, its infinities
// and -0 tagged as above; a category's values are strings and `null` (missing); decoding gives it back plain, and
// normalization types it. Any other column is tagged, `{ "$column": "int32", "values": […] }` (a category:
// `{ "$column": "category", "dict": […], "codes": […] }`), and decodes typed with no language. Parsing is strict
// (a duplicate key is an issue, not last-wins) and reports every issue it finds, each with its path. Written from meso*'s tree core; `custom` and
// columns are the delta.

import { type Codec, CodecIssues, type DecodeOptions, type Decoded } from '../../codec/codec.js';
import { type ColumnPlaces, declaredColumnPlaces, type Language } from '../../schema/index.js';
import {
  assertKeys,
  assertPrintable,
  CategoryColumn,
  type Column,
  type ColumnType,
  ColumnTypes,
  categoryFromCodes,
  columnType,
  IdHint,
  IdPattern,
  type Issue,
  isColumn,
  isFloatType,
  isNode,
  isReservedKey,
  issueOf,
  KindHint,
  KindPattern,
  MaxDepth,
  type Node,
  TreeIssues,
  type TreePath,
  toColumn,
  type Value,
  type ValueRecord,
} from '../../tree/index.js';

const NodeKeys = new Set(['kind', 'id', 'params', 'custom', 'children', 'args']);
const SpecialNumbers: Readonly<Record<string, number>> = {
  NaN: Number.NaN,
  Infinity: Number.POSITIVE_INFINITY,
  '-Infinity': Number.NEGATIVE_INFINITY,
  '-0': -0,
};

export interface JsonOptions {
  /** Indentation for people; compact by default. */
  readonly indent?: number;
  /**
   * The language whose declared column params are written as plain arrays; every other column is tagged
   * (`{ "$column": "int32", "values": […] }`). Without one, every column is tagged.
   */
  readonly language?: Language;
  /** The languages its slots name. */
  readonly languages?: (id: string) => Language | undefined;
}

/** The JSON-ready value of a node (for `JSON.stringify`). */
export function toJson(node: Node, options: JsonOptions = {}): unknown {
  const places = options.language
    ? declaredColumnPlaces(options.language, node, { languages: options.languages })
    : undefined;
  return nodeToJson(node, places);
}

function nodeToJson(node: Node, places: ColumnPlaces | undefined): unknown {
  assertPrintable(node);
  const out: Record<string, unknown> = { kind: node.kind };
  if (node.id !== undefined) out.id = node.id;
  if (node.params !== undefined) out.params = recordToJson(node.params, places);
  if (node.custom !== undefined) out.custom = recordToJson(node.custom, places);
  if (node.children !== undefined) out.children = node.children.map((c) => nodeToJson(c, places));
  if (node.args?.length) out.args = listToJson(node.args, places);
  return out;
}

function recordToJson(record: ValueRecord, places: ColumnPlaces | undefined): Record<string, unknown> {
  assertKeys(record);
  const declared = places?.get(record);
  const out: Record<string, unknown> = {};
  // an `undefined` value is absent, as in the other printer and the hash
  for (const key of Object.keys(record))
    if (record[key] !== undefined) out[key] = valueToJson(record[key], places, declared?.has(key) ?? false);
  return out;
}

function listToJson(list: readonly Value[], places: ColumnPlaces | undefined): unknown[] {
  const declared = places?.get(list);
  return list.map((v, i) => valueToJson(v, places, declared?.has(i) ?? false));
}

function numberToJson(value: number): unknown {
  if (Number.isFinite(value) && !Object.is(value, -0)) return value;
  return { $number: Object.is(value, -0) ? '-0' : String(value) };
}

/** A column's values: a float column's NaN `null` (format.md), a category's strings and `null`. */
function columnValuesToJson(column: Column): unknown[] {
  if (column instanceof CategoryColumn) {
    const out: (string | null)[] = new Array(column.length);
    for (let i = 0; i < column.length; i++) out[i] = column.get(i);
    return out;
  }
  const out: unknown[] = new Array(column.length);
  if (!isFloatType(columnType(column)!)) {
    for (let i = 0; i < column.length; i++) out[i] = column[i];
    return out;
  }
  for (let i = 0; i < column.length; i++) {
    const v = column[i];
    out[i] = Number.isNaN(v) ? null : numberToJson(v);
  }
  return out;
}

/** A column outside a declared column param, tagged so that it decodes typed with no language. */
function taggedColumn(column: Column): unknown {
  if (column instanceof CategoryColumn)
    return { $column: 'category', dict: column.dict, codes: Array.from(column.codes) };
  return { $column: columnType(column), values: columnValuesToJson(column) };
}

function valueToJson(value: Value, places: ColumnPlaces | undefined, declared: boolean): unknown {
  if (typeof value === 'number') return numberToJson(value);
  if (typeof value === 'bigint') return { $bigint: value.toString() };
  if (value === null || typeof value !== 'object') return value;
  if (isColumn(value)) return declared ? columnValuesToJson(value) : taggedColumn(value);
  if (Array.isArray(value)) return listToJson(value, places);
  if (isNode(value)) return nodeToJson(value, places);
  return recordToJson(value as ValueRecord, places);
}

/** The tree as JSON text. */
export function printJson(node: Node, options: JsonOptions = {}): string {
  return JSON.stringify(toJson(node, options), null, options.indent);
}

/** A tree from JSON text, strictly: a key given twice in an object is an issue; a byte order mark is skipped. */
export function parseJson(text: string): Decoded {
  const bom = text.charCodeAt(0) === 0xfeff ? 1 : 0;
  let value: unknown;
  try {
    value = JSON.parse(bom ? text.slice(1) : text);
  } catch (e) {
    // the engine's message quotes the text: only its position is kept
    const offset = /position (\d+)/.exec(String(e))?.[1];
    return {
      issues: [
        issueOf(CodecIssues, 'repr.codec.json.syntax', 'Not JSON', {
          offset: offset !== undefined ? +offset + bom : undefined,
        }),
      ],
    };
  }
  // JSON.parse keeps the last of a duplicate key: a second scan finds them (valid JSON only, so it can be simple)
  const duplicates = duplicateKeys(text, bom);
  if (duplicates.length) return { issues: duplicates };
  return fromJson(value);
}

/** Every key given twice in one object of valid JSON text, with its offset and path. */
function duplicateKeys(text: string, start: number): Issue[] {
  const issues: Issue[] = [];
  type Frame = { readonly object: boolean; keys?: Set<string>; key: string | number };
  const stack: Frame[] = [];
  const path = () => stack.map((f) => f.key);
  for (let i = start; i < text.length; i++) {
    const c = text.charCodeAt(i);
    if (c === 0x22 /* " */) {
      const from = i;
      let escaped = false;
      for (i++; text.charCodeAt(i) !== 0x22; i++) {
        if (text.charCodeAt(i) === 0x5c /* \ */) {
          escaped = true;
          i++;
        }
      }
      const top = stack[stack.length - 1];
      if (!top?.object) continue;
      let j = i + 1;
      while (text.charCodeAt(j) <= 0x20) j++;
      if (text.charCodeAt(j) !== 0x3a /* : */) continue;
      const key = escaped ? (JSON.parse(text.slice(from, i + 1)) as string) : text.slice(from + 1, i);
      top.key = key;
      if (top.keys!.has(key))
        issues.push(
          issueOf(CodecIssues, 'repr.codec.json.duplicate_key', 'A key given twice', { offset: from, path: path() }),
        );
      top.keys!.add(key);
    } else if (c === 0x7b /* { */) stack.push({ object: true, keys: new Set(), key: '' });
    else if (c === 0x5b /* [ */) stack.push({ object: false, key: 0 });
    else if (c === 0x7d || c === 0x5d) stack.pop();
    else if (c === 0x2c /* , */) {
      const top = stack[stack.length - 1];
      if (top && !top.object) top.key = (top.key as number) + 1;
    }
  }
  return issues;
}

/** A tree from a JSON value (as `JSON.parse` gives it), strictly. */
export function fromJson(value: unknown): Decoded {
  const issues: Issue[] = [];
  const node = decodeNode(value, [], issues, 1);
  return issues.length ? { issues } : { node, issues };
}

type Code = keyof typeof TreeIssues.codes;
function report(issues: Issue[], code: Code, path: TreePath, message: string, hint?: string) {
  issues.push(issueOf(TreeIssues, code, message, { path, hint }));
}

const isObject = (v: unknown): v is Record<string, unknown> => typeof v === 'object' && v !== null && !Array.isArray(v);

/** Depth counts nodes, lists and records, as the S-expression reader's brackets (MaxDepth). */
function tooDeep(depth: number, path: TreePath, issues: Issue[]) {
  if (depth <= MaxDepth) return false;
  if (!issues.some((i) => i.code === 'repr.tree.depth_exceeded'))
    report(issues, 'repr.tree.depth_exceeded', path, `Nested deeper than ${MaxDepth} levels`);
  return true;
}

function decodeNode(value: unknown, path: TreePath, issues: Issue[], depth: number): Node {
  if (tooDeep(depth, path, issues)) return { kind: '' };
  if (!isObject(value) || typeof value.kind !== 'string') {
    report(issues, 'repr.tree.node.shape', path, 'A node must be an object with a string `kind`');
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
    if (!isObject(v) || Object.hasOwn(v, 'kind'))
      report(
        issues,
        'repr.tree.node.shape',
        [...path, key],
        `\`${key}\` must be a record (an object without \`kind\`)`,
      );
    else node[key] = decodeRecord(v, [...path, key], issues, depth);
  }
  if (children !== undefined) {
    if (!Array.isArray(children))
      report(issues, 'repr.tree.node.shape', [...path, 'children'], '`children` must be a list of nodes');
    else node.children = children.map((c, i) => decodeNode(c, [...path, 'children', i], issues, depth + 1));
  }
  if (args !== undefined) {
    if (!Array.isArray(args) || node.id !== undefined || node.children !== undefined)
      report(
        issues,
        'repr.tree.node.shape',
        [...path, 'args'],
        '`args` must be a list of values, on a node without id or children',
      );
    else if (args.length) node.args = args.map((a, i) => decodeValue(a, [...path, 'args', i], issues, depth + 1));
  }
  return node;
}

function decodeRecord(value: Record<string, unknown>, path: TreePath, issues: Issue[], depth: number): ValueRecord {
  const out: Record<string, Value> = {};
  for (const key of Object.keys(value)) {
    if (isReservedKey(key))
      report(issues, 'repr.tree.value.reserved_key', [...path, key], 'A reserved key', 'kind, __proto__, $…');
    else out[key] = decodeValue(value[key], [...path, key], issues, depth);
  }
  return out;
}

function decodeValue(value: unknown, path: TreePath, issues: Issue[], depth: number): Value {
  switch (typeof value) {
    case 'string':
    case 'number':
    case 'boolean':
      return value;
  }
  if (value === null) return null;
  if ((Array.isArray(value) || isObject(value)) && tooDeep(depth + 1, path, issues)) return null;
  if (Array.isArray(value)) return value.map((v, i) => decodeValue(v, [...path, i], issues, depth + 1));
  if (isObject(value)) {
    if (Object.hasOwn(value, 'kind')) return decodeNode(value, path, issues, depth + 1);
    const keys = Object.keys(value);
    if (keys.length === 1 && (keys[0] === '$number' || keys[0] === '$bigint'))
      return decodeTag(keys[0], value, path, issues);
    if (Object.hasOwn(value, '$column')) return decodeColumn(value, path, issues);
    return decodeRecord(value, path, issues, depth + 1);
  }
  report(issues, 'repr.tree.value.shape', path, 'Not a value');
  return null;
}

function decodeTag(tag: '$number' | '$bigint', value: Record<string, unknown>, path: TreePath, issues: Issue[]): Value {
  const payload = value[tag];
  if (tag === '$number' && typeof payload === 'string' && Object.hasOwn(SpecialNumbers, payload))
    return SpecialNumbers[payload];
  if (tag === '$bigint' && typeof payload === 'string' && /^-?(0|[1-9][0-9]*)$/.test(payload)) return BigInt(payload);
  const hint = tag === '$number' ? 'NaN, Infinity, -Infinity or -0' : 'an integer in decimal, as a string';
  report(issues, 'repr.tree.value.tag', path, 'A malformed tagged atom', hint);
  return null;
}

/**
 * A tagged column, `{ "$column": type, "values": […] }` or `{ "$column": "category", "dict": […], "codes": […] }`,
 * typed here; one that does not read is `repr.tree.value.tag`.
 */
function decodeColumn(value: Record<string, unknown>, path: TreePath, issues: Issue[]): Value {
  const type = value.$column;
  const keys = Object.keys(value).sort().join(',');
  let column: Column | undefined;
  if (type === 'category') {
    const { dict, codes } = value;
    const fits =
      keys === '$column,codes,dict' &&
      Array.isArray(dict) &&
      dict.every((d) => typeof d === 'string') &&
      Array.isArray(codes) &&
      codes.every((c) => Number.isInteger(c) && c >= -1 && c < dict.length);
    if (fits) column = categoryFromCodes(dict as string[], codes as number[]);
  } else if (typeof type === 'string' && ColumnTypes.includes(type as ColumnType) && keys === '$column,values') {
    const values = value.values;
    const numbers = Array.isArray(values)
      ? values.map((v) =>
          isObject(v) && Object.keys(v).length === 1 && typeof v.$number === 'string'
            ? (SpecialNumbers[v.$number] ?? v)
            : v,
        )
      : undefined;
    column = toColumn(type as ColumnType, numbers);
  }
  if (column) return column;
  report(issues, 'repr.tree.value.tag', path, 'A malformed tagged column', '{ "$column": type, "values": […] }');
  return null;
}

/** The JSON codec: compact text (`indent` for people); `language` writes its declared columns as plain arrays. */
export function jsonCodec(options: JsonOptions = {}): Codec<string> {
  return {
    id: 'json',
    output: 'text',
    typedColumns: true,
    elides: false,
    encode: (node) => printJson(node, options),
    decode: (text: string, _options?: DecodeOptions) => parseJson(text),
  };
}

export const JsonCodec = jsonCodec();
