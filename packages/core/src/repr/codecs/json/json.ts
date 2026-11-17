/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The JSON encoding of the tree (_spec/README.md "Codecs"; format.md): the default interchange, readable by any tool.
// A node is `{ kind, id?, params?, custom?, children?, args? }`, a list an array, a record an object without `kind`
// (the keys `kind`, `__proto__` and `$…` are reserved); the atoms JSON lacks are tagged objects,
// `{ "$number": "NaN" }` (also `Infinity`, `-Infinity`, `-0`) and `{ "$bigint": "123" }`. A typed column is a plain
// array: a float column's NaN is `null`, its infinities and -0 are tagged as above; a category's values are strings
// and `null` (missing). Decoding is schema-free: a column comes back as a plain array, typed by normalization.
// Parsing is strict and reports every issue it finds, each with its path. Written from meso*'s tree core; `custom` and
// columns are the delta.

import type { Codec, DecodeOptions, Decoded } from '../../codec/codec.js';
import {
  assertKeys,
  assertPrintable,
  CategoryColumn,
  type Column,
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

/** The JSON-ready value of a node (for `JSON.stringify`). */
export function toJson(node: Node): unknown {
  assertPrintable(node);
  const out: Record<string, unknown> = { kind: node.kind };
  if (node.id !== undefined) out.id = node.id;
  if (node.params !== undefined) out.params = recordToJson(node.params);
  if (node.custom !== undefined) out.custom = recordToJson(node.custom);
  if (node.children !== undefined) out.children = node.children.map(toJson);
  if (node.args?.length) out.args = node.args.map(valueToJson);
  return out;
}

function recordToJson(record: ValueRecord): Record<string, unknown> {
  assertKeys(record);
  const out: Record<string, unknown> = {};
  // an `undefined` value is absent, as in the other printer and the hash
  for (const key of Object.keys(record)) if (record[key] !== undefined) out[key] = valueToJson(record[key]);
  return out;
}

function numberToJson(value: number): unknown {
  if (Number.isFinite(value) && !Object.is(value, -0)) return value;
  return { $number: Object.is(value, -0) ? '-0' : String(value) };
}

function columnToJson(column: Column): unknown[] {
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
  // a float column's missing value is NaN, written `null` (format.md)
  for (let i = 0; i < column.length; i++) {
    const v = column[i];
    out[i] = Number.isNaN(v) ? null : numberToJson(v);
  }
  return out;
}

function valueToJson(value: Value): unknown {
  if (typeof value === 'number') return numberToJson(value);
  if (typeof value === 'bigint') return { $bigint: value.toString() };
  if (value === null || typeof value !== 'object') return value;
  if (isColumn(value)) return columnToJson(value);
  if (Array.isArray(value)) return value.map(valueToJson);
  if (isNode(value)) return toJson(value);
  return recordToJson(value as ValueRecord);
}

/** The tree as JSON text. */
export function printJson(node: Node, indent?: number): string {
  return JSON.stringify(toJson(node), null, indent);
}

/** A tree from JSON text; a leading byte order mark is skipped (offsets count it). */
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
        {
          code: 'repr.codec.json.syntax',
          severity: 'error',
          message: 'Not JSON',
          ...(offset !== undefined ? { offset: +offset + bom } : {}),
        },
      ],
    };
  }
  return fromJson(value);
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

/** The JSON codec: compact text (`indent` for people). */
export function jsonCodec(options: { readonly indent?: number } = {}): Codec<string> {
  return {
    id: 'json',
    output: 'text',
    typedColumns: false,
    elides: false,
    encode: (node) => printJson(node, options.indent),
    decode: (text: string, _options?: DecodeOptions) => parseJson(text),
  };
}

export const JsonCodec = jsonCodec();
