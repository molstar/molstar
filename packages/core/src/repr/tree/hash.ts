/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The structural hash and equality (_spec/README.md "Identity"): a node's hash covers its kind, id, params, custom
// data and its children's hashes, over a canonical text (record keys sorted, every atom kind told apart), so it is
// stable across key order and across encodings: hashes are over the tree, never over an encoding. A typed column
// enters as its own hash (its element type and bytes, column.ts), computed once per column. Cached per immutable node,
// so an unchanged subtree costs nothing to hash again (Merkle). The hash is @molstar/core's MurmurHash3 x86 128-bit
// over the text's UTF-8 bytes. Absent `params`, `custom` and `children` hash as empty ones (absent or empty `args` as
// nothing): `{ kind }` and `{ kind, params: {} }` are the same node.

import { murmurHash3_128, murmurHash3_128_fromBytes } from '../../data/util/hash-functions.js';
import { type Column, columnsEqual, columnType, hashColumn, isColumn } from './column.js';
import { isNode, type Node, type Value, type ValueRecord } from './node.js';

const Seed = 0x6d65736f; // pinned: a stored hash depends on it (the pinned test)
const Cache = new WeakMap<Node, string>();

/** 32 hex digits; equal for equal trees, whatever the order of record keys. */
export function hashNode(node: Node): string {
  let hash = Cache.get(node);
  if (hash === undefined) {
    hash = murmurHash3_128(canonicalNode(node), Seed);
    Cache.set(node, hash);
  }
  return hash;
}

/** A content hash of bytes (32 hex digits; the same function and seed as `hashNode`). */
export function hashBytes(bytes: Uint8Array): string {
  return murmurHash3_128_fromBytes(bytes, Seed);
}

/** Whether two trees are equal: by their hashes (cached, so a second comparison of the same nodes is free). */
export function nodesEqual(a: Node, b: Node): boolean {
  return a === b || hashNode(a) === hashNode(b);
}

/**
 * Whether two values are equal: atoms by `Object.is` (`0` is not `-0`, NaN is NaN, `1` is not `1n`), columns by
 * element type and bytes, lists item by item, records key by key (order aside, `undefined` absent), nodes by hash.
 */
export function valuesEqual(a: Value | undefined, b: Value | undefined): boolean {
  if (Object.is(a, b)) return true;
  if (typeof a !== 'object' || typeof b !== 'object' || a === null || b === null) return false;
  if (isColumn(a) || isColumn(b)) return isColumn(a) && isColumn(b) && columnsEqual(a, b);
  if (Array.isArray(a) || Array.isArray(b)) {
    if (!Array.isArray(a) || !Array.isArray(b) || a.length !== b.length) return false;
    for (let i = 0; i < a.length; i++) if (!valuesEqual(a[i], b[i])) return false;
    return true;
  }
  if (isNode(a) || isNode(b)) return isNode(a) && isNode(b) && nodesEqual(a, b);
  const ra = a as ValueRecord;
  const rb = b as ValueRecord;
  const ka = Object.keys(ra).filter((k) => ra[k] !== undefined);
  if (ka.length !== Object.keys(rb).filter((k) => rb[k] !== undefined).length) return false;
  for (const k of ka) if (!Object.hasOwn(rb, k) || !valuesEqual(ra[k], rb[k])) return false;
  return true;
}

/** The node's canonical text, with its child nodes (in children and in params) as their hashes. */
function canonicalNode(node: Node): string {
  const out = ['(', JSON.stringify(node.kind), node.id === undefined ? '_' : JSON.stringify(node.id)];
  out.push(node.params ? canonicalRecord(node.params) : '{}', '[');
  if (node.children) for (const child of node.children) out.push(hashNode(child), ',');
  out.push(']');
  // a data form's args after its children, only when present: every other node's text is unchanged
  if (node.args?.length) out.push(canonicalValue(node.args));
  // custom data last, only when present: a tree without it hashes as it did before custom data existed
  if (node.custom && Object.keys(node.custom).some((k) => node.custom![k] !== undefined))
    out.push('$', canonicalRecord(node.custom));
  out.push(')');
  return out.join('');
}

function canonicalRecord(record: ValueRecord): string {
  // an `undefined` value is absent, as in the printers
  const keys = Object.keys(record)
    .filter((k) => record[k] !== undefined)
    .sort();
  let s = '{';
  for (const key of keys) s += `${JSON.stringify(key)}:${canonicalValue(record[key])},`;
  return `${s}}`;
}

function canonicalValue(value: Value): string {
  switch (typeof value) {
    case 'string':
      return JSON.stringify(value);
    case 'number':
      // bare numbers; String() spells NaN and the infinities, -0 needs its own
      return Object.is(value, -0) ? '-0' : String(value);
    case 'bigint':
      return `${value}n`;
    case 'boolean':
      return value ? 't' : 'f';
  }
  if (value === null) return 'n';
  if (isColumn(value)) return `%${columnType(value)}:${hashColumn(value as Column)}`;
  if (Array.isArray(value)) return `[${value.map(canonicalValue).join(',')}]`;
  if (isNode(value)) return `#${hashNode(value)}`;
  return canonicalRecord(value as ValueRecord);
}
