/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The node (_spec/README.md "The tree"): the one recursive form every language on `repr` is spelled in. A node is a
// `kind`, an optional local `id`, named `params`, extension data in `custom` and ordered `children`. A param value is
// an atom, a list, a record, a node (a slot: how a language embeds another) or a typed column. A *data form* (a query,
// a theme expression) holds positional values in `args` instead of an id and children; which kinds are data forms is
// the language's declaration (schema `KindSpec.data`), passed to the readers that need it.

import { type Column, isColumn } from './column.js';

/** Lossless atoms: float64 numbers including -0, NaN and the infinities, and bigints. */
export type Atom = string | number | bigint | boolean | null;

export type Value = Atom | Column | readonly Value[] | ValueRecord | Node;

/** A record: named values, no `kind` (an object with a `kind` is a node), no key starting with `$`. */
export interface ValueRecord {
  readonly [key: string]: Value;
}

export type Params = Readonly<Record<string, Value>>;

export interface Node {
  readonly kind: string;
  /** Local: unique among its siblings (IdPattern); the node's address is the path of ids from the root. */
  readonly id?: string;
  readonly params?: Params;
  /** Extension data: kept, hashed and encoded with the node, never validated against its kind. */
  readonly custom?: Params;
  readonly children?: readonly Node[];
  /** A data form's positional values, in order; a node with `args` has no id and no children. */
  readonly args?: readonly Value[];
}

/**
 * A kind is a word, optionally in a namespace (`ext:track`), so it is a symbol in the S-expression spelling. The
 * namespace is how an extension's kinds stay apart from the language's and from each other's.
 */
export const KindPattern = /^([a-z][a-z0-9_-]*:)?[a-z][a-z0-9_-]*$/;
export const KindHint = 'a lowercase letter, then a-z, 0-9, _ or -; an extension kind after its namespace and `:`';

/** The namespace of a kind (`ext` of `ext:track`), or undefined. */
export function kindNamespace(kind: string): string | undefined {
  const i = kind.indexOf(':');
  return i < 0 ? undefined : kind.slice(0, i);
}

/**
 * An id is local: unique among its siblings, and a node's address is the path of ids from the root (address.ts). So an
 * id is one address segment: letters, digits, `_`, `.`, `+`, `-`, opening with a letter, a digit or `_` (no separator,
 * no empty segment, nothing an address would have to quote).
 */
export const IdPattern = /^[\p{L}\p{N}_][\p{L}\p{N}_.+-]*$/u;
export const IdHint = 'a letter, digit or _, then letters, digits, _, ., + or -';

/**
 * Keys no param or record may have: `kind` (an object with a `kind` is a node), `$…` (the encodings' tagged atoms
 * and reserved spellings) and `__proto__` (an object's prototype in JavaScript, never an own key once assigned).
 */
export function isReservedKey(key: string): boolean {
  return key === 'kind' || key === '__proto__' || key.startsWith('$');
}

/** How deep a tree the readers build: nodes, lists and records each count one level (repr.tree.depth_exceeded). */
export const MaxDepth = 512;

/** `record[key]` when `key` is its own (never a prototype's: `toString`, `constructor`). */
export function own<T>(record: Readonly<Record<string, T>>, key: string): T | undefined {
  return Object.hasOwn(record, key) ? record[key] : undefined;
}

/** Printing a tree no reader would give back is a programmer error: this throws. */
export function assertPrintable(node: Node): void {
  if (!KindPattern.test(node.kind)) throw new Error('A node kind is not a word');
  if (node.id !== undefined && !IdPattern.test(node.id)) throw new Error('A node id is not an id');
  if (node.params) assertKeys(node.params);
  if (node.custom) assertKeys(node.custom);
  if (node.args !== undefined && (node.id !== undefined || node.children !== undefined))
    throw new Error('A data form (a node with args) has no id and no children');
}

export function assertKeys(record: ValueRecord): void {
  for (const key of Object.keys(record)) if (isReservedKey(key)) throw new Error('A param or record key is reserved');
}

/** Whether a value is an object that is neither a list nor a column (a node or a record). */
function isPlainObject(value: unknown): value is object {
  return typeof value === 'object' && value !== null && !Array.isArray(value) && !isColumn(value);
}

export function isNode(value: unknown): value is Node {
  return isPlainObject(value) && typeof (value as Node).kind === 'string';
}

export function isRecord(value: unknown): value is ValueRecord {
  return isPlainObject(value) && !('kind' in value);
}

/** Depth-first, parents before children; `visit` returning `false` skips the node's children. */
export function visitTree(root: Node, visit: (node: Node, parent: Node | undefined) => boolean | undefined) {
  const go = (node: Node, parent: Node | undefined) => {
    if (visit(node, parent) === false || !node.children) return;
    for (const child of node.children) go(child, node);
  };
  go(root, undefined);
}
