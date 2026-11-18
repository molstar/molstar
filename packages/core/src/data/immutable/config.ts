/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { defaultHash, lookupDefaultHash } from './hash.js';

/** A collection's hash and equality, given at creation (README §4). */
export interface CollectionConfig<K, V = unknown> {
  /** The key hash, an int32. Default `defaultHash`. */
  readonly hash?: (key: K) => number;
  /** Key equality; must agree with `hash`. Default `sameValueZero`. */
  readonly equals?: (a: K, b: K) => boolean;
  /** Value equality: no-op edits and `equals`. Default `Object.is`. */
  readonly valueEquals?: (a: V, b: V) => boolean;
  /** The value hash for `hashCode`; must agree with `valueEquals`. Default `defaultHash`, none when only `valueEquals` is set. */
  readonly valueHash?: (value: V) => number;
}

/** A reusable map cursor: `next()` fills `key` and `value`, and returns `false` at the end (README §5). */
export interface MapCursor<K, V> {
  next(): boolean;
  readonly key: K;
  readonly value: V;
}

/** A reusable cursor: `next()` fills `value`, and returns `false` at the end (README §5). */
export interface Cursor<T> {
  next(): boolean;
  readonly value: T;
}

/** JS `Map`'s key equality: 0 and −0 are one key, and NaN equals NaN. */
export function sameValueZero(a: unknown, b: unknown): boolean {
  return a === b || (Number.isNaN(a as number) && Number.isNaN(b as number));
}

/** A config with every function resolved; two configs are the same when their four functions are (README §4). */
export interface Resolved {
  readonly hash: (key: any) => number;
  /** `hash` for reads: NaN when the key cannot be present (an object `defaultHash` has not seen). */
  readonly lookupHash: (key: any) => number;
  readonly equals: (a: any, b: any) => boolean;
  /** `equals` is `sameValueZero`: inlined. */
  readonly defaultEquals: boolean;
  readonly valueEquals: (a: any, b: any) => boolean;
  /** Absent when `valueEquals` is given without it: `hashCode` covers keys and size only. */
  readonly valueHash: ((value: any) => number) | undefined;
}

export const DefaultConfig: Resolved = {
  hash: defaultHash,
  lookupHash: lookupDefaultHash,
  equals: sameValueZero,
  defaultEquals: true,
  valueEquals: Object.is,
  valueHash: defaultHash,
};

export function resolveConfig(config: CollectionConfig<any, any> | undefined): Resolved {
  if (!config) return DefaultConfig;
  const { hash, equals, valueEquals, valueHash } = config;
  if (!hash && !equals && !valueEquals && !valueHash) return DefaultConfig;
  const h = hash ?? defaultHash;
  const eq = equals ?? sameValueZero;
  const veq = valueEquals ?? Object.is;
  const r: Resolved = {
    hash: h,
    // a custom hash coerced as edits coerce it (champ.ts): NaN is absent only for the default's unseen objects
    lookupHash: h === defaultHash ? lookupDefaultHash : (key: any) => h(key) | 0,
    equals: eq,
    defaultEquals: eq === sameValueZero,
    valueEquals: veq,
    valueHash: valueHash ?? (valueEquals ? undefined : defaultHash),
  };
  return sameConfig(r, DefaultConfig) ? DefaultConfig : r;
}

export function sameConfig(a: Resolved, b: Resolved): boolean {
  return (
    a === b ||
    (a.hash === b.hash && a.equals === b.equals && a.valueEquals === b.valueEquals && a.valueHash === b.valueHash)
  );
}

/** Ownership of a transient's nodes (README §3): a node is written in place only under its owner's live token. */
export interface Owner {
  dead: boolean;
}

export function newOwner(): Owner {
  return { dead: false };
}

export function deadHandle(): never {
  throw new Error('This transient was frozen by asImmutable(): edit the collection it returned instead.');
}

/** The ordered fold of `hashCode` (vectors, ordered collections). */
export function hashOrdered(acc: number, h: number): number {
  return (Math.imul(31, acc) + h) | 0;
}

/** An entry's share of an unordered `hashCode`. */
export function hashEntry(kh: number, vh: number): number {
  return (kh ^ Math.imul(vh, 0x5bd1e995)) | 0;
}
