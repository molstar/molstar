/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The model-based runs (README §8): random edits on a collection and on a JS `Map` or `Set` side by side, checked
// after every step; every earlier version is checked again at the end (persistence), transients included.

import { expect } from 'vitest';
import type { CollectionConfig, MapCursor } from '../config.js';
import type { HashMap, MutableHashMap } from '../hash-map.js';
import type { HashSet, MutableHashSet } from '../hash-set.js';
import { keyPool, Rng, sorted } from './model.js';

export const norm = (k: unknown) => (Object.is(k, -0) ? 0 : k);

type AnyMap = HashMap<unknown, number> & { first?(): unknown; last?(): unknown };
type AnyMutableMap = MutableHashMap<unknown, number>;

export interface MapKind {
  name: string;
  ordered: boolean;
  empty(config?: CollectionConfig<unknown, number>): AnyMap;
  of(entries: Iterable<readonly [unknown, number]>, config?: CollectionConfig<unknown, number>): AnyMap;
}

export function checkMap(m: AnyMap | AnyMutableMap, model: Map<unknown, number>, pool: unknown[], ordered: boolean) {
  expect(m.size).toBe(model.size);
  for (const k of pool) {
    // compared first, asserted on a mismatch: an assertion a key is most of a long run
    if (m.has(k) !== model.has(k) || m.get(k) !== model.get(k)) {
      expect(m.has(k), `has ${String(k)}`).toBe(model.has(k));
      expect(m.get(k), `get ${String(k)}`).toBe(model.get(k));
    }
  }
  const entries = [...model].map(([k, v]) => [norm(k), v]);
  const arr = m.toArray().map(([k, v]) => [norm(k), v]);
  if (ordered) expect(arr).toEqual(entries);
  else expect(sorted(arr)).toEqual(sorted(entries));
  const seen: unknown[][] = [];
  m.forEach(function (this: unknown[][], v, k) {
    this.push([norm(k), v]);
  }, seen);
  expect(seen).toEqual(arr);
  const viaCursor: unknown[][] = [];
  const c: MapCursor<unknown, number> = m.cursor();
  while (c.next()) viaCursor.push([norm(c.key), c.value]);
  expect(viaCursor).toEqual(arr);
  expect(c.next()).toBe(false);
  expect([...m].map(([k, v]) => [norm(k), v])).toEqual(arr);
  expect([...m.keys()].map(norm)).toEqual(arr.map((e) => e[0]));
  expect([...m.values()]).toEqual(arr.map((e) => e[1]));
  if (ordered) {
    const first = model.size ? [...model.values()][0] : undefined;
    const last = model.size ? [...model.values()][model.size - 1] : undefined;
    expect((m as AnyMap).first!()).toBe(first);
    expect((m as AnyMap).last!()).toBe(last);
  }
}

/** Random sets, deletes and transient batches; each version checked at once and at the end. */
/**
 * Bad hashes for the model runs: every key alike, and one key in forty beside the rest, which leaves an inline entry
 * beside a collision node of more than the flat threshold (README §2: the node goes up to the root on deletion).
 */
export const BadHashes: Record<string, CollectionConfig<unknown, any>> = {
  zero: { hash: () => 0 },
  rare: { hash: (k: unknown) => (typeof k === 'number' && k % 40 === 0 ? 1 : 0) },
};

export function runMapModel(
  kind: MapKind,
  seed: number,
  steps: number,
  poolSize: number,
  config?: CollectionConfig<unknown, number>,
) {
  const rng = new Rng(seed);
  const pool = keyPool(rng, poolSize);
  let m = kind.empty(config);
  const model = new Map<unknown, number>();
  const versions: [AnyMap, Map<unknown, number>][] = [];
  for (let step = 0; step < steps; step++) {
    const op = rng.int(10);
    if (op < 5) {
      const k = rng.pick(pool);
      const v = rng.int(4);
      const prev = m;
      m = m.set(k, v) as AnyMap;
      if (model.has(k) && Object.is(model.get(k), v)) expect(m).toBe(prev);
      model.set(k, v);
    } else if (op < 8) {
      const k = rng.pick(pool);
      const prev = m;
      m = m.delete(k) as AnyMap;
      if (!model.has(k)) expect(m).toBe(prev);
      model.delete(k);
    } else {
      // a transient batch, read while it edits
      const t = m.asMutable() as AnyMutableMap;
      expect(t.asMutable()).toBe(t);
      const n = 1 + rng.int(20);
      for (let i = 0; i < n; i++) {
        const k = rng.pick(pool);
        if (rng.int(3) === 0) {
          expect(t.delete(k)).toBe(t);
          model.delete(k);
        } else {
          const v = rng.int(4);
          expect(t.set(k, v)).toBe(t);
          model.set(k, v);
        }
        if (i % 5 === 0) checkMap(t, model, pool, kind.ordered);
      }
      const next = t.asImmutable() as AnyMap;
      expect(next).not.toBe(t);
      expect(() => t.set(pool[0], 0)).toThrow();
      expect(() => t.asImmutable()).toThrow();
      m = next;
    }
    checkMap(m, model, pool, kind.ordered);
    versions.push([m, new Map(model)]);
    if (step % 16 === 0) {
      // equal contents built in another order are equal and hash alike
      const other = kind.of(rng.shuffle([...model]), config);
      if (kind.ordered) {
        const inOrder = kind.of([...model], config);
        expect(inOrder.equals(m)).toBe(true);
        expect(inOrder.hashCode()).toBe(m.hashCode());
      } else {
        expect(other.equals(m)).toBe(true);
        expect(m.equals(other)).toBe(true);
        expect(other.hashCode()).toBe(m.hashCode());
      }
    }
  }
  for (const [v, snapshot] of versions) checkMap(v, snapshot, pool, kind.ordered);
}

type AnySet = HashSet<unknown> & { first?(): unknown; last?(): unknown };

export interface SetKind {
  name: string;
  ordered: boolean;
  empty(config?: CollectionConfig<unknown>): AnySet;
  of(keys: Iterable<unknown>, config?: CollectionConfig<unknown>): AnySet;
}

export function checkSet(s: AnySet | MutableHashSet<unknown>, model: Set<unknown>, pool: unknown[], ordered: boolean) {
  expect(s.size).toBe(model.size);
  for (const k of pool) if (s.has(k) !== model.has(k)) expect(s.has(k), `has ${String(k)}`).toBe(model.has(k));
  const keys = [...model].map(norm);
  const arr = s.toArray().map(norm);
  if (ordered) expect(arr).toEqual(keys);
  else expect(sorted(arr)).toEqual(sorted(keys));
  const seen: unknown[] = [];
  s.forEach((v, k) => {
    expect(v).toBe(k);
    seen.push(norm(k));
  });
  expect(seen).toEqual(arr);
  const viaCursor: unknown[] = [];
  const c = s.cursor();
  while (c.next()) viaCursor.push(norm(c.value));
  expect(viaCursor).toEqual(arr);
  expect([...s].map(norm)).toEqual(arr);
  if (ordered) {
    expect(norm((s as AnySet).first!())).toBe(model.size ? keys[0] : undefined);
    expect(norm((s as AnySet).last!())).toBe(model.size ? keys[keys.length - 1] : undefined);
  }
}

export function runSetModel(
  kind: SetKind,
  seed: number,
  steps: number,
  poolSize: number,
  config?: CollectionConfig<unknown>,
) {
  const rng = new Rng(seed);
  const pool = keyPool(rng, poolSize);
  let s = kind.empty(config);
  const model = new Set<unknown>();
  const versions: [AnySet, Set<unknown>][] = [];
  for (let step = 0; step < steps; step++) {
    const op = rng.int(10);
    if (op < 5) {
      const k = rng.pick(pool);
      const prev = s;
      s = s.add(k) as AnySet;
      if (model.has(k)) expect(s).toBe(prev);
      model.add(k);
    } else if (op < 8) {
      const k = rng.pick(pool);
      s = s.delete(k) as AnySet;
      model.delete(k);
    } else {
      const t = s.asMutable();
      const n = 1 + rng.int(20);
      for (let i = 0; i < n; i++) {
        const k = rng.pick(pool);
        if (rng.int(3) === 0) {
          t.delete(k);
          model.delete(k);
        } else {
          t.add(k);
          model.add(k);
        }
        if (i % 5 === 0) checkSet(t, model, pool, kind.ordered);
      }
      s = t.asImmutable() as AnySet;
      expect(() => t.add(pool[0])).toThrow();
    }
    checkSet(s, model, pool, kind.ordered);
    versions.push([s, new Set(model)]);
    if (step % 16 === 0) {
      const other = kind.ordered ? kind.of([...model], config) : kind.of(rng.shuffle([...model]), config);
      expect(other.equals(s)).toBe(true);
      expect(other.hashCode()).toBe(s.hashCode());
    }
  }
  for (const [v, snapshot] of versions) checkSet(v, snapshot, pool, kind.ordered);
}
