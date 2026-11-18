/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { ChampNode, FlatNode, FlatThreshold, KindChamp, KindCollision, KindFlat } from '../champ.js';
import { HashMap, type HashMapImpl } from '../hash-map.js';
import { OrderedHashMap } from '../ordered-hash-map.js';
import { id22, keyPool, Rng } from './model.js';
import { BadHashes, checkMap, type MapKind, runMapModel } from './run-model.js';

const Plain: MapKind = { name: 'HashMap', ordered: false, empty: HashMap.empty as any, of: HashMap.of as any };
const Ordered: MapKind = {
  name: 'OrderedHashMap',
  ordered: true,
  empty: OrderedHashMap.empty as any,
  of: OrderedHashMap.of as any,
};

const zeroHash = { hash: () => 0 };
const fourHashes = { hash: (k: unknown) => (typeof k === 'number' ? k & 3 : 0) };

function impl<K, V>(m: unknown) {
  return m as HashMapImpl<K, V>;
}

/** The trie part's shape, for comparing shapes across edit orders: collision nodes and flat roots as sorted sizes. */
function shape(node: any): unknown {
  if (node.kind === KindFlat) return ['flat', node.entries.length];
  if (node.kind === KindCollision) return ['collision', node.hash, node.entries.length];
  const d = 2 * popcount(node.dataMap);
  const subs = [];
  for (let i = node.content.length - 1; i >= d; i--) subs.push(shape(node.content[i]));
  return [node.dataMap, node.nodeMap, node.content.slice(0, d).filter((_: unknown, i: number) => i % 2 === 0), subs];
}

function popcount(x: number) {
  let n = 0;
  for (; x; x &= x - 1) n++;
  return n;
}

describe('HashMap', () => {
  for (const kind of [Plain, Ordered]) {
    it(`${kind.name}: agrees with a JS Map over random edits (seeded)`, () => {
      for (const seed of [1, 2, 3]) runMapModel(kind, seed, 150, 60);
    });
    it(`${kind.name}: agrees with a JS Map under bad hashes (seeded)`, () => {
      for (const config of Object.values(BadHashes))
        for (const seed of [4, 5]) runMapModel(kind, seed, 150, 90, config);
    });
    it(`${kind.name}: a long random run`, { tags: ['heavy'] }, () => {
      for (const config of Object.values(BadHashes)) {
        for (let seed = 40; seed < 42; seed++) runMapModel(kind, seed, 1500, 300, config);
      }
      for (let seed = 10; seed < 20; seed++) runMapModel(kind, seed, 1500, 300);
    });
  }

  it('−0 and 0 are one key, NaN equals NaN (SameValueZero)', () => {
    const m = HashMap.empty<number, string>().set(-0, 'a').set(NaN, 'n');
    expect(m.get(0)).toBe('a');
    expect(m.set(0, 'b').size).toBe(2);
    expect(m.get(Number('x'))).toBe('n');
    const strict = HashMap.empty<number, string>({ equals: Object.is }).set(-0, 'a').set(0, 'b');
    expect(strict.size).toBe(2);
    expect(strict.get(-0)).toBe('a');
  });

  it('stores undefined: set and update keep it, only delete removes', () => {
    const m = HashMap.empty<string, number | undefined>().set('a', undefined);
    expect(m.has('a')).toBe(true);
    expect(m.size).toBe(1);
    const u = m.update('b', () => undefined);
    expect(u.has('b')).toBe(true);
    expect(u.update('b', (v) => (v ?? 0) + 1).get('b')).toBe(1);
  });

  it('a no-op edit returns the same map', () => {
    const m = HashMap.of([
      ['a', 1],
      ['b', 2],
    ]);
    expect(m.set('a', 1)).toBe(m);
    expect(m.delete('c')).toBe(m);
    expect(m.delete({} as any)).toBe(m);
    const structural = HashMap.of([['a', [1, 2]]], {
      valueEquals: (x, y) => x.length === y.length && x.every((v, i) => v === y[i]),
    });
    expect(structural.set('a', [1, 2])).toBe(structural);
  });

  it('flat up to the threshold, a trie above it, flat again on deletion at the threshold', () => {
    let m = HashMap.empty<number, number>();
    for (let i = 0; i < FlatThreshold; i++) m = m.set(i, i);
    expect(impl(m).root).toBeInstanceOf(FlatNode);
    m = m.set(FlatThreshold, 0);
    expect(impl(m).root).toBeInstanceOf(ChampNode);
    m = m.delete(3);
    expect(impl(m).root).toBeInstanceOf(FlatNode);
    expect(m.size).toBe(FlatThreshold);
    for (let i = 0; i <= FlatThreshold; i++) expect(m.get(i)).toBe(i === 3 ? undefined : i === FlatThreshold ? 0 : i);
  });

  it('the trie part has one shape for equal contents, whatever the edit order', () => {
    const rng = new Rng(7);
    const keys = keyPool(rng, 300);
    for (const config of [undefined, fourHashes, { hash: (k: unknown) => (typeof k === 'string' ? k.length : 1) }]) {
      const a = HashMap.of(
        keys.map((k) => [k, 1] as const),
        config,
      );
      // a superset, then deletes back down
      let b = HashMap.of(
        rng.shuffle([...keys, ...keyPool(new Rng(8), 200).slice(19)]).map((k) => [k, 1] as const),
        config,
      );
      for (const k of keyPool(new Rng(8), 200).slice(19)) if (!keys.includes(k)) b = b.delete(k);
      expect(b.size).toBe(a.size);
      expect(shape(impl(b).root)).toEqual(shape(impl(a).root));
      expect(b.equals(a)).toBe(true);
      expect(b.hashCode()).toBe(a.hashCode());
    }
  });

  it('a bad hash (every key 0) gives correct results', () => {
    for (const config of [zeroHash, fourHashes]) {
      const rng = new Rng(3);
      const pool = keyPool(rng, 40);
      let m = HashMap.empty<unknown, number>(config);
      const model = new Map<unknown, number>();
      for (let i = 0; i < 400; i++) {
        const k = rng.pick(pool);
        if (rng.int(3) === 0) {
          m = m.delete(k);
          model.delete(k);
        } else {
          m = m.set(k, i);
          model.set(k, i);
        }
        checkMap(m as any, model, pool, false);
      }
      expect(HashMap.of([...model], config).equals(m)).toBe(true);
    }
  });

  it('a collision node moved up to the root is filed under its own hash, not the deleted key’s', () => {
    // keys 0..17 collide (hash 0); 100 sits inline beside them (hash 1); deleting 100 moves the collision node up
    const m = HashMap.of(
      [...Array.from({ length: 18 }, (_, i) => i), 100].map((k) => [k, k] as const),
      { hash: (k: number) => (k < 100 ? 0 : 1) },
    );
    const d = m.delete(100);
    expect(d.size).toBe(18);
    for (let k = 0; k < 18; k++) expect(d.get(k)).toBe(k);
    expect(d.has(100)).toBe(false);
  });

  it('a collision node alone under a node goes up on deletion: the shape insertion makes', () => {
    // 1 and 33 collide; 2 shares their first fragment; the rest are elsewhere
    const hash = (k: number) => (k === 1 || k === 33 ? 1 : k === 2 ? 1 | (1 << 5) : k << 10);
    let m = HashMap.empty<number, number>({ hash });
    for (let i = 100; i < 100 + FlatThreshold; i++) m = m.set(i, i);
    const a = m.set(1, 1).set(33, 33);
    const b = m.set(2, 2).set(1, 1).set(33, 33).delete(2);
    expect(shape(impl(b).root)).toEqual(shape(impl(a).root));
    expect(b.equals(a)).toBe(true);
  });

  describe('transients', () => {
    it('read like the persistent map, edit in place, and die on asImmutable', () => {
      const base = HashMap.of(Array.from({ length: 50 }, (_, i) => [`k${i}`, i] as const));
      const t = base.asMutable();
      expect(t.asMutable()).toBe(t);
      t.set('x', 1).delete('k0');
      expect(t.get('x')).toBe(1);
      expect(t.has('k0')).toBe(false);
      expect(t.size).toBe(50);
      expect(base.has('k0')).toBe(true);
      expect(base.has('x')).toBe(false);
      const p = t.asImmutable();
      expect(p).not.toBe(t);
      expect(p.size).toBe(50);
      expect(() => t.set('y', 2)).toThrow(/asImmutable/);
      expect(() => t.delete('x')).toThrow();
      expect(() => t.asImmutable()).toThrow();
      expect(() => t.asMutable()).toThrow();
      // reads on the dead handle still answer
      expect(t.get('x')).toBe(1);
      // the next transient is isolated from p
      const t2 = p.asMutable();
      t2.set('x', 2);
      expect(p.get('x')).toBe(1);
      expect(t2.asImmutable().get('x')).toBe(2);
    });

    it('withMutations returns a persistent map', () => {
      const m = HashMap.empty<string, number>().withMutations((t) => {
        for (let i = 0; i < 100; i++) t.set(`${i}`, i);
      });
      expect(m.size).toBe(100);
      expect(m.set('0', 0)).toBe(m);
    });

    it('deleting from colliding keys writes no node another version shares', () => {
      // the regression StateTree's `renewOwner` works around in immutable.js (state/_spec/tree.md §2)
      for (const config of [zeroHash, fourHashes, undefined]) {
        const keys = Array.from({ length: 40 }, (_, i) => i * 4);
        const source = HashMap.of(
          keys.map((k) => [k, 'old'] as const),
          config,
        );
        const t = source.asMutable();
        t.delete(keys[0]);
        t.set(keys[1], 'new');
        t.delete(keys[2]);
        t.set(keys[3], 'new');
        expect(source.get(keys[1])).toBe('old');
        expect(source.get(keys[3])).toBe('old');
        expect(source.has(keys[0])).toBe(true);
        expect(source.size).toBe(keys.length);
        const done = t.asImmutable();
        expect(done.get(keys[1])).toBe('new');
        expect(done.size).toBe(keys.length - 2);
      }
    });
  });

  describe('equality and hashCode', () => {
    it('different configs are unequal; empty() and empty({}) are the same config', () => {
      const a = HashMap.of([['a', 1]]);
      expect(HashMap.empty().set('a', 1).equals(a)).toBe(true);
      expect(HashMap.empty<string, number>({}).set('a', 1).equals(a)).toBe(true);
      expect(HashMap.empty({}) as unknown).toBe(HashMap.empty());
      const other = HashMap.of([['a', 1]], { hash: (k: string) => k.length });
      expect(other.equals(a)).toBe(false);
      expect(a.equals(other)).toBe(false);
      expect(OrderedHashMap.of([['a', 1]]).equals(a as any)).toBe(false);
    });

    it('a structural valueEquals with its valueHash gives equal maps equal hashCodes', () => {
      const config = {
        valueEquals: (x: number[], y: number[]) => x.length === y.length && x.every((v, i) => v === y[i]),
        valueHash: (x: number[]) => x.reduce((h, v) => (h * 31 + v) | 0, 7),
      };
      const a = HashMap.of(
        Array.from({ length: 30 }, (_, i) => [i, [i, i + 1]] as const),
        config,
      );
      const b = HashMap.of(
        Array.from({ length: 30 }, (_, i) => [29 - i, [29 - i, 30 - i]] as const),
        config,
      );
      expect(a.equals(b)).toBe(true);
      expect(a.hashCode()).toBe(b.hashCode());
      expect(a.equals(b.set(3, [0]))).toBe(false);
      // valueEquals without valueHash: keys and size only
      const keysOnly = { valueEquals: config.valueEquals };
      const c = HashMap.of([[1, [1]]], keysOnly);
      expect(c.hashCode()).toBe(HashMap.of([[1, [2]]], keysOnly).hashCode());
    });

    it('values take part in equality and hashCode', () => {
      const ids = Array.from({ length: 20 }, () => id22(new Rng(5)));
      const a = HashMap.of(ids.map((k, i) => [k, i] as const));
      expect(a.equals(a.set(ids[0], -1))).toBe(false);
      expect(a.hashCode()).not.toBe(a.set(ids[0], -1).hashCode());
    });
  });

  it('a cursor is reused across maps and kinds', () => {
    const a = HashMap.of(Array.from({ length: 40 }, (_, i) => [i, i] as const));
    const b = OrderedHashMap.of([['x', 1]]);
    const c = a.cursor();
    let n = 0;
    while (c.next()) n++;
    expect(n).toBe(40);
    expect(b.cursor(c as any)).toBe(c);
    expect(c.next()).toBe(true);
    expect(c.key).toBe('x');
    expect(c.next()).toBe(false);
    expect(a.cursor(c)).toBe(c);
    expect(c.next()).toBe(true);
  });

  it('kinds of the root', () => {
    expect(impl(HashMap.empty()).root.kind).toBe(KindFlat);
    expect(
      impl(
        HashMap.of(
          Array.from({ length: 20 }, (_, i) => [i, i] as const),
          zeroHash,
        ),
      ).root.kind,
    ).toBe(KindChamp);
  });
});
