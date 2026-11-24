/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { CompactMin } from '../ordered-hash-map.js';
import { HashSet } from '../hash-set.js';
import { OrderedHashSet, type OrderedHashSetImpl } from '../ordered-hash-set.js';
import { OrderedHashMap, type OrderedHashMapImpl } from '../ordered-hash-map.js';
import { BadHashes, runSetModel, type SetKind } from './run-model.js';

const Plain: SetKind = { name: 'HashSet', ordered: false, empty: HashSet.empty as any, of: HashSet.of as any };
const Ordered: SetKind = {
  name: 'OrderedHashSet',
  ordered: true,
  empty: OrderedHashSet.empty as any,
  of: OrderedHashSet.of as any,
};

describe('HashSet and OrderedHashSet', () => {
  for (const kind of [Plain, Ordered]) {
    it(`${kind.name}: agrees with a JS Set over random edits (seeded)`, () => {
      for (const seed of [1, 2, 3]) runSetModel(kind, seed, 150, 60);
    });
    it(`${kind.name}: agrees with a JS Set under bad hashes (seeded)`, () => {
      for (const config of Object.values(BadHashes))
        for (const seed of [4, 5]) runSetModel(kind, seed, 150, 90, config);
    });
    it(`${kind.name}: a long random run`, { tags: ['heavy'] }, () => {
      for (let seed = 10; seed < 20; seed++) runSetModel(kind, seed, 1500, 300);
    });
  }

  it('an all-zero hash gives correct results', () => {
    const zero = { hash: () => 0 };
    const s = HashSet.of(
      Array.from({ length: 30 }, (_, i) => i),
      zero,
    );
    expect(s.size).toBe(30);
    for (let i = 0; i < 30; i++) expect(s.has(i)).toBe(true);
    expect(s.has(30)).toBe(false);
    let t = s;
    for (let i = 0; i < 30; i += 2) t = t.delete(i);
    expect(t.size).toBe(15);
    expect(t.equals(HashSet.of([29, 27, 25, 23, 21, 19, 17, 15, 13, 11, 9, 7, 5, 3, 1], zero))).toBe(true);
    // another function is another config
    expect(t.equals(HashSet.of(t.toArray(), { hash: () => 0 }))).toBe(false);
  });

  it('order is part of an ordered set’s equality, not of a set’s', () => {
    expect(HashSet.of([1, 2, 3]).equals(HashSet.of([3, 2, 1]))).toBe(true);
    expect(HashSet.of([1, 2, 3]).hashCode()).toBe(HashSet.of([3, 2, 1]).hashCode());
    expect(OrderedHashSet.of([1, 2, 3]).equals(OrderedHashSet.of([3, 2, 1]))).toBe(false);
    expect(OrderedHashSet.of([1, 2, 3]).equals(OrderedHashSet.of([1, 2, 3]))).toBe(true);
    expect(
      OrderedHashSet.of([1, 2, 3])
        .delete(2)
        .add(2)
        .equals(OrderedHashSet.of([1, 3, 2])),
    ).toBe(true);
  });

  it('adding a present key keeps its place and returns the same set', () => {
    const s = OrderedHashSet.of(['a', 'b', 'c']);
    expect(s.add('a')).toBe(s);
    expect(s.add('d').toArray()).toEqual(['a', 'b', 'c', 'd']);
    const m = OrderedHashMap.of([
      ['a', 1],
      ['b', 2],
    ]);
    expect(m.set('a', 3).toArray()).toEqual([
      ['a', 3],
      ['b', 2],
    ]);
  });

  it('holes compact past half once there are CompactMin slots; trailing holes are trimmed', () => {
    let s = OrderedHashSet.of(Array.from({ length: CompactMin }, (_, i) => i));
    for (let i = 0; i < CompactMin / 2; i++) s = s.delete(i);
    expect((s as unknown as OrderedHashSetImpl<number>).holes).toBe(CompactMin / 2);
    s = s.delete(CompactMin / 2);
    expect((s as unknown as OrderedHashSetImpl<number>).holes).toBe(0);
    expect((s as unknown as OrderedHashSetImpl<number>).slots.size).toBe(s.size);
    expect(s.first()).toBe(CompactMin / 2 + 1);
    const t = OrderedHashSet.of([1, 2, 3, 4]).delete(3).delete(4);
    expect((t as unknown as OrderedHashSetImpl<number>).slots.size).toBe(2);
    expect(t.last()).toBe(2);
    const u = OrderedHashSet.of([1, 2, 3, 4]).delete(3).delete(2).delete(4);
    expect((u as unknown as OrderedHashSetImpl<number>).slots.size).toBe(1);
    expect(u.last()).toBe(1);

    let m = OrderedHashMap.of(Array.from({ length: 40 }, (_, i) => [i, i] as const));
    for (let i = 0; i < 21; i++) m = m.delete(i);
    expect((m as unknown as OrderedHashMapImpl<number, number>).holes).toBe(0);
    expect(m.first()).toBe(21);
    expect(m.toArray().length).toBe(19);
  });

  it('a transient set: reads, edits in place, dies on asImmutable', () => {
    const base = OrderedHashSet.of(['a', 'b']);
    const t = base.asMutable();
    expect(t.asMutable()).toBe(t);
    expect(t.add('c').delete('a')).toBe(t);
    expect(t.toArray()).toEqual(['b', 'c']);
    expect(t.first()).toBe('b');
    expect(base.toArray()).toEqual(['a', 'b']);
    const p = t.asImmutable();
    expect(p.toArray()).toEqual(['b', 'c']);
    expect(() => t.add('x')).toThrow();
    // emptied in a transient
    const e = p.withMutations((m) => {
      m.delete('b').delete('c');
    });
    expect(e).toBe(OrderedHashSet.empty());
  });

  it('withMutations and of', () => {
    const s = HashSet.empty<number>().withMutations((m) => {
      for (let i = 0; i < 100; i++) m.add(i % 50);
    });
    expect(s.size).toBe(50);
    expect(s.equals(HashSet.of(Array.from({ length: 50 }, (_, i) => 49 - i)))).toBe(true);
  });
});
