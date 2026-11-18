/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Vector } from '../vector.js';
import { Rng } from './model.js';

type V = Vector<number>;

function check(v: { size: number; toArray(): number[] } & Omit<V, 'set'>, model: number[]) {
  expect(v.size).toBe(model.length);
  const arr = v.toArray();
  if (arr.length !== model.length || arr.some((x, i) => !Object.is(x, model[i]))) expect(arr).toEqual(model);
  for (let i = -1; i <= model.length; i++) if (!Object.is(v.get(i), model[i])) expect(v.get(i)).toBe(model[i]);
  expect(v.first()).toBe(model[0]);
  expect(v.last()).toBe(model[model.length - 1]);
  const seen: number[] = [];
  v.forEach(function (this: number[], x, i) {
    expect(i).toBe(this.length);
    this.push(x);
  }, seen);
  expect(seen).toEqual(model);
  const c = v.cursor();
  const viaCursor: number[] = [];
  while (c.next()) viaCursor.push(c.value);
  expect(viaCursor).toEqual(model);
}

/** Pushes and pops in runs, so the trie grows past a level (32 × 32 + 32) and shrinks back. */
function runVectorModel(seed: number, steps: number, maxRun: number) {
  const rng = new Rng(seed);
  let v: V = Vector.empty<number>();
  const model: number[] = [];
  const versions: [V, number[]][] = [];
  for (let step = 0; step < steps; step++) {
    const op = rng.int(12);
    if (op < 4) {
      const n = rng.int(maxRun);
      for (let i = 0; i < n; i++) {
        const x = rng.int(1000);
        v = v.push(x);
        model.push(x);
      }
    } else if (op < 6) {
      const n = rng.int(maxRun);
      for (let i = 0; i < n && model.length; i++) {
        v = v.pop();
        model.pop();
      }
    } else if (op < 8 && model.length) {
      const i = rng.int(model.length);
      const x = rng.int(4);
      const prev = v;
      v = v.set(i, x);
      if (Object.is(model[i], x)) expect(v).toBe(prev);
      model[i] = x;
    } else if (op < 9 && model.length) {
      if (rng.int(2)) {
        const i = rng.int(model.length);
        v = v.remove(i);
        model.splice(i, 1);
      } else {
        v = v.shift();
        model.shift();
      }
    } else {
      const t = v.asMutable();
      expect(t.asMutable()).toBe(t);
      const n = rng.int(maxRun);
      for (let i = 0; i < n; i++) {
        const r = rng.int(5);
        if (r < 3) {
          const x = rng.int(1000);
          expect(t.push(x)).toBe(t);
          model.push(x);
        } else if (r < 4 && model.length) {
          t.pop();
          model.pop();
        } else if (model.length) {
          const j = rng.int(model.length);
          t.set(j, -j);
          model[j] = -j;
        }
      }
      check(t, model);
      v = t.asImmutable();
      expect(() => t.push(0)).toThrow();
    }
    check(v, model);
    versions.push([v, model.slice()]);
  }
  for (const [x, m] of versions) check(x, m);
}

describe('Vector', () => {
  it('agrees with an Array over random pushes, pops and sets (seeded)', () => {
    for (const seed of [1, 2]) runVectorModel(seed, 60, 400);
  });
  it('a long random run', { tags: ['heavy'] }, () => {
    for (let seed = 10; seed < 16; seed++) runVectorModel(seed, 400, 3000);
  });

  it('remove and shift move the later elements', () => {
    const xs = Array.from({ length: 1100 }, (_, i) => i);
    const v = Vector.of(xs);
    for (const i of [0, 1, 31, 32, 33, 500, 1023, 1024, 1056, 1099]) {
      const r = v.remove(i);
      const m = xs.slice();
      m.splice(i, 1);
      check(r, m);
    }
    check(v.shift(), xs.slice(1));
    check(v, xs);
    const t = v.asMutable();
    t.shift().remove(10);
    const m = xs.slice(1);
    m.splice(10, 1);
    check(t, m);
    check(v, xs);
    expect(() => v.remove(1100)).toThrow(RangeError);
    expect(() => v.set(-1, 0)).toThrow(RangeError);
  });

  it('indexOf by valueEquals, find by predicate', () => {
    const v = Vector.of([1, NaN, -0, 3, 3]);
    expect(v.indexOf(3)).toBe(3);
    expect(v.indexOf(NaN)).toBe(1);
    expect(v.indexOf(0)).toBe(-1);
    expect(v.indexOf(-0)).toBe(2);
    expect(v.find((x, i) => i > 3 && x === 3)).toBe(3);
    expect(v.find((x) => x > 10)).toBe(undefined);
  });

  it('equality and hashCode are ordered; configs part of equality', () => {
    const a = Vector.of([1, 2, 3]);
    expect(a.equals(Vector.of([1, 2, 3]))).toBe(true);
    expect(a.hashCode()).toBe(Vector.of([1, 2, 3]).hashCode());
    expect(a.equals(Vector.of([3, 2, 1]))).toBe(false);
    expect(a.hashCode()).not.toBe(Vector.of([3, 2, 1]).hashCode());
    const big = Vector.of(Array.from({ length: 2000 }, (_, i) => i));
    expect(big.equals(big.set(1500, 0).set(1500, 1500))).toBe(true);
    expect(Vector.empty({}) as unknown).toBe(Vector.empty());
    const structural = { valueEquals: (x: number[], y: number[]) => x[0] === y[0], valueHash: (x: number[]) => x[0] };
    const s = Vector.of([[1], [2]], structural);
    expect(s.equals(Vector.of([[1], [2]], structural))).toBe(true);
    expect(s.equals(Vector.of([[1], [2]]) as any)).toBe(false);
    expect(s.set(0, [1])).toBe(s);
  });

  it('a transient pushing into a shared tail copies it first', () => {
    const base = Vector.of(Array.from({ length: 64 }, (_, i) => i));
    const t = base.asMutable();
    for (let i = 0; i < 40; i++) t.push(-i);
    t.set(10, 99).set(63, 99);
    check(
      base,
      Array.from({ length: 64 }, (_, i) => i),
    );
    expect(t.get(10)).toBe(99);
  });
});
