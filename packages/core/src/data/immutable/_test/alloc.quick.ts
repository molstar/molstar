/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The allocation targets (README §5): `get`, `has`, `forEach` with a context and a reused cursor allocate nothing warm,
// and a lookup of an unseen object assigns no id; a persistent edit copies a path, within `EditBound` bytes. By hand,
// never in CI: `npx vitest run -c packages/core/vitest.perf.config.ts` (`--reporter=verbose` prints the bytes a step).

import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import { lookupDefaultHash } from '../hash.js';
import { HashMap } from '../hash-map.js';
import { HashSet } from '../hash-set.js';
import { OrderedHashMap } from '../ordered-hash-map.js';
import { OrderedHashSet } from '../ordered-hash-set.js';
import { Vector } from '../vector.js';
import { id22, Rng } from './model.js';

/**
 * Bytes a persistent edit of a 1,000-entry map may allocate: the copied path (a node object and its content array a
 * level, about three levels at this size) and the new map object, with headroom.
 */
const EditBound = 1_200;

const rng = new Rng(1);
const ids = Array.from({ length: 1000 }, () => id22(rng));
const small = ids.slice(0, 16);

const zero = () => inject('perfThreshold');

function count(this: { n: number }) {
  this.n++;
}

describe('reads allocate nothing', () => {
  for (const [name, keys] of [
    ['flat (16)', small],
    ['trie (1,000)', ids],
  ] as const) {
    const m = HashMap.of(keys.map((k, i) => [k, i] as const));
    const s = HashSet.of(keys);
    const om = OrderedHashMap.of(keys.map((k, i) => [k, i] as const));
    const os = OrderedHashSet.of(keys);
    const ctx = { n: 0 };

    it(`HashMap ${name}: get, has`, ({ task }) => {
      const { r } = quickCase(task, `HashMap ${name} get`, (i) => {
        const k = keys[i % keys.length];
        if (m.get(k) === undefined || !m.has(k)) throw new Error('absent');
      });
      expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(zero());
    });

    it(`HashMap ${name}: forEach with a context, a reused cursor`, ({ task }) => {
      const c = m.cursor();
      const { r } = quickCase(task, `HashMap ${name} walk`, () => {
        m.forEach(count, ctx);
        m.cursor(c);
        while (c.next()) ctx.n += c.value;
      });
      expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(zero());
    });

    it(`HashSet ${name}: has, cursor`, ({ task }) => {
      const c = s.cursor();
      const { r } = quickCase(task, `HashSet ${name} has`, (i) => {
        if (!s.has(keys[i % keys.length])) throw new Error('absent');
        if (i % 64 === 0) {
          s.cursor(c);
          while (c.next()) ctx.n++;
        }
      });
      expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(zero());
    });

    it(`OrderedHashMap and OrderedHashSet ${name}: get, has, forEach, cursor`, ({ task }) => {
      const c = om.cursor();
      const d = os.cursor();
      const { r } = quickCase(task, `Ordered ${name} reads`, (i) => {
        const k = keys[i % keys.length];
        if (om.get(k) === undefined || !os.has(k)) throw new Error('absent');
        if (i % 64 === 0) {
          om.forEach(count, ctx);
          os.forEach(count, ctx);
          om.cursor(c);
          while (c.next()) ctx.n++;
          os.cursor(d);
          while (d.next()) ctx.n++;
        }
      });
      expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(zero());
    });
  }

  it('Vector: get, forEach with a context, cursor', ({ task }) => {
    const v = Vector.of(Array.from({ length: 1000 }, (_, i) => i));
    const c = v.cursor();
    const ctx = { n: 0 };
    const { r } = quickCase(task, 'Vector reads', (i) => {
      ctx.n += v.get(i % 1000)!;
      if (i % 64 === 0) {
        v.forEach(count, ctx);
        v.cursor(c);
        while (c.next()) ctx.n += c.value;
      }
    });
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(zero());
  });

  it('a lookup of an unseen object assigns no id', ({ task }) => {
    const objects = Array.from({ length: 100 }, (_, i) => ({ i }));
    const m = HashMap.of(objects.map((o, i) => [o, i] as const));
    const unseen = Array.from({ length: 100 }, (_, i) => ({ i }));
    const { r } = quickCase(task, 'HashMap unseen object get', (i) => {
      if (m.has(unseen[i % 100])) throw new Error('present');
    });
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(zero());
    for (const o of unseen) expect(Number.isNaN(lookupDefaultHash(o))).toBe(true);
  });
});

describe('an edit stays within its bound', () => {
  it('a persistent set and delete on a 1,000-entry map', ({ task }) => {
    const m = HashMap.of(ids.map((k, i) => [k, i] as const));
    const { r } = quickCase(
      task,
      'HashMap persistent set',
      (i) => {
        if (m.set(ids[i % 1000], -1).size !== 1000) throw new Error('size');
      },
      { control: false, why: 'an edit allocates by design: the control hides in it' },
    );
    console.info(`HashMap persistent set (1,000): ${r.bytesPerStep.toFixed(0)} bytes (${r.note})`);
    expect(r.bytesPerStep, r.note).toBeLessThan(EditBound);
    const { r: d } = quickCase(
      task,
      'HashMap persistent delete',
      (i) => {
        if (m.delete(ids[i % 1000]).size !== 999) throw new Error('size');
      },
      { control: false, why: 'an edit allocates by design: the control hides in it' },
    );
    console.info(`HashMap persistent delete (1,000): ${d.bytesPerStep.toFixed(0)} bytes (${d.note})`);
    expect(d.bytesPerStep, d.note).toBeLessThan(EditBound);
  });

  it('a transient set of present keys writes in place', ({ task }) => {
    const t = HashMap.of(ids.map((k, i) => [k, i] as const)).asMutable();
    // the first pass copies the paths under the transient's owner
    for (let i = 0; i < 1000; i++) t.set(ids[i], -1);
    const { r } = quickCase(task, 'HashMap transient set', (i) => {
      t.set(ids[i % 1000], i);
    });
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(zero());
  });
});
