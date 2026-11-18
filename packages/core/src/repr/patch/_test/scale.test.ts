/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A patch costs the paths it touches, not the tree: a one-key set on a large tree, with and without a language, and a
// diff of a wide parent reordered (a longest increasing subsequence, O(n log n)). Quick without the clock, heavy with
// it (itWithin).
// Each bound is about 4-8x what the case takes alone, the best of three runs (a parallel heavy run read one sample
// at 26 ms against 20, and 111 against 80).

import { bestOf, itWithin } from '@molstar/test-harness/timing';
import { describe, expect } from 'vitest';
import { applyPatch, diff } from '../../patch.js';
import { hashNode, type Node } from '../../tree/index.js';
import { PatchConformanceLanguage as PC } from '../conformance.js';

/** A doc of `boxes` boxes of `leaves` leaves each, every leaf with a ref to its box's first leaf. */
function wide(boxes: number, leaves: number): Node {
  return {
    kind: 'doc',
    children: Array.from({ length: boxes }, (_, i) => ({
      kind: 'box',
      id: `b${i}`,
      children: Array.from({ length: leaves }, (_, j) => ({
        kind: 'leaf',
        id: `l${j}`,
        params: { label: `${i}.${j}`, link: `/b${i}/l0` },
      })),
    })),
  };
}

describe('repr/patch: what a patch costs', () => {
  const big = wide(200, 100); // 20,000 leaves
  hashNode(big);

  itWithin('a one-key set on a tree of 20,000 nodes shares all but its path', (within) => {
    const op = { op: 'set', at: ['b100', 'l50'], key: 'label', value: 'x' } as const;
    applyPatch(big, { ops: [op] }); // warm
    const plain = within(5, () => applyPatch(big, { ops: [op] }));
    const typed = within(5, () => applyPatch(big, { ops: [op] }, { language: PC }));
    for (const r of [plain, typed]) {
      expect(r.ok).toBe(true);
      if (!r.ok) continue;
      expect(r.root.children![99]).toBe(big.children![99]);
      expect(r.root.children![100].children![49]).toBe(big.children![100].children![49]);
    }
  });

  itWithin('a diff of a parent of 20,000 children, reordered, keeps the longest run', (within) => {
    const a: Node = {
      kind: 'doc',
      children: big.children!.flatMap((b) => b.children!.map((l) => ({ ...l, id: `${b.id}-${l.id}` }))),
    };
    const kids = [...a.children!];
    kids.push(...kids.splice(10, 3)); // three moved to the end
    const b: Node = { kind: 'doc', children: kids };
    hashNode(a);
    hashNode(b);
    const p = bestOf(within, 40, () => diff(a, b));
    expect(p.ops).toHaveLength(3);
    const r = bestOf(within, 80, () => applyPatch(a, p, { maxOps: Number.POSITIVE_INFINITY }));
    expect(r.ok && hashNode(r.root)).toBe(hashNode(b));
  });

  itWithin('one move and one insert on a tree of 20,000 nodes and 20,000 refs, with a language', (within) => {
    const fresh = () => {
      const t = wide(200, 100);
      hashNode(t);
      return t;
    };
    const move = { op: 'move', at: ['b100'], parent: ['b101'], place: { last: true }, rename: 'm' } as const;
    const insert = {
      op: 'insert',
      parent: [],
      place: { last: true },
      node: { kind: 'box', id: 'n', params: { link: '/b5/l5' } },
    } as const;
    applyPatch(fresh(), { ops: [move] }, { language: PC }); // warm
    const moved = bestOf(within, 80, () => applyPatch(big, { ops: [move] }, { language: PC }));
    expect(moved.ok && moved.root.children![100].children![100].children![0].params?.link).toBe('/b101/m/l0');
    expect(moved.ok && moved.inverse?.ops).toHaveLength(1);
    const root = fresh();
    const inserted = bestOf(within, 80, () => applyPatch(root, { ops: [insert] }, { language: PC }));
    expect(inserted.ok).toBe(true);
  });

  itWithin('a diff of a parent of 10,000 children (moves, inserts, removes) applied', (within) => {
    const leaf = (id: string, label = id): Node => ({ kind: 'leaf', id, params: { label } });
    const a: Node = { kind: 'box', id: 'p', children: Array.from({ length: 10_000 }, (_, i) => leaf(`c${i}`)) };
    const kids = [...a.children!];
    for (let i = 0; i < 100; i++) kids.splice((i * 7919) % kids.length, 1); // removes
    for (let i = 0; i < 100; i++) kids.splice((i * 104729) % kids.length, 0, leaf(`n${i}`)); // inserts
    for (let i = 0; i < 100; i++) kids.push(...kids.splice((i * 15485863) % kids.length, 1)); // moves to the end
    const b: Node = { ...a, children: kids };
    hashNode(a);
    hashNode(b);
    const p = bestOf(within, 30, () => diff(a, b));
    expect(p.ops.length).toBeGreaterThan(250);
    const r = bestOf(within, 40, () => applyPatch(a, p, { maxOps: Number.POSITIVE_INFINITY }));
    expect(r.ok && hashNode(r.root)).toBe(hashNode(b));
  });
});
