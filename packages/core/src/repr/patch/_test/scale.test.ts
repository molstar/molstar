/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A patch costs the paths it touches, not the tree: a one-key set on a large tree, with and without a language, and a
// diff of a wide parent reordered (a longest increasing subsequence, O(n log n)). Quick without the clock, heavy with
// it (itWithin).

import { itWithin } from '@molstar/test-harness/timing';
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
    const p = within(20, () => diff(a, b));
    expect(p.ops).toHaveLength(3);
    const r = within(40, () => applyPatch(a, p, { maxOps: Number.POSITIVE_INFINITY }));
    expect(r.ok && hashNode(r.root)).toBe(hashNode(b));
  });
});
