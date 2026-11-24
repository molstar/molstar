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
import { applyPatch, diff, type Op, patchLanguage, patchToNode } from '../../patch.js';
import { hashNode, type Node } from '../../tree/index.js';
import { PatchConformanceLanguage as PC } from '../conformance.js';
import { refIndexOf } from '../refs.js';

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

  itWithin('edits of a p53-size spec with 30 stored patches cost no walk of the steps', (within) => {
    // data: 12 boxes of 25 leaves; views: 10 boxes of 15 leaves naming the data; a story of 30 steps of 4 ops each.
    // A stored `set` value is typed by its key alone, so no edit types the steps again (the old root's index is
    // derived, and a reorder builds none)
    const options = {
      language: PC,
      validate: { languages: (id: string) => (id === 'pc.patch' ? patchLanguage(PC) : undefined) },
    };
    const data = Array.from({ length: 12 }, (_, i) => ({
      kind: 'box',
      id: `d${i}`,
      children: Array.from({ length: 25 }, (_, j) => ({
        kind: 'leaf',
        id: `f${j}`,
        params: { link: `/data/d${i}/f0` },
      })),
    }));
    const views = Array.from({ length: 10 }, (_, i) => ({
      kind: 'box',
      id: `v${i}`,
      params: { link: `/data/d${i}`, links: [`/data/d${i}/f1`, `/data/d${(i + 1) % 12}`] },
      children: Array.from({ length: 15 }, (_, j) => ({
        kind: 'leaf',
        id: `t${j}`,
        params: { link: `/data/d${i}/f${j}` },
      })),
    }));
    const story = Array.from({ length: 30 }, (_, i) => {
      const v = `v${i % 10}`;
      const w = `v${(i + 3) % 10}`;
      const ops: Op[] = [
        { op: 'set', at: ['views', v, 't0'], key: 'link', value: `/data/d${i % 12}/f${i % 25}` },
        { op: 'move', at: ['views', v, 't1'], parent: ['views', v], place: { last: true }, rename: `r${i}` },
        {
          op: 'insert',
          parent: ['views', w],
          place: { last: true },
          node: { kind: 'leaf', id: `n${i}`, params: { link: `/data/d${(i + 5) % 12}/f3` } },
        },
        { op: 'set', at: ['views', w, `n${i}`], key: 'link', value: `/data/d${(i + 7) % 12}` },
      ];
      return { kind: 'box', id: `s${i}`, params: { step: patchToNode({ ops }) } };
    });
    const spec: Node = {
      kind: 'doc',
      children: [
        { kind: 'box', id: 'data', children: data },
        { kind: 'box', id: 'views', children: views },
        { kind: 'box', id: 'story', children: story },
      ],
    };
    hashNode(spec);
    const edits: Record<string, Op[]> = {
      rename: [{ op: 'move', at: ['data', 'd4'], parent: ['data'], place: { last: true }, rename: 'q4' }],
      insert: [
        {
          op: 'insert',
          parent: ['views'],
          place: { last: true },
          node: { kind: 'box', id: 'nv', params: { link: '/data/d1' } },
        },
      ],
      reorders: Array.from(
        { length: 10 },
        (_, i) =>
          ({ op: 'move', at: ['data', `d${i}`, 'f20'], parent: ['data', `d${i}`], place: { after: null } }) as Op,
      ),
      renames: Array.from(
        { length: 10 },
        (_, i) =>
          ({
            op: 'move',
            at: ['data', `d${i}`, 'f20'],
            parent: ['data', `d${i}`],
            place: { last: true },
            rename: `h${i}`,
          }) as Op,
      ),
    };
    // an editing session: the spec's index exists, each edit derives the next
    refIndexOf(spec, PC, options.validate);
    const opts = { ...options, allowDangling: true, allowRebind: true };
    for (const ops of Object.values(edits)) applyPatch(spec, { ops }, opts); // warm
    for (const [name, ops] of Object.entries(edits)) {
      const r = bestOf(within, 4, () => applyPatch(spec, { ops }, opts));
      expect(r.ok, name).toBe(true);
    }
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
