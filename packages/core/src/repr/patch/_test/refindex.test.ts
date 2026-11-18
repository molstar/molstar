/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ref index across revisions (_spec/README.md "Patches", Refs): a spec of 10,000 refs edited by moves and one-key
// sets, each applied to the previous result: what building the old root's index costs against deriving it, and what
// one move costs with the index derived. A move there rewrites 100 refs and takes about 5 ms, a set about 2 ms (the
// index: about 0.1 ms to derive, 0.2 ms a query); each bound is a sample at 6-10x that.

import { itWithin } from '@molstar/test-harness/timing';
import { describe, expect, it } from 'vitest';
import { applyPatch, patchLanguage, patchToNode } from '../../patch.js';
import { hashNode, type Node } from '../../tree/index.js';
import { PatchConformanceLanguage as PC } from '../conformance.js';
import { freshRefIndex, peekIndex, refIndexOf, type Site } from '../refs.js';

/** `targets` holds `n` boxes of one leaf; `refs` holds `refs` leaves, each linking to a target's leaf. */
function bookmarks(n: number, refs: number): Node {
  const root: Node = {
    kind: 'doc',
    children: [
      {
        kind: 'box',
        id: 'targets',
        children: Array.from({ length: n }, (_, i) => ({
          kind: 'box',
          id: `t${i}`,
          children: [{ kind: 'leaf', id: 'l', params: { label: `${i}` } }],
        })),
      },
      {
        kind: 'box',
        id: 'refs',
        children: Array.from({ length: refs }, (_, i) => ({
          kind: 'leaf',
          id: `r${i}`,
          params: { label: `${i}`, link: `/targets/t${i % n}/l` },
        })),
      },
    ],
  };
  hashNode(root);
  return root;
}

const move = (i: number) =>
  ({
    op: 'move',
    at: ['targets', `t${i}`],
    parent: ['targets'],
    place: { last: true },
    rename: `m${i}`,
  }) as const;
const setLabel = (i: number) => ({ op: 'set', at: ['refs', `r${i}`], key: 'label', value: `x${i}` }) as const;

const same = (a: readonly Site[], b: readonly Site[]) => {
  expect(a.length).toBe(b.length);
  a.forEach((s, i) => {
    expect(s.node).toBe(b[i].node);
    expect([s.value, s.path, s.at, s.anchor, s.tail]).toEqual([b[i].value, b[i].path, b[i].at, b[i].anchor, b[i].tail]);
  });
};

describe('repr/patch: the ref index across revisions', () => {
  it('100 moves and 100 one-key sets on 10,000 refs: building the old root index against deriving it', {
    tags: ['heavy'],
  }, () => {
    const n = 100;
    let root = bookmarks(100, 10_000);
    refIndexOf(root, PC); // the first revision's, built once
    const stats = { move: { build: 0, apply: 0, n: 0 }, set: { build: 0, apply: 0, n: 0 } };
    for (let i = 0; i < 2 * n; i++) {
      const isMove = i % 2 === 0;
      const op = isMove ? move(i / 2) : setLabel(i * 31);
      // what the apply cost before: the old root's index built from scratch
      const t0 = performance.now();
      freshRefIndex(root, PC);
      const t1 = performance.now();
      const r = applyPatch(root, { ops: [op] }, { language: PC });
      const t2 = performance.now();
      expect(r.ok).toBe(true);
      if (!r.ok) return;
      const s = stats[isMove ? 'move' : 'set'];
      s.build += t1 - t0;
      s.apply += t2 - t1;
      s.n++;
      expect(peekIndex(r.root, PC), `step ${i}`).toBeDefined();
      root = r.root;
    }
    const f = (x: number, k: number) => (x / k).toFixed(2);
    console.log(
      `ref index, per apply (ms): move: old root's index built ${f(stats.move.build, n)}, apply with it derived ${f(stats.move.apply, n)}; ` +
        `set: built ${f(stats.set.build, n)}, apply ${f(stats.set.apply, n)}`,
    );
    // the chain's last index equals a fresh build
    same(refIndexOf(root, PC).under([]) as Site[], freshRefIndex(root, PC).under([]) as Site[]);
  });

  itWithin('one move on 10,000 refs with the old root index derived, not built', (within) => {
    let root = bookmarks(100, 10_000);
    refIndexOf(root, PC);
    for (let i = 0; i < 3; i++) {
      // warm, and a chain: every index but the first derived
      const r = applyPatch(root, { ops: [setLabel(i)] }, { language: PC });
      expect(r.ok && r.root).toBeTruthy();
      if (r.ok) root = r.root;
    }
    for (let i = 0; i < 5; i++) {
      applyPatch(root, { ops: [move(i)] }, { language: PC });
      applyPatch(root, { ops: [setLabel(i + 9)] }, { language: PC });
    } // warm
    const moved = within(30, () => applyPatch(root, { ops: [move(50)] }, { language: PC }));
    expect(moved.ok && moved.root.children![1].children![50].params?.link).toBe('/targets/m50/l');
    expect(moved.ok && peekIndex(moved.root, PC)).toBeDefined();
    const set = within(20, () => applyPatch(root, { ops: [setLabel(77)] }, { language: PC }));
    expect(set.ok).toBe(true);
  });

  it('a stored patch whose set value is typed by a node another op replaces: the derived index follows', () => {
    const languages = (id: string) => (id === 'pc.patch' ? patchLanguage(PC) : undefined);
    const validate = { languages };
    const step = patchToNode({ ops: [{ op: 'set', at: ['b'], key: 'link', value: '/a' }] });
    const root: Node = {
      kind: 'doc',
      children: [
        { kind: 'box', id: 'a', params: { step } },
        { kind: 'box', id: 'b' },
      ],
    };
    hashNode(root);
    const first = refIndexOf(root, PC, validate);
    const before = first.under([]).length;
    // `b` becomes a crate: a `set` of `link` on it is typed by no param, so its ref is no longer one
    const r = applyPatch(
      root,
      { ops: [{ op: 'replace', at: ['b'], node: { kind: 'crate', id: 'b' } }] },
      { language: PC, validate },
    );
    expect(r.ok).toBe(true);
    if (!r.ok) return;
    const derived = peekIndex(r.root, PC, validate)!;
    expect(derived).toBeDefined();
    const fresh = freshRefIndex(r.root, PC, validate);
    expect(derived.under([]).length).toBe(fresh.under([]).length);
    expect(derived.under([]).length).toBeLessThan(before);
  });
});
