/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ref index across revisions (_spec/README.md "Patches", Refs): specs of 10,000 refs edited by moves, renames and
// one-key sets, each applied to the previous result: what building the old root's index costs against deriving it,
// what one edit costs with the index derived, with stored patches whose typed `set` reads the root, and with refs
// spread over many first ids. A move on the first spec rewrites 100 refs and takes about 5 ms, a set about 2 ms; each
// bound is a sample at 5-10x that, the best of three.

import { bestOf, itWithin } from '@molstar/test-harness/timing';
import { describe, expect, it } from 'vitest';
import { applyPatch, type Op, patchLanguage, patchToNode } from '../../patch.js';
import { defineLanguage } from '../../schema/check.js';
import type { Language } from '../../schema/language.js';
import { hashNode, type Node } from '../../tree/index.js';
import { PatchConformanceLanguage as PC } from '../conformance.js';
import { freshRefIndex, peekIndex, refIndexOf, type Site } from '../refs.js';

const languages = (id: string) => (id === 'pc.patch' ? patchLanguage(PC) : undefined);
const validate = { languages };

/** `targets` holds `n` boxes of one leaf; `refs` holds `refs` leaves, each linking to a target's leaf. */
function bookmarks(n: number, refs: number, stories = 0): Node {
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
      // a story: stored patches whose `set` is typed through the node it writes, in the wide `refs`
      {
        kind: 'box',
        id: 'story',
        children: Array.from({ length: stories }, (_, i) => ({
          kind: 'box',
          id: `s${i}`,
          params: {
            step: patchToNode({
              ops: [{ op: 'set', at: ['refs', `r${i}`], key: 'link', value: `/targets/t${i % n}` }],
            }),
          },
        })),
      },
    ],
  };
  hashNode(root);
  return root;
}

/** `count` top-level boxes of `each` leaves, each linking to its box's next leaf: the refs over many first ids. */
function groups(count: number, each: number): Node {
  const root: Node = {
    kind: 'doc',
    children: Array.from({ length: count }, (_, g) => ({
      kind: 'box',
      id: `g${g}`,
      children: Array.from({ length: each }, (_, j) => ({
        kind: 'leaf',
        id: `l${j}`,
        params: { link: `/g${g}/l${(j + 1) % each}` },
      })),
    })),
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
/** One patch of 50 renames of leaves of group `g`, each rewriting the ref that names it. */
const renames = (g: number): Op[] =>
  Array.from({ length: 50 }, (_, j) => ({
    op: 'move',
    at: [`g${g}`, `l${2 * j}`],
    parent: [`g${g}`],
    place: { last: true },
    rename: `m${j}`,
  }));
const allow = { language: PC, allowDangling: true, allowRebind: true } as const;
const setLabel = (i: number) => ({ op: 'set', at: ['refs', `r${i}`], key: 'label', value: `x${i}` }) as const;

const same = (a: readonly Site[], b: readonly Site[]) => {
  expect(a.length).toBe(b.length);
  a.forEach((s, i) => {
    expect(s.node).toBe(b[i].node);
    expect([s.value, s.path, s.at, s.anchor, s.tail, s.rel, s.runs]).toEqual([
      b[i].value,
      b[i].path,
      b[i].at,
      b[i].anchor,
      b[i].tail,
      b[i].rel,
      b[i].runs,
    ]);
  });
};

describe('repr/patch: the ref index across revisions', () => {
  it('100 moves and 100 one-key sets on 10,000 refs: building the old root index against deriving it', {
    tags: ['heavy'],
  }, () => {
    const n = 100;
    let root = bookmarks(100, 10_000);
    refIndexOf(root, PC); // the first revision's, built once
    const stats = { move: { build: 0, apply: 0, derive: 0 }, set: { build: 0, apply: 0, derive: 0 } };
    for (let i = 0; i < 2 * n; i++) {
      const isMove = i % 2 === 0;
      const op = isMove ? move(i / 2) : setLabel(i * 31);
      // what an apply cost before: the old root's index built from scratch
      const t0 = performance.now();
      freshRefIndex(root, PC);
      const t1 = performance.now();
      const r = applyPatch(root, { ops: [op] }, { language: PC });
      const t2 = performance.now();
      expect(r.ok).toBe(true);
      if (!r.ok) return;
      // the derive's own time: the new root's index read
      expect(peekIndex(r.root, PC), `step ${i}`).toBeDefined();
      const t3 = performance.now();
      const s = stats[isMove ? 'move' : 'set'];
      s.build += t1 - t0;
      s.apply += t2 - t1;
      s.derive += t3 - t2;
      root = r.root;
    }
    const f = (x: number) => (x / n).toFixed(2);
    console.log(
      `ref index, per apply (ms): move: old root's index built ${f(stats.move.build)}, apply ${f(stats.move.apply)}, derive ${f(stats.move.derive)}; ` +
        `set: built ${f(stats.set.build)}, apply ${f(stats.set.apply)}, derive ${f(stats.set.derive)}`,
    );
    // the chain's last index equals a fresh build
    same(refIndexOf(root, PC).under([]) as Site[], freshRefIndex(root, PC).under([]) as Site[]);
  });

  it('a patch of 50 renames over 100 groups of 100 refs: refs over many first ids', { tags: ['heavy'] }, () => {
    let root = groups(100, 100);
    refIndexOf(root, PC);
    let apply = 0;
    for (let i = 0; i < 50; i++) {
      const t0 = performance.now();
      const r = applyPatch(root, { ops: renames(i) }, allow);
      apply += performance.now() - t0;
      expect(r.ok).toBe(true);
      if (!r.ok) return;
      root = r.root;
    }
    console.log(`ref index, 100 groups of 100 refs: a patch of 50 renames ${(apply / 50).toFixed(2)} ms per apply`);
  });

  it('one-key sets on 10,000 refs and 30 stored patches with typed sets', { tags: ['heavy'] }, () => {
    let root = bookmarks(100, 10_000, 30);
    refIndexOf(root, PC, validate);
    let apply = 0;
    for (let i = 0; i < 50; i++) {
      const t0 = performance.now();
      const r = applyPatch(root, { ops: [setLabel(i)] }, { language: PC, validate });
      apply += performance.now() - t0;
      expect(r.ok).toBe(true);
      if (!r.ok) return;
      root = r.root;
    }
    console.log(`ref index, 10,000 refs and 30 typed-set steps: a one-key set ${(apply / 50).toFixed(2)} ms per apply`);
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
    const moved = bestOf(within, 30, () => applyPatch(root, { ops: [move(50)] }, { language: PC }));
    expect(moved.ok && moved.root.children![0].children!.at(-1)!.id).toBe('m50');
    expect(moved.ok && moved.root.children![1].children![50].params?.link).toBe('/targets/m50/l');
    expect(moved.ok && peekIndex(moved.root, PC)).toBeDefined();
    const set = bestOf(within, 20, () => applyPatch(root, { ops: [setLabel(77)] }, { language: PC }));
    expect(set.ok).toBe(true);
  });

  itWithin('one-key set on 10,000 refs and 30 stored patches with typed sets', (within) => {
    let root = bookmarks(100, 10_000, 30);
    refIndexOf(root, PC, validate);
    for (let i = 0; i < 8; i++) {
      // a chain: a set per apply, the index never read in between; warm
      const r = applyPatch(root, { ops: [setLabel(i)] }, { language: PC, validate });
      expect(r.ok).toBe(true);
      if (r.ok) root = r.root;
    }
    const set = bestOf(within, 12, () => applyPatch(root, { ops: [setLabel(77)] }, { language: PC, validate }));
    expect(set.ok).toBe(true);
  });

  itWithin('a patch of 50 renames over 100 groups of 100 refs, many first ids', (within) => {
    let root = groups(100, 100);
    refIndexOf(root, PC);
    for (let i = 0; i < 8; i++) {
      const r = applyPatch(root, { ops: renames(i) }, allow);
      expect(r.ok).toBe(true);
      if (r.ok) root = r.root;
    }
    const renamed = bestOf(within, 12, () => applyPatch(root, { ops: renames(50) }, allow));
    expect(renamed.ok && renamed.root.children![50].children!.at(-1)!.id).toBe('m49');
  });

  describe('the chunks reused only at the position they were built', () => {
    const rename = (at: string, to: string) =>
      ({ op: 'move', at: [at], parent: [], place: { after: null }, rename: to }) as const;
    /** The root's index (built), the patch applied, and the next root's index against one built on it. */
    const derives = (root: Node, op: Op, language: Language = PC) => {
      hashNode(root);
      refIndexOf(root, language);
      const r = applyPatch(root, { ops: [op] }, { language, allowDangling: true, allowRebind: true });
      expect(r.ok).toBe(true);
      if (!r.ok) throw new Error('refused');
      const derived = peekIndex(r.root, language);
      expect(derived).toBeDefined();
      const a = derived!.under([]) as Site[];
      same(a, freshRefIndex(r.root, language).under([]) as Site[]);
      return a;
    };

    it('a renamed parent in place: the same path, another address', () => {
      const sites = derives(
        {
          kind: 'doc',
          children: [
            { kind: 'box', id: 'a', children: [{ kind: 'leaf', id: 'x', params: { link: '/b' } }] },
            { kind: 'box', id: 'b' },
          ],
        },
        rename('a', 'c'),
      );
      expect(sites).toHaveLength(1);
      expect(sites[0].at).toEqual(['c', 'x']);
    });

    it('an id-less node under a renamed parent: the same child position, another anchor', () => {
      const L = defineLanguage({
        id: 'ids',
        version: '1.0.0',
        root: 'doc',
        doc: 'Optional ids.',
        kinds: {
          doc: { doc: 'The root.', parents: [] },
          box: { doc: 'A box.', parents: ['doc'], ids: 'required', params: { link: { type: 'ref', doc: 'A ref.' } } },
          note: { doc: 'A note.', parents: ['box'], ids: 'optional', params: { link: { type: 'ref', doc: 'A ref.' } } },
        },
      });
      const sites = derives(
        {
          kind: 'doc',
          children: [
            { kind: 'box', id: 'a', children: [{ kind: 'note', params: { link: '/b' } }] },
            { kind: 'box', id: 'b' },
          ],
        },
        rename('a', 'c'),
        L,
      );
      expect(sites).toHaveLength(1);
      expect(sites[0].anchor).toEqual(['c']);
      expect(sites[0].tail).toEqual([0]);
    });

    it('a node object shared at two positions', () => {
      const shared: Node = { kind: 'leaf', id: 'x', params: { link: '/a' } };
      const sites = derives(
        {
          kind: 'doc',
          children: [
            { kind: 'box', id: 'a', children: [shared] },
            { kind: 'box', id: 'b', children: [shared] },
          ],
        },
        { op: 'set', at: ['a'], key: 'label', value: 'L' },
      );
      expect(sites.map((s) => s.at)).toEqual([
        ['a', 'x'],
        ['b', 'x'],
      ]);
    });
  });

  it("a stored patch's set value stays a ref when its target changes kind: typed by key, not by the spec", () => {
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
    // `b` becomes a crate, which has no `link`: the step's value is still a ref, typed by the language's `link`
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
    expect(derived.under([]).length).toBe(before);
  });
});
