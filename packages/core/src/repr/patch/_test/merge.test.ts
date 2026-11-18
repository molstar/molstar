/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The three-way merge (_spec/README.md "Merge"): the laws as seeded property tests, every conflict kind both
// ways round, dependents left out with their conflict, placement, orders, `invalid` with `partial` ours, a rebase of a
// stale patch whose move rewrites refs, and the merge's result consumed by `applyPatch` and `diff`.

import { itWithin } from '@molstar/test-harness/timing';
import { describe, expect, it } from 'vitest';
import {
  type ApplyOptions,
  applyPatch,
  type Conflict,
  diff,
  type MergeOptions,
  type MergeResult,
  merge3,
  type Op,
  type Patch,
  type PatchValidator,
  patchFromNode,
  patchLanguage,
  patchToNode,
} from '../../patch.js';
import { type Address, hashNode, type Issue, type Node, resolveAddress, type Value } from '../../tree/index.js';
import { PatchConformanceLanguage as PC } from '../conformance.js';

const languages = (id: string) => (id === 'pc.patch' ? patchLanguage(PC) : undefined);
const lang: MergeOptions & ApplyOptions = { language: PC, validate: { languages } };

type P = Record<string, Value>;
const leaf = (id: string, params?: P): Node => ({ kind: 'leaf', id, ...(params ? { params } : {}) });
const box = (id: string, params: P | undefined, ...children: Node[]): Node => ({
  kind: 'box',
  id,
  ...(params ? { params } : {}),
  ...(children.length ? { children } : {}),
});
const doc = (...children: Node[]): Node => ({ kind: 'doc', params: { title: 'T' }, children });
const base = () =>
  doc(box('a', { label: 'A' }, leaf('x'), leaf('y')), box('b', { label: 'B' }, leaf('x')), box('c', undefined));

/** `root` with the node at `at` changed by `f` (removed when it gives undefined). */
function update(root: Node, at: Address, f: (n: Node) => Node | undefined): Node {
  if (!at.length) return f(root)!;
  const [id, ...rest] = at;
  const children: Node[] = [];
  for (const c of root.children ?? []) {
    if (c.id !== id) children.push(c);
    else {
      const d = rest.length ? update(c, rest, f) : f(c);
      if (d) children.push(d);
    }
  }
  return { ...root, children };
}
const setP = (root: Node, at: Address, params: P) =>
  update(root, at, (n) => ({ ...n, params: { ...n.params, ...params } }));
const remove = (root: Node, at: Address) => update(root, at, () => undefined);
const kids = (root: Node, at: Address, f: (c: readonly Node[]) => Node[]) =>
  update(root, at, (n) => ({ ...n, children: f(n.children ?? []) }));
const ids = (root: Node, at: Address) => (resolveAddress(root, at)!.children ?? []).map((c) => c.id);

function failed(r: MergeResult): { conflicts: readonly Conflict[]; partial: Node } {
  expect(r.ok).toBe(false);
  if (r.ok) throw new Error('merged');
  return r;
}
function merged(r: MergeResult): Node {
  if (!r.ok) expect(r.conflicts).toEqual([]);
  if (!r.ok) throw new Error('conflicts');
  return r.root;
}
const shape = (c: Conflict) => ({
  kind: c.kind,
  at: c.at,
  ours: c.ours.map((o) => o.op),
  theirs: c.theirs.map((o) => o.op),
});

describe('repr/patch merge3: what both sides changed', () => {
  it('applies changes to different nodes and to different keys of one node, and equal changes once', () => {
    const b = base();
    const x = setP(setP(b, ['a'], { label: 'X' }), ['b', 'x'], { label: 'bx' });
    const y = setP(setP(b, ['a'], { size: 2 }), ['b', 'x'], { label: 'bx' });
    const root = merged(merge3(b, x, y));
    expect(hashNode(root)).toBe(hashNode(setP(x, ['a'], { size: 2 })));
  });

  it('says without a language that the merge was not validated, and gives the net change from ours', () => {
    const b = base();
    const x = setP(b, ['a'], { label: 'X' });
    const y = remove(b, ['c']);
    const r = merge3(b, x, y);
    expect(r.ok && r.issues.map((i) => [i.code, i.severity, i.hint])).toEqual([
      ['repr.patch.unvalidated', 'info', 'Give the language to merge with'],
    ]);
    expect(r.ok && r.changes.removed).toEqual([['c']]);
    expect(r.ok && r.changes.params).toEqual([]);
    const typed = merge3(b, x, y, lang);
    expect(typed.ok && typed.issues).toEqual([]);
  });

  it('gives ours itself when nothing of theirs applies', () => {
    const b = base();
    const x = setP(b, ['a'], { label: 'X' });
    expect(merge3(b, x, b).ok && (merge3(b, x, b) as { root: Node }).root).toBe(x);
    expect(merge3(b, x, x, lang).ok && (merge3(b, x, x, lang) as { root: Node }).root).toBe(x);
  });
});

describe('repr/patch merge3: conflicts, both ways round', () => {
  const both = (b: Node, x: Node, y: Node) => [failed(merge3(b, x, y)), failed(merge3(b, y, x))] as const;

  it('set: one key set to two values, or set and unset', () => {
    const b = base();
    const x = setP(b, ['a'], { label: 'X' });
    const y = setP(b, ['a'], { label: 'Y' });
    const [r, s] = both(b, x, y);
    expect(r.conflicts.map(shape)).toEqual([{ kind: 'set', at: ['a'], ours: ['set'], theirs: ['set'] }]);
    expect(r.partial).toBe(x);
    expect(s.partial).toBe(y);
    const u = update(b, ['a'], (n) => ({ ...n, params: {} }));
    const [r2, s2] = both(b, x, u);
    expect(r2.conflicts.map(shape)).toEqual([{ kind: 'set', at: ['a'], ours: ['set'], theirs: ['unset'] }]);
    expect(s2.conflicts.map(shape)).toEqual([{ kind: 'set', at: ['a'], ours: ['unset'], theirs: ['set'] }]);
  });

  it('remove: any change at or under what the other side removed', () => {
    const b = base();
    const gone = remove(b, ['a']);
    const changes: Node[] = [
      setP(b, ['a'], { label: 'X' }),
      setP(b, ['a', 'x'], { label: 'X' }),
      kids(b, ['a'], (c) => [...c, leaf('z')]),
      kids(b, ['a'], (c) => [...c].reverse()),
      update(b, ['a', 'x'], () => ({ kind: 'box', id: 'x' })),
    ];
    for (const y of changes) {
      const [r, s] = both(b, gone, y);
      expect(r.conflicts.map((c) => [c.kind, c.at, c.ours.map((o) => o.op)])).toEqual([['remove', ['a'], ['remove']]]);
      expect(r.conflicts[0].theirs).toEqual(diff(b, y).ops);
      expect(s.conflicts.map((c) => [c.kind, c.at, c.theirs.map((o) => o.op)])).toEqual([
        ['remove', ['a'], ['remove']],
      ]);
      expect(r.partial).toBe(gone);
      expect(s.partial).toBe(y);
    }
  });

  it('remove: removed on both sides, or under what the other removed, is no conflict', () => {
    const b = base();
    const gone = remove(b, ['a']);
    expect(merged(merge3(b, gone, gone))).toBe(gone);
    expect(merged(merge3(b, gone, remove(b, ['a', 'x'])))).toBe(gone);
    expect(hashNode(merged(merge3(b, remove(b, ['a', 'x']), gone)))).toBe(hashNode(gone));
  });

  it('insert: one id under one parent with other content; the same content once', () => {
    const b = base();
    const x = kids(b, ['c'], () => [leaf('z', { label: 'X' })]);
    const y = kids(b, ['c'], () => [leaf('z', { label: 'Y' })]);
    const [r, s] = both(b, x, y);
    expect(r.conflicts.map(shape)).toEqual([{ kind: 'insert', at: ['c', 'z'], ours: ['insert'], theirs: ['insert'] }]);
    expect(s.partial).toBe(y);
    expect(
      merged(
        merge3(
          b,
          x,
          kids(b, ['c'], () => [leaf('z', { label: 'X' })]),
        ),
      ),
    ).toBe(x);
  });

  it('replace: any change under it, or another replace of it', () => {
    const b = base();
    const x = update(b, ['a'], () => ({ kind: 'crate', id: 'a' }));
    for (const y of [setP(b, ['a', 'y'], { label: 'Y' }), setP(b, ['a'], { size: 1 }), remove(b, ['a', 'y'])]) {
      const [r, s] = both(b, x, y);
      expect(r.conflicts.map((c) => [c.kind, c.at, c.ours.map((o) => o.op)])).toEqual([
        ['replace', ['a'], ['replace']],
      ]);
      expect(s.conflicts.map((c) => [c.kind, c.at, c.theirs.map((o) => o.op)])).toEqual([
        ['replace', ['a'], ['replace']],
      ]);
    }
    const other = update(b, ['a'], () => ({ kind: 'crate', id: 'a', params: { label: 'Z' } }));
    const [r] = both(b, x, other);
    expect(r.conflicts.map(shape)).toEqual([{ kind: 'replace', at: ['a'], ours: ['replace'], theirs: ['replace'] }]);
    expect(
      merged(
        merge3(
          b,
          x,
          update(b, ['a'], () => ({ kind: 'crate', id: 'a' })),
        ),
      ),
    ).toBe(x);
    // a remove of what the other replaced is the remove row
    const [g, h] = both(b, remove(b, ['a']), x);
    expect(g.conflicts.map((c) => c.kind)).toEqual(['remove']);
    expect(h.conflicts.map((c) => c.kind)).toEqual(['remove']);
  });

  it('order: compared as orders over the children kept on both sides', () => {
    const b = kids(base(), ['a'], () => ['p', 'q', 'r', 's'].map((i) => leaf(i)));
    const order = (...o: string[]) => kids(b, ['a'], (c) => o.map((i) => c.find((n) => n.id === i)!));
    // one side: taken, whichever child diff chose to move
    const x = order('s', 'p', 'q', 'r');
    expect(ids(merged(merge3(b, b, x)), ['a'])).toEqual(['s', 'p', 'q', 'r']);
    expect(ids(merged(merge3(b, setP(b, ['c'], { size: 1 }), x)), ['a'])).toEqual(['s', 'p', 'q', 'r']);
    expect(ids(merged(merge3(b, x, setP(b, ['c'], { size: 1 }))), ['a'])).toEqual(['s', 'p', 'q', 'r']);
    // the same new order on both sides, once
    expect(merged(merge3(b, x, order('s', 'p', 'q', 'r')))).toBe(x);
    // a new order on one side over what the other removed: no conflict
    const y = order('q', 'p', 'r', 's');
    const removed = kids(b, ['a'], (c) => c.filter((n) => n.id !== 'p'));
    expect(ids(merged(merge3(b, removed, y)), ['a'])).toEqual(['q', 'r', 's']);
    // two different new orders: one conflict at the parent, partial keeps ours
    const [r, s] = both(b, x, y);
    expect(r.conflicts.map((c) => [c.kind, c.at])).toEqual([['order', ['a']]]);
    expect(r.conflicts[0].ours).toEqual(diff(b, x).ops);
    expect(r.conflicts[0].theirs).toEqual(diff(b, y).ops);
    expect(r.partial).toBe(x);
    expect(s.partial).toBe(y);
    // with an independent change of theirs beside it
    const z = setP(y, ['c'], { size: 3 });
    const t = failed(merge3(b, x, z));
    expect(ids(t.partial, ['a'])).toEqual(['s', 'p', 'q', 'r']);
    expect(resolveAddress(t.partial, ['c'])!.params).toEqual({ size: 3 });
  });

  it('order: one node both inserted, placed two ways, is one order conflict', () => {
    const b = doc(box('b', undefined), box('d', undefined));
    const c = box('c', { label: 'C' });
    const x = kids(b, [], (k) => [k[0], c, k[1]]); // b c d
    const y = kids(b, [], (k) => [k[0], k[1], c]); // b d c
    const r = failed(merge3(b, x, y));
    const s = failed(merge3(b, y, x));
    expect(r.conflicts.map(shape)).toEqual([{ kind: 'order', at: [], ours: ['insert'], theirs: ['insert'] }]);
    expect(s.conflicts.map(shape)).toEqual([{ kind: 'order', at: [], ours: ['insert'], theirs: ['insert'] }]);
    expect(r.partial).toBe(x);
    expect(s.partial).toBe(y);
    // placed alike: once, whichever side reordered the rest
    expect(
      merged(
        merge3(
          b,
          x,
          kids(b, [], (k) => [k[0], c, k[1]]),
        ),
      ),
    ).toBe(x);
  });
});

describe('repr/patch merge3: placement', () => {
  it('puts a child after its nearest predecessor in theirs that still stands, else first', () => {
    const b = base(); // a: x, y
    const theirs = kids(b, ['a'], (c) => [c[0], leaf('z'), c[1], leaf('w')]); // x z y w
    expect(ids(merged(merge3(b, remove(b, ['a', 'x']), theirs)), ['a'])).toEqual(['z', 'y', 'w']);
    expect(ids(merged(merge3(b, remove(b, ['a', 'y']), theirs)), ['a'])).toEqual(['x', 'z', 'w']);
    const both = kids(b, ['a'], () => []);
    expect(ids(merged(merge3(b, both, theirs)), ['a'])).toEqual(['z', 'w']);
  });

  it('keeps what ours inserted after its own predecessor, beside a new order of theirs', () => {
    const b = base();
    const x = kids(b, ['a'], (c) => [c[0], leaf('o'), c[1]]); // x o y
    const y = kids(b, ['a'], (c) => [c[1], c[0], leaf('t')]); // y x t
    expect(ids(merged(merge3(b, x, y)), ['a'])).toEqual(['y', 'x', 'o', 't']);
  });
});

describe('repr/patch merge3: dependents and the merged root', () => {
  it('leaves out a ref theirs wrote to a node only theirs inserted, when ours inserted another with that id', () => {
    const b = base();
    const x = kids(b, ['c'], () => [leaf('n', { label: 'ours' })]);
    let y = kids(b, ['c'], () => [leaf('n', { label: 'theirs' }), leaf('m')]);
    y = setP(y, ['b'], { link: '/c/n', size: 2 });
    y = setP(y, ['a'], { label: 'Y' });
    const r = failed(merge3(b, x, y, lang));
    expect(r.conflicts.map((c) => [c.kind, c.at])).toEqual([['insert', ['c', 'n']]]);
    // the insert, the sibling placed after it, the ref to it; not the size beside the ref
    expect(r.conflicts[0].theirs).toEqual([
      { op: 'set', at: ['b'], key: 'link', value: '/c/n' },
      { op: 'insert', parent: ['c'], place: { after: null }, node: leaf('n', { label: 'theirs' }) },
      { op: 'insert', parent: ['c'], place: { after: 'n' }, node: leaf('m') },
    ]);
    expect(hashNode(r.partial)).toBe(hashNode(setP(setP(x, ['a'], { label: 'Y' }), ['b'], { size: 2 })));
  });

  it('leaves out a ref naming what theirs inserted through an address nested in its form', () => {
    const b = doc(box('a', undefined), box('b', undefined, leaf('x')), box('c', undefined));
    const x = kids(b, ['c'], () => [leaf('n', { label: 'ours' })]);
    const y = setP(
      kids(b, ['c'], () => [leaf('n', { label: 'theirs' })]),
      ['a'],
      { link: '/b/(leaf :near /c/n)' },
    );
    const r = failed(merge3(b, x, y, lang));
    expect(r.conflicts.map(shape)).toEqual([
      { kind: 'insert', at: ['c', 'n'], ours: ['insert'], theirs: ['set', 'insert'] },
    ]);
    expect(r.partial).toBe(x);
  });

  it("makes no dependency of a stored patch's refs: they name nodes of its own dry run", () => {
    const b = base();
    const x = kids(b, ['c'], () => [leaf('n', { label: 'ours' })]);
    const step = patchToNode({
      ops: [
        { op: 'remove', at: ['c', 'n'] },
        { op: 'insert', parent: ['c'], place: { last: true }, node: box('n', { label: 'step' }) },
      ],
    });
    const y = setP(
      kids(b, ['c'], () => [leaf('n', { label: 'theirs' })]),
      ['a'],
      { step },
    );
    const r = failed(merge3(b, x, y, lang));
    expect(r.conflicts.map(shape)).toEqual([{ kind: 'insert', at: ['c', 'n'], ours: ['insert'], theirs: ['insert'] }]);
    expect(hashNode(r.partial)).toBe(hashNode(setP(x, ['a'], { step })));
  });

  it('refuses a merged root its validators refuse: two links making a cycle', () => {
    const cycle: PatchValidator = (_before, after) => {
      const links = new Map<string, string>();
      for (const c of after.children ?? [])
        if (typeof c.params?.link === 'string') links.set(`/${c.id}`, c.params.link);
      const out: Issue[] = [];
      for (const [from] of links) {
        let at = links.get(from);
        for (let n = 0; at !== undefined && n <= links.size; n++, at = links.get(at))
          if (at === from) {
            out.push({ code: 'pc.cycle', severity: 'error', message: 'A cycle', address: from });
            break;
          }
      }
      return out;
    };
    const b = base();
    const x = setP(b, ['a'], { link: '/b' });
    const y = setP(setP(b, ['b'], { link: '/a' }), ['c'], { size: 1 });
    const r = failed(merge3(b, x, y, { ...lang, validators: [cycle] }));
    expect(r.conflicts.map((c) => [c.kind, c.at, c.ours, c.theirs])).toEqual([['invalid', ['a'], [], diff(b, y).ops]]);
    expect(r.conflicts[0].issues!.map((i) => i.code)).toEqual(['pc.cycle', 'pc.cycle']);
    expect(r.partial).toBe(x);
    // each side alone is fine
    expect(merge3(b, b, y, { ...lang, validators: [cycle] }).ok).toBe(true);
  });

  it('refuses a ref of ours left dangling by theirs, unless allowed; `at` is the first error', () => {
    const b = base();
    const x = setP(b, ['c'], { link: '/b/x' });
    const y = remove(b, ['b', 'x']);
    const r = failed(merge3(b, x, y, lang));
    expect(r.conflicts.map((c) => [c.kind, c.at])).toEqual([['invalid', ['c']]]);
    expect(r.conflicts[0].issues!.map((i) => i.code)).toEqual(['repr.patch.dangling']);
    expect(r.partial).toBe(x);
    const allowed = merge3(b, x, y, { ...lang, allowDangling: true });
    expect(allowed.ok && allowed.issues.map((i) => [i.code, i.severity])).toEqual([['repr.patch.dangling', 'warning']]);
  });

  it("keeps partial ours when theirs' independent changes are refused, with an invalid conflict beside", () => {
    const b = base();
    const x = setP(setP(b, ['c'], { link: '/b/x' }), ['a'], { label: 'X' });
    const y = setP(remove(b, ['b', 'x']), ['a'], { label: 'Y' });
    const r = failed(merge3(b, x, y, lang));
    expect(r.conflicts.map((c) => c.kind)).toEqual(['set', 'invalid']);
    expect(r.conflicts[1].theirs).toEqual([{ op: 'remove', at: ['b', 'x'] }]);
    expect(r.partial).toBe(x);
  });
});

describe('repr/patch merge3: rebase', () => {
  it("rebases a stale patch, its move's ref rewrites following", () => {
    const b = doc(box('a', undefined, leaf('x', { label: 'X' })), box('b', undefined), box('c', { link: '/a/x' }));
    const p: Patch = {
      ops: [{ op: 'move', at: ['a', 'x'], parent: ['b'], place: { last: true }, rename: 'w' }],
      base: hashNode(b),
    };
    const current = setP(b, ['a'], { label: 'A' });
    expect(applyPatch(current, p, lang).ok).toBe(false); // stale
    const theirs = applyPatch(b, p, lang);
    expect(theirs.ok).toBe(true);
    if (!theirs.ok) return;
    const root = merged(merge3(b, current, theirs.root, lang));
    expect(resolveAddress(root, ['c'])!.params).toEqual({ link: '/b/w' });
    expect(resolveAddress(root, ['a'])!.params).toEqual({ label: 'A' });
    const fresh = applyPatch(current, { ops: p.ops }, lang);
    expect(fresh.ok && hashNode(fresh.root)).toBe(hashNode(root));
  });
});

describe('repr/patch merge3: its result consumed', () => {
  it('by applyPatch (its changes and root) and by diff (the rebased patch read back)', () => {
    const b = base();
    const x = setP(
      kids(b, ['a'], (c) => [c[1], c[0]]),
      ['c'],
      { link: '/a/y' },
    );
    const y = kids(setP(remove(b, ['b', 'x']), ['a', 'x'], { link: '/c' }), ['c'], () => [leaf('n')]);
    for (const options of [{}, lang]) {
      const r = merge3(b, x, y, options);
      const root = merged(r);
      if (!r.ok) return;
      const again = applyPatch(x, diff(x, root), { ...options, maxOps: Number.POSITIVE_INFINITY });
      expect(again.ok && hashNode(again.root)).toBe(hashNode(root));
      expect(again.ok && again.changes).toEqual(r.changes);
      const read = patchFromNode(patchToNode(diff(x, root)));
      expect(read.issues).toEqual([]);
      const back = applyPatch(x, read.patch!, options);
      expect(back.ok && hashNode(back.root)).toBe(hashNode(root));
    }
  });
});

// the laws, seeded

function rng(seed: number) {
  let s = seed >>> 0;
  return () => {
    s = (s + 0x6d2b79f5) >>> 0;
    let t = s;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function addresses(root: Node): Address[] {
  const out: Address[] = [[]];
  const go = (n: Node, at: Address) => {
    for (const c of n.children ?? []) {
      out.push([...at, c.id!]);
      go(c, [...at, c.id!]);
    }
  };
  go(root, []);
  return out;
}

const Pool = ['a', 'b', 'c', 'd', 'e'];
const Labels = ['p', 'q', 'r'];

/** One random edit, valid in the conformance language: labels, links, inserts, removes, reorders. */
function editOnce(next: () => number, root: Node): Node {
  const pick = <T>(xs: readonly T[]) => xs[Math.floor(next() * xs.length)];
  const all = addresses(root);
  const at = pick(all);
  const node = resolveAddress(root, at)!;
  switch (Math.floor(next() * 6)) {
    case 0:
      return at.length ? setP(root, at, { label: pick(Labels) }) : root;
    case 1:
      return at.length
        ? update(root, at, (n) => {
            const { label: _, ...rest } = n.params ?? {};
            return { ...n, params: rest };
          })
        : root;
    case 2: {
      if (node.kind === 'leaf' || at.length > 2) return root;
      const taken = new Set((node.children ?? []).map((c) => c.id));
      const free = Pool.filter((i) => !taken.has(i));
      if (!free.length) return root;
      const id = pick(free);
      const child = at.length && next() < 0.5 ? leaf(id, { label: pick(Labels) }) : box(id, { label: pick(Labels) });
      const i = Math.floor(next() * ((node.children?.length ?? 0) + 1));
      return kids(root, at, (c) => [...c.slice(0, i), child, ...c.slice(i)]);
    }
    case 3:
      return at.length ? remove(root, at) : root;
    case 4: {
      const n = node.children?.length ?? 0;
      if (n < 2) return root;
      const i = Math.floor(next() * n);
      const j = Math.floor(next() * n);
      return kids(root, at, (c) => {
        const d = [...c];
        [d[i], d[j]] = [d[j], d[i]];
        return d;
      });
    }
    default:
      return at.length ? setP(root, at, { link: `/${pick(all).join('/')}` }) : root;
  }
}

function edits(next: () => number, root: Node, n: number) {
  let r = root;
  for (let i = 0; i < n; i++) r = editOnce(next, r);
  return r;
}

/**
 * `n` (at `at`) with the children only one side inserted moved last, sorted by id, and the rest kept in order:
 * equality up to the order of siblings one side alone inserted.
 */
function oneSideLast(n: Node, b: Node, x: Node, y: Node, at: Address = []): Node {
  if (!n.children) return n;
  const has = (root: Node, id: string) => !!resolveAddress(root, [...at, id]);
  const alone = (c: Node) => !has(b, c.id!) && !(has(x, c.id!) && has(y, c.id!));
  const go = (c: Node) => oneSideLast(c, b, x, y, [...at, c.id!]);
  const last = n.children.filter(alone).sort((p, q) => (p.id! < q.id! ? -1 : 1));
  return { ...n, children: [...n.children.filter((c) => !alone(c)), ...last].map(go) };
}

/** `x` and `y` with one node inserted on both, the same content, each at a place of its own, under a box both have. */
function insertBoth(next: () => number, x: Node, y: Node): [Node, Node] {
  const boxes = addresses(x).filter(
    (a) => resolveAddress(x, a)?.kind !== 'leaf' && resolveAddress(y, a)?.kind !== 'leaf',
  );
  const at = boxes[Math.floor(next() * boxes.length)];
  const node = at.length ? leaf('q', { label: 'both' }) : box('q', { label: 'both' });
  const put = (root: Node) =>
    resolveAddress(root, [...at, 'q'])
      ? root
      : kids(root, at, (c) => {
          const i = Math.floor(next() * (c.length + 1));
          return [...c.slice(0, i), node, ...c.slice(i)];
        });
  return [put(x), put(y)];
}

const Seeds = 150;

describe('repr/patch merge3: the laws, seeded', () => {
  const start = (next: () => number) => edits(next, doc(box('a', undefined), box('b', undefined)), 12);

  it('merge3(b, x, b) is x, merge3(b, b, y) is y, merge3(b, x, x) is x', () => {
    for (let seed = 1; seed <= Seeds; seed++) {
      const next = rng(seed);
      const b = start(next);
      const x = edits(next, b, 1 + Math.floor(next() * 5));
      expect(merged(merge3(b, x, b)), `seed ${seed}`).toBe(x);
      expect(hashNode(merged(merge3(b, b, x))), `seed ${seed}`).toBe(hashNode(x));
      expect(merged(merge3(b, x, x)), `seed ${seed}`).toBe(x);
    }
  });

  it('is symmetric without conflicts, up to the order of siblings one side inserted', () => {
    let clean = 0;
    for (let seed = 1; seed <= Seeds; seed++) {
      const next = rng(seed);
      const b = start(next);
      let x = edits(next, b, 1 + Math.floor(next() * 4));
      let y = edits(next, b, 1 + Math.floor(next() * 4));
      if (seed % 2 === 0) [x, y] = insertBoth(next, x, y);
      const r = merge3(b, x, y);
      const s = merge3(b, y, x);
      expect(r.ok, `seed ${seed}`).toBe(s.ok);
      if (!r.ok || !s.ok) {
        const kinds = (m: MergeResult) => (m.ok ? [] : m.conflicts.map((c) => c.kind).sort());
        expect(kinds(r), `seed ${seed}`).toEqual(kinds(s));
        continue;
      }
      clean++;
      expect(hashNode(oneSideLast(r.root, b, x, y)), `seed ${seed}`).toBe(hashNode(oneSideLast(s.root, b, x, y)));
    }
    expect(clean).toBeGreaterThan(Seeds / 5);
  });

  it('gives a partial (and a root) valid under the options, with a language', () => {
    const options = { ...lang, allowDangling: true, allowRebind: true };
    let conflicted = 0;
    for (let seed = 1; seed <= Seeds; seed++) {
      const next = rng(seed);
      const b = start(next);
      const x = edits(next, b, 1 + Math.floor(next() * 4));
      const y = edits(next, b, 1 + Math.floor(next() * 4));
      const r = merge3(b, x, y, options);
      const out = r.ok ? r.root : r.partial;
      if (!r.ok) conflicted++;
      const again = applyPatch(x, diff(x, out), { ...options, maxOps: Number.POSITIVE_INFINITY });
      expect(again.ok, `seed ${seed}`).toBe(true);
      if (!r.ok && r.conflicts.some((c) => c.kind === 'invalid')) expect(r.partial).toBe(x);
    }
    expect(conflicted).toBeGreaterThan(0);
  });
});

// what a merge costs: two edits of a tree of 20,000 nodes. Bounds are one sample, a few times what a case takes alone.

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

describe('repr/patch merge3: what a merge costs', () => {
  itWithin('two edits of a tree of 20,000 nodes and 20,000 refs, with a language', (within) => {
    const b = wide(200, 100);
    hashNode(b);
    const ops = (k: number): Op[] => [
      { op: 'set', at: [`b${k}`, 'l5'], key: 'label', value: 'x' },
      { op: 'remove', at: [`b${k + 1}`, 'l7'] },
      { op: 'insert', parent: [`b${k + 2}`], place: { after: 'l3' }, node: leaf('n', { link: `/b${k}/l1` }) },
      { op: 'move', at: [`b${k + 3}`, 'l9'], parent: [`b${k + 3}`], place: { after: null } },
    ];
    const x = applyPatch(b, { ops: ops(10) }, lang);
    const y = applyPatch(b, { ops: ops(100) }, lang);
    if (!x.ok || !y.ok) throw new Error('edits');
    merge3(b, x.root, y.root, lang); // warm
    // About 0.5 ms warm (the diffs skip equal subtrees by hash); one sample, so 5 ms.
    const r = within(5, () => merge3(b, x.root, y.root, lang));
    expect(r.ok).toBe(true);
    if (!r.ok) return;
    expect(r.changes.made).toEqual([['b102', 'n']]);
    expect(r.root.children![50]).toBe(b.children![50]);
  });
});
