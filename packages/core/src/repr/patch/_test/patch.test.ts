/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { parseJson, printJson } from '../../codecs/json/index.js';
import { parseSexpr, printSexpr } from '../../codecs/sexpr/index.js';
import {
  AnyPatchLanguage,
  type ApplyOptions,
  applyPatch,
  diff,
  type Patch,
  patchFromNode,
  patchLanguage,
  patchToNode,
} from '../../patch.js';
import { validate } from '../../schema/index.js';
import { hashNode, type Issue, type Node, type Value } from '../../tree/index.js';
import { PatchConformanceLanguage as PC } from '../conformance.js';
import { longestIncreasing } from '../diff.js';

const leaf = (id: string, params?: Node['params']): Node => ({ kind: 'leaf', id, ...(params ? { params } : {}) });
const box = (id: string, params: Node['params'] | undefined, ...children: Node[]): Node => ({
  kind: 'box',
  id,
  ...(params ? { params } : {}),
  ...(children.length ? { children } : {}),
});
const doc = (...children: Node[]): Node => ({ kind: 'doc', children });
const languages = (id: string) => (id === 'pc.patch' ? patchLanguage(PC) : undefined);
const lang: ApplyOptions = { language: PC, validate: { languages } };

function applied(root: Node, patch: Patch, options?: ApplyOptions) {
  const r = applyPatch(root, patch, options);
  if (!r.ok) throw new Error(r.issues.map((i) => i.code).join());
  return r;
}

describe('repr/patch: Changes are net', () => {
  const root = doc(box('a', { label: 'A' }, leaf('x'), leaf('y')), box('b', undefined, leaf('z')));

  it('an insert then a move is one made node, at its last address', () => {
    const r = applied(root, {
      ops: [
        { op: 'insert', parent: ['a'], place: { last: true }, node: leaf('n') },
        { op: 'move', at: ['a', 'n'], parent: ['b'], place: { after: null }, rename: 'm' },
      ],
    });
    expect(r.changes).toEqual({ made: [['b', 'm']], params: [], children: [['b']], removed: [], moved: [] });
  });

  it('a set then a remove is only the remove', () => {
    const r = applied(root, {
      ops: [
        { op: 'set', at: ['a', 'x'], key: 'label', value: 'X' },
        { op: 'remove', at: ['a', 'x'] },
      ],
    });
    expect(r.changes).toEqual({ made: [], params: [], children: [['a']], removed: [['a', 'x']], moved: [] });
  });

  it('a swap by three renames is two moves and one changed child list', () => {
    const tree = doc(box('a', undefined, leaf('x', { label: 'X' }), leaf('y')));
    const r = applied(tree, {
      ops: [
        { op: 'move', at: ['a', 'x'], parent: ['a'], place: { after: 'y' }, rename: 't' },
        { op: 'move', at: ['a', 'y'], parent: ['a'], place: { after: null }, rename: 'x' },
        { op: 'move', at: ['a', 't'], parent: ['a'], place: { last: true }, rename: 'y' },
      ],
    });
    expect(r.changes.made).toEqual([]);
    expect(r.changes.children).toEqual([['a']]);
    expect(r.changes.moved).toEqual([
      { from: ['a', 'y'], to: ['a', 'x'] },
      { from: ['a', 'x'], to: ['a', 'y'] },
    ]);
    expect(r.root.children![0].children).toEqual([leaf('x'), leaf('y', { label: 'X' })]);
  });

  it('a set back to the old value is no change at all', () => {
    const r = applied(root, {
      ops: [
        { op: 'set', at: ['a'], key: 'label', value: 'B' },
        { op: 'set', at: ['a'], key: 'label', value: 'A' },
      ],
    });
    expect(r.changed).toBe(false);
    expect(r.root).toBe(root);
    expect(r.inverse).toBeUndefined();
  });

  it('a reorder within a parent is in children only; a ref rewrite is in params', () => {
    const tree = doc(box('a', undefined, leaf('x')), box('b', { link: '/a/x' }));
    const r = applied(
      tree,
      {
        ops: [
          { op: 'move', at: ['b'], parent: [], place: { after: null } },
          { op: 'move', at: ['a', 'x'], parent: ['a'], place: { last: true }, rename: 'w' },
        ],
      },
      lang,
    );
    expect(r.changes.children).toEqual([[], ['a']]);
    expect(r.changes.params).toEqual([{ at: ['b'], keys: ['link'], where: 'params' }]);
    expect(r.changes.moved).toEqual([{ from: ['a', 'x'], to: ['a', 'w'] }]);
  });

  it('a node moved into a subtree that is then removed is removed at its old address', () => {
    const r = applied(root, {
      ops: [
        { op: 'move', at: ['b', 'z'], parent: ['a'], place: { last: true } },
        { op: 'remove', at: ['a'] },
      ],
    });
    expect(r.changes.removed).toEqual([['a'], ['b', 'z']]);
    expect(r.changes.children).toEqual([[], ['b']]);
  });
});

describe('repr/patch: the inverse, with what it needs only', () => {
  const back = (root: Node, patch: Patch, options: ApplyOptions) => {
    const r = applied(root, patch, options);
    const undo = applyPatch(r.root, r.inverse!, { ...lang, ...r.inverseNeeds });
    expect(undo.ok && hashNode(undo.root)).toBe(hashNode(root));
    return r;
  };

  it('a stored set made before the node it targets is typed again, so a later move takes its ref along', () => {
    // the inverse makes `b` (its step setting /a/a's link to /a) before it makes `/a`, then moves `/a` to `/b/d`
    const root: Node = {
      kind: 'doc',
      children: [
        {
          kind: 'box',
          id: 'b',
          params: {
            link: '/b/b',
            step: patchToNode({ ops: [{ op: 'set', at: ['b', 'd', 'a'], key: 'link', value: '/b/d' }] }),
          },
          children: [leaf('b', { link: '/b/d' }), { kind: 'box', id: 'd', children: [leaf('a')] }],
        },
      ],
    };
    const patch: Patch = {
      ops: [
        { op: 'set', at: ['b', 'b'], key: 'link', value: '/b/b' },
        { op: 'move', at: ['b', 'd'], parent: [], place: { last: true }, rename: 'a' },
        { op: 'remove', at: ['a'] },
        { op: 'remove', at: ['b'] },
      ],
    };
    back(root, patch, { ...lang, allowDangling: true, allowRebind: true });
  });

  it("a stored set of the old root's that a later op types is typed again, so a later move takes its ref along", () => {
    // the step sets /x/e's link to /s: typed while /x/e is there. The patch moves /s to /t (the ref follows), then
    // removes /x/e; the inverse makes /x/e again, and its move back must take the step's ref along
    const step = patchToNode({ ops: [{ op: 'set', at: ['x', 'e'], key: 'link', value: '/s' }] });
    const root = doc(box('x', undefined, leaf('e')), box('s', { step }));
    const patch: Patch = {
      ops: [
        { op: 'move', at: ['s'], parent: [], place: { last: true }, rename: 't' },
        { op: 'remove', at: ['x', 'e'] },
      ],
    };
    back(root, patch, { ...lang, allowDangling: true, allowRebind: true });
  });

  it('undoes a remove whose content held a ref already dangling', () => {
    back(doc(box('a', { link: '/zz' })), { ops: [{ op: 'remove', at: ['a'] }] }, lang);
    back(
      doc(box('a', { link: '/zz' }, leaf('x'))),
      { ops: [{ op: 'replace', at: ['a'], node: box('a', undefined) }] },
      lang,
    );
  });

  it('a replace that keeps a ref already dangling is not refused for it', () => {
    const r = applyPatch(
      doc(box('a', { link: '/zz' })),
      { ops: [{ op: 'replace', at: ['a'], node: box('a', { link: '/zz', label: 'x' }) }] },
      lang,
    );
    expect(r.ok).toBe(true);
  });

  it('a ref whose node was removed and whose address now names another node is rebound', () => {
    const tree = doc(
      box('a', undefined, leaf('y')),
      box('b', undefined, leaf('y')),
      box('r', { link: '/b/y', links: ['/a/y'] }),
    );
    const patch: Patch = {
      ops: [
        { op: 'remove', at: ['b'] },
        { op: 'move', at: ['a'], parent: [], place: { after: null }, rename: 'b' },
      ],
    };
    const refused = applyPatch(tree, patch, lang);
    expect(refused.ok).toBe(false);
    expect(refused.issues.map((i) => [i.code, i.path])).toEqual([
      ['repr.patch.rebound', ['children', 1, 'params', 'link']],
    ]);
    const r = back(tree, patch, { ...lang, allowRebind: true });
    expect(r.root.children![1].params).toEqual({ link: '/b/y', links: ['/b/y'] });
  });

  it("the inverse's moves never drag a ref the forward patch did not rewrite", () => {
    // /c/x dangles before and after; the move back of /m (now /c) would take it along
    const tree = doc(box('a', { link: '/c/x' }), box('m', undefined, leaf('y')));
    back(
      tree,
      { ops: [{ op: 'move', at: ['m'], parent: [], place: { last: true }, rename: 'c' }] },
      { ...lang, allowRebind: true },
    );
    // a ref an op set under the destination, and content removed after the move
    const two = doc(box('a', { link: '/q' }), box('m', undefined, leaf('y')), box('n', { link: '/c/z' }));
    back(
      two,
      {
        ops: [
          { op: 'move', at: ['m'], parent: [], place: { last: true }, rename: 'c' },
          { op: 'set', at: ['a'], key: 'link', value: '/c/y' },
          { op: 'remove', at: ['n'] },
        ],
      },
      { ...lang, allowRebind: true },
    );
  });
});

describe('repr/patch: structural sharing', () => {
  const root = doc(box('a', { label: 'A' }, leaf('x')), box('b', undefined, leaf('z')));

  it('an untouched subtree is the same object after a patch', () => {
    const r = applied(root, { ops: [{ op: 'set', at: ['a', 'x'], key: 'label', value: 'X' }] });
    expect(r.root.children![1]).toBe(root.children![1]);
    expect(r.root.children![0]).not.toBe(root.children![0]);
  });

  it('an equal set or replace keeps the node, and the root', () => {
    expect(applied(root, { ops: [{ op: 'set', at: ['a'], key: 'label', value: 'A' }] }).root).toBe(root);
    const copy = structuredClone(root.children![1]);
    expect(applied(root, { ops: [{ op: 'replace', at: ['b'], node: copy }] }).root).toBe(root);
  });

  it('a refusal leaves the root as it was, by hash and by object', () => {
    const before = hashNode(structuredClone(root));
    const r = applyPatch(root, {
      ops: [
        { op: 'set', at: ['a'], key: 'label', value: 'Q' },
        { op: 'insert', parent: ['a'], place: { last: true }, node: leaf('x') },
      ],
    });
    expect(r.ok).toBe(false);
    expect(hashNode(structuredClone(root))).toBe(before);
    expect(root.children![0].params!.label).toBe('A');
  });
});

describe('repr/patch: issues', () => {
  it('an op issue has the patch tree path', () => {
    const r = applyPatch(doc(), {
      ops: [
        { op: 'set', at: [], key: 'title', value: 'x' },
        { op: 'remove', at: ['q'] },
      ],
    });
    expect(r.ok).toBe(false);
    expect(r.issues[0]).toMatchObject({ code: 'repr.patch.target', path: ['children', 1, 'params', 'at'] });
  });

  it("names the patch tree's laid-out params: after, not place; kind, not op", () => {
    const place = applyPatch(doc(box('a', undefined)), {
      ops: [{ op: 'insert', parent: ['a'], place: { after: 'zz' }, node: leaf('n') }],
    });
    expect(place.issues.map((i) => i.path)).toEqual([['children', 0, 'params', 'after']]);
    const kind = applyPatch(doc(), { ops: [{ op: 'swap' } as never] });
    expect(kind.issues.map((i) => i.path)).toEqual([['children', 0, 'kind']]);
  });

  it('result is checked before base', () => {
    const root = doc();
    expect(applyPatch(root, { ops: [], base: '0'.repeat(32), result: hashNode(root) }).ok).toBe(true);
  });

  it('without a language, refs are not rewritten, and the result says so', () => {
    const tree = doc(box('a', undefined, leaf('x')), box('b', { link: '/a/x' }));
    const r = applied(tree, { ops: [{ op: 'move', at: ['a'], parent: [], place: { last: true }, rename: 'c' }] });
    expect(r.refs).toBe('not-rewritten');
    expect(r.root.children![0].params!.link).toBe('/a/x');
  });

  it("a dangling ref's issue has the spec's path and address; allowed, a warning, and its inverse needs a rebinding", () => {
    const tree = doc(box('a', undefined, leaf('x')), box('b', { link: '/a/x' }));
    const refused = applyPatch(tree, { ops: [{ op: 'remove', at: ['a'] }] }, lang);
    expect(refused.ok).toBe(false);
    expect(refused.issues[0]).toMatchObject({
      code: 'repr.patch.dangling',
      severity: 'error',
      path: ['children', 0, 'params', 'link'], // in the new root
      address: '/b',
    });
    const allowed = applied(tree, { ops: [{ op: 'remove', at: ['a'] }] }, { ...lang, allowDangling: true });
    expect(allowed.issues.map((i) => i.severity)).toEqual(['warning']);
    expect(allowed.inverseNeeds).toEqual({ allowRebind: true });
    expect(applyPatch(allowed.root, allowed.inverse!, lang).ok).toBe(false);
  });

  it("a language's validators run last, on what changed; their errors refuse", () => {
    const seen: string[] = [];
    const v = (before: Node, after: Node): Issue[] => {
      seen.push(`${hashNode(before) !== hashNode(after)}`);
      return [{ code: 'test.no', severity: 'error', message: 'no' }];
    };
    const r = applyPatch(
      doc(),
      { ops: [{ op: 'set', at: [], key: 'title', value: 'x' }] },
      { ...lang, validators: [v] },
    );
    expect(r.ok).toBe(false);
    expect(r.issues).toEqual([{ code: 'test.no', severity: 'error', message: 'no' }]);
    expect(seen).toEqual(['true']);
  });

  it('validation reports only what the patch introduced', () => {
    const tree = doc(box('a', { size: -1 }), box('b', { size: -2 }));
    const r = applyPatch(tree, { ops: [{ op: 'set', at: ['b'], key: 'size', value: -3 }] }, lang);
    expect(r.ok).toBe(false);
    expect(r.issues.map((i) => [i.code, i.address])).toEqual([['repr.schema.param.range', '/b']]);
    const unset = applyPatch(
      doc({ kind: 'crate', id: 'c' }),
      { ops: [{ op: 'set', at: ['c'], key: 'label', value: 'x' }] },
      lang,
    );
    expect(unset.ok).toBe(false);
    expect(unset.issues.map((i) => i.code)).toEqual(['repr.schema.param.unknown']);
  });
});

describe('repr/patch: the patch as a tree', () => {
  const patch: Patch = {
    base: 'a'.repeat(32),
    ops: [
      { op: 'set', at: ['views', 'main', 't1'], key: 'height', value: 40 },
      { op: 'unset', at: [], key: 'title', where: 'custom' },
      { op: 'insert', parent: ['a'], place: { after: null }, node: leaf('n', { label: 'n' }) },
      { op: 'insert', parent: ['a'], place: { last: true }, node: leaf('m') },
      { op: 'remove', at: ['a', 'x'] },
      { op: 'move', at: ['a', 'y'], parent: ['b'], place: { after: 'z' }, rename: 'w' },
      { op: 'replace', at: ['b'], node: box('b', undefined) },
    ],
  };

  it('round-trips through patchToNode, the codecs and patchFromNode', () => {
    const node = patchToNode(patch);
    const text = printSexpr(node);
    expect(text).toContain('(set :at /views/main/t1 :key height :value 40)');
    expect(patchFromNode(parseSexpr(text).node!).patch).toEqual(patch);
    expect(patchFromNode(parseJson(printJson(node)).node!).patch).toEqual(patch);
  });

  it('is valid in its patch language; one language object per language, its id with .patch', () => {
    const pl = patchLanguage(PC);
    expect(pl).toBe(patchLanguage(PC));
    expect(pl.id).toBe('pc.patch');
    const ops = patch.ops.filter((o) => o.op !== 'replace' && o.op !== 'insert');
    expect(validate(pl, patchToNode({ ...patch, ops }), { languages: (id) => (id === 'pc' ? PC : undefined) })).toEqual(
      [],
    );
    expect(validate(pl, patchToNode(patch), { languages: (id) => (id === 'pc' ? PC : undefined) })).toEqual([]);
    expect(validate(AnyPatchLanguage, patchToNode(patch))).toEqual([]);
  });

  it('refuses a malformed op tree with repr.patch.op at its path', () => {
    const r = patchFromNode({
      kind: 'patch',
      children: [
        { kind: 'set', params: { at: 'nope', key: 'k', value: 1 } },
        { kind: 'insert', params: { parent: '/a', first: true, last: true, node: leaf('q') } },
        { kind: 'swap' },
      ],
    });
    expect(r.patch).toBeUndefined();
    expect(r.issues.map((i) => [i.code, i.path])).toEqual([
      ['repr.patch.op', ['children', 0, 'params', 'at']],
      ['repr.patch.op', ['children', 1, 'params']],
      ['repr.patch.op', ['children', 2, 'kind']],
    ]);
  });
});

describe('repr/patch: diff', () => {
  it('emits removes, then placements in b order, then descents', () => {
    const a = doc(box('p', undefined, leaf('a'), leaf('b'), leaf('c'), leaf('d', { label: 'x' })));
    const b = doc(box('p', undefined, leaf('d', { label: 'y' }), leaf('a'), leaf('n'), leaf('c')));
    expect(diff(a, b).ops).toEqual([
      { op: 'remove', at: ['p', 'b'] },
      { op: 'move', at: ['p', 'd'], parent: ['p'], place: { after: null } },
      { op: 'insert', parent: ['p'], place: { after: 'a' }, node: leaf('n') },
      { op: 'set', at: ['p', 'd'], key: 'label', value: 'y' },
    ]);
  });

  it('falls back to a replace: a kind change, children without ids, a root of another kind', () => {
    expect(diff(doc(leaf('a')), doc(box('a', undefined))).ops).toEqual([
      { op: 'replace', at: ['a'], node: box('a', undefined) },
    ]);
    const bag = (n: number): Node =>
      doc({ kind: 'bag', id: 'g', children: Array.from({ length: n }, () => ({ kind: 'item' })) });
    expect(diff(bag(1), bag(2)).ops).toEqual([{ op: 'replace', at: ['g'], node: bag(2).children![0] }]);
    expect(diff(doc(), leaf('r')).ops).toEqual([{ op: 'replace', at: [], node: leaf('r') }]);
  });

  it('keeps the longest run in order', () => {
    expect(longestIncreasing([3, 0, 1, 4, 2, 5])).toEqual([1, 2, 4, 5]);
    expect(longestIncreasing([])).toEqual([]);
  });
});

/** A seeded generator (mulberry32). */
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

const Ids = ['a', 'b', 'c', 'd', 'e', 'f'];

/** A random tree in the conformance language: boxes of boxes and leaves, labels, sizes, columns and refs. */
function randomTree(r: () => number, refs: string[]): Node {
  const pick = <T>(xs: readonly T[]) => xs[Math.floor(r() * xs.length)];
  const params = (kind: 'box' | 'leaf') => {
    const p: Record<string, Value> = {};
    if (r() < 0.5) p.label = pick(['x', 'y', 'z']);
    if (r() < 0.3) p.link = pick(refs);
    if (kind === 'box' && r() < 0.3) p.size = pick([1, 2, 3]);
    if (kind === 'box' && r() < 0.2) p.cols = new Float32Array([r() < 0.5 ? 1 : 2, 3]);
    if (kind === 'box' && r() < 0.2) p.links = [pick(refs), pick(refs)];
    return Object.keys(p).length ? p : undefined;
  };
  const kids = (depth: number): Node[] => {
    const ids = Ids.filter(() => r() < 0.5);
    return ids.map((id) => {
      if (depth < 3 && r() < 0.5) return { ...box(id, params('box')), ...nonEmpty(kids(depth + 1)) };
      return leaf(id, params('leaf'));
    });
  };
  return { kind: 'doc', ...(r() < 0.5 ? { params: { title: pick(['t', 'u']) } } : {}), children: kids(0).map(asBox) };
}
const nonEmpty = (children: Node[]) => (children.length ? { children } : {});
/** The root's children are boxes (a leaf sits under a box only). */
const asBox = (n: Node): Node => (n.kind === 'box' ? n : box(n.id!, n.params));

/** A tree with nodes moved, removed, renamed, changed in kind, with params set and unset. */
function mutate(r: () => number, root: Node, refs: string[]): Node {
  const other = randomTree(r, refs);
  const mix = (a: Node, b: Node, depth: number): Node => {
    const roll = r();
    if (roll < 0.15 && depth > 0) return b;
    const kids = [...(a.children ?? [])];
    if (r() < 0.3) kids.reverse();
    if (r() < 0.3 && kids.length) kids.splice(Math.floor(r() * kids.length), 1);
    if (r() < 0.3) for (const c of b.children ?? []) if (!kids.some((k) => k.id === c.id)) kids.push(c);
    const params = r() < 0.4 ? b.params : a.params;
    const out = kids.map((k) => {
      const twin = b.children?.find((c) => c.id === k.id && c.kind === k.kind);
      return twin && r() < 0.6 ? mix(k, twin, depth + 1) : k;
    });
    return {
      kind: a.kind,
      ...(a.id !== undefined ? { id: a.id } : {}),
      ...(params ? { params } : {}),
      ...nonEmpty(out),
    };
  };
  return mix(root, other, 0);
}

describe('repr/patch: the laws of diff (seeded)', () => {
  const refs = ['/a', '/b/c', '/a/b', '/c/d/e', '/f'];
  const pairs = (n: number, seed: number) => {
    const r = rng(seed);
    return Array.from({ length: n }, () => {
      const a = randomTree(r, refs);
      return [a, mutate(r, a, refs)] as const;
    });
  };

  it('without a language: apply(a, diff(a, b)) is b; diff(a, a) is empty; the inverse gives a back', () => {
    for (const [a, b] of pairs(200, 7)) {
      expect(diff(a, a).ops).toEqual([]);
      const p = diff(a, b);
      const r = applied(a, p, { maxOps: Number.POSITIVE_INFINITY });
      expect(hashNode(r.root)).toBe(hashNode(b));
      if (!r.inverse) continue;
      const back = applied(r.root, r.inverse, { ...r.inverseNeeds, maxOps: Number.POSITIVE_INFINITY });
      expect(hashNode(back.root)).toBe(hashNode(a));
    }
  });

  it('with a language, on valid trees, with dangling and rebinding allowed', () => {
    const options: ApplyOptions = { ...lang, allowDangling: true, allowRebind: true, maxOps: Number.POSITIVE_INFINITY };
    for (const [a, b] of pairs(200, 11)) {
      expect(validate(PC, a).filter((i) => i.severity === 'error')).toEqual([]);
      expect(validate(PC, b).filter((i) => i.severity === 'error')).toEqual([]);
      const r = applied(a, diff(a, b), options);
      expect(hashNode(r.root)).toBe(hashNode(b));
      if (!r.inverse) continue;
      const back = applied(r.root, r.inverse, { ...options, ...r.inverseNeeds });
      expect(hashNode(back.root)).toBe(hashNode(a));
    }
  });

  it('a patch and its inverse, with a move that rewrites refs', () => {
    const r = rng(3);
    for (let i = 0; i < 100; i++) {
      const a = randomTree(r, refs);
      const kids = a.children ?? [];
      if (kids.length < 2) continue;
      const from = kids[0].id!;
      const to = kids[1].id!;
      const rename = Ids.find((id) => !(kids[1].children ?? []).some((c) => c.id === id))!;
      const p: Patch = { ops: [{ op: 'move', at: [from], parent: [to], place: { last: true }, rename }] };
      const res = applyPatch(a, p, { ...lang, allowDangling: true, allowRebind: true });
      if (!res.ok) continue; // a leaf under the root's box: the parent rule may refuse
      const back = applied(res.root, res.inverse!, { ...lang, ...res.inverseNeeds });
      expect(hashNode(back.root)).toBe(hashNode(a));
    }
  });
});
