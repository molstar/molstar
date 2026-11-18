/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Law 3 over seeded multi-op patches with refs and a language: whatever a patch does (sets of refs, inserts, removes,
// moves and renames that rewrite refs, replaces), its inverse applied with what it needs only gives the root back.
// Short in the quick lane, long in the heavy one.

import { describe, expect, it } from 'vitest';
import { type ApplyOptions, applyPatch, type Op, type Patch, patchLanguage } from '../../patch.js';
import { type Address, hashNode, type Node, resolveAddress, type Value } from '../../tree/index.js';
import { PatchConformanceLanguage as PC } from '../conformance.js';
import { freshRefIndex, peekIndex, type RefIndex, refIndexOf, type Site } from '../refs.js';

const languages = (id: string) => (id === 'pc.patch' ? patchLanguage(PC) : undefined);
const lang: ApplyOptions = { language: PC, validate: { languages } };

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

const Ids = ['a', 'b', 'c', 'd'];

/** Every box's address in a tree (boxes hold boxes and leaves; the root's children are boxes). */
function boxes(root: Node): Address[] {
  const out: Address[] = [];
  const go = (n: Node, at: Address) => {
    for (const c of n.children ?? []) {
      if (c.kind !== 'box' || c.id === undefined) continue;
      out.push([...at, c.id]);
      go(c, [...at, c.id]);
    }
  };
  go(root, []);
  return out;
}

function nodes(root: Node): Address[] {
  const out: Address[] = [];
  const go = (n: Node, at: Address) => {
    for (const c of n.children ?? []) {
      if (c.id === undefined) continue;
      out.push([...at, c.id]);
      go(c, [...at, c.id]);
    }
  };
  go(root, []);
  return out;
}

const text = (a: Address) => `/${a.join('/')}`;

/** A site without its node's identity (compared apart), steps by what they point at. */
function view(s: Site) {
  const { node: _node, step, ...rest } = s;
  return { ...rest, ...(step ? { step: { ...step, slot: hashNode(step.slot) } } : {}) };
}

/** A derived index against one built from scratch on the same root: every site, `dangling`, and `under` by prefix. */
function expectSameIndex(derived: RefIndex, root: Node) {
  const fresh = freshRefIndex(root, PC, lang.validate);
  const a = derived.under([]) as readonly Site[];
  const b = fresh.under([]) as readonly Site[];
  expect(a.map(view)).toEqual(b.map(view));
  a.forEach((s, i) => expect(s.node).toBe(b[i].node));
  expect((derived.dangling() as readonly Site[]).map(view)).toEqual((fresh.dangling() as readonly Site[]).map(view));
  for (const at of nodes(root).slice(0, 12))
    expect((derived.under(at) as readonly Site[]).map(view)).toEqual((fresh.under(at) as readonly Site[]).map(view));
}

/**
 * `derived`: the old root has its index, so the apply derives the new root's, and the inverse is applied over it and
 * over a copy of the root with none (the two consumers: link integrity and the index itself).
 */
function fuzz(cases: number, seed: number, derived = false) {
  const r = rng(seed);
  const pick = <T>(xs: readonly T[]): T => xs[Math.floor(r() * xs.length)];
  const refTo = (root: Node) => {
    const all = nodes(root);
    // mostly a node there is, sometimes one there is not, sometimes below one
    const roll = r();
    if (roll < 0.6 && all.length) return text(pick(all));
    if (roll < 0.8 && all.length) return `${text(pick(all))}/${pick(Ids)}`;
    return `/${pick(Ids)}/${pick(Ids)}`;
  };
  const params = (root: Node, kind: 'box' | 'leaf') => {
    const p: Record<string, Value> = {};
    if (r() < 0.5) p.link = refTo(root);
    if (kind === 'box' && r() < 0.3) p.links = [refTo(root), refTo(root)];
    if (r() < 0.3) p.label = pick(['x', 'y']);
    return Object.keys(p).length ? p : undefined;
  };
  const tree = (): Node => {
    const root: { kind: string; children: Node[] } = { kind: 'doc', children: [] };
    const fill = (n: { children?: Node[] }, depth: number) => {
      const kids: Node[] = [];
      for (const id of Ids) {
        if (r() < 0.5) continue;
        const isBox = depth === 0 || r() < 0.5;
        const c: { kind: string; id: string; params?: Record<string, Value>; children?: Node[] } = {
          kind: isBox ? 'box' : 'leaf',
          id,
        };
        if (isBox && depth < 2) fill(c, depth + 1);
        kids.push(c);
      }
      if (kids.length) n.children = kids;
    };
    fill(root, 0);
    // refs, once the nodes they may name are there
    const withRefs = (n: Node): Node => ({
      ...n,
      ...(n.kind === 'doc' ? {} : { params: params(root, n.kind as 'box' | 'leaf') }),
      ...(n.children ? { children: n.children.map(withRefs) } : {}),
    });
    return withRefs(root);
  };
  const op = (root: Node): Op | undefined => {
    const all = nodes(root);
    const bx = boxes(root);
    const roll = r();
    if (roll < 0.25 && all.length) return { op: 'set', at: pick(all), key: 'link', value: refTo(root) };
    if (roll < 0.35 && all.length) return { op: 'unset', at: pick(all), key: 'link' };
    if (roll < 0.5) {
      const parent = r() < 0.3 || !bx.length ? [] : pick(bx);
      const node: Node = parent.length
        ? { kind: 'leaf', id: pick(Ids), params: { link: refTo(root) } }
        : { kind: 'box', id: pick(Ids), params: { link: refTo(root) } };
      return { op: 'insert', parent, place: r() < 0.5 ? { last: true } : { after: null }, node };
    }
    if (roll < 0.65 && all.length) return { op: 'remove', at: pick(all) };
    if (roll < 0.9 && all.length) {
      const at = pick(all);
      const n = resolveAddress(root, at)!;
      const parents = n.kind === 'box' ? [[] as Address, ...bx] : bx;
      return {
        op: 'move',
        at,
        parent: pick(parents),
        place: { last: true },
        ...(r() < 0.6 ? { rename: pick(Ids) } : {}),
      };
    }
    if (all.length) {
      const at = pick(all);
      const n = resolveAddress(root, at)!;
      return { op: 'replace', at, node: { ...n, params: { link: refTo(root) } } };
    }
    return undefined;
  };

  let applied = 0;
  for (let c = 0; c < cases; c++) {
    const root = tree();
    // ops made one by one against the tree as the earlier ones leave it
    const ops: Op[] = [];
    let state = root;
    const count = 1 + Math.floor(r() * 6);
    for (let i = 0; i < count; i++) {
      const o = op(state);
      if (!o) break;
      const step = applyPatch(state, { ops: [o] });
      if (!step.ok) continue;
      ops.push(o);
      state = step.root;
    }
    const patch: Patch = { ops };
    const allow = { allowDangling: r() < 0.5 || undefined, allowRebind: r() < 0.5 || undefined };
    const before = hashNode(structuredClone(root));
    if (derived) refIndexOf(root, PC, lang.validate);
    const res = applyPatch(root, patch, { ...lang, ...allow });
    expect(hashNode(structuredClone(root))).toBe(before);
    if (!res.ok || !res.inverse) continue;
    applied++;
    const back = applyPatch(res.root, res.inverse, { ...lang, ...res.inverseNeeds });
    if (derived) {
      const index = peekIndex(res.root, PC, lang.validate);
      expect(index, `case ${c}: the new root has its index`).toBeDefined();
      expectSameIndex(index!, res.root);
      if (back.ok && back.changed) {
        const next = peekIndex(back.root, PC, lang.validate);
        expect(next, `case ${c}: the root after the inverse has its index`).toBeDefined();
        expectSameIndex(next!, back.root);
      }
      // link integrity over the derived index and over none
      const alone = applyPatch(structuredClone(res.root), res.inverse, { ...lang, ...res.inverseNeeds });
      expect(alone.ok, `case ${c}`).toBe(back.ok);
      if (alone.ok && back.ok) expect(hashNode(alone.root)).toBe(hashNode(back.root));
      else if (!alone.ok && !back.ok) expect(alone.issues.map((i) => i.code)).toEqual(back.issues.map((i) => i.code));
    }
    if (!back.ok || hashNode(back.root) !== hashNode(root)) {
      throw new Error(
        `case ${c}: ${back.ok ? 'not the root' : back.issues.map((i) => i.code).join()}\n${JSON.stringify({ root, patch, allow })}`,
      );
    }
  }
  return applied;
}

describe('repr/patch: law 3 over seeded multi-op patches with refs and a language', () => {
  it('the inverse with what it needs gives the root back (short)', () => {
    expect(fuzz(150, 1)).toBeGreaterThan(50);
  });

  it('the inverse with what it needs gives the root back (long)', { tags: ['heavy'] }, () => {
    expect(fuzz(10000, 2)).toBeGreaterThan(4000);
  });

  it('a derived ref index equals one built on the same root (short)', () => {
    expect(fuzz(150, 3, true)).toBeGreaterThan(50);
  });

  it('a derived ref index equals one built on the same root (long)', { tags: ['heavy'] }, () => {
    expect(fuzz(3000, 4, true)).toBeGreaterThan(1000);
  });
});
