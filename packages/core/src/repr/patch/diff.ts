/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The patch between two revisions (_spec/README.md "Patches", Diff). Merkle descent: equal hashes are skipped.
// Children are matched by id within a parent; moves within a parent are the fewest that keep the longest run in order
// (a longest increasing subsequence). Emission, so every place is valid when applied in order: removes, then `b`'s
// children left to right, each not in the kept run placed after its predecessor in `b`, then the descents into the
// kept children. Changes of parent and renames come out as a remove and an insert.

import { type Address, hashNode, type Node, own, type Params, valuesEqual } from '../tree/index.js';
import type { Op, Patch, Where } from './patch.js';

/** A patch taking `a` to `b`: `base` is `a`'s revision, `result` `b`'s. */
export function diff(a: Node, b: Node): Patch {
  const ops: Op[] = [];
  const base = hashNode(a);
  const result = hashNode(b);
  if (base !== result) {
    if (a.kind !== b.kind || a.id !== b.id || !argsEqual(a, b)) ops.push({ op: 'replace', at: [], node: b });
    else descend(a, b, [], ops);
  }
  return { ops, base, result };
}

const argsEqual = (a: Node, b: Node) => valuesEqual(a.args ?? [], b.args ?? []);

/** Whether every child has an id and no two share one: else the parent changes by `replace`. */
function addressable(children: readonly Node[] | undefined): boolean {
  if (!children) return true;
  const seen = new Set<string>();
  for (const c of children) {
    if (c.id === undefined || seen.has(c.id)) return false;
    seen.add(c.id);
  }
  return true;
}

function childrenEqual(a: readonly Node[] | undefined, b: readonly Node[] | undefined): boolean {
  const x = a ?? [];
  const y = b ?? [];
  if (x.length !== y.length) return false;
  for (let i = 0; i < x.length; i++) if (x[i] !== y[i] && hashNode(x[i]) !== hashNode(y[i])) return false;
  return true;
}

function keys(a: Params | undefined, b: Params | undefined, at: Address, where: Where, ops: Op[]) {
  const all = new Set([...Object.keys(a ?? {}), ...Object.keys(b ?? {})]);
  for (const key of [...all].sort()) {
    const va = a ? own(a, key) : undefined;
    const vb = b ? own(b, key) : undefined;
    if (vb === undefined) {
      if (va !== undefined) ops.push(where === 'params' ? { op: 'unset', at, key } : { op: 'unset', at, key, where });
    } else if (!valuesEqual(va, vb))
      ops.push(where === 'params' ? { op: 'set', at, key, value: vb } : { op: 'set', at, key, value: vb, where });
  }
}

function descend(a: Node, b: Node, at: Address, ops: Op[]) {
  const sameChildren = childrenEqual(a.children, b.children);
  if (!sameChildren && (!addressable(a.children) || !addressable(b.children))) {
    ops.push({ op: 'replace', at, node: b });
    return;
  }
  keys(a.params, b.params, at, 'params', ops);
  keys(a.custom, b.custom, at, 'custom', ops);
  if (sameChildren) return;
  const ac = a.children ?? [];
  const bc = b.children ?? [];
  const inA = new Map<string, number>();
  ac.forEach((c, i) => inA.set(c.id!, i));
  const inB = new Set(bc.map((c) => c.id!));
  for (const c of ac) if (!inB.has(c.id!)) ops.push({ op: 'remove', at: [...at, c.id!] });
  // the kept children in b's order, by their positions in a: the longest increasing run stays where it is
  const kept: number[] = [];
  for (let j = 0; j < bc.length; j++) if (inA.has(bc[j].id!)) kept.push(j);
  const stay = new Set(longestIncreasing(kept.map((j) => inA.get(bc[j].id!)!)).map((k) => kept[k]));
  for (let j = 0; j < bc.length; j++) {
    if (stay.has(j)) continue;
    const place = { after: j > 0 ? bc[j - 1].id! : null };
    if (inA.has(bc[j].id!)) ops.push({ op: 'move', at: [...at, bc[j].id!], parent: at, place });
    else ops.push({ op: 'insert', parent: at, place, node: bc[j] });
  }
  for (const j of kept) {
    const x = ac[inA.get(bc[j].id!)!];
    const y = bc[j];
    if (x === y || hashNode(x) === hashNode(y)) continue;
    const there = [...at, y.id!];
    if (x.kind !== y.kind || !argsEqual(x, y)) ops.push({ op: 'replace', at: there, node: y });
    else descend(x, y, there, ops);
  }
}

/** The indices of a longest strictly increasing subsequence of `xs` (patience sorting, O(n log n)). */
export function longestIncreasing(xs: readonly number[]): number[] {
  const tails: number[] = []; // index into xs of the smallest tail of each length
  const prev = new Int32Array(xs.length).fill(-1);
  for (let i = 0; i < xs.length; i++) {
    let lo = 0;
    let hi = tails.length;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (xs[tails[mid]] < xs[i]) lo = mid + 1;
      else hi = mid;
    }
    if (lo > 0) prev[i] = tails[lo - 1];
    tails[lo] = i;
  }
  const out: number[] = new Array(tails.length);
  for (let k = tails.length - 1, i = tails[tails.length - 1]; k >= 0; k--, i = prev[i]) out[k] = i;
  return out;
}
