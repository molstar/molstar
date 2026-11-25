/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The combinators (_spec/README.md "Queries", Language): `merge`, `intersect`
// and `except` over results of one type (renamed from `or`, `and`, `not` by the author, 2026-11-20, "Rename", MolQL's
// names; `and`, `or`, `not` are predicate logic now). `except` takes its first operand less the others; with one
// operand it is the complement within what the prefix names (its systems' positions, its annotations' items, its
// nodes' children). A location less a location is exact for rects too: rects subtract, a rect has no complement.

import { type Named, type Node, printAddress, type Selected, type Value } from '@molstar/core/repr/tree';
import { evaluate } from './evaluate.js';
import {
  annotationsOf,
  checkLimit,
  checkParams,
  type Ctx,
  fail,
  prefixNodes,
  type Resolve,
  required,
  stage,
  systemsOf,
} from './form-kit.js';
import { EmptyKeys, intersectKeys, keySet, subtractKeys, unionKeys } from './keys.js';
import { intersectLocations, location, subtractLocations, unionLocations } from './location.js';
import type { QueryResult } from './result.js';

const nodeKey = (s: Selected) => printAddress(s.address);

export function combineNodes(
  a: readonly Selected[],
  b: readonly Selected[],
  keep: (inA: boolean, inB: boolean) => boolean,
): Selected[] {
  const inA = new Set(a.map(nodeKey));
  const inB = new Set(b.map(nodeKey));
  const seen = new Set<string>();
  const out: Selected[] = [];
  for (const s of [...a, ...b]) {
    const k = nodeKey(s);
    if (seen.has(k) || !keep(inA.has(k), inB.has(k))) continue;
    seen.add(k);
    out.push(s);
  }
  return out;
}

type Op = 'merge' | 'intersect' | 'except';

/** Two results of one type combined; `seq.query.form.type` when they differ. */
function pair(kind: Op, a: QueryResult, b: QueryResult): QueryResult {
  const keep =
    kind === 'merge'
      ? (x: boolean, y: boolean) => x || y
      : kind === 'intersect'
        ? (x: boolean, y: boolean) => x && y
        : (x: boolean, y: boolean) => x && !y;
  if (a.type === 'nodes' && b.type === 'nodes') return { type: 'nodes', nodes: combineNodes(a.nodes, b.nodes, keep) };
  if (a.type === 'keys' && b.type === 'keys')
    return {
      type: 'keys',
      keys: (kind === 'merge' ? unionKeys : kind === 'intersect' ? intersectKeys : subtractKeys)(a.keys, b.keys),
    };
  if (a.type === 'location' && b.type === 'location')
    return {
      type: 'location',
      location: (kind === 'merge' ? unionLocations : kind === 'intersect' ? intersectLocations : subtractLocations)(
        a.location,
        b.location,
      ),
    };
  return fail('seq.query.form.type', `A form of kind ${kind} mixes result types (${a.type}, ${b.type})`);
}

function limited(ctx: Ctx, r: QueryResult): QueryResult {
  if (r.type === 'nodes') checkLimit(ctx, r.nodes.length);
  if (r.type === 'keys') checkLimit(ctx, r.keys.keys.length);
  return r;
}

/** The operands of `:of`: a non-empty list (`except` also takes one value alone). */
function operandValues(form: Node, lone: boolean): readonly Value[] {
  checkParams(form, ['of']);
  const of = required(form, 'of');
  if (Array.isArray(of) && of.length) return of as readonly Value[];
  if (lone && !Array.isArray(of)) return [of];
  return fail('seq.query.form.param', `A form of kind ${form.kind} takes a list of forms`, {
    hint: `(${form.kind} :of [(range :from 1 :to 9) (at :label 20)])`,
  });
}

const fold =
  (kind: 'merge' | 'intersect'): Resolve =>
  (ctx, form, prefix) => {
    const rs = operandValues(form, false).map((v) => evaluate(ctx, v, prefix));
    let acc = rs[0];
    for (let i = 1; i < rs.length; i++) acc = pair(kind, acc, rs[i]);
    return limited(ctx, acc);
  };

/** The complement of a result within what the prefix names. */
function complement(ctx: Ctx, r: QueryResult, prefix: Named): QueryResult {
  if (r.type === 'location') {
    if (r.location.rects.length)
      return fail('seq.query.form.unsupported', 'A rect has no complement that is a location', {
        hint: '(except :of [(rect :a … :b …) (rect :a … :b …)])',
      });
    const all = location({
      parts: systemsOf(ctx, prefix, 'except').map((s) => ({ system: s.address, ranges: [0, s.length] })),
    });
    return { type: 'location', location: subtractLocations(all, r.location) };
  }
  if (r.type === 'keys') {
    let all = EmptyKeys;
    for (const a of annotationsOf(ctx, prefix, 'except')) {
      const read = a.items();
      stage(ctx, read.pending);
      all = unionKeys(all, keySet(read.items.map((it) => ({ owner: a.address, key: it.key }))));
    }
    return limited(ctx, { type: 'keys', keys: subtractKeys(all, r.keys) });
  }
  const children: Selected[] = [];
  for (const m of prefixNodes(prefix, 'except'))
    for (const c of m.node.children ?? [])
      if (c.id !== undefined) children.push({ node: c, address: [...m.address, c.id] });
  return limited(ctx, { type: 'nodes', nodes: combineNodes(children, r.nodes, (x, y) => x && !y) });
}

const except: Resolve = (ctx, form, prefix) => {
  const rs = operandValues(form, true).map((v) => evaluate(ctx, v, prefix));
  if (rs.length === 1) return complement(ctx, rs[0], prefix);
  let acc = rs[0];
  for (let i = 1; i < rs.length; i++) acc = pair('except', acc, rs[i]);
  return limited(ctx, acc);
};

/** The combinators, by kind. */
export const Combinators = { merge: fold('merge'), intersect: fold('intersect'), except } as const;
