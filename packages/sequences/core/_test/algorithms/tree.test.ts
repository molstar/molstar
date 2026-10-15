/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import {
  alignmentDistances,
  type GuideTree,
  MAX_TREE_ROWS,
  midpointRoot,
  neighborJoining,
  TreeSizeError,
  upgma,
} from '../../src/algorithms/index.js';

/** Row-major square matrix from a triangle written row by row (`[[b], [c, d]...]` as upper rows). */
function square(n: number, upper: number[][]): Float64Array {
  const d = new Float64Array(n * n);
  upper.forEach((row, i) => {
    row.forEach((v, k) => {
      d[i * n + (i + 1 + k)] = d[(i + 1 + k) * n + i] = v;
    });
  });
  return d;
}

/** The ancestors of a node, itself first, with the cumulative branch length to each. */
function up(t: GuideTree, v: number): { node: number; at: number }[] {
  const out = [];
  let at = 0;
  for (let u = v; u >= 0; u = t.parent[u]) {
    out.push({ node: u, at });
    at += t.branchLength[u] || 0;
  }
  return out;
}
/** The leaf-to-leaf path length. */
function path(t: GuideTree, a: number, b: number): number {
  const ua = up(t, a);
  const ub = new Map(up(t, b).map((x) => [x.node, x.at]));
  for (const x of ua) if (ub.has(x.node)) return x.at + (ub.get(x.node) as number);
  throw new Error('no common ancestor');
}
/** Each non-root inner node's leaves, as sorted strings `0,1`, sorted. */
function clades(t: GuideTree): string[] {
  const leaves = new Map<number, number[]>();
  for (let v = 0; v < t.leaves; v++)
    for (const { node } of up(t, v)) leaves.set(node, [...(leaves.get(node) ?? []), v]);
  const out: string[] = [];
  for (const [node, ls] of leaves)
    if (node >= t.leaves && t.parent[node] >= 0) out.push(ls.sort((x, y) => x - y).join(','));
  return out.sort();
}
/** The unrooted tree's splits: each inner clade as the side without leaf 0 (the root's place does not matter). */
function splits(t: GuideTree): string[] {
  return clades(t)
    .map((c) => c.split(',').map(Number))
    .map((ls) => (ls.includes(0) ? Array.from({ length: t.leaves }, (_, k) => k).filter((k) => !ls.includes(k)) : ls))
    .map((ls) => ls.join(','))
    .sort();
}
const total = (t: GuideTree) => t.branchLength.reduce((s, v) => s + (Number.isNaN(v) ? 0 : v), 0);
const near = (a: number, b: number) => Math.abs(a - b) < 1e-9;

describe('alignmentDistances', () => {
  const enc = (rows: string[]) => Uint8Array.from(rows.join('').replace(/-/g, '\0'), (c) => c.charCodeAt(0));
  it('is the share of mismatches over the columns where both rows have a letter', async () => {
    const rows = ['ACDEF', 'ACDEG', 'A-DEG', '-----'];
    const d = await alignmentDistances(enc(rows), 4, 5);
    expect(d[0 * 4 + 1]).toBeCloseTo(1 / 5);
    expect(d[1 * 4 + 0]).toBeCloseTo(1 / 5);
    expect(d[0 * 4 + 2]).toBeCloseTo(1 / 4); // the gap column is out
    expect(d[1 * 4 + 2]).toBe(0);
    expect(d[0 * 4 + 0]).toBe(0);
    expect(d[0 * 4 + 3]).toBe(1); // nothing in common: the largest
    expect(d[3 * 4 + 3]).toBe(0);
  });
  it('refuses more rows than the limit before allocating', async () => {
    await expect(alignmentDistances(new Uint8Array(0), MAX_TREE_ROWS + 1, 0)).rejects.toBeInstanceOf(TreeSizeError);
  });
});

describe('neighborJoining', () => {
  // The Wikipedia example (Saitou and Nei's): additive, so the tree is the one with a:2 b:3 c:4 d:2 e:1 and the inner
  // edges 3 ((a, b) to c) and 2 (c to (d, e)); the root sits where the last three nodes join.
  const D = square(5, [[5, 9, 9, 8], [10, 10, 9], [8, 7], [3]]);
  it('recovers an additive tree: its clades, leaf edges and every path length', async () => {
    const t = await neighborJoining(D, 5);
    expect(splits(t)).toEqual(['2,3,4', '3,4']); // (a, b) | (c, d, e) and (d, e) | (a, b, c)
    expect([0, 1, 2, 3, 4].map((v) => t.branchLength[v])).toEqual([2, 3, 4, 2, 1]);
    for (let i = 0; i < 5; i++) for (let j = i + 1; j < 5; j++) expect(near(path(t, i, j), D[i * 5 + j])).toBe(true);
    expect(near(total(t), 17)).toBe(true);
    expect(t.parent.filter((p) => p < 0)).toHaveLength(1);
    expect(Number.isNaN(t.branchLength[t.parent.indexOf(-1)])).toBe(true);
  });
  it('does not change the matrix and is deterministic, ties by row order', async () => {
    const copy = Float64Array.from(D);
    const a = await neighborJoining(D, 5);
    expect(D).toEqual(copy);
    expect(await neighborJoining(D, 5)).toEqual(a);
    const flat = square(4, [[1, 1, 1], [1, 1], [1]]);
    const t = await neighborJoining(flat, 4);
    expect(t.parent[0]).toBe(t.parent[1]); // the first pair in row order is joined first
    expect(t.parent[2]).not.toBe(t.parent[0]);
  });
  it('recovers random additive trees', async () => {
    let seed = 7;
    const rnd = () => {
      seed = (Math.imul(seed, 1103515245) + 12345) & 0x7fffffff;
      return seed / 0x7fffffff;
    };
    for (const n of [4, 9, 23]) {
      // A random tree: each leaf attached to the middle of a random edge, with random positive lengths.
      const parent: number[] = [-1, 0, 0];
      const len: number[] = [0, 0.1 + rnd(), 0.1 + rnd()];
      for (let leaves = 2; leaves < n; leaves++) {
        const e = 1 + Math.floor(rnd() * (parent.length - 1));
        const mid = parent.length;
        parent.push(parent[e], mid);
        len.push(0.1 + rnd(), 0.1 + rnd());
        parent[e] = mid;
      }
      const hasChild = new Set(parent);
      const leafNodes = parent.map((_, i) => i).filter((i) => !hasChild.has(i));
      const ancestors = (v: number) => {
        const out: [number, number][] = [];
        let at = 0;
        for (let u = v; u >= 0; u = parent[u]) {
          out.push([u, at]);
          at += len[u];
        }
        return out;
      };
      const rows = leafNodes.length;
      const D2 = new Float64Array(rows * rows);
      for (let i = 0; i < rows; i++)
        for (let j = i + 1; j < rows; j++) {
          const bj = new Map(ancestors(leafNodes[j]));
          const first = ancestors(leafNodes[i]).find(([u]) => bj.has(u)) as [number, number];
          D2[i * rows + j] = D2[j * rows + i] = first[1] + (bj.get(first[0]) as number);
        }
      const t = await neighborJoining(D2, rows);
      for (let i = 0; i < rows; i++)
        for (let j = i + 1; j < rows; j++) expect(near(path(t, i, j), D2[i * rows + j])).toBe(true);
    }
  });
  it('handles one, two and three rows', async () => {
    expect((await neighborJoining(new Float64Array(1), 1)).parent).toEqual(Int32Array.of(-1));
    const two = await neighborJoining(square(2, [[3]]), 2);
    expect(Array.from(two.parent)).toEqual([2, 2, -1]);
    expect([two.branchLength[0], two.branchLength[1]]).toEqual([1.5, 1.5]);
    const three = await neighborJoining(square(3, [[3, 4], [5]]), 3);
    expect(Array.from(three.branchLength.slice(0, 3))).toEqual([1, 2, 3]);
    expect(three.parent.slice(0, 3)).toEqual(Int32Array.of(3, 3, 3));
  });
  it('sets a negative branch length to 0', async () => {
    // Not additive: d(a, b) is too long for the other distances.
    const t = await neighborJoining(square(4, [[10, 2, 2], [2, 2], [1]]), 4);
    expect(t.branchLength.every((v) => Number.isNaN(v) || v >= 0)).toBe(true);
  });
});

describe('midpointRoot', () => {
  const D = square(5, [[5, 9, 9, 8], [10, 10, 9], [8, 7], [3]]);
  /** The farthest leaf below each of the root's two children, from the root. */
  function reach(t: GuideTree): number[] {
    const root = t.parent.indexOf(-1);
    const kids = Array.from(t.parent.keys()).filter((v) => t.parent[v] === root);
    const far = kids.map(() => 0);
    for (let v = 0; v < t.leaves; v++) {
      const chain = up(t, v);
      const above = new Set(chain.map((c) => c.node));
      const i = kids.findIndex((k) => above.has(k));
      far[i] = Math.max(far[i], chain[chain.length - 1].at);
    }
    return far;
  }
  it('puts a binary root at the middle of the longest path and keeps every distance', async () => {
    const t = await neighborJoining(D, 5);
    const r = midpointRoot(t);
    expect(r.leaves).toBe(5);
    expect(r.parent).toHaveLength(t.parent.length + 1);
    expect(r.parent.filter((p) => p === r.parent.length - 1)).toHaveLength(2);
    for (let i = 0; i < 5; i++) for (let j = i + 1; j < 5; j++) expect(near(path(r, i, j), D[i * 5 + j])).toBe(true);
    const [x, y] = reach(r);
    expect(near(x, y)).toBe(true);
    expect(near(x, 5)).toBe(true); // half of the longest path, 10
    expect(near(total(r), 17)).toBe(true);
    expect([...new Set(splits(r))]).toEqual(splits(t)); // the root's two sides are one split
  });
  it("keeps the leaves' rows and leaves small trees alone", async () => {
    const two = await neighborJoining(square(2, [[3]]), 2);
    expect(midpointRoot(two)).toBe(two);
    const big = square(6, [[2, 4, 6, 6, 6], [4, 6, 6, 6], [6, 6, 6], [2, 4], [4]]);
    const r = midpointRoot(await neighborJoining(big, 6));
    for (let i = 0; i < 6; i++) for (let j = i + 1; j < 6; j++) expect(near(path(r, i, j), big[i * 6 + j])).toBe(true);
  });
  it('leaves a tree of identical rows as it is: no leaf gets children, no branch is NaN (CB-02)', async () => {
    const t = await neighborJoining(new Float64Array(16), 4);
    const r = midpointRoot(t);
    expect(r).toBe(t);
    for (let v = 0; v < r.leaves; v++) expect(r.parent.includes(v)).toBe(false);
    for (let v = 0; v < r.parent.length; v++) if (r.parent[v] >= 0) expect(Number.isNaN(r.branchLength[v])).toBe(false);
  });
});

describe('upgma', () => {
  // The Wikipedia example: ((a:8.5, b:8.5):2.5, e:11):5.5 and (c:14, d:14):2.5 under one root, at 16.5.
  const D = square(5, [[17, 21, 31, 23], [30, 34, 21], [28, 39], [43]]);
  it('builds the ultrametric tree of the textbook example', async () => {
    const t = await upgma(D, 5);
    expect(clades(t)).toEqual(['0,1', '0,1,4', '2,3']);
    expect([0, 1, 2, 3, 4].map((v) => t.branchLength[v])).toEqual([8.5, 8.5, 14, 14, 11]);
    const ab = t.parent[0];
    expect(t.branchLength[ab]).toBe(2.5);
    expect(t.branchLength[t.parent[ab]]).toBe(5.5);
    expect(t.branchLength[t.parent[2]]).toBe(2.5);
    for (let v = 0; v < 5; v++) expect(near(up(t, v).at(-1)?.at as number, 16.5)).toBe(true); // every leaf at one height
    expect(t.parent.filter((p) => p < 0)).toHaveLength(1);
  });
  it('is deterministic, ties by row order, and handles one and two rows', async () => {
    const flat = square(4, [[1, 1, 1], [1, 1], [1]]);
    const t = await upgma(flat, 4);
    expect(t.parent[0]).toBe(t.parent[1]);
    expect(await upgma(flat, 4)).toEqual(t);
    expect((await upgma(new Float64Array(1), 1)).parent).toEqual(Int32Array.of(-1));
    const two = await upgma(square(2, [[3]]), 2);
    expect([two.branchLength[0], two.branchLength[1]]).toEqual([1.5, 1.5]);
  });
});

describe('abort and size', () => {
  const big = (n: number) => {
    const d = new Float64Array(n * n);
    let seed = 1;
    for (let i = 0; i < n; i++)
      for (let j = i + 1; j < n; j++) {
        seed = (seed * 48271) % 2147483647;
        d[i * n + j] = d[j * n + i] = (seed % 1000) + 1;
      }
    return d;
  };
  it("rejects with the signal's reason, before starting or while running", async () => {
    const pre = AbortSignal.abort(new Error('newer request'));
    await expect(neighborJoining(big(8), 8, pre)).rejects.toThrow('newer request');
    await expect(upgma(big(8), 8, pre)).rejects.toThrow('newer request');
    await expect(alignmentDistances(new Uint8Array(8 * 8), 8, 8, pre)).rejects.toThrow('newer request');
    for (const build of [neighborJoining, upgma]) {
      const c = new AbortController();
      const running = build(big(700), 700, c.signal);
      setTimeout(() => c.abort(new Error('stop')), 1);
      await expect(running).rejects.toThrow('stop');
    }
  });
  it('lets other tasks run during a long build', async () => {
    let ticks = 0;
    const timer = setInterval(() => ticks++, 1);
    await neighborJoining(big(500), 500);
    clearInterval(timer);
    expect(ticks).toBeGreaterThan(0);
  });
  it('refuses too many rows', async () => {
    await expect(neighborJoining(new Float64Array(0), MAX_TREE_ROWS + 1)).rejects.toBeInstanceOf(TreeSizeError);
    await expect(upgma(new Float64Array(0), MAX_TREE_ROWS + 1)).rejects.toBeInstanceOf(TreeSizeError);
  });
});
