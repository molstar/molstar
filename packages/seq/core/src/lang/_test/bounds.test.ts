/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The bounds on what a spec sizes (_spec/README.md, Language, Bounds): each refused with its coded issue at the bound
// and one past it, nothing read past a bound, and the nesting depth refused the same way by the whole spec's check and
// by a patch's incremental one.

import type { Changes } from '@molstar/core/repr/patch';
import { categoryColumn, type Node } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import { checkSeqSpec } from '../check.js';
import { checkSeqChanges } from '../check-changes.js';
import { checkSpecBounds, DefaultSpecBounds, measureNode, queryDepth, specBounds } from '../bounds.js';
import { SeqLanguageVersion } from '../language.js';
import { p53Spec } from './p53.js';

/** A chain of `n` nodes, one under the other (ids `c`), as a child of the root. */
function chain(n: number): Node {
  let node: Node = { kind: 'ext:note', id: 'c' };
  for (let i = 1; i < n; i++) node = { kind: 'ext:note', id: 'c', children: [node] };
  return node;
}

const specOf = (...children: Node[]): Node => ({ kind: 'spec', params: { version: SeqLanguageVersion }, children });
const codes = (spec: Node, bounds = DefaultSpecBounds) => checkSpecBounds(spec, bounds).map((i) => i.code);

describe('the defaults', () => {
  it('are the security table’s', () => {
    expect(DefaultSpecBounds).toEqual({
      maxBytes: 64 * 1024 * 1024,
      maxDepth: 256,
      maxQueryDepth: 64,
      maxColumnLength: 2 ** 31 - 1,
      maxNodeIssues: 100,
      maxPlaces: 10,
    });
    expect(specBounds({ maxDepth: 4 }).maxDepth).toBe(4);
    expect(specBounds({ maxDepth: 4 }).maxBytes).toBe(64 * 1024 * 1024);
  });
});

describe('the nesting depth', () => {
  it('passes 256 levels and refuses the 257th', () => {
    expect(codes(specOf(chain(256)))).toEqual([]);
    const over = checkSpecBounds(specOf(chain(257)));
    expect(over.map((i) => i.code)).toEqual(['seq.lang.bound.depth']);
    expect(over[0].severity).toBe('error');
    expect(over[0].path?.length).toBe(2 * 257);
  });

  it('refuses a 20,000-deep spec with its issue alone, no stack spent', () => {
    const issues = checkSeqSpec(specOf(chain(20_000)));
    expect(issues.map((i) => i.code)).toEqual(['seq.lang.bound.depth']);
  });

  it('counts the values nested in params', () => {
    let v: unknown = 1;
    for (let i = 0; i < 300; i++) v = [v];
    const spec = specOf({ kind: 'ext:note', id: 'n', params: { v } as never });
    expect(codes(spec)).toEqual(['seq.lang.bound.depth']);
    expect(codes(spec, specBounds({ maxDepth: 400 }))).toEqual([]);
  });

  it('is a setting', () => {
    expect(codes(specOf(chain(5)), specBounds({ maxDepth: 5 }))).toEqual([]);
    expect(codes(specOf(chain(5)), specBounds({ maxDepth: 4 }))).toEqual(['seq.lang.bound.depth']);
  });
});

describe('the same depth in a patch’s incremental check', () => {
  const noChanges: Changes = { made: [], moved: [], removed: [], params: [], children: [] };
  const before = specOf({ kind: 'data', id: 'data' });
  const after = (n: number) => specOf({ kind: 'data', id: 'data', children: [chain(n)] });

  it('refuses what the whole check refuses, and passes what it passes', () => {
    const bounds = specBounds({ maxDepth: 8 });
    for (const n of [5, 6, 7, 8]) {
      const root = after(n);
      const changes: Changes = { ...noChanges, made: [['data', 'c']] };
      // `data` is a level, the chain starts at 2: the deepest is at 1 + n
      const whole = checkSpecBounds(root, bounds).map((i) => i.code);
      const incremental = checkSeqChanges(before, root, changes, bounds).map((i) => i.code);
      expect(incremental.filter((c) => c.startsWith('seq.lang.bound'))).toEqual(whole);
    }
    const root = after(7);
    expect(checkSeqChanges(before, root, { ...noChanges, made: [['data', 'c']] }, bounds).map((i) => i.code)).toEqual(
      [],
    );
    const deep = after(8);
    expect(checkSeqChanges(before, deep, { ...noChanges, made: [['data', 'c']] }, bounds).map((i) => i.code)).toEqual([
      'seq.lang.bound.depth',
    ]);
  });

  it('refuses a moved subtree too, and a 20,000-deep one without recursing', () => {
    const root = after(20_000);
    const issues = checkSeqChanges(before, root, { ...noChanges, moved: [{ from: ['x'], to: ['data', 'c'] }] });
    expect(issues.map((i) => i.code)).toEqual(['seq.lang.bound.depth']);
  });

  it('reads the changed node’s own params for a query that nests too deep', () => {
    const q = (n: number) => `/data/${'(a :b '.repeat(n)}x${')'.repeat(n)}`;
    const doc = (n: number) =>
      specOf({ kind: 'data', id: 'data', children: [{ kind: 'ext:note', id: 'c', params: { to: q(n) } as never }] });
    const bounds = specBounds({ maxQueryDepth: 3 });
    const changes: Changes = { ...noChanges, params: [{ at: ['data', 'c'], keys: ['to'], where: 'params' as const }] };
    expect(checkSeqChanges(before, doc(3), changes, bounds)).toEqual([]);
    expect(checkSeqChanges(before, doc(4), changes, bounds).map((i) => i.code)).toEqual(['seq.lang.bound.query_depth']);
  });
});

describe('the query depth', () => {
  it('counts the nesting of forms outside quotes', () => {
    expect(queryDepth('/a/b')).toBe(0);
    expect(queryDepth('/a/(x :k (y :k (z)))')).toBe(3);
    expect(queryDepth('/a/(x :k "((((" )')).toBe(1);
    expect(queryDepth('/a/(x :k "\\"(((" (y))')).toBe(2);
    expect(queryDepth(')))((')).toBe(2);
  });

  it('passes 64 nested forms and refuses the 65th, in a spec', () => {
    const q = (n: number) => `/data/${'(a :b '.repeat(n)}x${')'.repeat(n)}`;
    const spec = (n: number) => specOf({ kind: 'ext:note', id: 'n', params: { to: q(n) } as never });
    expect(codes(spec(64))).toEqual([]);
    const over = checkSpecBounds(spec(65));
    expect(over.map((i) => i.code)).toEqual(['seq.lang.bound.query_depth']);
    expect(over[0].path).toEqual(['children', 0, 'params', 'to']);
  });
});

describe('the spec size', () => {
  // the walk counts a node (8 and its kind), its params record (8), a key (8 and its length) and its value: 26 and the string
  const sized = (n: number): Node => ({ kind: 'k', params: { a: 'x'.repeat(n) } });

  it('passes the bound and refuses one byte past it', () => {
    const bounds = specBounds({ maxBytes: 1026 });
    expect(codes(sized(1000), bounds)).toEqual([]);
    const over = checkSpecBounds(sized(1001), bounds);
    expect(over.map((i) => i.code)).toEqual(['seq.lang.bound.size']);
  });

  it('counts columns by their bytes', () => {
    const bounds = specBounds({ maxBytes: 4000 });
    const col = (n: number): Node => ({ kind: 'k', params: { a: new Float64Array(n) } as never });
    expect(codes(col(496), bounds)).toEqual([]);
    expect(codes(col(497), bounds)).toEqual(['seq.lang.bound.size']);
  });

  it('stops reading at the bound: what is after it is never visited', () => {
    let visited = 0;
    const later: Node = {
      kind: 'k',
      get params() {
        visited++;
        return {};
      },
    } as Node;
    const spec: Node = { kind: 'k', children: [sized(100), later] };
    expect(codes(spec, specBounds({ maxBytes: 50 }))).toEqual(['seq.lang.bound.size']);
    expect(visited).toBe(0);
    expect(codes(spec, specBounds({ maxBytes: 5000 }))).toEqual([]);
    expect(visited).toBe(1);
  });

  it('is refused alone by the whole check: the description never reads it', () => {
    const issues = checkSeqSpec(
      specOf({ kind: 'ext:note', id: 'n', params: { s: 'x'.repeat(2000) } as never }),
      specBounds({ maxBytes: 1000 }),
    );
    expect(issues.map((i) => i.code)).toEqual(['seq.lang.bound.size']);
  });
});

describe('a column’s length', () => {
  it('passes the bound and refuses one past it', () => {
    const bounds = specBounds({ maxColumnLength: 10 });
    const col = (n: number): Node => ({ kind: 'k', params: { a: new Uint8Array(n) } as never });
    expect(codes(col(10), bounds)).toEqual([]);
    const over = checkSpecBounds(col(11), bounds);
    expect(over.map((i) => i.code)).toEqual(['seq.lang.bound.column']);
    expect(over[0].path).toEqual(['params', 'a']);
  });
});

describe('the issues kept', () => {
  const notes = (n: number) =>
    specOf({ kind: 'data', id: 'data', children: Array.from({ length: n }, () => ({ kind: 'ext:note' }) as Node) });

  it('lists the first 10 places of a code and counts the rest', () => {
    const issues = checkSeqSpec(notes(10));
    expect(issues.filter((i) => i.code === 'seq.lang.id.required')).toHaveLength(10);
    expect(issues.some((i) => i.code === 'seq.lang.bound.issues')).toBe(false);
    const more = checkSeqSpec(notes(35));
    expect(more.filter((i) => i.code === 'seq.lang.id.required')).toHaveLength(10);
    const summary = more.filter((i) => i.code === 'seq.lang.bound.issues');
    expect(summary).toHaveLength(1);
    expect(summary[0].severity).toBe('warning');
    expect(summary[0].message).toContain('25 more');
  });

  it('keeps at most a bound’s worth a node, counting the rest', () => {
    // one ranges node with three problems: an interval, a bound past the letters, a strand
    const spec = specOf({
      kind: 'data',
      id: 'data',
      children: [
        { kind: 'sequence', id: 'p', params: { alphabet: 'protein', letters: 'MEEP' } },
        {
          kind: 'ranges',
          id: 'r',
          params: {
            target: '/data/p',
            rows: { start: new Int32Array([5]), end: new Int32Array([3]), strand: categoryColumn(['x']) },
          },
        },
      ],
    });
    const of = (b: Partial<typeof DefaultSpecBounds>) =>
      checkSeqSpec(spec, specBounds(b))
        .map((i) => i.code)
        .sort();
    expect(of({})).toEqual(['seq.lang.ranges.bounds', 'seq.lang.ranges.interval', 'seq.lang.ranges.strand']);
    expect(of({ maxNodeIssues: 2 })).toEqual([
      'seq.lang.bound.issues',
      'seq.lang.ranges.bounds',
      'seq.lang.ranges.interval',
    ]);
    expect(of({ maxNodeIssues: 0 })).toEqual([
      'seq.lang.bound.issues',
      'seq.lang.bound.issues',
      'seq.lang.bound.issues',
    ]);
  });
});

describe('a spec at the bounds’ edges', () => {
  const wide = Array.from({ length: 20_000 }, (_, i) => ({ kind: 'ext:note', id: `n${i}` }) as Node);
  // 256 levels deep, 10 places of a code, 20,000 nodes
  const edge = () =>
    specOf(
      chain(256),
      { kind: 'data', id: 'data', children: wide },
      ...Array.from({ length: 10 }, () => ({ kind: 'ext:note' }) as Node),
    );

  it('passes the walk at the edges, and the whole check keeps its places', () => {
    const spec = edge();
    expect(checkSpecBounds(spec)).toEqual([]);
    expect(checkSeqSpec(spec).filter((i) => i.code === 'seq.lang.id.required')).toHaveLength(10);
    expect(checkSpecBounds(p53Spec())).toEqual([]);
  });

  it('is walked in a fraction of the whole check’s time', { tags: ['heavy'] }, () => {
    const spec = edge();
    const time = (f: () => unknown) => {
      const t = performance.now();
      f();
      return performance.now() - t;
    };
    checkSpecBounds(spec);
    checkSeqSpec(spec);
    const walkMs = Math.min(...[1, 2, 3].map(() => time(() => checkSpecBounds(spec))));
    const wholeMs = Math.min(...[1, 2, 3].map(() => time(() => checkSeqSpec(spec))));
    expect(walkMs / wholeMs).toBeLessThan(0.5);
  });
});

describe('custom data and args, walked as values', () => {
  const deep = (n: number) => {
    let v: unknown = 1;
    for (let i = 0; i < n; i++) v = { a: v };
    return v;
  };
  const withCustom = (custom: unknown): Node => specOf({ kind: 'ext:note', id: 'n', custom } as never);
  const withArgs = (args: unknown): Node => specOf({ kind: 'ext:note', args } as never);

  it('count a string in custom or args in the size', () => {
    const bounds = specBounds({ maxBytes: 1000 });
    expect(codes(withCustom({ blob: 'x'.repeat(5000) }), bounds)).toEqual(['seq.lang.bound.size']);
    expect(codes(withArgs(['x'.repeat(5000)]), bounds)).toEqual(['seq.lang.bound.size']);
    expect(codes(withCustom({ blob: 'x'.repeat(100) }), bounds)).toEqual([]);
    expect(checkSeqSpec(withCustom({ blob: 'x'.repeat(5000) }), bounds).map((i) => i.code)).toEqual([
      'seq.lang.bound.size',
    ]);
  });

  it('count their nesting in the depth, 20,000 deep without a stack', () => {
    expect(codes(withCustom(deep(20_000)))).toEqual(['seq.lang.bound.depth']);
    expect(codes(withArgs([deep(20_000)]))).toEqual(['seq.lang.bound.depth']);
    expect(checkSeqSpec(withCustom(deep(20_000))).map((i) => i.code)).toEqual(['seq.lang.bound.depth']);
    expect(codes(withCustom(deep(100)))).toEqual([]);
  });

  it('count a query in args, not in custom', () => {
    const q = `/data/${'(a :b '.repeat(70)}x${')'.repeat(70)}`;
    expect(codes(withArgs([q]))).toEqual(['seq.lang.bound.query_depth']);
    expect(codes(withCustom({ q }))).toEqual([]);
  });

  it('count a column in custom against the column bound', () => {
    expect(codes(withCustom({ c: new Uint8Array(11) }), specBounds({ maxColumnLength: 10 }))).toEqual([
      'seq.lang.bound.column',
    ]);
  });

  it('are walked by the patch’s check too', () => {
    const noChanges: Changes = { made: [], moved: [], removed: [], params: [], children: [] };
    const mk = (custom: unknown) =>
      specOf({ kind: 'data', id: 'data', children: [{ kind: 'ext:note', id: 'c', custom } as never] });
    const bounds = specBounds({ maxDepth: 20 });
    const made: Changes = { ...noChanges, made: [['data', 'c']] };
    expect(checkSeqChanges(specOf(), mk(deep(30)), made, bounds).map((i) => i.code)).toEqual(['seq.lang.bound.depth']);
    // only the keys a change set are walked
    const other: Changes = { ...noChanges, params: [{ at: ['data', 'c'], keys: ['x'], where: 'custom' }] };
    expect(checkSeqChanges(specOf(), mk({ a: deep(30) }), other, bounds).map((i) => i.code)).toEqual([]);
    const set: Changes = { ...noChanges, params: [{ at: ['data', 'c'], keys: ['a'], where: 'custom' }] };
    expect(checkSeqChanges(specOf(), mk({ a: deep(30) }), set, bounds).map((i) => i.code)).toEqual([
      'seq.lang.bound.depth',
    ]);
    const blob = mk({ s: 'x'.repeat(500) });
    expect(checkSeqChanges(specOf(), blob, made, specBounds({ maxBytes: 100 })).map((i) => i.code)).toEqual([
      'seq.lang.bound.size',
    ]);
  });
});

describe('the size estimate charges structure', () => {
  const bounds = specBounds({ maxBytes: 100 });

  it('refuses a million nodes and a million empty records before pushing them', () => {
    const kid: Node = { kind: 'k' };
    const kids = new Array<Node>(2_000_000).fill(kid);
    const before = process.memoryUsage().heapUsed;
    expect(codes({ kind: 'k', children: kids }, bounds)).toEqual(['seq.lang.bound.size']);
    const records = new Array(2_000_000).fill({});
    expect(codes({ kind: 'k', params: { a: records } as never }, bounds)).toEqual(['seq.lang.bound.size']);
    // 2M frames would be hundreds of MB: the walk refused on the count
    expect(process.memoryUsage().heapUsed - before).toBeLessThan(80 * 1024 * 1024);
  });

  it('charges a fixed overhead for each key and element, and a bigint by its digits', () => {
    const one = (params: unknown): number => measureNode({ kind: 'k', params }, specBounds({ maxBytes: 1e9 })).bytes;
    const empty = one({});
    expect(one({ a: 1 }) - empty).toBe(8 + 1 + 8);
    expect(one({ a: [1, 2, 3] }) - empty).toBe(8 + 1 + (3 + 1) * 8 + 24);
    expect(one({ a: 12345678901234567890n }) - empty).toBe(8 + 1 + 20);
    expect(one({ a: null }) - empty).toBe(8 + 1 + 8);
  });
});

describe('a node built in code', () => {
  const shape = (spec: unknown) => codes(spec as Node);
  it('is refused, never thrown on', () => {
    expect(shape({ kind: 'k', children: [null] })).toEqual(['seq.lang.bound.shape']);
    expect(shape({ kind: 'k', children: ['text'] })).toEqual(['seq.lang.bound.shape']);
    expect(shape({ kind: 'k', children: [{ id: 'x' }] })).toEqual(['seq.lang.bound.shape']);
    expect(shape({ kind: 'k', children: { length: 1e9 } })).toEqual(['seq.lang.bound.shape']);
    expect(shape('not a node')).toEqual(['seq.lang.bound.shape']);
    expect(checkSeqSpec({ kind: 'spec', children: [null] } as never).map((i) => i.code)).toEqual([
      'seq.lang.bound.shape',
    ]);
  });
});

describe('the patch check keeps its issues as the whole check does', () => {
  it('lists 10 places of a code and counts the rest, for 20,000 made nodes with a dangling source', () => {
    const n = 20_000;
    const kids = Array.from({ length: n }, (_, i) => ({
      kind: 'sequence',
      id: `s${i}`,
      params: { alphabet: 'protein', source: '/data/missing' },
    })) as Node[];
    const after = specOf({ kind: 'data', id: 'data', children: kids });
    const changes: Changes = {
      made: kids.map((k) => ['data', k.id!]),
      moved: [],
      removed: [],
      params: [],
      children: [],
    };
    const issues = checkSeqChanges(specOf({ kind: 'data', id: 'data' }), after, changes);
    expect(issues.filter((i) => i.code === 'seq.lang.ref.target')).toHaveLength(10);
    const summary = issues.filter((i) => i.code === 'seq.lang.bound.issues');
    expect(summary).toHaveLength(1);
    expect(summary[0].message).toContain(`${n - 10} more`);
  });
});
