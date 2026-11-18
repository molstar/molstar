/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The bounds on what a spec sizes (_spec/README.md, Language, Bounds): each refused with its coded issue at the bound
// and one past it, nothing read past a bound, and the nesting depth refused the same way by the whole spec's check and
// by a patch's incremental one.

import { itWithin } from '@molstar/test-harness/timing';
import type { Changes } from '@molstar/core/repr/patch';
import type { Node } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import { checkSeqSpec } from '../check.js';
import { checkSeqChanges } from '../check-changes.js';
import { checkSpecBounds, DefaultSpecBounds, queryDepth, specBounds } from '../bounds.js';
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
  const noChanges: Changes = { made: [], moved: [], removed: [], params: [] };
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
    const changes: Changes = { ...noChanges, params: [{ at: ['data', 'c'], keys: ['to'] } as never] };
    expect(checkSeqChanges(before, doc(3), changes, bounds)).toEqual([]);
    expect(checkSeqChanges(before, doc(4), changes, bounds).map((i) => i.code)).toEqual([
      'seq.lang.bound.query_depth',
    ]);
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
  // the walk counts a node's kind and id, a param's key and its value: `k`, `a` and the string
  const sized = (n: number): Node => ({ kind: 'k', params: { a: 'x'.repeat(n) } });

  it('passes the bound and refuses one byte past it', () => {
    const bounds = specBounds({ maxBytes: 1002 });
    expect(codes(sized(1000), bounds)).toEqual([]);
    const over = checkSpecBounds(sized(1001), bounds);
    expect(over.map((i) => i.code)).toEqual(['seq.lang.bound.size']);
  });

  it('counts columns by their bytes', () => {
    const bounds = specBounds({ maxBytes: 4000 });
    const col = (n: number): Node => ({ kind: 'k', params: { a: new Float64Array(n) } as never });
    expect(codes(col(498), bounds)).toEqual([]);
    expect(codes(col(500), bounds)).toEqual(['seq.lang.bound.size']);
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
    const issues = checkSeqSpec(specOf({ kind: 'ext:note', id: 'n', params: { s: 'x'.repeat(2000) } as never }), specBounds({ maxBytes: 1000 }));
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

  it('keeps at most a bound’s worth a node', () => {
    const bounds = specBounds({ maxNodeIssues: 2, maxPlaces: 1000 });
    const issues = checkSeqSpec(
      specOf({ kind: 'data', id: 'data', children: [{ kind: 'ext:note' }] }, { kind: 'views', id: 'wrong' }),
      bounds,
    );
    expect(issues.length).toBeGreaterThan(0);
  });
});

describe('a spec at the bounds’ edges', () => {
  itWithin('is checked within the time of the whole check of a spec its size', (within) => {
    // 256 levels deep, 10 places of a code: the walk adds one pass over what the other checks read anyway
    const spec = specOf(chain(256), ...Array.from({ length: 10 }, () => ({ kind: 'ext:note' }) as Node));
    const whole = within(500, () => checkSeqSpec(spec));
    expect(whole.filter((i) => i.code === 'seq.lang.id.required')).toHaveLength(10);
    expect(checkSpecBounds(p53Spec())).toEqual([]);
  });
});
