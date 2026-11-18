/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import {
  type AddressPattern,
  exactAddress,
  type FormResolver,
  type FormResolverInput,
  matchesAddress,
  type Named,
  type Node,
  parseAddress,
  type Segment,
  type Selected,
  printAddress,
  resolveAddress,
  selectAddress,
} from '../index.js';
import { ChildIndexThreshold, childIndex } from '../address.js';

const print = (p: AddressPattern) => printAddress(p.segments, p.subtree);

const Spec: Node = {
  kind: 'spec',
  children: [
    {
      kind: 'data',
      id: 'data',
      children: [
        { kind: 'sequence', id: 'p53', params: { alphabet: 'protein', tags: ['a', 'b'] } },
        { kind: 'sequence', id: 'tp53', params: { alphabet: 'dna' } },
        { kind: 'ranges', id: 'domains', params: { target: '/data/p53' }, children: [{ kind: 'x', id: 'x1' }] },
      ],
    },
    {
      kind: 'views',
      id: 'views',
      children: [
        {
          kind: 'view',
          id: 'main',
          children: [
            { kind: 'section', id: 'protein', params: { of: '/data/p53' } },
            { kind: 'section', id: 'gene', params: { of: '/data/tp53', follow: '/views/main/protein' } },
          ],
        },
      ],
    },
  ],
};

describe('repr/tree: address grammar', () => {
  it('reads and prints canonically: parse(print(a)) = a', () => {
    const texts = [
      '/',
      '/*',
      '/data/p53',
      '/data/*',
      '/data/p53/(range :from 94 :numbering uniprot :to 292)',
      '/data/1tup/(chain :id A)/(at :label 94A :numbering auth)',
      '/data/globins/(row :id "HBA_HUMAN/1-141")',
      '/data/globins/(rect :columns (range :from 30 :to 60) :rows (rows :from "HBA_HUMAN/1-141" :to "HBB_HUMAN/1-146"))',
      '/data/meth/(rect :positions (range :from 7610000 :to 7620000) :states (states :in [s1 s2]))',
      '/views/main/(section :of /data/p53)',
      '/data/genes/(node :key T1)',
      '/data/x/(f :a -1 :b 1e+21 :c nan :d nil :e true :f "a b" :g 5n :h [] :i ext:kind)',
      '/data/(ext:item :"odd key" 1)',
    ];
    for (const text of texts) {
      const { pattern, issues } = parseAddress(text);
      expect(issues, text).toEqual([]);
      expect(print(pattern!)).toBe(text);
      expect(parseAddress(print(pattern!)).pattern).toEqual(pattern);
    }
  });

  it('reads atoms as numbers when they parse as JavaScript numbers, else as symbols', () => {
    const { pattern } = parseAddress('/a/(f :n 1e5 :m -1 :l 94A :s HBA_HUMAN :q "94" :addr /b/(c :x 1))');
    expect((pattern!.segments[1] as Node).params).toEqual({
      n: 1e5,
      m: -1,
      l: '94A',
      s: 'HBA_HUMAN',
      q: '94',
      addr: '/b/(c :x 1)',
    });
  });

  it('canonicalizes spacing, key order and an address value', () => {
    const { pattern } = parseAddress('/a/( f   :y /b/( c :z 1 :a 2 )  :x  [ 1  2 ] )');
    expect(print(pattern!)).toBe('/a/(f :x [1 2] :y /b/(c :a 2 :z 1))');
  });

  it('reports what is not an address, by code and offset', () => {
    const cases: [string, string][] = [
      ['a', 'repr.tree.address.syntax'],
      ['/a/', 'repr.tree.address.syntax'],
      ['/a//b', 'repr.tree.address.syntax'],
      ['/*/a', 'repr.tree.address.syntax'],
      ['/a/(f :x 1', 'repr.tree.address.syntax'],
      ['/a/(f x)', 'repr.tree.address.form'],
      ['/a/(F)', 'repr.tree.node.kind'],
      ['/a/(f :x)', 'repr.tree.address.form'],
      ['/a/(f :x 1 :x 2)', 'repr.tree.address.form'],
      ['/a/(f :kind 1)', 'repr.tree.value.reserved_key'],
      ['/a/(f :x (:y 1))', 'repr.tree.address.form'],
      ['/a/(f :x /b//c)', 'repr.tree.address.syntax'],
    ];
    for (const [text, code] of cases) {
      const { pattern, issues } = parseAddress(text);
      expect(pattern, text).toBeUndefined();
      expect(
        issues.map((i) => i.code),
        text,
      ).toEqual([code]);
      expect(issues[0].severity).toBe('error');
      expect(typeof issues[0].offset).toBe('number');
    }
  });

  it('names one node by an exact address', () => {
    expect(exactAddress('/data/p53')).toBe('/data/p53');
    expect(exactAddress('/data/*')).toBeUndefined();
    expect(exactAddress('/data/(sequence)')).toBeUndefined();
    expect(exactAddress('/')).toBeUndefined();
  });
});

describe('repr/tree: selection by child selection (the default resolver)', () => {
  const select = (a: string) => selectAddress(Spec, a, { limit: 10 });
  const addresses = (a: string) => select(a).matches.map((m) => printAddress(m.address));

  it('steps by ids and selects by forms, in document order', () => {
    expect(addresses('/data/p53')).toEqual(['/data/p53']);
    expect(addresses('/data/(sequence)')).toEqual(['/data/p53', '/data/tp53']);
    expect(addresses('/data/(sequence :alphabet dna)')).toEqual(['/data/tp53']);
    expect(addresses('/data/(sequence :tags a)')).toEqual(['/data/p53']);
    expect(addresses('/data/(sequence :tags [a b])')).toEqual(['/data/p53']);
    expect(addresses('/views/main/(section :of /data/p53)')).toEqual(['/views/main/protein']);
    expect(addresses('/data/*')).toEqual(['/data', '/data/p53', '/data/tp53', '/data/domains', '/data/domains/x1']);
    expect(addresses('/data/nope')).toEqual([]);
  });

  it('is bounded: more matches than the limit is an issue, never a truncation', () => {
    const r = selectAddress(Spec, '/data/*', { limit: 2 });
    expect(r.matches).toEqual([]);
    expect(r.issues.map((i) => i.code)).toEqual(['repr.tree.address.limit_exceeded']);
    expect(selectAddress(Spec, '/data', { limit: -1 }).issues[0].code).toBe('repr.tree.address.option');
    expect(selectAddress(Spec, '/data/*', { limit: 10, maxDepth: 1 }).issues[0].code).toBe(
      'repr.tree.address.depth_exceeded',
    );
  });

  it("checks forms against the caller's kinds and params", () => {
    const kinds = new Set(['sequence']);
    expect(selectAddress(Spec, '/data/(ranges)', { limit: 9, kinds }).issues[0].code).toBe(
      'repr.tree.address.unknown_kind',
    );
    const params = new Map([['sequence', new Set(['alphabet'])]]);
    expect(selectAddress(Spec, '/data/(sequence :x 1)', { limit: 9, params }).issues[0].code).toBe(
      'repr.tree.address.unknown_param',
    );
  });

  it('stops on an aborted signal', () => {
    const c = new AbortController();
    c.abort();
    expect(selectAddress(Spec, '/data', { limit: 9, signal: c.signal }).issues[0].code).toBe(
      'repr.tree.address.cancelled',
    );
  });

  it('matches a pattern against an address', () => {
    const p = parseAddress('/data/(sequence :alphabet protein)/*').pattern!;
    expect(matchesAddress(p, ['data', 'p53'], Spec)).toBe(true);
    expect(matchesAddress(p, ['data', 'tp53'], Spec)).toBe(false);
    expect(() => matchesAddress(p, ['data', 'p53'])).toThrow();
    expect(matchesAddress(parseAddress('/data/*').pattern!, ['data', 'x'])).toBe(true);
  });
});

describe('repr/tree: the form-resolver table', () => {
  const nodesOf = (n: Named) => (n.kind === 'nodes' ? n.matches : []);
  /** A range over what the prefix names: a value the tree does not know; a warning when it is empty. */
  const range: FormResolver = ({ form, prefix }) => ({
    result: {
      kind: 'other',
      type: 'location',
      value: nodesOf(prefix).map((m) => ({ system: printAddress(m.address), from: form.params?.from })),
    },
    ...(form.params?.from === 0 ? { issues: [{ code: 'test.empty', severity: 'warning', message: 'Empty' }] } : {}),
    ...(form.params?.from === -1 ? { complete: false, pending: ['region 1'] } : {}),
  });
  /** Children by kind at any depth below: nodes, so ids and `/*` may follow. */
  const deep: FormResolver = ({ form, prefix }) => {
    const out: Selected[] = [];
    const go = (n: Node, a: readonly Segment[]) => {
      for (const c of n.children ?? []) {
        if (c.id === undefined) continue;
        if (c.kind === form.params?.of) out.push({ node: c, address: [...a, c.id] });
        go(c, [...a, c.id]);
      }
    };
    for (const m of nodesOf(prefix)) go(m.node, m.address);
    return { result: { kind: 'nodes', matches: out } };
  };
  /** A dynamic node no id names, addressed by its form: `/data/genes/(node :key T1)`. */
  const dynamic: FormResolver = ({ form, prefix }) => ({
    result: {
      kind: 'nodes',
      matches: nodesOf(prefix).map((m) => ({
        node: { kind: 'gene', params: { key: form.params?.key ?? null } },
        address: [...m.address, form],
      })),
    },
  });
  /** What it is given, to check the bounds reach it. */
  let seen: FormResolverInput | undefined;
  const spy: FormResolver = (input) => ((seen = input), { result: input.prefix });
  const failing: FormResolver = () => ({ issues: [{ code: 'test.fail', severity: 'error', message: 'No' }] });
  const resolvers = { range, deep, node: dynamic, spy, failing };

  it('hands what the prefix names to the resolver, and returns its result', () => {
    const r = selectAddress(Spec, '/data/p53/(range :from 94 :to 292)', { limit: 9, resolvers });
    expect(r.issues).toEqual([]);
    expect(r.matches).toEqual([]);
    expect(r.result).toEqual({ kind: 'other', type: 'location', value: [{ system: '/data/p53', from: 94 }] });
  });

  it('passes the bounds to the resolver', () => {
    const signal = new AbortController().signal;
    selectAddress(Spec, '/data/(spy)', { limit: 3, maxDepth: 7, resolvers, signal, deadline: 1e15 });
    expect(seen).toMatchObject({ limit: 3, maxDepth: 7, signal, deadline: 1e15, root: Spec });
  });

  it("carries a resolver's warnings and its completeness; only its errors fail", () => {
    const w = selectAddress(Spec, '/data/p53/(range :from 0)', { limit: 9, resolvers });
    expect(w.result?.type).toBe('location');
    expect(w.issues.map((i) => [i.code, i.severity])).toEqual([['test.empty', 'warning']]);
    const p = selectAddress(Spec, '/data/p53/(range :from -1)', { limit: 9, resolvers });
    expect([p.complete, p.pending]).toEqual([false, ['region 1']]);
    const e = selectAddress(Spec, '/data/(failing)', { limit: 9, resolvers });
    expect([e.matches, e.result, e.issues.map((i) => i.code)]).toEqual([[], undefined, ['test.fail']]);
  });

  it('names a node by a form where no id does', () => {
    const r = selectAddress(Spec, '/data/(node :key T1)', { limit: 9, resolvers });
    expect(r.matches.map((m) => printAddress(m.address))).toEqual(['/data/(node :key T1)']);
  });

  it('lets ids, child selection and /* follow a resolver of nodes; the limit bounds the final matches', () => {
    const r = selectAddress(Spec, '/(deep :of section)', { limit: 9, resolvers });
    expect(r.matches.map((m) => printAddress(m.address))).toEqual(['/views/main/protein', '/views/main/gene']);
    const s = selectAddress(Spec, '/(deep :of ranges)/*', { limit: 9, resolvers });
    expect(s.matches.map((m) => printAddress(m.address))).toEqual(['/data/domains', '/data/domains/x1']);
    // two sections in a step, one in the end: within a limit of one
    const first: FormResolver = ({ prefix }) => ({ result: { kind: 'nodes', matches: nodesOf(prefix).slice(-1) } });
    const one = selectAddress(Spec, '/(deep :of section)/(first)', { limit: 1, resolvers: { ...resolvers, first } });
    expect(one.matches.map((m) => printAddress(m.address))).toEqual(['/views/main/gene']);
    expect(selectAddress(Spec, '/(deep :of section)', { limit: 1, resolvers }).issues[0].code).toBe(
      'repr.tree.address.limit_exceeded',
    );
    // a resolver's form kind is the resolver's to check, not the language's
    expect(selectAddress(Spec, '/(deep :of view)', { limit: 9, resolvers, kinds: new Set() }).issues).toEqual([]);
  });

  it('stops on an aborted signal between steps', () => {
    const c = new AbortController();
    const abort: FormResolver = ({ prefix }) => (c.abort(), { result: prefix });
    const r = selectAddress(Spec, '/data/(abort)/p53', { limit: 9, resolvers: { abort }, signal: c.signal });
    expect(r.issues.map((i) => i.code)).toEqual(['repr.tree.address.cancelled']);
  });

  it('refuses an id, a form without a resolver or /* after a result that is not nodes', () => {
    for (const a of [
      '/data/p53/(range :from 1)/x',
      '/data/p53/(range :from 1)/(item)',
      '/data/p53/(range :from 1)/*',
    ]) {
      const r = selectAddress(Spec, a, { limit: 9, resolvers });
      expect(
        r.issues.map((i) => i.code),
        a,
      ).toEqual(['repr.tree.address.resolver']);
    }
  });
});

describe('repr/tree: a numeric-looking key', () => {
  it('matches a string id only when quoted', () => {
    const tree: Node = { kind: 'r', children: [{ kind: 'chain', id: 'c', params: { id: '1' } }] };
    const ids = (a: string) => selectAddress(tree, a, { limit: 9 }).matches.length;
    expect(ids('/(chain :id "1")')).toBe(1);
    expect(ids('/(chain :id 1)')).toBe(0);
  });
});

describe('repr/tree: the child index', () => {
  const wide = (n: number): Node => ({
    kind: 'p',
    children: Array.from({ length: n }, (_, i) => ({ kind: 'c', id: `c${i}`, params: { i } })),
  });

  it('finds a child by id below and above the flat threshold, the first of two with one id', () => {
    for (const n of [ChildIndexThreshold, ChildIndexThreshold + 1, 1000]) {
      const tree = wide(n);
      expect(childIndex(tree.children, `c${n - 1}`)).toBe(n - 1);
      expect(resolveAddress(tree, [`c${n - 1}`])?.params?.i).toBe(n - 1);
      expect(childIndex(tree.children, 'nope')).toBe(-1);
      const twice = [...tree.children!, { kind: 'c', id: 'c0', params: { i: -1 } }];
      expect(childIndex(twice, 'c0')).toBe(0);
    }
    expect(childIndex(undefined, 'a')).toBe(-1);
  });

  it('is per children array: a path copy is indexed again', () => {
    const tree = wide(20);
    expect(childIndex(tree.children, 'c5')).toBe(5);
    const copy = [{ kind: 'c', id: 'new' }, ...tree.children!];
    expect(childIndex(copy, 'c5')).toBe(6);
    expect(childIndex(tree.children, 'c5')).toBe(5);
  });
});
