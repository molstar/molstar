/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `seq.lang.*` checks of a change: on every patch here, what `checkSeqChanges` reports is exactly what the whole
// check finds after the patch and did not find before it (by code, node and place in the node, through moves); an
// error already present is not reported again; and as the patch validator, a patch is refused for what it introduced.

import { applyPatch, type Changes, type Patch } from '@molstar/core/repr/patch';
import type { Issue, Node } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import { createSeqSpec } from '../builder.js';
import { checkSeqLang } from '../check.js';
import { checkSeqChanges, SeqPatchValidator } from '../check-changes.js';
import { SeqLanguage, SeqValidateOptions } from '../language.js';
import { moleculesSpec } from './molecules.js';
import { P53_LETTERS, p53Spec } from './p53.js';

const options = { language: SeqLanguage, validate: SeqValidateOptions, allowDangling: true, allowRebind: true };

const key = (i: Issue, address = i.address ?? '/') => {
  const path = (i.path ?? []) as readonly (string | number)[];
  const at = path.indexOf('params');
  return `${i.code}|${address}|${at < 0 ? '' : path.slice(at + 1).join('/')}`;
};

function forward(address: string, changes: Changes): string {
  let ids = address.split('/').filter(Boolean);
  for (const m of changes.moved)
    if (m.from.every((s, i) => ids[i] === s) && ids.length >= m.from.length)
      ids = [...m.to, ...ids.slice(m.from.length)];
  return `/${ids.join('/')}`;
}

/** The patch applied (without the validator), and both answers: the change's, and the whole check's difference. */
function both(root: Node, patch: Patch) {
  const r = applyPatch(root, patch, options);
  if (!r.ok) throw new Error(r.issues.map((i) => i.code).join());
  const was = new Set(checkSeqLang(root).map((i) => key(i, i.address && forward(i.address, r.changes))));
  const whole = checkSeqLang(r.root)
    .map((i) => key(i))
    .filter((k) => !was.has(k))
    .sort();
  const change = checkSeqChanges(root, r.root, r.changes)
    .map((i) => key(i))
    .sort();
  return { whole, change, root: r.root };
}

/** p53 with a user ref the domains layer draws through, and a domain past the protein's end. */
function withRefAndError(): Node {
  const b = createSeqSpec({ title: 'p53' });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53_LETTERS });
  const domains = b.data.ranges('domains', {
    target: p53,
    rows: { start: new Int32Array([0, 17]), end: new Int32Array([14, 30]) },
  });
  b.data.ranges('broken', { target: p53, rows: { start: new Int32Array([0]), end: new Int32Array([99]) } });
  const r = b.refs.ref('d', { to: domains });
  b.views.view('main', (v) => {
    const protein = v.section('protein', { of: p53 });
    protein.track('domains').layer('ranges', { type: 'ranges', data: String(r) });
    v.section('again', { of: p53, follow: protein });
  });
  return b.build();
}

const Cases: { name: string; root: () => Node; patch: Patch; introduces: readonly string[] }[] = [
  {
    name: 'a layer hidden',
    root: p53Spec,
    patch: { ops: [{ op: 'set', at: ['views', 'main', 'protein', 'letters', 'letters'], key: 'hidden', value: true }] },
    introduces: [],
  },
  {
    name: 'the protein shortened: the domains and the values past it',
    root: p53Spec,
    patch: { ops: [{ op: 'set', at: ['data', 'p53'], key: 'letters', value: 'MEEPQ' }] },
    introduces: ['seq.lang.ranges.bounds|/data/domains|rows/end/0', 'seq.lang.values.length|/data/conservation|values'],
  },
  {
    name: 'a follow cycle closed',
    root: p53Spec,
    patch: { ops: [{ op: 'set', at: ['views', 'main', 'protein'], key: 'follow', value: '/views/main/gene' }] },
    introduces: [
      'seq.lang.section.follow|/views/main/gene|follow',
      'seq.lang.section.follow|/views/main/protein|follow',
    ],
  },
  {
    name: 'a section range past its system',
    root: p53Spec,
    patch: { ops: [{ op: 'set', at: ['views', 'main', 'protein'], key: 'range', value: [0, 300] }] },
    introduces: ['seq.lang.section.range|/views/main/protein|range'],
  },
  {
    name: 'a ref a layer draws through retargeted to a sequence: the layer, through the ref',
    root: withRefAndError,
    patch: { ops: [{ op: 'set', at: ['refs', 'd'], key: 'to', value: '/data/p53' }] },
    introduces: ['seq.lang.layer.data|/views/main/protein/domains/ranges|data'],
  },
  {
    name: 'an edit beside an error already there: nothing new',
    root: withRefAndError,
    patch: { ops: [{ op: 'set', at: ['data', 'broken'], key: 'target', value: '/data/p53' }] },
    introduces: [],
  },
  {
    name: 'a section renamed: its follower follows, its old issues are not new',
    root: withRefAndError,
    patch: {
      ops: [
        { op: 'move', at: ['views', 'main', 'protein'], parent: ['views', 'main'], place: { last: true }, rename: 'p' },
      ],
    },
    introduces: [],
  },
  {
    name: 'a template made plain: its layer has no data',
    root: () => {
      const b = createSeqSpec();
      const src = b.data.source('genes', { format: 'bed', uri: 'genes.bed' });
      const g = b.data.sequence('g', { alphabet: 'dna', letters: 'ACGT' });
      b.views.view('main', (v) =>
        v
          .section('s', { of: g })
          .track('t', { for_each: { source: src, node_kind: 'ranges' } })
          .layer('l', { type: 'ranges' }),
      );
      return b.build();
    },
    patch: { ops: [{ op: 'unset', at: ['views', 'main', 's', 't'], key: 'for_each' }] },
    introduces: ['seq.lang.layer.data|/views/main/s/t/l|'],
  },
  {
    name: 'a map replaced with params that do not fit its method',
    root: p53Spec,
    patch: {
      ops: [
        {
          op: 'replace',
          at: ['data', 'p53-on-gene'],
          node: {
            kind: 'map',
            id: 'p53-on-gene',
            params: { from: '/data/p53', to: '/data/tp53', method: 'offset' },
          },
        },
      ],
    },
    introduces: ['seq.lang.map.params|/data/p53-on-gene|'],
  },
  {
    name: 'values no layer draws removed',
    root: p53Spec,
    patch: { ops: [{ op: 'remove', at: ['data', 'conservation'] }] },
    introduces: [],
  },
  {
    name: 'a glycan’s residues cut: its own links past them',
    root: moleculesSpec,
    patch: { ops: [{ op: 'set', at: ['data', 'glycan'], key: 'residues', value: ['NAG'] }] },
    introduces: [
      'seq.lang.link.bounds|/data/glycan/bonds|rows/a_start/0',
      'seq.lang.link.bounds|/data/glycan/bonds|rows/b_start/1',
    ],
  },
  {
    name: 'a unit’s range cut: its numbering, the complex’s links and the ranges on the unit past it',
    root: moleculesSpec,
    patch: { ops: [{ op: 'set', at: ['data', 'c', 'A'], key: 'range', value: [0, 8] }] },
    introduces: [
      'seq.lang.link.bounds|/data/c/links|rows/a_start/0',
      'seq.lang.numbering.length|/data/c/A/auth|labels',
      'seq.lang.ranges.bounds|/data/sites|rows/end/0',
      'seq.lang.unit.observed|/data/c/A|observed/start/0',
    ],
  },
  {
    name: 'a unit removed: a link naming it',
    root: moleculesSpec,
    patch: { ops: [{ op: 'remove', at: ['data', 'c', 'B'] }] },
    introduces: ['seq.lang.link.unit|/data/c/links|rows/b_unit/1'],
  },
  {
    name: 'a unit renamed: a link naming its old id',
    root: moleculesSpec,
    patch: { ops: [{ op: 'move', at: ['data', 'c', 'G'], parent: ['data', 'c'], place: { last: true }, rename: 'X' }] },
    introduces: ['seq.lang.link.unit|/data/c/links|rows/b_unit/0'],
  },
  {
    name: 'a sequence shortened: the units over it',
    root: moleculesSpec,
    patch: { ops: [{ op: 'set', at: ['data', 'fc'], key: 'letters', value: 'MKTAYIAKQR' }] },
    introduces: ['seq.lang.unit.range|/data/c/A|range', 'seq.lang.unit.range|/data/c/B|range'],
  },
  {
    name: 'a molecule’s atoms dropped: a block of a unit over it drawn as atoms',
    root: () => {
      const r = applyPatch(
        moleculesSpec(),
        { ops: [{ op: 'set', at: ['views', 'v', 'blocks', 'atp'], key: 'of', value: '/data/c/L' }] },
        options,
      );
      if (!r.ok) throw new Error('setup');
      return r.root;
    },
    patch: { ops: [{ op: 'unset', at: ['data', 'atp'], key: 'atoms' }] },
    introduces: ['seq.lang.molecule.form|/views/v/blocks/atp|form'],
  },
];

/** A ref through a ref with segments after it: `/refs/q` names `/refs/p/domains`, `/refs/p` the protein's section. */
function refThroughRef(): Node {
  const b = createSeqSpec({ title: 't' });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53_LETTERS });
  b.views.view('main', (v) => v.section('protein', { of: p53 }).track('domains').layer('l', { type: 'letters' }));
  b.refs.ref('p', { to: '/views/main/protein' });
  b.refs.ref('q', { to: '/refs/p/domains' });
  return b.build();
}

/** A layer drawing through two refs. */
function twoRefs(): Node {
  const b = createSeqSpec({ title: 't' });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53_LETTERS });
  const dom = b.data.ranges('domains', { target: p53, rows: { start: new Int32Array([0]), end: new Int32Array([4]) } });
  const d = b.refs.ref('d', { to: dom });
  const e = b.refs.ref('e', { to: String(d) });
  b.views.view('main', (v) =>
    v
      .section('s', { of: p53 })
      .track('t')
      .layer('l', { type: 'ranges', data: String(e) }),
  );
  return b.build();
}

/** A chain of eight refs, one more than... and a ref beside it. */
function longChain(): Node {
  const b = createSeqSpec({ title: 't' });
  let prev = String(b.data.sequence('p53', { alphabet: 'protein', letters: P53_LETTERS }));
  for (let i = 0; i < 8; i++) prev = String(b.refs.ref(`r${i}`, { to: prev }));
  b.refs.ref('x', { to: '/data/p53' });
  return b.build();
}

const Through: { name: string; root: () => Node; patch: Patch }[] = [
  {
    name: 'the node a ref names through another ref removed',
    root: refThroughRef,
    patch: { ops: [{ op: 'remove', at: ['views', 'main', 'protein', 'domains'] }] },
  },
  {
    name: 'the node a ref names through another ref renamed',
    root: refThroughRef,
    patch: {
      ops: [
        {
          op: 'move',
          at: ['views', 'main', 'protein', 'domains'],
          parent: ['views', 'main', 'protein'],
          place: { last: true },
          rename: 'dd',
        },
      ],
    },
  },
  {
    name: 'the data a layer draws through two refs replaced by another kind',
    root: twoRefs,
    patch: {
      ops: [
        {
          op: 'replace',
          at: ['data', 'domains'],
          node: { kind: 'sequence', id: 'domains', params: { alphabet: 'protein', letters: 'MEE' } },
        },
      ],
    },
  },
  {
    name: 'a chain of refs lengthened past the bound',
    root: longChain,
    patch: { ops: [{ op: 'set', at: ['refs', 'r0'], key: 'to', value: '/refs/x' }] },
  },
  {
    name: 'a sequence renamed',
    root: p53Spec,
    patch: { ops: [{ op: 'move', at: ['data', 'p53'], parent: ['data'], place: { last: true }, rename: 'q' }] },
  },
];

/** A spec with refs to nodes, refs through refs with segments, and layers drawing through refs. */
function chains(): Node {
  const b = createSeqSpec({ title: 'chains' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: P53_LETTERS });
  const r = (i: number) =>
    b.data.ranges(`r${i}`, { target: p, rows: { start: new Int32Array([0]), end: new Int32Array([4 + i]) } });
  const ranges = [r(0), r(1), r(2)];
  b.views.view('main', (v) => {
    const s = v.section('s', { of: p });
    ranges.forEach((x, i) => s.track(`t${i}`).layer('l', { type: 'ranges', data: i === 2 ? '/refs/d2' : x }));
    v.section('f', { of: p, follow: s });
  });
  b.refs.ref('sec', { to: '/views/main/s' });
  b.refs.ref('t1', { to: '/refs/sec/t1' });
  b.refs.ref('d2', { to: String(ranges[2]) });
  b.refs.ref('dd', { to: '/refs/d2' });
  b.refs.ref('set', { to: ['/refs/t1', '/data/p'] });
  return b.build();
}

/** A seeded generator (a 32-bit LCG). */
function seeded(seed: number): () => number {
  let s = seed >>> 0;
  return () => {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
    return s / 0x100000000;
  };
}

/** A random edit of the chains spec: removes, renames, retargets, letters, follows. */
function randomPatch(root: Node, rand: () => number): Patch {
  const pick = <T>(xs: readonly T[]) => xs[Math.floor(rand() * xs.length)];
  const kids = (at: string[]) => {
    let n: Node | undefined = root;
    for (const id of at) n = n?.children?.find((c) => c.id === id);
    return (n?.children ?? []).flatMap((c) => (c.id === undefined ? [] : [c.id]));
  };
  const node = () => {
    const level = pick([['data'], ['refs'], ['views', 'main'], ['views', 'main', 's']]);
    const id = pick(kids(level));
    return id === undefined ? undefined : [...level, id];
  };
  const at = node();
  if (!at) return { ops: [] };
  const anyAddress = () => {
    const a = node();
    return a ? `/${a.join('/')}` : '/data/p';
  };
  switch (pick(['remove', 'rename', 'to', 'letters', 'follow', 'data'])) {
    case 'remove':
      return { ops: [{ op: 'remove', at }] };
    case 'rename':
      return { ops: [{ op: 'move', at, parent: at.slice(0, -1), place: { last: true }, rename: `${at.at(-1)}x` }] };
    case 'to':
      return pick(kids(['refs'])) === undefined
        ? { ops: [] }
        : { ops: [{ op: 'set', at: ['refs', pick(kids(['refs']))], key: 'to', value: anyAddress() }] };
    case 'letters':
      return { ops: [{ op: 'set', at: ['data', 'p'], key: 'letters', value: P53_LETTERS.slice(0, pick([2, 5, 30])) }] };
    case 'follow':
      return { ops: [{ op: 'set', at: ['views', 'main', 's'], key: 'follow', value: '/views/main/f' }] };
    default: {
      const t = pick(kids(['views', 'main', 's']));
      return t === undefined
        ? { ops: [] }
        : { ops: [{ op: 'set', at: ['views', 'main', 's', t, 'l'], key: 'data', value: anyAddress() }] };
    }
  }
}

describe('checkSeqChanges through refs', () => {
  for (const c of Through)
    it(`reports what the whole check reports: ${c.name}`, () => {
      const { whole, change } = both(c.root(), c.patch);
      expect(change).toEqual(whole);
    });

  it('a ref through a removed ref prefix is reported unresolved', () => {
    const { change } = both(refThroughRef(), Through[0].patch);
    expect(change).toEqual(['seq.lang.ref.unresolved|/refs/q|to']);
  });

  it('holds over seeded random patches of refs in chains, along the revisions', () => {
    const rand = seeded(7);
    for (let run = 0; run < 40; run++) {
      let root = chains();
      for (let step = 0; step < 12; step++) {
        const p = randomPatch(root, rand);
        const r = applyPatch(root, p, options);
        if (!r.ok) continue;
        const { whole, change } = both(root, p);
        expect(change, JSON.stringify(p.ops)).toEqual(whole);
        root = r.root;
      }
    }
  });
});

describe('checkSeqChanges', () => {
  for (const c of Cases)
    it(`reports what the change introduced: ${c.name}`, () => {
      const { whole, change } = both(c.root(), c.patch);
      expect(change).toEqual(whole);
      expect(change).toEqual([...c.introduces].sort());
    });

  it('as the patch validator: refuses what a patch introduced, never an error already there', () => {
    const v = { ...options, allowDangling: false, allowRebind: false, validators: [SeqPatchValidator] };
    const shorten = applyPatch(
      p53Spec(),
      { ops: [{ op: 'set', at: ['data', 'p53'], key: 'letters', value: 'ME' }] },
      v,
    );
    expect(shorten.ok).toBe(false);
    expect(shorten.issues.map((i) => i.code)).toContain('seq.lang.ranges.bounds');
    // the broken ranges are there before: an edit beside them passes
    const beside = applyPatch(withRefAndError(), { ops: [{ op: 'set', at: [], key: 'title', value: 'x' }] }, v);
    expect(beside.ok).toBe(true);
  });

  it('keeps its index per revision: a chain of edits answers as the whole check does', () => {
    let root = p53Spec();
    const steps: Patch[] = [
      { ops: [{ op: 'set', at: ['data', 'p53'], key: 'letters', value: 'MEEPQSDPSV' }] },
      {
        ops: [
          {
            op: 'insert',
            parent: ['data'],
            place: { last: true },
            node: { kind: 'ranges', id: 'more', params: { target: '/data/p53', rows: { start: [2], end: [40] } } },
          },
        ],
      },
      { ops: [{ op: 'set', at: ['data', 'p53'], key: 'letters', value: P53_LETTERS.repeat(2) }] },
      { ops: [{ op: 'move', at: ['data', 'p53'], parent: ['data'], place: { after: null }, rename: 'tp53p' }] },
      { ops: [{ op: 'set', at: ['data', 'tp53p'], key: 'letters', value: 'M' }] },
    ];
    for (const p of steps) {
      // the next step's `before` is this one's `after`: its index was derived along the changes
      const r = both(root, p);
      expect(r.change).toEqual(r.whole);
      root = r.root;
    }
  });
});

/** A spec of `n` ranges on one protein, each drawn by a layer of its own track. */
function wide(n: number): Node {
  const b = createSeqSpec({ title: 'wide' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSV' });
  const rs = Array.from({ length: n }, (_, i) =>
    b.data.ranges(`r${i}`, { target: p, rows: { start: new Int32Array([0]), end: new Int32Array([4]) } }),
  );
  b.views.view('main', (v) => {
    const s = v.section('s', { of: p });
    rs.forEach((r, i) => s.track(`t${i}`).layer('l', { type: 'ranges', data: r }));
  });
  return b.build();
}

/** The least of a few timings (ms), each what `f` returns. */
function least(f: () => number, tries = 3): number {
  let best = Number.POSITIVE_INFINITY;
  for (let i = 0; i < tries; i++) best = Math.min(best, f());
  return best;
}
const clock = (f: () => unknown) => () => {
  const t0 = performance.now();
  f();
  return performance.now() - t0;
};
const hideOne = { ops: [{ op: 'set', at: ['views', 'main', 's', 't5', 'l'], key: 'hidden', value: true }] } as Patch;

describe('checkSeqChanges: its cost follows the change', () => {
  const measure = (n: number) => {
    const root = wide(n);
    const whole = least(clock(() => checkSeqLang(root)));
    // the first check after a load builds the index: a fresh tree each try
    const first = least(() => {
      const fresh = wide(n);
      const h = applyPatch(fresh, hideOne, options);
      if (!h.ok) throw new Error('refused');
      return clock(() => checkSeqChanges(fresh, h.root, h.changes))();
    });
    // the protein's letters: every ranges node names it (a fan-out)
    const r = applyPatch(
      root,
      { ops: [{ op: 'set', at: ['data', 'p'], key: 'letters', value: 'MEEPQSDPSVA' }] },
      options,
    );
    if (!r.ok) throw new Error('refused');
    checkSeqChanges(root, r.root, r.changes);
    const fanOut = least(clock(() => checkSeqChanges(root, r.root, r.changes)));
    // one layer hidden, the index kept
    const h = applyPatch(root, hideOne, options);
    if (!h.ok) throw new Error('refused');
    const one = least(clock(() => checkSeqChanges(root, h.root, h.changes)));
    return { whole, first, fanOut, one };
  };

  it('scales with the spec at first, with the change after: 1k and 16k layers', {
    tags: ['heavy'],
    timeout: 120_000,
  }, () => {
    const small = measure(1000);
    const large = measure(16000);
    // no worse than the whole check, at first and on a fan-out; a one-node edit a small part of it
    expect(large.first).toBeLessThan(2 * large.whole + 5);
    expect(large.fanOut).toBeLessThan(2 * large.whole + 5);
    expect(large.first / small.first).toBeLessThan(2 * Math.max(16, large.whole / small.whole));
    expect(large.one).toBeLessThan(large.whole / 10);
  });
});
