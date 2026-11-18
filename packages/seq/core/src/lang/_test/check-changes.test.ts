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
];

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
