/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { applyPatch } from '@molstar/core/repr/patch';
import { parseSexpr, printSexpr } from '@molstar/core/repr/codecs/sexpr';
import { normalize } from '@molstar/core/repr/schema';
import { hashNode, type Node } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import { createSeqSpec } from '../builder.js';
import { checkSeqSpec } from '../check.js';
import { SeqLanguage } from '../language.js';
import { followRef, refAddresses } from '../refs.js';
import { P53_LETTERS } from './p53.js';

// User-defined refs: the `refs` container and the `ref` kind (one address or a set), checked by the language (a ref
// naming nothing or leading back to itself is a warning: the rest of the spec opens), followed as a rewrite of an
// address's leading ids, and rewritten by `repr/patch` when what they name moves.

function spec(): Node {
  const b = createSeqSpec({ title: 'refs' });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53_LETTERS });
  const domains = b.data.ranges('domains', {
    target: p53,
    rows: { start: new Int32Array([0, 17]), end: new Int32Array([14, 30]), id: ['TAD1', 'TAD2'] },
  });
  b.refs.ref('protein', { to: p53, label: 'p53' });
  b.refs.ref('tad', { to: `${domains}/(item :id TAD1)`, color: '#ff8800' });
  b.refs.ref('start', { to: '/refs/protein/(range :from 1 :to 10)' });
  b.refs.ref('both', { to: ['/refs/tad', `${p53}/(at :label 20)`] });
  b.refs.ref('alias', { to: '/refs/start' });
  b.views.view('main', (v) => {
    const s = v.section('protein', { of: p53 });
    s.track('domains').layer('ranges', { type: 'ranges', data: '/refs/domains-ref' });
  });
  b.refs.ref('domains-ref', { to: domains });
  return b.build();
}

describe('the refs container', () => {
  it('builds, prints and reads back, and checks without issues', () => {
    const s = spec();
    expect(s.children?.map((c) => c.id)).toEqual(['data', 'views', 'refs']);
    expect(checkSeqSpec(s)).toEqual([]);
    const text = printSexpr(s, { language: SeqLanguage });
    expect(text).toContain('(ref both :to [/refs/tad /data/p53/(at :label 20)])');
    const back = parseSexpr(text);
    expect(back.issues).toEqual([]);
    expect(hashNode(normalize(SeqLanguage, back.node as Node))).toBe(hashNode(s));
  });

  it('reports a ref naming nothing and a chain that comes back, as warnings', () => {
    const b = createSeqSpec();
    b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT' });
    b.refs.ref('gone', { to: ['/data/s', '/data/x/(range :from 1 :to 2)'] });
    b.refs.ref('a', { to: '/refs/b/(range :from 1 :to 2)' });
    b.refs.ref('b', { to: '/refs/a' });
    const issues = checkSeqSpec(b.build());
    expect(issues.map((i) => [i.code, i.severity, i.address, i.path?.slice(-2)])).toEqual([
      ['seq.lang.ref.unresolved', 'warning', '/refs/gone', ['to', 1]],
      ['seq.lang.ref.cycle', 'warning', '/refs/a', ['params', 'to']],
      ['seq.lang.ref.cycle', 'warning', '/refs/b', ['params', 'to']],
    ]);
  });

  it('refuses a container under another id, and a ref that is no address', () => {
    const bad: Node = {
      kind: 'spec',
      params: { version: SeqLanguage.version },
      children: [{ kind: 'refs', id: 'bookmarks' }],
    };
    expect(checkSeqSpec(bad).map((i) => i.code)).toEqual(['seq.lang.container.id']);
    const b = createSeqSpec();
    b.refs.ref('x', { to: 'data/s' });
    expect(checkSeqSpec(b.build()).map((i) => i.code)).toEqual(['repr.schema.ref.address']);
  });

  it('lets a layer draw through a ref, of a kind its type draws', () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT' });
    b.refs.ref('seq', { to: s });
    b.views.view('v', (v) => {
      v.section('sec', { of: s }).track('t').layer('l', { type: 'ranges', data: '/refs/seq' });
    });
    expect(checkSeqSpec(b.build()).map((i) => i.code)).toEqual(['seq.lang.layer.data']);
  });
});

describe('following refs', () => {
  const s = spec();

  it('rewrites the leading ids, through a chain, keeping the segments after them', () => {
    expect(followRef(s, '/refs/alias')).toEqual({
      addresses: ['/data/p53/(range :from 1 :to 10)'],
      via: ['/refs/alias', '/refs/start', '/refs/protein'],
    });
    expect(followRef(s, '/refs/protein/(at :label 3)')).toEqual({
      addresses: ['/data/p53/(at :label 3)'],
      via: ['/refs/protein'],
    });
    expect(followRef(s, '/data/p53')).toEqual({ addresses: ['/data/p53'], via: [] });
  });

  it('follows a set member by member', () => {
    expect(followRef(s, '/refs/both')).toEqual({
      addresses: ['/data/domains/(item :id TAD1)', '/data/p53/(at :label 20)'],
      via: ['/refs/both', '/refs/tad'],
    });
  });

  it('refuses a set followed by more segments, a cycle and a chain past the bound', () => {
    expect(followRef(s, '/refs/both/(range :from 1 :to 2)')).toMatchObject({ refused: 'set' });
    expect(followRef(s, '/refs/alias', 2)).toMatchObject({ refused: 'depth' });
    const b = createSeqSpec();
    b.refs.ref('a', { to: '/refs/b' });
    b.refs.ref('b', { to: '/refs/a' });
    expect(followRef(b.build(), '/refs/a')).toEqual({ refused: 'cycle', via: ['/refs/a', '/refs/b'] });
  });

  it('reads one address or a set', () => {
    expect(refAddresses({ kind: 'ref', params: { to: '/a' } })).toEqual(['/a']);
    expect(refAddresses({ kind: 'ref', params: { to: ['/a', '/b'] } })).toEqual(['/a', '/b']);
    expect(refAddresses({ kind: 'ref' })).toEqual([]);
  });
});

describe('refs under patches', () => {
  it('follow a rename of what they name, a form after its ids kept', () => {
    const s = spec();
    const r = applyPatch(
      s,
      { ops: [{ op: 'move', at: ['data', 'p53'], parent: ['data'], place: { last: true }, rename: 'tp53' }] },
      { language: SeqLanguage },
    );
    expect(r.ok).toBe(true);
    if (!r.ok) return;
    const refs = r.root.children!.find((c) => c.id === 'refs')!;
    expect(refs.children!.map((c) => c.params?.to)).toEqual([
      '/data/tp53',
      '/data/domains/(item :id TAD1)',
      '/refs/protein/(range :from 1 :to 10)',
      ['/refs/tad', '/data/tp53/(at :label 20)'],
      '/refs/start',
      '/data/domains',
    ]);
    expect(checkSeqSpec(r.root)).toEqual([]);
  });

  it('refuse a removal that would leave a ref dangling, unless allowed', () => {
    const s = spec();
    const remove = { ops: [{ op: 'remove' as const, at: ['data', 'domains'] }] };
    const refused = applyPatch(s, remove, { language: SeqLanguage });
    expect(refused.ok).toBe(false);
    expect(refused.issues.map((i) => i.code)).toContain('repr.patch.dangling');
  });
});
