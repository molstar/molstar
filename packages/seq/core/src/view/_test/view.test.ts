/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The resolved view of the p53 example and of small specs: sections and axes, tracks, layers with their data, paths,
// themes and fallbacks, follow links with their scales, and what does not resolve.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { OpenLength, resolveSpec } from '@molstar/seqx-core/resolve/index';
import { describe, expect, it } from 'vitest';
import { EXONS, P53_LETTERS, p53Spec } from '../../lang/_test/p53.js';
import { followRange, resolveView, resolveViews } from '../view.js';

const model = resolveSpec(p53Spec());
const main = resolveView(model, '/views/main')!;
const codes = (issues: readonly { code: string }[]) => issues.map((i) => i.code);

describe('resolved view: the p53 example', () => {
  it('has its sections, each an axis on a system', () => {
    expect(main.address).toBe('/views/main');
    expect(main.layout).toBe('tracks');
    expect(main.issues).toEqual([]);
    expect(main.sections.map((s) => [s.address, s.axis])).toEqual([
      ['/views/main/protein', { system: '/data/p53', length: 30, cyclic: false, ranges: [0, 30] }],
      ['/views/main/gene', { system: '/data/tp53', cyclic: false, ranges: [0, OpenLength] }],
    ]);
    expect(resolveViews(model)).toEqual([main]);
    expect(resolveView(model, '/views/main')).toBe(main);
    expect(resolveView(model, '/views/none')).toBeUndefined();
  });

  it('has its tracks and layers, each layer’s data on its section’s system', () => {
    const [letters, domains] = main.sections[0].tracks;
    expect([letters.address, letters.label, domains.label]).toEqual([
      '/views/main/protein/letters',
      'letters',
      'Domains',
    ]);
    expect(letters.layers).toEqual([
      {
        index: 0,
        type: 'letters',
        fallbacks: ['letters'],
        params: {},
        theme: {},
        data: { kind: 'address', address: '/data/p53', system: '/data/p53', shape: 'letters' },
        placed: true,
      },
    ]);
    expect(domains.layers[0].data).toEqual({
      kind: 'address',
      address: '/data/domains',
      system: '/data/p53',
      shape: 'ranges',
    });
    expect(domains.layers[0].path).toBeUndefined();
    expect(main.sections[1].tracks).toEqual([]);
  });

  it('resolves the gene’s follow to the codon path, three nucleotides a residue', () => {
    const follow = main.sections[1].follow!;
    expect([follow.leader, follow.scale, follow.path.from, follow.path.to, follow.path.step]).toEqual([
      '/views/main/protein',
      3,
      '/data/p53',
      '/data/tp53',
      1,
    ]);
    expect(followRange(main.sections[1], 0, 1)).toEqual({ start: 197, end: 200 });
    expect(followRange(main.sections[1], 0, 30)).toEqual({ start: 100, end: 200 });
    expect(followRange(main.sections[0], 0, 1)).toBeUndefined();
  });
});

describe('resolved view: layers on other systems, queries, themes', () => {
  /** p53 and its gene inline; a view with the protein's domains drawn on the gene and on a lone sequence. */
  function spec(): Node {
    const b = createSeqSpec();
    const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53_LETTERS });
    const gene = b.data.sequence('tp53', { alphabet: 'dna', letters: 'A'.repeat(200) });
    const lone = b.data.sequence('lone', { alphabet: 'protein', letters: 'MKW' });
    b.data.map('p53-on-gene', { from: p53, to: gene, method: 'codon', exons: EXONS, strand: '+' });
    const domains = b.data.ranges('domains', {
      target: p53,
      rows: { start: new Int32Array([0]), end: new Int32Array([14]) },
    });
    b.views.view('v', (v) => {
      const g = v.section('gene', { of: gene, range: [90, 200] });
      g.track('t').layer('d', { type: 'ranges', data: domains, color: { kind: 'uniform', value: 'red' } });
      g.track('q').layer({
        type: 'ranges',
        data: '/data/domains/(items :overlapping /data/p53/(range :from 1 :to 5))',
      });
      const l = v.section('lone', { of: lone });
      l.track('t').layer({ type: 'ranges', data: domains });
      v.section('follower', { of: lone, follow: g });
    });
    return b.build();
  }
  const m = resolveSpec(spec());
  const v = resolveView(m, '/views/v')!;

  it('carries a layer’s data along the static path to its section, with its theme as a value', () => {
    const d = v.sections[0].tracks[0].layers[0];
    expect([d.address, d.placed, d.path?.step, d.fallbacks]).toEqual(['/views/v/gene/t/d', true, 1, ['ranges']]);
    expect(d.theme).toEqual({ color: { kind: 'uniform', value: 'red' } });
    expect(d.params).toEqual({ color: { kind: 'uniform', value: 'red' } });
  });

  it('keeps a query as a query, on the system its first node is on', () => {
    const q = v.sections[0].tracks[1].layers[0];
    expect(q.data).toEqual({
      kind: 'query',
      address: '/data/domains/(items :overlapping /data/p53/(range :from 1 :to 5))',
      system: '/data/p53',
      shape: 'ranges',
    });
    expect(q.placed && q.path?.to).toBe('/data/tp53');
  });

  it('reports data that does not reach its section and a follow with no path, and draws neither', () => {
    const lone = v.sections[1].tracks[0].layers[0];
    expect([lone.placed, lone.path]).toEqual([false, undefined]);
    expect(v.sections[2].follow).toBeUndefined();
    expect(codes(v.issues)).toEqual(['seq.view.layer.unplaced', 'seq.view.follow.unplaced']);
    expect(v.issues.map((i) => i.address)).toEqual(['/views/v/lone/t', '/views/v/follower']);
    expect(v.issues[0].hint).toContain('No path of maps joins the systems');
  });

  it('keeps a section’s range', () => {
    expect(v.sections[0].axis).toEqual({ system: '/data/tp53', length: 200, cyclic: false, ranges: [90, 200] });
  });

  it('reports a layer whose data names nothing (an unchecked spec)', () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'protein', letters: 'MK' });
    b.views.view('v', (view) => view.section('s', { of: s }).track('t').layer({ type: 'ranges', data: '/data/none' }));
    const r = resolveView(resolveSpec(b.build()), '/views/v')!;
    expect(codes(r.issues)).toEqual(['seq.view.layer.data']);
    expect(r.sections[0].tracks[0].layers[0].placed).toBe(false);
  });
});

describe('resolved view: cyclic sections', () => {
  const b = createSeqSpec();
  const a = b.data.sequence('a', { alphabet: 'dna', letters: 'ACGTAC' });
  const pl = b.data.sequence('plasmid', { alphabet: 'dna', letters: 'ACGTACGTAC', cyclic: true });
  b.data.map('m', { from: a, to: pl, method: 'offset', by: 7 });
  b.views.view('v', (v) => {
    const lead = v.section('a', { of: a });
    v.section('pl', { of: pl, range: [8, 2], follow: lead });
  });
  const v = resolveView(resolveSpec(b.build()), '/views/v')!;

  it('a range that wraps the origin is two ranges (EC-9, EC-10)', () => {
    expect(v.sections[1].axis.ranges).toEqual([0, 2, 8, 10]);
  });

  it('a follower on a cyclic system takes the shortest arc (EC-12)', () => {
    expect(v.sections[1].follow?.scale).toBe(1);
    expect(followRange(v.sections[1], 1, 5)).toEqual({ start: 8, end: 2 });
    expect(followRange(v.sections[1], 0, 2)).toEqual({ start: 7, end: 9 });
  });
});
