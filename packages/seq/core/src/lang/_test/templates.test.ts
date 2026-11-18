/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { createSeqSpec } from '../builder.js';

// Track templates and item maps in the language: a template's layer without data draws each instance's node; the
// template's source names a source; a ranges node's `map` makes each item a coding sequence.

describe('track templates', () => {
  function spec(source: string, data = false) {
    const b = createSeqSpec();
    const file = b.data.source('genes', { format: 'seq', uri: 'genes.seq' });
    const gene = b.data.sequence('gene', { alphabet: 'dna', source: file });
    b.data.ranges('cds', { target: gene, rows: { start: [0, 20], end: [9, 29], id: ['T1', 'T1'] }, map: 'codon' });
    b.views.view('main', (v) => {
      const s = v.section('gene', { of: gene });
      s.track('transcripts', {
        for_each: { source, node_kind: 'ranges', depth: 0, order: '-key', filter: { map: 'codon' } },
      }).layer({ type: 'ranges', ...(data ? { data: '/data/cds' } : {}) });
      s.track('plain').layer({ type: 'ranges' });
    });
    return b;
  }

  it('a template layer draws its node without data; a plain layer still needs data', () => {
    // the plain track's layer alone: the second track of the section
    expect(
      spec('/data/genes')
        .check()
        .map((i) => [i.code, i.path?.slice(0, 8)]),
    ).toEqual([['seq.lang.layer.data', ['children', 1, 'children', 0, 'children', 0, 'children', 1]]]);
    expect(spec('/data/genes', true).check()).toHaveLength(1);
  });

  it('the template’s source names a source', () => {
    expect(
      spec('/data/gene')
        .check()
        .map((i) => [i.code, i.path?.slice(-2)]),
    ).toContainEqual(['seq.lang.ref.target', ['for_each', 'source']]);
  });

  it('a node kind, an order and a filter are checked by the description', () => {
    const b = createSeqSpec();
    b.data.source('genes', { format: 'seq', uri: 'genes.seq' });
    b.views.view('main', (v) =>
      // @ts-expect-error: not a kind a template makes tracks for
      v.section('s', { of: '/data/x' }).track('t', { for_each: { source: '/data/genes', node_kind: 'view' } }),
    );
    expect(b.check().map((i) => i.code)).toContain('repr.schema.param.enum');
  });

  it('expands annotations only: a template over source nodes is refused by the description', () => {
    const b = createSeqSpec();
    b.data.source('genes', { format: 'seq', uri: 'genes.seq' });
    b.views.view('main', (v) =>
      // @ts-expect-error: a source node has nothing a layer draws
      v.section('s', { of: '/data/x' }).track('t', { for_each: { source: '/data/genes', node_kind: 'source' } }),
    );
    expect(b.check().map((i) => i.code)).toContain('repr.schema.param.enum');
  });
});
