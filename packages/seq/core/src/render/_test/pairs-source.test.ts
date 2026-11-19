/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { inlineSeq, SeqRuntime, seqToValue } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { type PairsLay, recordDrawLists, renderSvg, viewLegends } from '../index.js';
import { settle } from './helpers.js';
import { letters, pairRows } from './pairs-fixtures.js';

// Pairs over a sequence whose letters arrive from a source (_spec/README.md "Render", Pairs): the resolved axis has no
// pieces until the source gave its length, the layout places it from the coverage, and the arcs' ends are read on that
// axis, in the canvas's lists and in the export alike; the arcs' legend lists the pair kinds, titled by their track.

const N = 40;
const Pairs = [
  [0, 39],
  [1, 38],
  [5, 20],
] as const;

function spec() {
  const b = createSeqSpec();
  const src = b.data.source('file', {
    format: 'seq',
    value: seqToValue(
      inlineSeq(
        { name: 's', length: N, alphabet: 'rna' },
        { letters: { form: 'letters', start: 0, letters: letters(N) } },
      ),
    ),
  });
  const s = b.data.sequence('s', { alphabet: 'rna', source: src });
  const p = b.data.pairs('p', { a: s, rows: pairRows(Pairs, ['canonical', 'canonical', 'wobble']) });
  b.views.view('v', (v) => {
    const sec = v.section('s', { of: s });
    sec.track('pairs', { label: 'Base pairs' }).layer('arcs', { type: 'pairs', data: p });
    sec.track('letters').layer('letters', { type: 'letters' });
  });
  return b.build();
}

const paths = (lists: Parameters<typeof recordDrawLists>[0]) =>
  recordDrawLists(lists, ['data'])
    .split('\n')
    .filter((l) => l.startsWith('path'));

describe('pairs over a sequence read from a source', () => {
  it('the canvas places every pair, and the export draws the same arcs without an issue', async () => {
    const s = await settle(spec(), [SeqRuntime]);
    const f = await s.frame({ size: { width: 640, height: 200 }, dpr: 1, settled: true });
    const lay = f.layout.sections[0].sectionLays[0].lay as PairsLay;
    expect(lay.n).toBe(3);
    expect(paths(f.lists)).toHaveLength(3);
    const r = await renderSvg(spec(), { size: { width: 640, height: 200 }, runtimes: [SeqRuntime] });
    expect(r.issues).toEqual([]);
    expect(paths(r.lists)).toEqual(paths(f.lists));
  });

  it('the arcs’ legend lists the pair kinds, titled by the track that holds them', async () => {
    const s = await settle(spec(), [SeqRuntime]);
    const f = await s.frame({ size: { width: 640, height: 200 }, dpr: 1, settled: true });
    const legend = viewLegends(f.layout).legends.find((l) => l.layers.includes('/views/v/s/pairs/arcs'));
    expect(legend?.title).toBe('Base pairs');
    expect(legend?.entries.kind === 'swatches' && legend.entries.items.map((i) => i.label)).toEqual([
      'canonical',
      'wobble',
    ]);
  });
});
