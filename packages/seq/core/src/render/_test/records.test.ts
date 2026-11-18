/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The committed draw-list records (the recording backend's format), compared as text: p53 as tracks, p53's protein on
// its gene with the gene's section following, and values below a pixel a position (a bar a device column, from a
// coarse level). A change to one is a change to look at: `SEQ_RECORDS=update` writes
// them again.

import { existsSync, readFileSync, writeFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';
import { recordDrawLists, renderSvg } from '../index.js';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { p53GeneSpec, p53TracksSpec, seeded } from './fixtures.js';

/** 50,000 values (a walk, some missing, negative and positive) on a system of their own: the store adds levels. */
function columnsSpec() {
  const b = createSeqSpec({ title: 'values in columns' });
  const n = 50_000;
  const r = seeded(11);
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(n / 4) });
  let x = 0;
  const values = Float32Array.from({ length: n }, (_, i) => {
    x += r() - 0.5;
    return i % 5_000 < 40 ? Number.NaN : x;
  });
  const v = b.data.values('v', { target: s, values });
  b.views.view('main', (view) =>
    view.section('s', { of: s }).track('v', { height: 24 }).layer({ type: 'values', data: v }),
  );
  return b.build();
}

function compare(name: string, record: string) {
  const file = fileURLToPath(new URL(`./records/${name}.txt`, import.meta.url));
  if (process.env.SEQ_RECORDS === 'update') writeFileSync(file, record);
  // a missing record is a failure: it is written on purpose, never by accident
  expect(existsSync(file), `no record ${name}: run with SEQ_RECORDS=update`).toBe(true);
  expect(record).toBe(readFileSync(file, 'utf8'));
}

describe('draw-list records', () => {
  it('p53 as tracks: letters, domains in lanes, conservation with missing values', async () => {
    const r = await renderSvg(p53TracksSpec(), {
      size: { width: 640, height: 116 },
      viewport: { start: 84, end: 112 },
      emphasis: {
        hover: { layer: '/views/main/protein/letters:0', section: 'protein', position: 90, system: '/data/p53' },
      },
    });
    expect(r.settled).toBe(true);
    expect(r.issues).toEqual([]);
    compare('p53-tracks', recordDrawLists(r.lists));
  });

  it('p53 on its gene: the gene’s section follows across an intron, the domain cut into its exons', async () => {
    const r = await renderSvg(p53GeneSpec(), { size: { width: 640, height: 144 }, viewport: { start: 150, end: 160 } });
    expect(r.settled).toBe(true);
    expect(r.issues).toEqual([]);
    compare('p53-gene', recordDrawLists(r.lists));
  });

  it('values in columns: a bar a device column, its extent from a coarse level, missing bins left out', async () => {
    // 164 device columns over 30,000 positions: 183 a column, so, at coarser 2, a level of 256 (up to two columns' positions) is read
    const r = await renderSvg(columnsSpec(), {
      size: { width: 300, height: 40 },
      viewport: { start: 10_000, end: 40_000 },
      config: { values: { coarser: 2 } },
    });
    expect(r.settled).toBe(true);
    expect(r.issues).toEqual([]);
    const record = recordDrawLists(r.lists, ['data']);
    // a level coarser than a column's positions: not what the default (coarser 1, a pixel's level) gives
    const fine = await renderSvg(columnsSpec(), {
      size: { width: 300, height: 40 },
      viewport: { start: 10_000, end: 40_000 },
    });
    expect(recordDrawLists(fine.lists, ['data'])).not.toBe(record);
    compare('values-columns', record);
  });
});
