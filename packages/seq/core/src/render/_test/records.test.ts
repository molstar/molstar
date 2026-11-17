/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The committed draw-list records (the recording backend's format), compared as text: p53 as tracks, and p53's protein
// on its gene with the gene's section following. A change to one is a change to look at: `SEQ_RECORDS=update` writes
// them again.

import { existsSync, readFileSync, writeFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';
import { recordDrawLists, renderSvg } from '../index.js';
import { p53GeneSpec, p53TracksSpec } from './fixtures.js';

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
});
