/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A canonical transcript on the gene line (_spec/README.md "Render", Layer types, The gene line): an item whose
// `canonical` field has a value (GFF3's `Ensembl_canonical` and `MANE_Select` tags) is packed first, so it takes its
// gene's first lane, and its label is bold; with no tag the order and the labels are as they were.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { GffRuntime } from '../../runtimes/gff.js';
import type { RangesLay } from '../index.js';
import { recordDrawLists } from '../index.js';
import { settle } from './helpers.js';
import { layerLay } from './small-helpers.js';

/** One gene of three overlapping transcripts, the third the tagged one (when `tag`). */
function gff(tag: string): string {
  const t = (id: string, end: number, extra = '') =>
    ['c', '.', 'mRNA', '101', String(end), '.', '+', '.', `ID=${id};Name=${id.toUpperCase()}${extra}`].join('\t');
  const exon = (id: string, end: number) =>
    ['c', '.', 'exon', '101', String(end), '.', '+', '.', `Parent=${id}`].join('\t');
  return [
    '##gff-version 3',
    '##sequence-region c 1 2000',
    t('t1', 900),
    exon('t1', 900),
    t('t2', 800),
    exon('t2', 800),
    t('t3', 700, tag),
    exon('t3', 700),
  ].join('\n');
}

async function frame(text: string) {
  const b = createSeqSpec({ title: 'canonical' });
  const file = b.data.source('gff', { format: 'gff3', value: text as never });
  const chr = b.data.sequence('c', { alphabet: 'dna', source: file, name: 'c' });
  const tx = b.data.ranges('tx', { target: chr, source: file });
  b.views.view('main', (v) => {
    v.section('c', { of: chr })
      .track('tx', { height: 80 } as never)
      .layer('tx', { type: 'ranges', data: tx, line: 'gene' } as never);
  });
  const s = await settle(b.build(), [GffRuntime]);
  const f = await s.frame({ size: { width: 2400, height: 120 } });
  return { f, lay: layerLay<RangesLay>(f, 0) };
}

const laneOf = (lay: RangesLay, key: string) => {
  for (let j = 0; j < lay.m; j++) if (String(lay.keys[lay.mitem[j]]) === key) return lay.mlane[j];
  return -1;
};

/** The label texts the labels pass draws, each with its font. */
function labelFonts(f: Awaited<ReturnType<typeof frame>>['f']) {
  return recordDrawLists(f.lists, ['labels'])
    .split('\n')
    .filter((l) => l.startsWith('text ') && !l.startsWith('text "tx"'))
    .map((l) => /^text "([^"]*)" [\d.-]+ [\d.-]+ \w+ "([^"]*)"/.exec(l)!)
    .map((m) => [m[1], m[2]] as const);
}

describe('a canonical transcript on the gene line', () => {
  it('with no tag nothing is marked: the lanes follow the start, every label plain', async () => {
    const { f, lay } = await frame(gff(''));
    expect([laneOf(lay, 't1'), laneOf(lay, 't2'), laneOf(lay, 't3')]).toEqual([0, 1, 2]);
    expect(labelFonts(f).length).toBe(3);
    expect(labelFonts(f).every(([, font]) => !font.startsWith('bold'))).toBe(true);
  });

  it('a tagged transcript takes the first lane of its gene’s stack and its label is bold, the others plain', async () => {
    const { f, lay } = await frame(gff(';tag=basic,Ensembl_canonical'));
    expect(laneOf(lay, 't3')).toBe(0);
    expect([laneOf(lay, 't1'), laneOf(lay, 't2')]).toEqual([1, 2]);
    const fonts = new Map(labelFonts(f));
    expect(fonts.get('T3')).toMatch(/^bold /);
    expect(fonts.get('T1')).not.toMatch(/^bold /);
    expect(fonts.get('T2')).not.toMatch(/^bold /);
  });
});
