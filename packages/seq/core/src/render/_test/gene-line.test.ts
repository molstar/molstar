/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ranges layer's gene line (`line: 'gene'`; _spec/README.md "Render", Layer types, The gene line): a transcript's
// pieces in one lane, laid by its whole extent; boxes full between `thick_start` and `thick_end`, thin outside; a
// hairline across the extent with chevrons by the strand (turned by a reversed axis); a hit on the hairline is the
// transcript's; a view inside an intron draws the hairline from the span rows alone (EC-229); rendered to SVG too.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { GffRuntime } from '../../runtimes/gff.js';
import { hitView, type RangesLay, recordDrawLists } from '../index.js';
import { renderSvg } from '../render-svg.js';
import { settle } from './helpers.js';
import { layerLay, readPaths, readRects } from './small-helpers.js';

const SIZE = { width: 1000, height: 120 } as const;

/** Two transcripts of one gene (+) and one on the other strand, 1-based inclusive as GFF3 writes them. */
export const GeneModelsGff = [
  '##gff-version 3',
  '##sequence-region c 1 2000',
  ...[
    ['c', '.', 'gene', '101', '900', '.', '+', '.', 'ID=g1;Name=ALPHA'],
    ['c', '.', 'mRNA', '101', '900', '.', '+', '.', 'ID=t1;Parent=g1;Name=ALPHA-201'],
    ['c', '.', 'exon', '101', '200', '.', '+', '.', 'Parent=t1'],
    ['c', '.', 'exon', '401', '500', '.', '+', '.', 'Parent=t1'],
    ['c', '.', 'exon', '801', '900', '.', '+', '.', 'Parent=t1'],
    ['c', '.', 'CDS', '151', '200', '.', '+', '0', 'Parent=t1'],
    ['c', '.', 'CDS', '401', '500', '.', '+', '1', 'Parent=t1'],
    ['c', '.', 'CDS', '801', '850', '.', '+', '0', 'Parent=t1'],
    ['c', '.', 'lnc_RNA', '101', '500', '.', '+', '.', 'ID=t2;Parent=g1;Name=ALPHA-202'],
    ['c', '.', 'exon', '101', '200', '.', '+', '.', 'Parent=t2'],
    ['c', '.', 'exon', '401', '500', '.', '+', '.', 'Parent=t2'],
    ['c', '.', 'mRNA', '1201', '1800', '.', '-', '.', 'ID=t3;Name=BETA-201'],
    ['c', '.', 'exon', '1201', '1300', '.', '-', '.', 'Parent=t3'],
    ['c', '.', 'exon', '1701', '1800', '.', '-', '.', 'Parent=t3'],
    ['c', '.', 'CDS', '1251', '1300', '.', '-', '0', 'Parent=t3'],
    ['c', '.', 'CDS', '1701', '1750', '.', '-', '0', 'Parent=t3'],
  ].map((l) => l.join('\t')),
].join('\n');

function spec(options: { line?: string; reversed?: boolean } = {}) {
  const b = createSeqSpec({ title: 'gene line' });
  const file = b.data.source('gff', { format: 'gff3', value: GeneModelsGff as never });
  const chr = b.data.sequence('c', { alphabet: 'dna', source: file, name: 'c' });
  const tx = b.data.ranges('tx', { target: chr, source: file });
  b.views.view('main', (v) => {
    v.section('c', { of: chr, ...(options.reversed ? { direction: 'reversed' } : {}) } as never)
      .track('tx', { height: 80 } as never)
      .layer('tx', { type: 'ranges', data: tx, line: options.line ?? 'gene' } as never);
  });
  return b.build();
}

async function frameOf(
  options: Parameters<typeof spec>[0] = {},
  viewport?: { start: number; end: number },
  size: { width: number; height: number } = SIZE,
) {
  const s = await settle(spec(options), [GffRuntime]);
  const f = await s.frame({ size, ...(viewport ? { viewport } : {}) });
  return { s, f, lay: layerLay<RangesLay>(f, 0) };
}

/** Each item key's lanes and marks (x0, x1, thin). */
function byItem(lay: RangesLay) {
  const out = new Map<string, { lanes: Set<number>; marks: { x0: number; x1: number; thin: number }[] }>();
  for (let j = 0; j < lay.m; j++) {
    const key = String(lay.keys[lay.mitem[j]]);
    let g = out.get(key);
    if (!g) out.set(key, (g = { lanes: new Set(), marks: [] }));
    g.lanes.add(lay.mlane[j]);
    g.marks.push({ x0: lay.mx0[j], x1: lay.mx1[j], thin: lay.mthin[j] });
  }
  for (const g of out.values()) g.marks.sort((a, b) => a.x0 - b.x0);
  return out;
}

describe('the gene line', () => {
  it('a transcript’s pieces in one lane; thick between thick_start and thick_end, thin outside; a non-coding one full', async () => {
    const { lay } = await frameOf();
    const items = byItem(lay);
    expect([...items.keys()].sort()).toEqual(['t1', 't2', 't3']);
    for (const g of items.values()) expect(g.lanes.size).toBe(1);
    // t1 and t2 overlap: two lanes; t3 beside t1 shares t1's
    expect(items.get('t1')!.lanes).not.toEqual(items.get('t2')!.lanes);
    const x = (p: number) => lay.axis.x(p);
    const t1 = items.get('t1')!.marks;
    // exon 1: a thin UTR 100–150, then thick 150–200; exon 2 thick; exon 3 thick 800–850, thin 850–900
    expect(t1.map((m) => m.thin)).toEqual([1, 0, 0, 0, 1]);
    expect(t1[0].x0).toBeCloseTo(x(100), 0);
    expect(t1[1].x0).toBeCloseTo(x(150), 0);
    expect(t1[4].x1).toBeCloseTo(x(900), 0);
    expect(items.get('t2')!.marks.map((m) => m.thin)).toEqual([0, 0]);
  });

  it('a hairline across each transcript’s extent, chevrons by its strand; a hit on it is the transcript’s', async () => {
    const { f, lay } = await frameOf();
    const lines = Array.from({ length: lay.cn }, (_, j) => ({
      key: String(lay.keys[lay.citem[j]]),
      x0: lay.cx0[j],
      x1: lay.cx1[j],
      dir: lay.cdir[j],
      lane: lay.clane[j],
    }));
    expect(lines.map((l) => [l.key, l.dir]).sort()).toEqual([
      ['t1', 1],
      ['t2', 1],
      ['t3', 2],
    ]);
    const t1 = lines.find((l) => l.key === 't1')!;
    expect(t1.x0).toBeCloseTo(lay.axis.x(100), 0);
    expect(t1.x1).toBeCloseTo(lay.axis.x(900), 0);
    // in the intron 200–400
    const y = lay.laneY + t1.lane * lay.laneHeight + lay.laneHeight / 2;
    expect(hitView(f.layout, lay.axis.x(300), y)?.element?.key).toBe('t1');
    // the chevrons are paths in the data pass, one every 12 px or so of the hairlines
    const width = lines.reduce((w, l) => w + l.x1 - l.x0, 0);
    expect(readPaths(f.lists).length).toBeGreaterThan(width / 12 - 2 * lines.length);
  });

  it('a reversed axis turns the chevrons', async () => {
    const { lay } = await frameOf({ reversed: true });
    const dirs = new Map(Array.from({ length: lay.cn }, (_, j) => [String(lay.keys[lay.citem[j]]), lay.cdir[j]]));
    expect(dirs.get('t1')).toBe(2);
    expect(dirs.get('t3')).toBe(1);
  });

  it('a view inside an intron draws each transcript’s hairline from its span rows, no box (EC-229)', async () => {
    const { lay } = await frameOf({}, { start: 250, end: 350 });
    expect(lay.m).toBe(0);
    const keys = Array.from({ length: lay.cn }, (_, j) => String(lay.keys[lay.citem[j]])).sort();
    expect(keys).toEqual(['t1', 't2']);
    for (let j = 0; j < lay.cn; j++) {
      expect(lay.cx0[j]).toBeCloseTo(lay.box.x, 0);
      expect(lay.cx1[j]).toBeCloseTo(lay.box.x + lay.box.width, 0);
    }
  });

  it('as blocks (no `line`) a span row draws no mark', async () => {
    const { lay } = await frameOf({ line: 'blocks' }, { start: 250, end: 350 });
    expect(lay.m).toBe(0);
  });

  it('rendered to SVG: boxes of two heights, the hairlines and their chevrons; no issue', async () => {
    const r = await renderSvg(spec(), { size: SIZE, runtimes: [GffRuntime] });
    expect(r.settled).toBe(true);
    expect(r.issues).toEqual([]);
    expect(r.svg).toContain('<path');
    const heights = new Set(
      readRects(r.lists)
        .filter((x) => x.w > 3 && x.h > 2)
        .map((x) => Math.round(x.h)),
    );
    expect(heights.size).toBe(2);
    expect(readPaths(r.lists).length).toBeGreaterThan(10);
  });

  it('inline items: the track’s lanes counted an item a lane by its extent, its pieces in it', async () => {
    const lanes: number[] = [];
    for (const line of ['blocks', 'gene']) {
      const b = createSeqSpec({ title: 'inline' });
      const g = b.data.sequence('g', { alphabet: 'dna', letters: 'ACGT'.repeat(25) });
      // t1's and t2's pieces never overlap, their extents do
      const tx = b.data.ranges('tx', {
        target: g,
        rows: { start: [0, 50, 20, 70], end: [10, 60, 30, 80], id: ['t1', 't1', 't2', 't2'] },
      });
      b.views.view('main', (v) => {
        v.section('g', { of: g })
          .track('tx')
          .layer('tx', { type: 'ranges', data: tx, ...(line === 'gene' ? { line } : {}) } as never);
      });
      const s = await settle(b.build());
      lanes.push(layerLay<RangesLay>(await s.frame({ size: SIZE }), 0).laneCount);
    }
    expect(lanes).toEqual([1, 2]);
  });

  it('a label is plain text left of its model, whole, never a box over its pieces: each transcript with room has one', async () => {
    // wide enough that the first model's 100 positions before it hold its label
    const { f, lay } = await frameOf({}, undefined, { width: 2400, height: 120 });
    const texts = recordDrawLists(f.lists, ['labels'])
      .split('\n')
      .filter((l) => l.startsWith('text ') && !l.startsWith('text "tx"'))
      .map((l) => {
        const m = /^text "([^"]*)" ([\d.-]+) ([\d.-]+) (\w+) .* width ([\d.]+)/.exec(l)!;
        return { text: m[1], x: Number(m[2]), anchor: m[4], width: Number(m[5]) };
      });
    expect(texts.map((t) => t.text).sort()).toEqual(['ALPHA-201', 'ALPHA-202', 'BETA-201']);
    for (const t of texts) {
      expect(t.anchor).toBe('start');
      expect(t.text.includes('…')).toBe(false);
      // it ends before its model starts on the screen
      const start = lay.axis.x(t.text.startsWith('BETA') ? 1200 : 100);
      expect(t.x + t.width).toBeLessThanOrEqual(start + 0.5);
      expect(t.x).toBeGreaterThanOrEqual(lay.box.x);
    }
    // no filled box in the labels pass (a label is text alone)
    expect(readRects(f.lists, 'labels')).toEqual([]);
    // a label is its item's
    const g = lay.gene;
    const k = Array.from({ length: g.ln }, (_, i) => i).find((i) => lay.keys[g.litem[i]] === 't3')!;
    expect(hitView(f.layout, g.lx[k] + 2, g.ly[k] - 3)?.element?.key).toBe('t3');
  });
});
