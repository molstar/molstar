/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The connector of a ranges layer (_spec/README.md "Render",
// Layer types and Configuration): an item projected into pieces with shown positions between them is joined by a
// hairline in the item's colour; across a collapsed break, or in one run, nothing; `connector: none` draws none, the
// preset gives its width, ink and alpha; a hit on it is the item's and the hover outlines it.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { hitView, type RangesLay } from '../index.js';
import { settle } from './helpers.js';
import { layerLay, readRects } from './small-helpers.js';

const SIZE = { width: 800, height: 60 } as const;

/** A 120 nt gene with exons 12–46 and 76–105 coding a 20-residue protein; one domain item on the protein. */
function spec(
  rows: { start: number[]; end: number[]; id: string[] },
  options: {
    pieces?: 'whole' | 'coding';
    connector?: 'line' | 'none';
    direction?: string;
    strand?: '+' | '-';
    pad?: number;
  } = {},
): Node {
  const b = createSeqSpec({ title: 'connector' });
  const gene = b.data.sequence('gene', { alphabet: 'dna', letters: 'ACGT'.repeat(30) });
  const protein = b.data.sequence('protein', { alphabet: 'protein', letters: 'M'.repeat(20) });
  b.data.map('map', {
    from: protein,
    to: gene,
    method: 'codon',
    exons: { start: [12, 76], end: [46, 105] },
    strand: '+',
  });
  const domain = b.data.ranges('domain', { target: protein, rows });
  b.views.view('main', (v) => {
    v.section('gene', {
      of: gene,
      direction: (options.direction ?? 'genomic') as never,
      pieces: options.pieces ?? 'whole',
      ...(options.pad ? { pad: options.pad } : {}),
    })
      .track('domain')
      .layer('domain', {
        type: 'ranges',
        data: domain,
        ...(options.connector ? { connector: options.connector } : {}),
      } as never);
  });
  return b.build();
}

const SPLIT = { start: [4], end: [14], id: ['D'] };

/** The domain layer's marks and connectors of a frame. */
async function frameOf(node: Node, config = {}) {
  const s = await settle(node, [], config);
  const f = await s.frame({ size: SIZE });
  const lay = layerLay<RangesLay>(f, 0);
  return { s, f, lay };
}

const marksOf = (lay: RangesLay) =>
  Array.from({ length: lay.m }, (_, j) => ({ x0: lay.mx0[j], x1: lay.mx1[j] })).sort((a, b) => a.x0 - b.x0);

describe('the ranges connector', () => {
  it('joins the two pieces of a domain split by a shown intron: a hairline from one piece’s end to the next’s start', async () => {
    const { f, lay } = await frameOf(spec(SPLIT));
    const marks = marksOf(lay);
    expect(marks.length).toBe(2);
    expect(lay.cn).toBe(1);
    expect(lay.cx0[0]).toBeCloseTo(marks[0].x1, 0);
    expect(lay.cx1[0]).toBeCloseTo(marks[1].x0, 0);
    // drawn in the data pass as a line in the item's colour at 0.6, centred in the lane, between the two marks
    const line = readRects(f.lists).filter((r) => r.h === 1 && r.w > 10);
    expect(line.length).toBe(1);
    expect(line[0].x).toBeCloseTo(lay.cx0[0], 3);
    expect(line[0].w).toBeCloseTo(lay.cx1[0] - lay.cx0[0], 3);
    expect(Math.abs(line[0].y + 0.5 - (lay.laneY + lay.laneHeight / 2))).toBeLessThanOrEqual(0.5);
    expect(line[0].fill & 0xff).toBe(Math.round(255 * 0.6));
    expect(line[0].fill >>> 8).toBe(lay.fills[0] >>> 8);
  });

  it('the pieces end where the domain’s residues’ codons end (residues 5–14: coding 12–42)', async () => {
    const { lay } = await frameOf(spec(SPLIT));
    const marks = marksOf(lay);
    const x = (p: number) => lay.axis.x(p);
    // 12 + 12 = 24 to the exon's end 46; the next exon's start 76 to 76 + (42 − 34) = 84
    expect(marks[0].x0).toBeCloseTo(x(24), 0);
    expect(marks[0].x1).toBeCloseTo(x(46), 0);
    expect(marks[1].x0).toBeCloseTo(x(76), 0);
    expect(marks[1].x1).toBeCloseTo(x(84), 0);
  });

  it('a hit between the pieces is the item’s, and its hover outlines the connector with it', async () => {
    const { s, f, lay } = await frameOf(spec(SPLIT));
    const mid = (lay.cx0[0] + lay.cx1[0]) / 2;
    const y = lay.laneY + lay.laneHeight / 2;
    const hit = hitView(f.layout, mid, y);
    expect(hit?.element?.key).toBe('D');
    // not on an empty lane elsewhere
    expect(hitView(f.layout, lay.box.x + 2, y)?.element).toBeUndefined();
    const hovered = s.renderer.frame({ size: SIZE, emphasis: { hover: hit } });
    const line = readRects(hovered.lists, 'hover').filter((r) => r.h <= 2 && r.w > 10);
    expect(line.length).toBe(1);
    expect(line[0].x).toBeCloseTo(lay.cx0[0], 3);
  });

  it('draws nothing across a collapsed break (zero units wide), and nothing for an item in one run', async () => {
    const collapsed = await frameOf(spec(SPLIT, { pieces: 'coding' }));
    expect(marksOf(collapsed.lay).length).toBe(2);
    expect(collapsed.lay.cn).toBe(0);
    expect(readRects(collapsed.f.lists).filter((r) => r.h === 1 && r.w > 10)).toEqual([]);
    const one = await frameOf(spec({ start: [0], end: [3], id: ['A'] }));
    expect(marksOf(one.lay).length).toBe(1);
    expect(one.lay.cn).toBe(0);
  });

  it('`connector: none` draws none, and a hit between the pieces is the lane’s alone', async () => {
    const { f, lay } = await frameOf(spec(SPLIT, { connector: 'none' }));
    expect(marksOf(lay).length).toBe(2);
    expect(lay.cn).toBe(0);
    expect(readRects(f.lists).filter((r) => r.h === 1 && r.w > 10)).toEqual([]);
    const marks = marksOf(lay);
    expect(hitView(f.layout, (marks[0].x1 + marks[1].x0) / 2, lay.laneY + lay.laneHeight / 2)?.element).toBeUndefined();
  });

  it('the preset gives its width, ink and alpha', async () => {
    const { f, lay } = await frameOf(spec(SPLIT), { ranges: { connector: { width: 2, ink: 'muted', alpha: 1 } } });
    const line = readRects(f.lists).filter((r) => r.h === 2 && r.w > 10);
    expect(line.length).toBe(1);
    expect(line[0].fill).toBe(lay.tokens.colors['msui-text-muted']);
    // alpha 0 draws none
    const none = await frameOf(spec(SPLIT), { ranges: { connector: { alpha: 0 } } });
    expect(readRects(none.f.lists).filter((r) => r.h === 1 && r.w > 10)).toEqual([]);
  });

  it('on a reversed axis and over a minus-strand map the connector still joins the gap between the pieces', async () => {
    for (const options of [
      { direction: 'reversed' },
      { strand: '-' as const },
      { direction: 'reversed', strand: '-' as const },
    ]) {
      const { lay } = await frameOf(spec(SPLIT, options));
      const marks = marksOf(lay);
      expect(marks.length, JSON.stringify(options)).toBe(2);
      expect(lay.cn, JSON.stringify(options)).toBe(1);
      expect(lay.cx0[0]).toBeCloseTo(marks[0].x1, 0);
      expect(lay.cx1[0]).toBeCloseTo(marks[1].x0, 0);
    }
  });

  it('with `pad` the connector runs through the padded gap across the collapsed break', async () => {
    const { lay } = await frameOf(spec(SPLIT, { pieces: 'coding', pad: 4 }));
    const marks = marksOf(lay);
    expect(marks.length).toBe(2);
    expect(lay.cn).toBe(1);
    // from the first piece's end through its padding, the break and the next piece's padding to its start
    expect(lay.cx0[0]).toBeCloseTo(marks[0].x1, 0);
    expect(lay.cx1[0]).toBeCloseTo(marks[1].x0, 0);
    expect(lay.cx1[0] - lay.cx0[0]).toBeGreaterThan(8 * lay.axis.scale - 1);
  });

  it('a selected item’s connector is drawn in the accent, an unselected one’s is not', async () => {
    const { s, lay } = await frameOf(spec(SPLIT));
    const selected = s.renderer.frame({
      size: SIZE,
      emphasis: { elements: { keys: [{ owner: lay.owner, key: 'D' }] } },
    });
    const accent = lay.tokens.colors['msui-accent'];
    const line = (f: typeof selected) => readRects(f.lists, 'selection').filter((r) => r.h <= 2 && r.w > 10);
    expect(line(selected).length).toBe(1);
    expect(line(selected)[0].fill).toBe(accent);
    expect(line(selected)[0].x).toBeCloseTo(lay.cx0[0], 3);
    const other = s.renderer.frame({
      size: SIZE,
      emphasis: { elements: { keys: [{ owner: lay.owner, key: 'other' }] } },
    });
    expect(line(other)).toEqual([]);
  });
});
