/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Ranges drawn as lines and with arrowheads (_spec/README.md "Render", Layer types; EC-128, EC-232): a style anchored to
// the content (the dash phase moves with the data), an invalid style drawn solid and said once, an arrowhead on the piece
// holding the end it points to and never on a cut end.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { hitView, type RangesLay } from '../index.js';
import { settle } from './helpers.js';
import { layerLay, readPaths, readRects, rng } from './small-helpers.js';

const SIZE = { width: 800, height: 40 } as const;

interface Rows {
  start: number[];
  end: number[];
  id: string[];
  category?: string[];
}

function spec(
  rows: Rows,
  layer: { line?: string; arrow?: string },
  options: {
    length?: number;
    cyclic?: boolean;
    section?: Record<string, unknown>;
    pieces?: { start: number[]; end: number[] };
  } = {},
): Node {
  const length = options.length ?? 400;
  const b = createSeqSpec({ title: 'ranges lines' });
  const s = b.data.sequence('s', {
    alphabet: 'dna',
    letters: 'ACGT'.repeat(Math.ceil(length / 4)).slice(0, length),
    ...(options.cyclic ? { cyclic: true } : {}),
  });
  const items = b.data.ranges('items', { target: s, rows: rows as never });
  const pieces = options.pieces ? b.data.ranges('pieces', { target: s, rows: options.pieces }) : undefined;
  b.views.view('main', (v) => {
    v.section('s', { of: s, ...(pieces ? { pieces } : {}), ...options.section })
      .track('items')
      .layer('items', { type: 'ranges', data: items, ...(layer as object) } as never);
  });
  return b.build();
}

const ONE: Rows = { start: [40], end: [300], id: ['a'] };

/** Each mark's item key and tips (bit 1: on the left, bit 2: on the right). */
function marks(lay: RangesLay) {
  const out: { key: string | number; tips: number; x0: number; x1: number }[] = [];
  for (let j = 0; j < lay.m; j++)
    out.push({ key: lay.keys[lay.mitem[j]], tips: lay.mtip[j], x0: lay.mx0[j], x1: lay.mx1[j] });
  return out;
}

describe('ranges: line styles', () => {
  it('each style draws its mark: a rect, a dashed path, a wave, two rects, an outline, a translucent block', async () => {
    const draw = async (line: string) => {
      const f = await (await settle(spec(ONE, { line }))).frame({ size: SIZE });
      return { rects: readRects(f.lists), paths: readPaths(f.lists), lay: layerLay<RangesLay>(f, 0), f };
    };
    const solid = await draw('solid');
    expect(solid.rects.length).toBe(1);
    expect(solid.rects[0].h).toBe(2);
    expect(solid.rects[0].y + 1).toBeCloseTo(solid.lay.box.y + solid.lay.laneHeight / 2, 3);
    const dashed = await draw('dashed');
    expect(dashed.paths.map((p) => [p.dash, p.strokeWidth])).toEqual([[6, 2]]);
    const dotted = await draw('dotted');
    expect(dotted.paths.map((p) => p.dash)).toEqual([2]);
    const wavy = await draw('wavy');
    expect(wavy.paths.length).toBe(1);
    const w = wavy.paths[0].xy;
    expect(w.length / 2).toBeGreaterThan(40);
    // a corner every 3 px, above and below the lane's middle in turn
    const cy = solid.rects[0].y + 1;
    for (let i = 0; i < 20; i++) expect(w[2 * i + 1]).toBeCloseTo(cy + (i % 2 ? 2 : -2), 3);
    const double = await draw('double');
    expect(double.rects.map((r) => r.h)).toEqual([1, 1]);
    const box = await draw('box');
    expect(box.rects.length).toBe(1);
    expect(box.rects[0].stroke).not.toBe(0);
    expect(box.rects[0].fill).toBe(0);
    const highlight = await draw('highlight');
    expect(highlight.rects[0].fill & 0xff).toBe(89);
    // none of them is a block: the lane's height is the hit box all the same
    const h = hitView(solid.f.layout, 200, solid.lay.box.y + solid.lay.laneHeight / 2);
    expect(h?.element?.key).toBe('a');
  });

  it('an invalid style is drawn solid and said once a layout', async () => {
    const s = await settle(spec(ONE, { line: 'zigzag' }));
    for (const w of [30, 60, 90]) {
      const f = await s.frame({ size: { width: 800, height: 40 }, viewport: { start: 0, end: 400 - w } });
      expect(f.layout.issues.filter((i) => i.code === 'seq.render.line-style').length).toBe(1);
      const rects = readRects(f.lists);
      expect(rects.length).toBe(1);
      expect(rects[0].h).toBe(2);
    }
    // a style given as the empty string is invalid too; none given: blocks, nothing said
    expect((await (await settle(spec(ONE, { line: '' }))).frame({ size: SIZE })).layout.issues.length).toBe(1);
    const blocks = await (await settle(spec(ONE, {}))).frame({ size: SIZE });
    expect(blocks.layout.issues).toEqual([]);
    expect(readRects(blocks.lists)[0].h).toBeGreaterThan(2);
  });

  for (const direction of ['genomic', 'reversed']) {
    it(`the dash and the wave run from the content’s start (${direction}): their phase moves with it over any pan`, async () => {
      const r = rng(7);
      for (const [line, period] of [
        ['dashed', 12],
        ['dotted', 4],
        ['wavy', 6],
      ] as const) {
        const s = await settle(spec(ONE, { line }, { section: { direction } }));
        let off = 0;
        for (let i = 0; i < 25; i++) {
          // views that start inside, before and past the item's start, and end inside or after it
          const w = 10 + r() * 150;
          const start = r() * (400 - w);
          const f = await s.frame({ size: SIZE, viewport: { start, end: start + w } });
          const lay = layerLay<RangesLay>(f, 0);
          const p = readPaths(f.lists)[0];
          if (!p) continue;
          // the content's start: the item's first position's edge (the right one reversed), unclipped
          const anchor = lay.axis.x(40);
          const xs = p.xy[0];
          const turns = (direction === 'genomic' ? xs - anchor : anchor - xs) / period;
          expect(Math.abs(turns - Math.round(turns))).toBeLessThan(1e-3);
          expect(Math.round(turns)).toBeGreaterThanOrEqual(0);
          if (direction === 'genomic' ? anchor < lay.box.x - period : anchor > lay.box.x + lay.box.width + period)
            off++;
        }
        // some views had the start out of view: the phase was kept from there
        expect(off).toBeGreaterThan(3);
      }
    });
  }
});

describe('ranges: labels', () => {
  it('inside a mark narrower than a label’s room only when the label fits whole (a residue’s code)', async () => {
    // 40 positions over the data's width: about 17 px a position, room for a one-letter code, not a long label
    const rows = { start: [5, 20], end: [6, 21], id: ['a', 'b'], label: ['I', 'a long label'] };
    const lay = layerLay<RangesLay>(
      await (await settle(spec(rows as never, {}, { length: 40 }))).frame({ size: SIZE }),
      0,
    );
    const w = lay.mx1[0] - lay.mx0[0];
    expect(w).toBeLessThan(24);
    const labelled = Array.from({ length: lay.m }, (_, j) => (lay.mlabel[j] > 0 ? lay.keys[lay.mitem[j]] : undefined));
    expect(labelled).toEqual(['a', undefined]);
    expect(lay.mlabel[0] + 6).toBeLessThanOrEqual(w);
  });
});

describe('ranges: arrowheads', () => {
  const MANY: Rows = { start: [10, 100, 250], end: [60, 220, 390], id: ['a', 'b', 'c'] };

  it('`end` is on the mark holding the item’s end, on the right going forward and the left reversed', async () => {
    for (const [direction, tip] of [
      ['genomic', 2],
      ['reversed', 1],
    ] as const) {
      const s = await settle(spec(MANY, { arrow: 'end' }, { section: { direction } }));
      const f = await s.frame({ size: SIZE });
      const m = marks(layerLay<RangesLay>(f, 0));
      expect(m.map((x) => x.tips)).toEqual([tip, tip, tip]);
      expect(f.lists.passes.data.paths.n).toBe(3);
    }
  });

  it('`start` is at the other end, `both` at both, `none` and no param at neither', async () => {
    const tipsOf = async (arrow?: string) =>
      marks(layerLay<RangesLay>(await (await settle(spec(MANY, arrow ? { arrow } : {}))).frame({ size: SIZE }), 0)).map(
        (x) => x.tips,
      );
    expect(await tipsOf('start')).toEqual([1, 1, 1]);
    expect(await tipsOf('both')).toEqual([3, 3, 3]);
    expect(await tipsOf('none')).toEqual([0, 0, 0]);
    expect(await tipsOf()).toEqual([0, 0, 0]);
  });

  it('is never on a cut end: the view’s, a piece’s edge that is not the item’s end', async () => {
    const s = await settle(spec(ONE, { arrow: 'both' }));
    const tips = async (start: number, end: number) =>
      marks(layerLay<RangesLay>(await s.frame({ size: SIZE, viewport: { start, end } }), 0))[0].tips;
    expect(await tips(0, 400)).toBe(3);
    expect(await tips(100, 200)).toBe(0);
    expect(await tips(100, 310)).toBe(2);
    expect(await tips(30, 200)).toBe(1);
    // collapsed: the item [5, 25) over two pieces: the start on the first, the end on the second
    const pieces = { start: [0, 20], end: [10, 30] };
    const c = await settle(spec({ start: [5], end: [25], id: ['c'] }, { arrow: 'both' }, { length: 40, pieces }));
    const m = marks(layerLay<RangesLay>(await c.frame({ size: SIZE }), 0));
    expect(m.map((x) => x.tips)).toEqual([1, 2]);
  });

  it('across the origin of a cyclic system: the end on the row after it, the start on the row before', async () => {
    const wrap: Rows = { start: [90, 0], end: [100, 5], id: ['w', 'w'] };
    for (const [direction, start, end] of [
      ['genomic', 1, 2],
      ['reversed', 2, 1],
    ] as const) {
      const s = await settle(spec(wrap, { arrow: 'both' }, { length: 100, cyclic: true, section: { direction } }));
      const lay = layerLay<RangesLay>(await s.frame({ size: SIZE }), 0);
      for (let j = 0; j < lay.m; j++) {
        const row = lay.starts[lay.mitem[j]];
        // the row from 90 holds the start, the row after the origin the end
        expect(lay.mtip[j]).toBe(row === 90 ? start : end);
      }
    }
  });

  it('is decided per item by its width, not by the layer’s step: at the bars step a wide item keeps its arrowhead', async () => {
    // 8,000 positions over 800 px: a tenth of a pixel a position, the bars step
    const rows: Rows = {
      start: [500, 3000, 3020, 6000],
      end: [2500, 3010, 3030, 7000],
      id: ['wide', 'n1', 'n2', 'w2'],
    };
    const s = await settle(spec(rows, { arrow: 'end' }, { length: 8000 }));
    const lay = layerLay<RangesLay>(await s.frame({ size: SIZE }), 0);
    expect(lay.mode).toBe('bars');
    const tips = new Map(marks(lay).map((x) => [x.key, x.tips]));
    expect(tips.get('wide')).toBe(2);
    expect(tips.get('w2')).toBe(2);
    // a narrow item has none; narrow ones close together still merge into one run
    expect([...tips.keys()].filter((k) => k !== 'wide' && k !== 'w2').every((k) => !tips.get(k))).toBe(true);
    expect(lay.m).toBe(3);
  });

  it('a line gets a triangle at its end and stops where it begins; hits are the same', async () => {
    const s = await settle(spec(ONE, { line: 'solid', arrow: 'end' }));
    const f = await s.frame({ size: SIZE });
    const rects = readRects(f.lists);
    const triangle = readPaths(f.lists);
    expect(triangle.length).toBe(1);
    expect(triangle[0].xy.length / 2).toBe(7);
    expect(rects.length).toBe(1);
    const lay = layerLay<RangesLay>(f, 0);
    // the line ends where the head begins, the head's tip is at the item's end
    const tip = Math.max(...triangle[0].xy.filter((_, i) => i % 2 === 0));
    expect(tip).toBeCloseTo(lay.axis.snap(lay.axis.x(300)), 3);
    expect(rects[0].x + rects[0].w).toBeLessThan(tip);
    expect(hitView(f.layout, tip - 1, lay.box.y + lay.laneHeight / 2)?.element?.key).toBe('a');
  });
});
