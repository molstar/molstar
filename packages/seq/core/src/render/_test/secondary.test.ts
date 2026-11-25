/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The secondary layer (_spec/README.md "Render", Layer types): blocks drawn against the data (a mark for every item in
// view, every mark hit), the strand's arrowhead on the piece holding its C-terminal end whatever the direction, the
// origin or the path, cut ends square, the copies of a chain drawn once (EC-130, EC-67), the fallback to ranges.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { hitView, recordDrawLists, type FrameInput, type SecondaryLay } from '../index.js';
import { BigRuntime, GENE_EXONS, geneLetters, P53, SecondaryRows } from './fixtures.js';
import { settle } from './helpers.js';
import { frameWithout, layerLay } from './small-helpers.js';

type Rows = {
  start: number[] | Int32Array;
  end: number[] | Int32Array;
  id?: string[];
  category: string[];
};

/** Item keys with their tip (0 none, 1 right, 2 left) a mark, in mark order. */
function tips(lay: SecondaryLay): { key: string | number; tip: number; x0: number; x1: number }[] {
  const out = [];
  for (let j = 0; j < lay.m; j++)
    out.push({ key: lay.rows.keys[lay.mitem[j]], tip: lay.mtip[j], x0: lay.mx0[j], x1: lay.mx1[j] });
  return out;
}

/** A sequence of `length` letters with the rows as secondary structure, in a section with the options. */
function spec(
  rows: Rows,
  options: {
    length?: number;
    cyclic?: boolean;
    section?: Record<string, unknown>;
    pieces?: { start: number[]; end: number[] };
  } = {},
): Node {
  const length = options.length ?? 100;
  const b = createSeqSpec({ title: 'secondary' });
  const s = b.data.sequence('s', {
    alphabet: 'dna',
    letters: 'ACGT'.repeat(Math.ceil(length / 4)).slice(0, length),
    ...(options.cyclic ? { cyclic: true } : {}),
  });
  const ss = b.data.ranges('ss', { target: s, rows: { ...rows, id: rows.id } as never });
  const pieces = options.pieces ? b.data.ranges('pieces', { target: s, rows: options.pieces }) : undefined;
  b.views.view('main', (v) => {
    v.section('s', { of: s, ...(pieces ? { pieces } : {}), ...(options.section as object) })
      .track('ss')
      .layer('ss', { type: 'secondary', data: ss });
  });
  return b.build();
}

const SIZE = { width: 800, height: 40 } as const;

describe('secondary: drawn against the data', () => {
  it('p53 at the fit view: a mark for each item, each mark hit as its item', async () => {
    const b = createSeqSpec({ title: 'p53 secondary' });
    const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
    const ss = b.data.ranges('ss', { target: p53, rows: SecondaryRows });
    b.views.view('main', (v) => {
      v.section('protein', { of: p53 }).track('ss').layer('ss', { type: 'secondary', data: ss });
    });
    const s = await settle(b.build());
    const f = await s.frame({ size: SIZE });
    const lay = layerLay<SecondaryLay>(f, 0);
    expect(lay.m).toBe(SecondaryRows.id.length);
    const box = f.layout.sections[0].tracks[0].box;
    for (const m of tips(lay)) {
      const h = hitView(f.layout, (m.x0 + m.x1) / 2, box.y + box.height / 2);
      expect(h?.element?.key).toBe(m.key);
      expect(h?.hull).toBeDefined();
    }
    expect(f.layout.issues).toEqual([]);
  });

  it('across a tile edge: an item over the edge is one mark, however the view meets the edge', async () => {
    const b = createSeqSpec({ title: 'tiles' });
    const src = b.data.source('synthetic', { format: 'synthetic', uri: 'synthetic:chr' });
    const chr = b.data.sequence('chr', { alphabet: 'dna', source: src });
    const items = b.data.ranges('items', { target: chr, source: src });
    b.views.view('main', (v) => {
      v.section('chr', { of: chr }).track('ss').layer('ss', { type: 'secondary', data: items });
    });
    const s = await settle(b.build(), [BigRuntime]);
    // items every 50, 70 long: one starts at 1000 and runs over the tile edge at 1024
    for (const [start, end] of [
      [990, 1060],
      [1010, 1030],
      [1000, 1100],
    ]) {
      const f = await s.frame({ size: SIZE, viewport: { start, end } });
      const lay = layerLay<SecondaryLay>(f, 0);
      const keys = tips(lay).map((m) => m.key);
      expect(new Set(keys).size).toBe(keys.length);
      // the items overlapping the view: starts p with p + 70 > start and p < end, p a multiple of 50
      const want: string[] = [];
      for (let p = 0; p < end; p += 50) if (p + 70 > start) want.push(`r${p}`);
      expect(keys.sort()).toEqual(want.sort());
    }
  });

  it('far out the layer is one strip of the items’ density (the source answers bins)', async () => {
    const b = createSeqSpec({ title: 'density' });
    const src = b.data.source('synthetic', { format: 'synthetic', uri: 'synthetic:chr' });
    const chr = b.data.sequence('chr', { alphabet: 'dna', source: src });
    const items = b.data.ranges('items', { target: chr, source: src });
    b.views.view('main', (v) => {
      v.section('chr', { of: chr }).track('ss').layer('ss', { type: 'secondary', data: items });
    });
    const s = await settle(b.build(), [BigRuntime]);
    const f = await s.frame({ size: { width: 1440, height: 40 }, viewport: { start: 0, end: 3_000_000 } });
    const lay = layerLay<SecondaryLay>(f, 0);
    expect(lay.density).toBeDefined();
    expect(lay.m).toBe(0);
    expect(f.lists.passes.data.images.n).toBe(1);
  });
});

describe('secondary: the arrowhead (EC-130)', () => {
  const strand: Rows = { start: [40, 70], end: [55, 90], id: ['s1', 's2'], category: ['strand', 'helix'] };

  it('is on the strand’s C-terminal end, pointing the way the residues run on the screen', async () => {
    for (const [section, want] of [
      [{}, 1],
      [{ direction: 'reversed' }, 2],
    ] as const) {
      const s = await settle(spec(strand, { section }));
      const f = await s.frame({ size: SIZE });
      const m = tips(layerLay<SecondaryLay>(f, 0));
      expect(m.map((x) => [x.key, x.tip])).toEqual(
        expect.arrayContaining([
          ['s1', want],
          ['s2', 0],
        ]),
      );
      expect(m.filter((x) => x.tip).length).toBe(1);
    }
  });

  it('is not drawn where the view clips the end, drawn when only the N end is clipped', async () => {
    const s = await settle(spec(strand));
    const clippedC = tips(layerLay<SecondaryLay>(await s.frame({ size: SIZE, viewport: { start: 30, end: 50 } }), 0));
    expect(clippedC.find((x) => x.key === 's1')?.tip).toBe(0);
    const clippedN = tips(layerLay<SecondaryLay>(await s.frame({ size: SIZE, viewport: { start: 45, end: 80 } }), 0));
    expect(clippedN.find((x) => x.key === 's1')?.tip).toBe(1);
    // the end a position short of the view's edge is shown
    const edge = tips(layerLay<SecondaryLay>(await s.frame({ size: SIZE, viewport: { start: 30, end: 55 } }), 0));
    expect(edge.find((x) => x.key === 's1')?.tip).toBe(1);
  });

  it('across the origin of a cyclic system it is on the piece after the origin, in either direction', async () => {
    // the strand `x` is two rows: 90 to the origin and the origin to 5; read 90, 91, …, 99, 0, …, 4
    const wrap: Rows = {
      start: [90, 0, 10],
      end: [100, 5, 20],
      id: ['x', 'x', 'y'],
      category: ['strand', 'strand', 'strand'],
    };
    for (const [section, side] of [
      [{}, 1],
      [{ direction: 'reversed' }, 2],
      [{ range: [80, 30] }, 1],
      [{ range: [80, 30], direction: 'reversed' }, 2],
    ] as const) {
      const s = await settle(spec(wrap, { cyclic: true, section }));
      const f = await s.frame({ size: SIZE });
      const lay = layerLay<SecondaryLay>(f, 0);
      const x = tips(lay).filter((m) => m.key === 'x');
      expect(x.length).toBe(2);
      // the row after the origin holds the tip: the one whose item rows start at 0
      const tipped = x.filter((m) => m.tip);
      expect(tipped.length).toBe(1);
      expect(tipped[0].tip).toBe(side);
      const at = lay.axis.positionAt((tipped[0].x0 + tipped[0].x1) / 2);
      expect(at).toBeLessThan(5);
      // the lone strand `y` has its own
      expect(tips(lay).find((m) => m.key === 'y')?.tip).toBe(side);
    }
  });

  it('on a collapsed axis it is on the piece holding the end, whatever the direction', async () => {
    const pieces = { start: [0, 20], end: [10, 30] };
    const crossing: Rows = { start: [5], end: [25], id: ['c'], category: ['strand'] };
    for (const [direction, side] of [
      ['genomic', 1],
      ['reversed', 2],
    ] as const) {
      const s = await settle(spec(crossing, { length: 40, pieces, section: { direction } }));
      const f = await s.frame({ size: SIZE });
      const lay = layerLay<SecondaryLay>(f, 0);
      const m = tips(lay);
      expect(m.length).toBe(2);
      expect(m.filter((x) => x.tip).length).toBe(1);
      const tipped = m.find((x) => x.tip)!;
      expect(tipped.tip).toBe(side);
      // the piece with the C end is [20, 30): its cells are at the right going forward, at the left reversed
      expect(lay.axis.positionAt((tipped.x0 + tipped.x1) / 2)).toBeGreaterThanOrEqual(20);
    }
  });

  it('through a codon map: the same side on the protein and on the gene that reads along it, the other where it does not', async () => {
    const b = createSeqSpec({ title: 'secondary on the gene' });
    const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
    const gene = b.data.sequence('tp53', { alphabet: 'dna', letters: geneLetters() });
    b.data.map('p53-on-gene', { from: p53, to: gene, method: 'codon', exons: GENE_EXONS, strand: '-' });
    const ss = b.data.ranges('ss', { target: p53, rows: SecondaryRows });
    b.views.view('main', (v) => {
      const protein = v.section('protein', { of: p53 });
      protein.track('ss').layer('ss', { type: 'secondary', data: ss });
      v.section('along', { of: gene, follow: protein }).track('ss').layer('ss', { type: 'secondary', data: ss });
      v.section('genomic', { of: gene, follow: protein, direction: 'genomic' })
        .track('ss')
        .layer('ss', { type: 'secondary', data: ss });
    });
    const s = await settle(b.build());
    const f = await s.frame({ size: { width: 1200, height: 200 } });
    const strands = SecondaryRows.id.filter((_, i) => SecondaryRows.category[i] === 'strand');
    for (const [section, side] of [
      [0, 1],
      [1, 1],
      [2, 2],
    ] as const) {
      const m = tips(layerLay<SecondaryLay>(f, 0, 0, section));
      const tipped = m.filter((x) => x.tip);
      // one tip a strand, on the same side
      expect(tipped.map((x) => x.key).sort()).toEqual([...strands].sort());
      for (const x of tipped) expect(x.tip).toBe(side);
    }
  });
});

describe('secondary: cut ends are square', () => {
  it('a mark cut by the view has no round end on the cut side; its real end is round', async () => {
    const rows: Rows = { start: [20], end: [60], id: ['h'], category: ['helix'] };
    const s = await settle(spec(rows));
    const lay = layerLay<SecondaryLay>(await s.frame({ size: SIZE, viewport: { start: 40, end: 80 } }), 0);
    // cut at its left: the left end square, the right (position 60) round
    expect(lay.mround[0] & 1).toBe(0);
    expect(lay.mround[0] & 2).toBe(2);
    const whole = layerLay<SecondaryLay>(await s.frame({ size: SIZE }), 0);
    expect(whole.mround[0]).toBe(3);
    const both = layerLay<SecondaryLay>(await s.frame({ size: SIZE, viewport: { start: 30, end: 50 } }), 0);
    expect(both.mround[0]).toBe(0);
  });

  it('across a piece edge: a piece’s edge cuts the mark, the item’s own end does not', async () => {
    const s = await settle(
      spec(
        { start: [5], end: [25], id: ['h'], category: ['helix'] },
        { length: 40, pieces: { start: [0, 20], end: [10, 30] } },
      ),
    );
    const lay = layerLay<SecondaryLay>(await s.frame({ size: SIZE }), 0);
    // [5, 10) is cut at its right by the piece; [20, 25) is cut at its left and ends at 25: round only at its right
    expect(lay.mround[0]).toBe(1);
    expect(lay.mround[1]).toBe(2);
  });
});

describe('secondary: copies (EC-67)', () => {
  /** The chain's structure `copies` times over, one id each. */
  const copies = (n: number): Rows => {
    const rows: Rows = { start: [], end: [], id: [], category: [] };
    const start = [...SecondaryRows.start];
    const end = [...SecondaryRows.end];
    for (let c = 0; c < n; c++)
      for (let i = 0; i < start.length; i++) {
        (rows.start as number[]).push(start[i]);
        (rows.end as number[]).push(end[i]);
        rows.id!.push(`${SecondaryRows.id[i]}-copy${c}`);
        rows.category.push(SecondaryRows.category[i]);
      }
    return rows;
  };

  it('200 identical copies draw what one does, in the lists and in the work', async () => {
    const one = await settle(spec(copies(1), { length: 393 }));
    const many = await settle(spec(copies(200), { length: 393 }));
    const input: FrameInput = { size: { width: 800, height: 40 } };
    const a = await one.frame(input);
    const lay = layerLay<SecondaryLay>(a, 0);
    const marks = lay.m;
    const record = recordDrawLists(a.lists, ['data']);
    const b = await many.frame(input);
    expect(layerLay<SecondaryLay>(b, 0).m).toBe(marks);
    expect(recordDrawLists(b.lists, ['data'])).toBe(record);
    // and a hit names an item of the group, the first
    const box = b.layout.sections[0].tracks[0].box;
    const m = tips(layerLay<SecondaryLay>(b, 0))[0];
    const h = hitView(b.layout, (m.x0 + m.x1) / 2, box.y + box.height / 2);
    expect(String(h?.element?.key)).toMatch(/-copy0$/);
  });

  it('copies that differ in kind or extent are all drawn', async () => {
    const rows: Rows = {
      start: [10, 10, 10],
      end: [30, 30, 31],
      id: ['a', 'b', 'c'],
      category: ['helix', 'strand', 'helix'],
    };
    const s = await settle(spec(rows));
    expect(layerLay<SecondaryLay>(await s.frame({ size: SIZE }), 0).m).toBe(3);
  });
});

describe('secondary: the fallback', () => {
  it('a target without the type draws ranges: the same items as marks of the ranges layer', async () => {
    const rows = SecondaryRows;
    const sp = (() => {
      const b = createSeqSpec({ title: 'fallback' });
      const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
      const ss = b.data.ranges('ss', { target: p53, rows });
      b.views.view('main', (v) => {
        v.section('protein', { of: p53 }).track('ss').layer('ss', { type: 'secondary', data: ss });
      });
      return b.build();
    })();
    const f = await frameWithout(sp, ['secondary'], { size: SIZE });
    expect(f.layout.issues).toEqual([]);
    const layer = f.layout.sections[0].tracks[0].layers[0];
    expect(layer.type.name).toBe('ranges');
    expect(layer.lay).toHaveProperty('keys');
    const none = await frameWithout(sp, ['secondary', 'ranges'], { size: SIZE });
    expect(none.layout.issues.map((i) => i.code)).toEqual(['seq.render.layer.type']);
  });
});
