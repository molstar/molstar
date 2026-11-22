/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The layout's edges: items that only touch the view draw and hit nothing; a follower placed above its leader; a
// cyclic view wrapping a short system reads each tile once; items past the lanes collapse into a density row (or the
// layer); values show a column's extent; the text metrics' version in the key; an export's deadline covers opening.

import { categoryColumn, type Node } from '@molstar/core/repr/tree';
import type { FormatRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { FakeClock } from '../../data/_test/synthetic.js';
import { hitView, type RangesLay, recordDrawLists } from '../index.js';
import { renderSvg } from '../render-svg.js';
import { bigSpec, p53GeneSpec } from './fixtures.js';
import { blocks, settle } from './helpers.js';

const Size = { width: 728, height: 120 };

function rangesSpec(rows: { start: number[]; end: number[]; id: string[] }, letters = 'M'.repeat(40)): Node {
  const b = createSeqSpec({ title: 'ranges' });
  const s = b.data.sequence('s', { alphabet: 'protein', letters });
  const r = b.data.ranges('r', {
    target: s,
    rows: { start: Int32Array.from(rows.start), end: Int32Array.from(rows.end), id: categoryColumn(rows.id) },
  });
  b.views.view('v', (v) => v.section('s', { of: s }).track('r').layer('ranges', { type: 'ranges', data: r }));
  return b.build();
}

describe('the layout’s edges', () => {
  it('items that only touch the view’s edges draw nothing and hit nothing', async () => {
    const s = await settle(rangesSpec({ start: [0, 20], end: [10, 30], id: ['left', 'right'] }));
    const f = await s.frame({ size: Size, viewport: { start: 10, end: 20 } });
    const lay = f.layout.sections[0].tracks[0].layers[0].lay as RangesLay;
    expect(lay.m).toBe(0);
    expect(blocks(f.lists)[0].rects).toBe(0);
    const { box } = f.layout.sections[0].tracks[0];
    for (const x of [box.x, box.x + 0.5, box.x + box.width - 0.5])
      expect(hitView(f.layout, x, box.y + 2)?.element).toBeUndefined();
  });

  it('pieces of one id in different tiles are all drawn; a range over a tile edge is drawn once', async () => {
    // a transcript's two exons share its id, read in two tiles; `edge` crosses a tile's edge and is in both
    const letters = 'M'.repeat(8_000);
    const rows = { start: [100, 6_000, 1_000], end: [200, 6_200, 1_100], id: ['exons', 'exons', 'edge'] };
    const s = await settle(rangesSpec(rows, letters));
    const f = await s.frame({ size: Size, viewport: { start: 0, end: 8_000 } });
    const lay = f.layout.sections[0].tracks[0].layers[0].lay as RangesLay;
    const items = Array.from({ length: lay.items }, (_, j) => [lay.keys[j], lay.starts[j]]);
    expect(items.sort((a, b) => (a[1] as number) - (b[1] as number))).toEqual([
      ['exons', 100],
      ['edge', 1_000],
      ['exons', 6_000],
    ]);
  });

  it('a hidden layer is kept in the spec and not drawn', async () => {
    const spec = (hidden: boolean) => {
      const b = createSeqSpec({ title: 'ranges' });
      const s = b.data.sequence('s', { alphabet: 'protein', letters: 'M'.repeat(40) });
      const r = b.data.ranges('r', { target: s, rows: { start: [0, 10], end: [5, 20] } });
      b.views.view('v', (v) =>
        v.section('s', { of: s }).track('r').layer('ranges', { type: 'ranges', data: r, hidden }),
      );
      return b.build();
    };
    const shown = await renderSvg(spec(false), { size: Size });
    const hidden = await renderSvg(spec(true), { size: Size });
    expect(blocks(shown.lists)[0].rects).toBeGreaterThan(0);
    expect(blocks(hidden.lists)).toEqual([]);
  });

  it('ranges at bars are squircles with round ends where a run really ends, as their emphasis is', async () => {
    const b = createSeqSpec({ title: 'bars' });
    const s = b.data.sequence('s', { alphabet: 'protein', letters: 'M'.repeat(4_000) });
    // two touching ranges merge into one run; the third stands alone
    const r = b.data.ranges('r', { target: s, rows: { start: [0, 1_000, 2_500], end: [1_000, 2_000, 3_500] } });
    b.views.view('v', (v) => v.section('s', { of: s }).track('r').layer('ranges', { type: 'ranges', data: r }));
    const f = await renderSvg(b.build(), { size: Size });
    const rects = recordDrawLists(f.lists, ['data'])
      .split('\n')
      .filter((l) => l.startsWith('rect '));
    expect(rects.length).toBe(2);
    for (const l of rects) expect(l).toMatch(/ squircle [\d.]+$/);
  });

  it('a follower placed above its leader follows it still', async () => {
    const reordered = p53GeneSpec();
    const view = reordered.children!.find((c) => c.kind === 'views')!.children![0];
    (view.children as Node[]).reverse();
    const viewport = { start: 150, end: 160 };
    const a = await settle(p53GeneSpec());
    const b = await settle(reordered);
    const fa = await a.frame({ size: Size, viewport });
    const fb = await b.frame({ size: Size, viewport });
    expect(fb.layout.issues).toEqual([]);
    const gene = (f: typeof fa) => f.layout.sections.find((s) => s.section.address.endsWith('/gene'))!.axis!;
    expect(fb.layout.sections[0].section.address).toBe('/views/main/gene');
    expect([gene(fb).pieces, gene(fb).u0, gene(fb).u1]).toEqual([gene(fa).pieces, gene(fa).u0, gene(fa).u1]);
    // continuous: the leader's window edges carried by `point`, fractional within the follower's pieces
    expect(gene(fa).pieces[1] - gene(fa).pieces[0]).toBeGreaterThanOrEqual(30);
  });

  it('a section following one outside its view, or a cycle, says so and shows its own range', async () => {
    const b = createSeqSpec({ title: 'cycle' });
    const s = b.data.sequence('s', { alphabet: 'protein', letters: 'M'.repeat(20) });
    b.views.view('v', (v) => {
      const one = v.section('one', { of: s, follow: '/views/v/two' });
      one.track('l').layer('letters', { type: 'letters' });
      v.section('two', { of: s, follow: one }).track('l').layer('letters', { type: 'letters' });
    });
    const st = await settle(b.build());
    const f = await st.frame({ size: Size });
    expect(f.layout.issues.map((i) => i.code)).toEqual(['seq.render.follow.cycle']);
    expect(f.layout.sections.every((x) => x.axis)).toBe(true);
  });

  it('a cyclic view wrapping a short system reads each tile once', async () => {
    const b = createSeqSpec({ title: 'cyclic' });
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(5), cyclic: true });
    const v = b.data.values('v', { target: s, values: Float32Array.from({ length: 20 }, (_, i) => i + 1) });
    b.views.view('main', (view) =>
      view
        .section('s', { of: s, range: [15, 5] })
        .track('v')
        .layer('values', { type: 'values', data: v }),
    );
    const r = await renderSvg(b.build(), { size: Size });
    expect(r.issues).toEqual([]);
    // ten bars and the baseline's hairline
    expect(blocks(r.lists)[0].rects).toBe(11);
  });

  it('items past the lanes collapse into a density row at the bottom; or the whole layer, as configured', async () => {
    // seven over one another: past the five lanes a track capped at 78 px holds
    const rows = {
      start: [0, 1, 2, 3, 4, 5, 6],
      end: [30, 30, 30, 30, 30, 30, 30],
      id: ['a', 'b', 'c', 'd', 'e', 'f', 'g'],
    };
    const cap = { track: { maxPx: 78, overflow: 'density' as const }, export: { lanes: 'view' as const } };
    const extra = await renderSvg(rangesSpec(rows), { size: Size, config: cap });
    const [e] = blocks(extra.lists);
    expect(e).toMatchObject({ rects: 4, images: 1 });
    // the density row is the fifth lane, below the padding (4 px): 8 + 4 + 4 × 14
    expect(recordDrawLists(extra.lists, ['data'])).toMatch(/\nimage 96 68 624 14 624x1 /);
    const layer = await renderSvg(rangesSpec(rows), { size: Size, config: { ...cap, ranges: { overflow: 'layer' } } });
    expect(blocks(layer.lists)[0]).toMatchObject({ rects: 0, images: 1 });
    expect(recordDrawLists(layer.lists, ['data'])).toMatch(/\nimage 96 8 624 78 624x1 /);
  });

  it('values below a pixel a position show a column’s extent (negative extremes), or its maximum', async () => {
    const b = createSeqSpec({ title: 'signal' });
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'A'.repeat(2000) });
    const v = b.data.values('v', {
      target: s,
      values: Float32Array.from({ length: 2000 }, (_, i) => (i % 2 ? 2 : -1)),
    });
    b.views.view('main', (view) =>
      view.section('s', { of: s }).track('v').layer('values', { type: 'values', data: v }),
    );
    const spec = b.build();
    const minmax = await renderSvg(spec, { size: Size });
    // the track is 32 px from y 8: every column spans the domain [-1, 2], the negative extreme below the baseline
    expect(minmax.svg).toContain('data-seq-domains="/views/main/s/v/values=-1,2"');
    // a profile (the area and its line), and under it the columns' extents below the baseline: a block a column
    const rects = (r: { lists: Parameters<typeof recordDrawLists>[0] }) =>
      recordDrawLists(r.lists, ['data'])
        .split('\n')
        .filter((l) => l.startsWith('rect ')).length;
    expect(recordDrawLists(minmax.lists, ['data'])).toContain('\npath ');
    // the domain's 3 units over 90% of the 32 px: the baseline (zero) at 40 - 9.6; a column's block from there (the
    // profile's area covers the baseline to the maximum) down to the track's floor (-1), one column wide
    expect(recordDrawLists(minmax.lists, ['data'])).toContain('\nrect 96 30.4 1 9.6 #5b7c9959\n');
    // the maximum alone hides it: the domain starts at 0
    const max = await renderSvg(spec, { size: Size, config: { values: { aggregate: 'max' } } });
    expect(max.svg).toContain('data-seq-domains="/views/main/s/v/values=0,2"');
    expect(rects(minmax)).toBeGreaterThan(rects(max) + 100);
  });
});

describe('keys and deadlines', () => {
  it('clearing the text metrics changes the key', async () => {
    const s = await settle(rangesSpec({ start: [0], end: [10], id: ['a'] }));
    const first = await s.frame({ size: Size });
    expect(s.renderer.frame({ size: Size }).rebuilt).toBe(false);
    s.renderer.options.text.clear();
    const after = s.renderer.frame({ size: Size });
    expect(after.rebuilt).toBe(true);
    expect(after.key).not.toBe(first.key);
  });

  it('an export’s deadline covers opening its sources: one never opening times out, and is aborted', async () => {
    const clock = new FakeClock();
    let signal: AbortSignal | undefined;
    const never: FormatRuntime = {
      name: 'never',
      formats: ['synthetic'],
      open: (_s, ctx) => {
        signal = ctx.signal;
        return new Promise(() => {});
      },
    };
    const out = renderSvg(bigSpec(), { size: Size, runtimes: [never], clock, timeoutMs: 1000 });
    await clock.advance(1000);
    const r = await out;
    expect(r.settled).toBe(false);
    expect(r.issues.map((i) => i.code)).toEqual(['seq.render.axis.pending', 'seq.render.unsettled']);
    expect(signal?.aborted).toBe(true);
  });
});
