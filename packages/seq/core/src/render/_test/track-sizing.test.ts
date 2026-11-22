/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Track sizing (_spec/README.md "Render", Track heights): a track of lanes is as tall as its lanes as known (a source's
// grown as they settle, one loading lane while pending) up to `track.maxPx`, past which its lanes scroll inside it
// with the header and the tracks around in place; an export draws every lane; a ring compacts; a gene model's label
// with no room floats at the data's left edge; a rows band with no height takes its rows up to `rows.maxPx`.

import type { Node } from '@molstar/core/repr/tree';
import { RecordsRuntime } from '@molstar/seqx-core/data/records';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { GffRuntime } from '../../runtimes/gff.js';
import { DefaultRenderConfig, type RangesLay, recordDrawLists } from '../index.js';
import { renderSvg } from '../render-svg.js';
import { lanesHeight } from '../sizing.js';
import { settle } from './helpers.js';
import { layerLay } from './small-helpers.js';

const Size = { width: 800, height: 900 };
const Section = '/views/main/s';

/** `n` rows each over the one before (one lane each). */
function stacked(n: number, at = 0) {
  return {
    start: Array.from({ length: n }, (_, i) => at + i),
    end: Array.from({ length: n }, () => at + 150),
    id: Array.from({ length: n }, (_, i) => `r${at}-${i}`),
  };
}

/** One section of ranges tracks (inline rows), each with a lane count; the first track's own height when given. */
function spec(lanes: number[], options: { height?: number; layout?: string } = {}): Node {
  const b = createSeqSpec({ title: 'sizing' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(100) });
  b.views.view('main', (v) => {
    const sec = v.section('s', { of: s, ...(options.layout ? { layout: options.layout } : {}) } as never);
    lanes.forEach((n, k) => {
      const r = b.data.ranges(`r${k}`, { target: s, rows: stacked(n) });
      sec
        .track(`t${k}`, k === 0 && options.height !== undefined ? { height: options.height } : {})
        .layer('l', { type: 'ranges', data: r });
    });
  });
  return b.build();
}

/** A GFF3 of `n` transcripts over one another (one lane each), the first starting at the region's start. */
function gff(n: number, at = 0): string {
  const lines = ['##gff-version 3', '##sequence-region c 1 4000'];
  for (let t = 0; t < n; t++) {
    const s = at + 1 + (t % 3);
    lines.push(
      ['c', '.', 'mRNA', `${s}`, `${at + 900}`, '.', '+', '.', `ID=x${at}_${t};Name=T${t}`].join('\t'),
      ['c', '.', 'exon', `${s}`, `${at + 300}`, '.', '+', '.', `Parent=x${at}_${t}`].join('\t'),
      ['c', '.', 'exon', `${at + 700}`, `${at + 900}`, '.', '+', '.', `Parent=x${at}_${t}`].join('\t'),
    );
  }
  return lines.join('\n');
}

function geneSpec(text: string, height?: number): Node {
  const b = createSeqSpec({ title: 'genes' });
  const file = b.data.source('gff', { format: 'gff3', value: text as never });
  const chr = b.data.sequence('c', { alphabet: 'dna', source: file, name: 'c' });
  const tx = b.data.ranges('tx', { target: chr, source: file });
  b.views.view('main', (v) => {
    v.section('c', { of: chr })
      .track('tx', height !== undefined ? ({ height } as never) : {})
      .layer('tx', { type: 'ranges', data: tx, line: 'gene' } as never);
  });
  return b.build();
}

const bandOf = (track: string, bandTop: number) => ({ [track]: { bandTop } });

describe('a track of lanes is as tall as its lanes, to its cap, then its lanes scroll inside it', () => {
  it('the GFF transcripts are sized to their lanes as they settle: one loading lane while pending, no blank band, no rail', async () => {
    const s = await settle(geneSpec(gff(9)), [GffRuntime]);
    const pending = await s.frame({ size: Size, settled: false });
    expect(pending.layout.sections[0].tracks[0].box.height).toBe(lanesHeight(1, DefaultRenderConfig));
    const f = await s.frame({ size: Size, settled: true });
    const t = f.layout.sections[0].tracks[0];
    // the design's picture: nine lanes, 134 px, tracks below moved once
    expect(t.box.height).toBe(134);
    const lay = layerLay<RangesLay>(f, 0);
    expect(lay.laneCount).toBe(9);
    expect(lay.scrollMax).toBe(0);
    // a pan never shrinks it
    const pan = await s.frame({ size: Size, settled: true, viewport: { start: 3000, end: 3100 } });
    expect(pan.layout.sections[0].tracks[0].box.height).toBe(134);
  });

  it('a 40-lane track is capped at 240 px and its lanes scroll inside it; the header and the tracks around stay', async () => {
    const s = await settle(spec([3, 40, 3]));
    const t1 = '/views/main/s/t1';
    const at = async (bandTop: number) =>
      s.frame({
        size: Size,
        settled: true,
        viewports: { [Section]: { start: 0, end: 200, ...(bandTop ? { bands: bandOf(t1, bandTop) } : {}) } },
      });
    const top = await at(0);
    const tracks = top.layout.sections[0].tracks;
    // the 40-lane track is far from the short ones: it takes the cap, they fit their own lanes
    const short = lanesHeight(3, DefaultRenderConfig);
    expect(tracks.map((t) => t.box.height)).toEqual([short, 240, short]);
    const lay = layerLay<RangesLay>(top, 1);
    expect(lay.scrolls).toBe(true);
    // 240 px hold (240 - 8) / 14 = 16 lanes: 24 more to scroll to
    expect(lay.room).toBe(16);
    expect(lay.scrollMax).toBe(40 - 16);
    expect(lay.top).toBe(0);
    // only the lanes in view are marks
    expect(lay.m).toBeLessThanOrEqual(18);
    const y0 = lay.laneY;
    const scrolled = await at(10);
    const lay2 = layerLay<RangesLay>(scrolled, 1);
    expect(lay2.top).toBe(10);
    expect(lay2.laneY).toBe(y0 - 10 * 14);
    expect(lay2.m).toBeLessThanOrEqual(18);
    // the tracks around: the same boxes, the header's label at the track's top
    expect(scrolled.layout.sections[0].tracks.map((t) => [t.box.y, t.box.height])).toEqual(
      tracks.map((t) => [t.box.y, t.box.height]),
    );
    // the band's top is kept inside its lanes
    const past = await at(99);
    expect(layerLay<RangesLay>(past, 1).top).toBe(24);
    // the rail and the cut edges are drawn only while lanes overflow
    const drawn = recordDrawLists(scrolled.lists, ['labels']);
    expect(drawn).toContain('rect');
  });

  it('a track’s own height is its cap: lanes past it scroll', async () => {
    const s = await settle(spec([20], { height: 120 }));
    const f = await s.frame({ size: Size, settled: true });
    expect(f.layout.sections[0].tracks[0].box.height).toBe(120);
    const lay = layerLay<RangesLay>(f, 0);
    expect(lay.scrolls).toBe(true);
    expect(lay.room).toBe(8);
    expect(lay.scrollMax).toBe(12);
  });

  it('`track.overflow` `compact` squeezes the lanes into the band, `density` puts the rest in a density row', async () => {
    const compact = await settle(spec([20], { height: 120 }), [], { track: { maxPx: 240, overflow: 'compact' } });
    const c = layerLay<RangesLay>(await compact.frame({ size: Size, settled: true }), 0);
    expect(c.scrolls).toBe(false);
    expect(c.laneHeight).toBeLessThan(14);
    const density = await settle(spec([20], { height: 120 }), [], { track: { maxPx: 240, overflow: 'density' } });
    const d = layerLay<RangesLay>(await density.frame({ size: Size, settled: true }), 0);
    expect(d.scrolls).toBe(false);
    expect(d.overflow).toBeGreaterThan(0);
  });

  it('EC-245: a pan to a region of more lanes keeps the band’s top and the track’s height (the reachable region’s lanes are known from the first frame)', async () => {
    // a region of 24 lanes, then one of 30 at the other end of the same source
    const text = [gff(24, 0), gff(30, 2000).split('\n').slice(2).join('\n')].join('\n');
    const s = await settle(geneSpec(text, 240), [GffRuntime]);
    const track = '/views/main/c/tx';
    const view = (start: number, end: number) => ({ '/views/main/c': { start, end, bands: bandOf(track, 5) } });
    const a = await s.frame({ size: Size, settled: true, viewports: view(0, 1000) });
    const la = layerLay<RangesLay>(a, 0);
    expect(la.top).toBe(5);
    // the lanes the whole region packs into (30), not those of the stretch in view (24)
    expect(la.scrollMax).toBeGreaterThan(la.laneCount - la.room);
    const before = la.scrollMax;
    const b = await s.frame({ size: Size, settled: true, viewports: view(2000, 3000) });
    const lb = layerLay<RangesLay>(b, 0);
    // the band's top stays; the lanes and the track as they were (a pan never resizes)
    expect(lb.top).toBe(5);
    expect(lb.scrollMax).toBe(before);
    expect(b.layout.sections[0].tracks[0].box.height).toBe(a.layout.sections[0].tracks[0].box.height);
    const c = await s.frame({ size: Size, settled: true, viewports: view(0, 1000) });
    expect(layerLay<RangesLay>(c, 0).scrollMax).toBe(before);
    expect(c.layout.sections[0].tracks[0].box.height).toBe(b.layout.sections[0].tracks[0].box.height);
  });

  it('EC-246: a filter leaving fewer lanes shrinks the track and clamps the band’s top (the data’s identity changed)', async () => {
    const track = '/views/main/s/t0';
    const many = await settle(spec([40]));
    const f1 = await many.frame({
      size: Size,
      settled: true,
      viewports: { [Section]: { start: 0, end: 200, bands: bandOf(track, 20) } },
    });
    expect(layerLay<RangesLay>(f1, 0).top).toBe(20);
    expect(f1.layout.sections[0].tracks[0].box.height).toBe(240);
    // the filter leaves ten: the same band top, now past the lanes
    const few = await settle(spec([10]));
    const f2 = await few.frame({
      size: Size,
      settled: true,
      viewports: { [Section]: { start: 0, end: 200, bands: bandOf(track, 20) } },
    });
    expect(f2.layout.sections[0].tracks[0].box.height).toBe(lanesHeight(10, DefaultRenderConfig));
    const lay = layerLay<RangesLay>(f2, 0);
    expect(lay.top).toBe(0);
    expect(lay.scrolls).toBe(false);
  });

  it('EC-248: two tracks of one kind, one overflowing: one height, each scrolls on its own', async () => {
    const s = await settle(spec([40, 20]));
    const t0 = '/views/main/s/t0';
    const f = await s.frame({
      size: Size,
      settled: true,
      viewports: { [Section]: { start: 0, end: 200, bands: bandOf(t0, 10) } },
    });
    const [a, b] = f.layout.sections[0].tracks;
    expect(a.box.height).toBe(240);
    expect(b.box.height).toBe(240);
    expect(layerLay<RangesLay>(f, 0).top).toBe(10);
    // the second track's lanes did not move with the first's
    const lb = layerLay<RangesLay>(f, 1);
    expect(lb.top).toBe(0);
    expect(lb.scrollMax).toBe(20 - 16);
  });

  it('EC-249: a track’s height under one lane draws one lane, the rest scrolled (never a band under a lane)', async () => {
    const s = await settle(spec([3], { height: 10 }));
    const f = await s.frame({ size: Size, settled: true });
    const t = f.layout.sections[0].tracks[0];
    expect(t.box.height).toBe(10);
    const lay = layerLay<RangesLay>(f, 0);
    expect(lay.scrolls).toBe(true);
    expect(lay.room).toBe(1);
    expect(lay.m).toBeGreaterThan(0);
    // a band under a lane shows the lane's top: anchored at the band's top, cut at its foot only
    expect(lay.laneY).toBe(t.box.y);
    // and its rail has a positive length, however short the band
    const short = await settle(spec([3], { height: 3 }));
    const g = await short.frame({ size: Size, settled: true });
    expect(layerLay<RangesLay>(g, 0).laneY).toBe(g.layout.sections[0].tracks[0].box.y);
    for (const r of g.layout.sections[0].lanesRails ?? []) expect(r.width).toBeGreaterThan(0);
    for (let j = 0; j < lay.m; j++) expect(lay.mlane[j]).toBe(0);
  });
});

describe('a gene model’s label with no room floats at the data’s left edge', () => {
  it('the first model of a region starts at its edge: its label is pinned there over a backdrop, never dropped', async () => {
    const s = await settle(geneSpec(gff(3)), [GffRuntime]);
    const f = await s.frame({ size: Size, settled: true, viewport: { start: 0, end: 1000 } });
    const lay = layerLay<RangesLay>(f, 0);
    const g = lay.gene;
    // every shown lane has its label
    expect(g.ln).toBe(3);
    let floating = 0;
    for (let k = 0; k < g.ln; k++) if (g.lfloat[k]) floating++;
    expect(floating).toBeGreaterThan(0);
    for (let k = 0; k < g.ln; k++) expect(g.lx[k]).toBeGreaterThanOrEqual(lay.box.x);
  });
});

describe('an export draws every lane; a ring compacts', () => {
  it('`renderSvg` draws a 40-lane track at its content height; `export.lanes: view` the band as on screen', async () => {
    const all = await renderSvg(spec([40]), { size: { width: 800, height: 'auto' }, runtimes: [RecordsRuntime] });
    const view = await renderSvg(spec([40]), {
      size: { width: 800, height: 'auto' },
      runtimes: [RecordsRuntime],
      config: { export: { lanes: 'view' } },
    });
    const rects = (r: typeof all) =>
      recordDrawLists(r.lists, ['data'])
        .split('\n')
        .filter((l) => l.startsWith('rect ')).length;
    expect(rects(all)).toBeGreaterThanOrEqual(40);
    expect(rects(view)).toBeLessThan(rects(all));
    expect(all.svg).toContain('<svg');
    // bounded: past `export.maxTrackPx` the track is drawn to it
    const bounded = await renderSvg(spec([40]), {
      size: { width: 800, height: 'auto' },
      runtimes: [RecordsRuntime],
      config: { export: { maxTrackPx: 300 } },
    });
    expect(rects(bounded)).toBeLessThan(rects(all));
    expect(rects(bounded)).toBeGreaterThan(rects(view));
  });

  it('a ring’s lanes track takes at most `rings.trackMaxPx` radially, no inner scroll: its lanes compact', async () => {
    const s = await settle(spec([20], { layout: 'rings' }));
    const f = await s.frame({ size: { width: 700, height: 700 }, settled: true });
    const lay = layerLay<RangesLay>(f, 0);
    expect(lay.scrolls).toBe(false);
    expect(lay.box.height).toBeLessThanOrEqual(DefaultRenderConfig.rings.trackMaxPx);
    // every lane stays inside its annulus
    for (let j = 0; j < lay.m; j++) {
      const y = lay.laneY + lay.mlane[j] * lay.laneHeight;
      expect(y).toBeGreaterThanOrEqual(lay.box.y - 0.01);
      expect(y + lay.laneHeight).toBeLessThanOrEqual(lay.box.y + lay.box.height + 0.01);
    }
  });
});

describe('the rows band with no `height` takes its rows up to `rows.maxPx`', () => {
  it('the config’s caps are the design’s', () => {
    expect(DefaultRenderConfig.track).toEqual({ maxPx: 240, overflow: 'scroll', empty: 'note', emptyPx: 18 });
    expect(DefaultRenderConfig.rows.maxPx).toBe(480);
    expect(DefaultRenderConfig.rings.trackMaxPx).toBe(48);
    expect(DefaultRenderConfig.export).toEqual({ lanes: 'all', maxTrackPx: 4000 });
    expect(DefaultRenderConfig.ranges.maxLanes).toBe(1000);
  });

  it('an alignment of 60 rows in a tall frame is capped at 480 px; with `height: fill` it takes the room', async () => {
    const rows = 60;
    const letters = Array.from({ length: rows }, (_, k) => `ACGT`.repeat(5).slice(0, 20 - (k % 3)));
    const alignment = (height?: string) => {
      const b = createSeqSpec({ title: 'rows' });
      const a = b.data.alignment('a', {
        rows: letters.map((l, k) => ({ name: `R${k}`, letters: l.padEnd(20, '-') })),
        alphabet: 'dna',
      } as never);
      b.views.view('main', (v) =>
        v
          .section('s', { of: a } as never)
          .track('rows')
          .layer('r', { type: 'rows', data: a, ...(height ? { height } : {}) } as never),
      );
      return b.build();
    };
    const tall = { width: 800, height: 1400 };
    const plain = await settle(alignment());
    const fp = await plain.frame({ size: tall, settled: true });
    const hp = fp.layout.sections[0].tracks.at(-1)!.box.height;
    expect(hp).toBeLessThanOrEqual(DefaultRenderConfig.rows.maxPx);
    expect(hp).toBeGreaterThan(240);
    const filled = await settle(alignment('fill'));
    const ff = await filled.frame({ size: tall, settled: true });
    expect(ff.layout.sections[0].tracks.at(-1)!.box.height).toBeGreaterThan(hp);
  });
});

describe('tracks of one kind share only when close', () => {
  it('a track’s own height still overrides, and a far short neighbour fits its own lanes', async () => {
    const s = await settle(spec([40, 1], { height: 100 }));
    const f = await s.frame({ size: Size, settled: true, viewports: { [Section]: { start: 0, end: 200 } } });
    expect(f.layout.sections[0].tracks.map((t) => t.box.height)).toEqual([100, lanesHeight(1, DefaultRenderConfig)]);
  });
});
