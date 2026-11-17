/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Sticky lanes while data streams: an item keeps its lane when a neighboring tile arrives, though a fresh packing
// would move it; lanes are repacked only once the view has settled and a packing needs fewer lanes.

import { categoryColumn } from '@molstar/core/repr/tree';
import { DataStore, type OpenedSource, type RegionAnswer, type RegionRequest } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { resolveViews } from '@molstar/seqx-core/view';
import { describe, expect, it } from 'vitest';
import { flush } from '../../data/_test/synthetic.js';
import { bundledTextMetrics, LightTokens, type RangesLay, Renderer, StickyLanes } from '../index.js';

/** A source whose answers the test releases tile by tile. */
class Deferred implements OpenedSource {
  readonly coverage = {
    systems: [{ name: 'chr', length: 4096 }],
    shapes: [{ system: 'chr', shape: 'ranges' as const }],
    levels: [1],
  };
  private readonly waiting = new Map<number, (a: RegionAnswer) => void>();

  constructor(private readonly items: Record<number, { start: number; end: number; id: string }[]>) {}

  read(r: RegionRequest): Promise<RegionAnswer> {
    return new Promise((resolve) => this.waiting.set(r.start, resolve));
  }

  release(start: number) {
    const rows = this.items[start] ?? [];
    this.waiting.get(start)?.({
      ok: true,
      tables: {
        ranges: {
          form: 'list',
          start: Int32Array.from(rows, (r) => r.start),
          end: Int32Array.from(rows, (r) => r.end),
          id: categoryColumn(rows.map((r) => r.id)),
        },
      },
    });
  }

  dispose() {}
}

function spec() {
  const b = createSeqSpec({ title: 'streaming ranges' });
  const src = b.data.source('src', { format: 'deferred', uri: 'deferred:chr' });
  const chr = b.data.sequence('chr', { alphabet: 'dna', source: src });
  const items = b.data.ranges('items', { target: chr, source: src });
  b.views.view('main', (v) => v.section('chr', { of: chr }).track('items').layer({ type: 'ranges', data: items }));
  return b.build();
}

function lanesOf(lay: RangesLay): Record<string, number> {
  const out: Record<string, number> = {};
  for (let i = 0; i < lay.keys.length; i++) out[String(lay.keys[i])] = lay.lanes[i];
  return out;
}

describe('sticky lanes', () => {
  const tile0 = [
    { start: 950, end: 990, id: 'a' },
    { start: 1000, end: 1060, id: 'long' },
  ];
  const tile1 = [
    { start: 1030, end: 1050, id: 'b' },
    { start: 1100, end: 1150, id: 'c' },
  ];

  it('keep an item’s lane when a neighboring tile arrives', async () => {
    const model = resolveSpec(spec());
    const view = resolveViews(model)[0];
    const store = new DataStore();
    const source = new Deferred({ 0: tile0, 1024: tile1 });
    store.attach('/data/src', source);
    const renderer = new Renderer({ model, view, data: store, tokens: LightTokens, text: bundledTextMetrics() });
    const input = { size: { width: 800, height: 80 }, viewport: { start: 900, end: 1200 } };
    store.want(renderer.frame(input).wants);
    await flush();

    source.release(1024);
    await flush();
    const first = renderer.frame(input);
    const lay = first.layout.sections[0].tracks[0].layers[0].lay as RangesLay;
    expect(lanesOf(lay)).toEqual({ b: 0, c: 0 });

    // the left neighbor arrives: `long` starts first, yet `b` keeps lane 0, settled or not
    source.release(0);
    await flush();
    renderer.frame(input);
    expect(lanesOf(lay)).toEqual({ a: 0, long: 1, b: 0, c: 0 });
    renderer.frame({ ...input, settled: true });
    expect(lanesOf(lay)).toEqual({ a: 0, long: 1, b: 0, c: 0 });

    // a renderer with no history packs them fresh: `long` first
    const fresh = new Renderer({ model, view, data: store, tokens: LightTokens, text: bundledTextMetrics() });
    const f = fresh.frame(input);
    expect(lanesOf(f.layout.sections[0].tracks[0].layers[0].lay as RangesLay)).toEqual({ a: 0, long: 0, b: 1, c: 0 });
    store.retire();
  });

  it('repack only when settled and sparse', () => {
    const lanes = new StickyLanes();
    const out = new Int32Array(2);
    // two overlapping items take two lanes
    expect(lanes.assign(['x', 'y'], [0, 5], [10, 15], [0, 1], 3, false, out)).toBe(2);
    expect([...out]).toEqual([0, 1]);
    // `x` leaves the view: `y` keeps lane 1 while the view moves
    const one = new Int32Array(1);
    expect(lanes.assign(['y'], [5], [15], [0], 3, false, one)).toBe(2);
    expect(one[0]).toBe(1);
    // settled: the lanes are sparse (one item, two lanes), so they are packed again
    expect(lanes.assign(['y'], [5], [15], [0], 3, true, one)).toBe(1);
    expect(one[0]).toBe(0);
    expect(lanes.laneOf('y')).toBe(0);
  });

  it('an item that fits no lane overflows (−1)', () => {
    const lanes = new StickyLanes();
    const out = new Int32Array(3);
    expect(lanes.assign(['a', 'b', 'c'], [0, 1, 2], [10, 10, 10], [0, 1, 2], 2, false, out)).toBe(2);
    expect([...out]).toEqual([0, 1, -1]);
  });
});
