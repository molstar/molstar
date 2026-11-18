/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The values layer's placement (_spec/README.md "Render", Layer types and Configuration): data on the section's own
// system placed directly gives the lists the general mapping gives (`values.direct`), over an axis of two pieces (a
// cyclic section wrapping its origin), a map with gaps (a protein's values on its gene's exons), any view and level;
// and the level a layer reads is at most `values.coarser` device pixels' positions a bin.

import { describe, expect, it } from 'vitest';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import type { Node } from '@molstar/core/repr/tree';
import { AxisTransform, DefaultRenderConfig, type FrameInput, recordDrawLists } from '../index.js';
import { layerLevel, type LayerSource } from '../layer-data.js';
import type { StoreReader } from '../layer-type.js';
import { BIG, BigLevels, BigRuntime, bigSpec, conservation, GENE_EXONS, geneLetters, P53, seeded } from './fixtures.js';
import { blocks, settle } from './helpers.js';

const CYCLIC = 6000;

/** A cyclic system's values (a few missing) on a section wrapping its origin and on a whole one; p53's on its gene. */
function spec(): Node {
  const b = createSeqSpec({ title: 'values placement' });
  const r = seeded(7);
  const letters = Array.from({ length: CYCLIC }, () => 'ACGT'[Math.floor(r() * 4)]).join('');
  const ring = b.data.sequence('ring', { alphabet: 'dna', letters, cyclic: true });
  const signal = Float32Array.from({ length: CYCLIC }, (_, i) => (i % 97 === 0 ? Number.NaN : r() * 2 - 0.5));
  const v = b.data.values('signal', { target: ring, values: signal });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  const gene = b.data.sequence('tp53', { alphabet: 'dna', letters: geneLetters() });
  b.data.map('p53-on-gene', { from: p53, to: gene, method: 'codon', exons: GENE_EXONS, strand: '-' });
  const c = b.data.values('conservation', { target: p53, values: conservation() });
  b.views.view('main', (view) => {
    view
      .section('wrap', { of: ring, range: [5000, 1000] })
      .track('v')
      .layer({ type: 'values', data: v });
    view.section('whole', { of: ring }).track('v').layer({ type: 'values', data: v });
    const protein = view.section('protein', { of: p53 });
    protein.track('c').layer({ type: 'values', data: c });
    view.section('gene', { of: gene, follow: protein }).track('c').layer({ type: 'values', data: c });
  });
  return b.build();
}

/** Seeded views: a window anywhere (fractional edges), a width, a device pixel ratio. */
function views(seed: number, count: number): FrameInput[] {
  const r = seeded(seed);
  const out: FrameInput[] = [];
  for (let i = 0; i < count; i++) {
    const w = r() < 0.3 ? 1 + r() * 60 : 50 + r() * 6000;
    const start = r() * 2000 - w / 4;
    out.push({
      size: { width: 300 + Math.floor(r() * 1300), height: 200 },
      dpr: [1, 1.5, 2][Math.floor(r() * 3)],
      viewport: { start, end: start + w },
      settled: r() < 0.5,
    });
  }
  return out;
}

describe('values: placed directly or through the mapping, the same lists', () => {
  for (const coarser of [1, 2]) {
    it(`over pieces, a map's gaps, views and levels (values.coarser ${coarser})`, async () => {
      const s = spec();
      const direct = await settle(s, [], { values: { coarser, direct: true } });
      const general = await settle(s, [], { values: { coarser, direct: false } });
      let bars = 0;
      let crossing = 0;
      for (const input of views(coarser, 16)) {
        const f = await direct.frame(input);
        const a = recordDrawLists(f.lists, ['data']);
        const b = recordDrawLists((await general.frame(input)).lists, ['data']);
        expect(a).toBe(b);
        // the values' rects (a block a track: each track holds one values layer)
        for (const block of blocks(f.lists)) bars += block.rects;
        // a view across the wrap: the section's two pieces (its origin at display unit 1,000) both in view
        const axis = f.layout.sections[0].axis!;
        if (axis.u0 < axis.offsets[1] && axis.u1 > axis.offsets[1]) crossing++;
      }
      expect(bars).toBeGreaterThan(16 * 8);
      expect(crossing).toBeGreaterThanOrEqual(4);
    });
  }

  it('at the edges of the view: a window starting and ending inside a position', async () => {
    const s = spec();
    const direct = await settle(s, [], { values: { direct: true } });
    const general = await settle(s, [], { values: { direct: false } });
    for (const [start, end] of [
      [0.5, 40.25],
      [999.5, 1003.5],
      [1990.75, 2000],
      [-3, 7.5],
      [0, 2000],
    ]) {
      const input: FrameInput = { size: { width: 700, height: 200 }, dpr: 2, viewport: { start, end } };
      const a = recordDrawLists((await direct.frame(input)).lists, ['data']);
      expect(a).toBe(recordDrawLists((await general.frame(input)).lists, ['data']));
    }
  });
});

describe('values.coarser: the level a values layer reads', () => {
  const reader = (levels: number[]): StoreReader => ({
    coverage: () => ({ systems: [], shapes: [], levels }),
    keysOf: () => [],
    get: () => undefined,
  });
  const src: LayerSource = { source: 's', system: 'chr', dataSystem: 'chr', shape: 'values', ratio: 1 };
  // 1,304 CSS px at dpr 2: 2,608 device columns
  const axisOver = (positions: number) =>
    new AxisTransform('chr', [0, 4 * BIG], 0, 1304, 2, { start: 0, end: positions });

  it('1: the coarsest bin no wider than a device pixel’s positions; 2: up to two pixels’', () => {
    const data = reader(BigLevels);
    // 100,000 positions: 38 a device pixel
    expect(layerLevel(src, data, axisOver(100_000), false, 1)).toBe(1);
    expect(layerLevel(src, data, axisOver(100_000), false, 2)).toBe(64);
    // 160,000: 61 a pixel, 122 two
    expect(layerLevel(src, data, axisOver(160_000), false, 1)).toBe(1);
    expect(layerLevel(src, data, axisOver(160_000), false, 2)).toBe(64);
    // 200,000: 77 a pixel
    expect(layerLevel(src, data, axisOver(200_000), false, 1)).toBe(64);
    expect(layerLevel(src, data, axisOver(200_000), false, 2)).toBe(64);
    // 6,000,000 positions in view: 2,300 a pixel, 4,600 two
    expect(layerLevel(src, data, axisOver(6_000_000), false, 1)).toBe(64);
    expect(layerLevel(src, data, axisOver(6_000_000), false, 2)).toBe(4096);
    // the finest whatever the factor
    expect(layerLevel(src, data, axisOver(200_000), true, 2)).toBe(1);
  });

  it('the default is 1: a values layer over 100,000 positions wants positions, at 2 it wants the 64 level', async () => {
    expect(DefaultRenderConfig.values.coarser).toBe(1);
    const input: FrameInput = {
      size: { width: 1440, height: 400 },
      dpr: 2,
      viewport: { start: BIG / 2 - 50_000, end: BIG / 2 + 50_000 },
    };
    const levels = async (coarser?: number) => {
      const s = await settle(bigSpec(), [BigRuntime], coarser === undefined ? {} : { values: { coarser } });
      const f = s.renderer.frame(input);
      return new Set(f.wants.filter((w) => w.priority === 'view').map((w) => w.level));
    };
    // the ranges layer reads a pixel's level (positions) in both
    expect([...(await levels())]).toEqual([1]);
    expect([...(await levels(2))].sort((a, b) => a - b)).toEqual([1, 64]);
  });
});
