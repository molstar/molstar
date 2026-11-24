/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Coverage, DataStore, DefaultStoreSettings, type OpenedSource, type Want } from '@molstar/seqx-core/data';
import { describe, expect, it } from 'vitest';
import { DefaultRenderConfig } from '../config.js';
import { gridLevel, itemsGate, type LayerSource, pairsLevels, readGridTiles } from '../layer-data.js';

// The level rules by shape for pairs and grids read by region (_spec/README.md "Render", Data): a grid's finest level
// at a device pixel, pairs' bins at a pixel's level and items at a tile size, the items gate's hysteresis, the 2D read
// with its margin ring and its pending rects (EC-169, EC-170, EC-171).

const reader = { settings: DefaultStoreSettings };

describe('pairs levels', () => {
  it('read items over a 10 Mb view within maxTilesPerWant tiles, bins at a couple of pixels (EC-170)', () => {
    const span = 10_000_000;
    // 1,000 device px across: 10,000 positions a pixel
    const { bins, items } = pairsLevels([1, 1024, 65_536], span / 1000, span, reader, DefaultRenderConfig.pairs.binPx);
    expect(items).toBe(1024);
    expect(Math.ceil(span / (DefaultStoreSettings.tileBins * items)) + 1).toBeLessThanOrEqual(
      DefaultStoreSettings.maxTilesPerWant,
    );
    expect(Math.ceil(span / (DefaultStoreSettings.tileBins * items))).toBe(10);
    expect(bins).toBe(1024);
  });

  it('keep items across a pan at the threshold, never flipping (EC-171: bins drawn instead past it)', () => {
    const { maxArcs, itemsHysteresis } = DefaultRenderConfig.pairs;
    // a pan's frames: bin pairs in view around the gate, the items loaded once asked
    const counts = [4000, 4090, 4100, 4200, 4150, 4096, 4300, 5000, 5200, 4096, 4000];
    let showing = false;
    const shown: boolean[] = [];
    for (const n of counts) {
      showing = itemsGate(showing, n, showing ? n : undefined, maxArcs, itemsHysteresis);
      shown.push(showing);
    }
    // asked at 4,000, kept to 5,120 (1.25 × 4,096), let go past it, asked again under the gate
    expect(shown).toEqual([true, true, true, true, true, true, true, true, false, true, true]);
    // never asked past the gate from bins: no "load anyway"
    expect(itemsGate(false, 4097, undefined, maxArcs, itemsHysteresis)).toBe(false);
  });
});

describe('grid levels', () => {
  it('take the finest level whose bin is at least a device pixel on both axes', () => {
    const levels = [1, 2, 4, 8, 16, 32, 64];
    expect(gridLevel(levels, { x: 3, y: 1 }, { x: 3000, y: 1000 }, reader)).toBe(4);
    expect(gridLevel(levels, { x: 0.5, y: 0.5 }, { x: 500, y: 500 }, reader)).toBe(1);
    // past the store's bound on x tiles × y tiles: coarser
    const small = { settings: { ...DefaultStoreSettings, maxTilesPerWant: 4 } };
    expect(gridLevel(levels, { x: 1, y: 1 }, { x: 1024, y: 1024 }, small)).toBe(4);
  });
});

/** A one-system symmetric grid source of `length` positions answering any tile with no cells. */
function gridSource(length: number): OpenedSource {
  const coverage: Coverage = {
    systems: [{ name: 'chr1', length }],
    shapes: [{ system: 'chr1', shape: 'grid', symmetric: true, aggregate: 'sum' }],
    levels: [1, 2, 4, 8],
  };
  return { coverage, read: () => new Promise(() => {}), dispose() {} };
}

const Src: LayerSource = { source: '/data/hic', system: 'chr1', dataSystem: '/data/chr1', shape: 'grid', ratio: 1 };

describe('a region grid’s read', () => {
  it('wants the view and a ring of tiles round it, the view’s tiles pending until they land (EC-169)', () => {
    const store = new DataStore({ settings: { summaries: false } });
    store.attach('/data/hic', gridSource(1_000_000));
    const wants: Want[] = [];
    const read = readGridTiles(
      { want: (w) => wants.push(w) },
      store,
      Src,
      [1000, 3000],
      { system: 'chr1', start: 600, end: 1100 },
      2,
    );
    expect(wants.map((w) => w.priority)).toEqual(['view', 'margin']);
    // x tiles 1000–3000 at level 2 (512 a tile): 1..5; y 600–1100: 1..2; nothing asked yet: every rect pending
    expect(read.entries).toEqual([]);
    expect(read.pending.length / 4).toBe(10);
    // across the diagonal: the tiles below it are their mirrors
    const below = readGridTiles(
      { want: () => {} },
      store,
      Src,
      [0, 1024],
      { system: 'chr1', start: 1024, end: 2048 },
      4,
      0,
    );
    const above = readGridTiles(
      { want: () => {} },
      store,
      Src,
      [1024, 2048],
      { system: 'chr1', start: 0, end: 1024 },
      4,
      0,
    );
    expect(below.pending).toEqual([0, 1024, 1024, 2048]);
    expect(above.pending).toEqual([1024, 2048, 0, 1024]);
    expect(
      store.tilesOf({ ...wants[0], start: 1024, end: 2048, y: { system: 'chr1', start: 0, end: 1024 } })[0].transposed,
    ).toBe(true);
  });

  it('holds the reference view’s 2D bytes with its margin ring under 44 MB (counted)', () => {
    // a 1440 × 900 plot at DPR 2, a cell a device pixel: 2,880 × 1,800 cells at level 1, tiles of 256² dense floats
    const store = new DataStore({ settings: { summaries: false } });
    store.attach('/data/hic', gridSource(10_000_000));
    const wants: Want[] = [];
    readGridTiles(
      { want: (w) => wants.push(w) },
      store,
      Src,
      [100_100, 102_980],
      { system: 'chr1', start: 500_100, end: 501_900 },
      1,
    );
    const keys = new Set<string>();
    for (const w of wants) for (const k of store.keysOf(w)) keys.add(k);
    const view = store.keysOf(wants[0]).length;
    expect(view).toBeLessThanOrEqual(13 * 9);
    const bytes = keys.size * 256 * 256 * 4;
    expect(bytes).toBeLessThan(44 * 1024 * 1024);
    expect(bytes).toBeLessThan(store.settings.maxBytes);
  });
});
