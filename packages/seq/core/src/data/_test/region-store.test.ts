/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { Coverage, GridDense, RegionAnswer, RegionRequest } from '../protocol.js';
import { openSeq } from '../seq-runtime.js';
import { DataStore, type StoreSettings, type Want } from '../store.js';
import { withSummaries } from '../summaries.js';
import type { SliceClock } from '../summaries-region.js';
import { FakeClock, SyntheticSource } from './synthetic.js';
import { contactCells } from './region-fixtures.js';

// The store's 2D wants and its bounds on a layer's wants together (_spec/README.md "Data", Store): x tiles × y tiles,
// a symmetric grid's folded, the product bounded; a layer's tiles and read cost bounded; an answer's bytes bounded;
// items tiles and sided systems keyed apart; 2D wants ordered and cancelled as 1D (EC-172, EC-173).

const S = '/data/hic';
const Now: SliceClock = { now: () => 0, yield: () => Promise.resolve() };

/** A grid on chr1 (1,000,000) with itself and chr2 (500,000), symmetric; each read answered with an empty tile. */
function gridCoverage(extra: Partial<Coverage> = {}): Coverage {
  return {
    systems: [
      { name: 'chr1', length: 1_000_000 },
      { name: 'chr2', length: 500_000 },
    ],
    shapes: [
      { system: 'chr1', shape: 'grid', y: '*', symmetric: true, aggregate: 'sum' },
      { system: 'chr2', shape: 'grid', y: '*', symmetric: true, aggregate: 'sum' },
      { system: 'chr1', shape: 'pairs' },
    ],
    levels: [1, 1000],
    ...extra,
  };
}

function emptyGrid(r: RegionRequest): RegionAnswer {
  if (!r.y) return { ok: true, tables: {} };
  const L = r.level;
  const x = Math.floor(r.start / L);
  const y = Math.floor(r.y!.start / L);
  const grid: GridDense = { form: 'dense', x, y, nx: 0, ny: 0, value: new Float32Array(0) };
  return { ok: true, tables: r.y ? { grid } : {} };
}

function setup(settings: Partial<StoreSettings> = {}, source?: (clock: FakeClock) => SyntheticSource) {
  const clock = new FakeClock();
  const src = source?.(clock) ?? new SyntheticSource(gridCoverage(), clock, { answer: emptyGrid, delayMs: 10 });
  const store = new DataStore({ clock, settings: { summaries: false, ...settings } });
  store.attach(S, src);
  return { clock, source: src, store };
}

const grid = (x: [number, number], y: [number, number], extra: Partial<Want> = {}): Want => ({
  source: S,
  system: 'chr1',
  start: x[0],
  end: x[1],
  level: 1000,
  y: { system: 'chr1', start: y[0], end: y[1] },
  ...extra,
});

describe('2D wants', () => {
  it('cover their x tiles × y tiles of 256 bins a side, a symmetric grid folded to the upper triangle', async () => {
    const { clock, source, store } = setup();
    // 3 × 3 tiles of 256,000 positions at level 1,000; the 3 below the diagonal are their mirrors' keys
    const w = grid([0, 768_000], [0, 768_000]);
    const tiles = store.tilesOf(w);
    expect(tiles.length).toBe(9);
    expect(tiles.filter((t) => t.transposed).length).toBe(3);
    expect(store.keysOf(w).length).toBe(6);
    store.want([w]);
    await clock.advance(20);
    expect(source.reads.length).toBe(6);
    for (const r of source.reads) {
      expect(r.request.shapes).toEqual(['grid']);
      expect(r.request.start).toBeLessThanOrEqual(r.request.y!.start);
    }
  });

  it('clip a tile to its system: a chromosome shorter than a tile is one small tile', () => {
    const { store } = setup();
    const t = store.tilesOf(grid([0, 1_000_000], [0, 500_000], { y: { system: 'chr2', start: 0, end: 500_000 } }));
    // x: 4 tiles (the last clipped), y: 2 (the last clipped at 500,000)
    expect(t.length).toBe(8);
    const last = t[t.length - 1].region;
    expect(last.end).toBe(1_000_000);
    expect(last.y!.end).toBe(500_000);
  });

  it('are bounded by the product of their tiles (maxTilesPerWant)', () => {
    const { store } = setup({ maxTilesPerWant: 8 });
    const w = grid([0, 768_000], [0, 768_000]);
    expect(store.tooMany(w)).toBe(true);
    expect(store.keysOf(w)).toEqual([]);
    store.want([w]);
    expect(store.issues.value.map((i) => i.code)).toEqual(['seq.store.too_many_tiles']);
  });

  it('never share a key with a 1D tile, and ask a 1D tile for no grid', async () => {
    const { clock, source, store } = setup();
    const one: Want = { source: S, system: 'chr1', start: 0, end: 256_000, level: 1000 };
    expect(store.keysOf(one)).not.toContain(store.keysOf(grid([0, 256_000], [0, 256_000]))[0]);
    store.want([one]);
    await clock.advance(20);
    expect(source.reads[0].request.shapes).toEqual(['pairs']);
  });

  it('are ordered and cancelled as 1D wants (EC-172, EC-173)', async () => {
    const { clock, source, store } = setup({ maxInFlight: 2 });
    store.want([
      grid([0, 256_000], [0, 256_000], { priority: 'margin' }),
      grid([256_000, 512_000], [256_000, 512_000]),
    ]);
    // the view's tile goes first
    expect(source.reads[0].request.start).toBe(256_000);
    store.want([]);
    await clock.advance(400);
    expect(source.reads.every((r) => r.aborted || r.answered)).toBe(true);
    expect(store.settled.value).toBe(true);
  });
});

describe("a layer's wants together", () => {
  it('past maxTilesPerLayer are not asked, and say so', () => {
    const { store } = setup({ maxTilesPerLayer: 3 });
    const wants = [
      grid([0, 256_000], [0, 512_000], { layer: 'hic' }),
      grid([256_000, 768_000], [0, 256_000], { layer: 'hic' }),
    ];
    // the second want's first tile is the first's mirror: three tiles together
    expect(store.layerLoad(wants)).toMatchObject({ tiles: 3, fits: true });
    const more = [...wants, grid([768_000, 1_000_000], [0, 256_000], { layer: 'hic' })];
    expect(store.layerLoad(more).fits).toBe(false);
    store.want(more);
    expect(store.keysOf(more[0]).length).toBeGreaterThan(0);
    expect(store.issues.value[0].message).toMatch(/layer's wanted regions cover 4 tiles/);
  });

  it('past maxBytesPerLayer of read cost are not asked: the layer reads coarser', async () => {
    const clock = new FakeClock();
    const src = new SyntheticSource(gridCoverage(), clock, { answer: emptyGrid, delayMs: 10 });
    // a tile costs what its index says: 6 MiB at level 1,000, 1 MiB at 10,000
    const costly = Object.assign(src, { readCost: (r: RegionRequest) => (r.level === 1000 ? 6 : 1) * 1024 * 1024 });
    const store = new DataStore({ clock, settings: { summaries: false } });
    store.attach(S, costly);
    const fine = [grid([0, 768_000], [0, 256_000], { layer: 'vcf' })];
    expect(store.layerLoad(fine)).toMatchObject({ tiles: 3, cost: 18 * 1024 * 1024, fits: false });
    store.want(fine);
    await clock.advance(20);
    expect(src.reads.length).toBe(0);
    expect(store.issues.value[0].message).toMatch(/cost/);
  });
});

describe('bytes and keys', () => {
  it('refuse an answer past maxTileBytes as too large, never kept', async () => {
    const big = (r: RegionRequest): RegionAnswer => {
      const n = 256;
      const grid: GridDense = {
        form: 'dense',
        x: Math.floor(r.start / 1000),
        y: Math.floor(r.y!.start / 1000),
        nx: n,
        ny: n,
        value: new Float32Array(n * n),
      };
      return { ok: true, tables: { grid } };
    };
    const { clock, store } = setup(
      { maxTileBytes: 100_000 },
      (c) => new SyntheticSource(gridCoverage(), c, { answer: big, delayMs: 10 }),
    );
    const w = grid([0, 256_000], [0, 256_000]);
    store.want([w]);
    await clock.advance(20);
    const e = store.get(store.keysOf(w)[0]);
    expect(e?.status === 'failed' && e.kind).toBe('too-large');
    expect(store.bytes).toBe(0);
  });

  it('read a sided system on its side, keyed apart', async () => {
    const clock = new FakeClock();
    const cov: Coverage = {
      systems: [
        { name: 'chr1', length: 10_000, side: 'a' },
        { name: 'chr1', length: 20_000, side: 'b' },
      ],
      shapes: [{ system: 'chr1', shape: 'pairs' }],
      levels: [1],
    };
    const src = new SyntheticSource(cov, clock, { answer: () => ({ ok: true, tables: {} }), delayMs: 1 });
    const store = new DataStore({ clock, settings: { summaries: false } });
    store.attach(S, src);
    const a: Want = { source: S, system: 'chr1', side: 'a', start: 0, end: 50_000, level: 1 };
    const b: Want = { ...a, side: 'b' };
    // each side's length clips its want
    expect(store.tilesOf(a).length).toBe(10);
    expect(store.tilesOf(b).length).toBe(20);
    expect(store.keysOf(a)[0]).not.toBe(store.keysOf(b)[0]);
    store.want([a]);
    await clock.advance(5);
    expect(src.reads[0].request.side).toBe('a');
  });

  it('keep an items tile apart from a bins tile at one level', () => {
    const { store } = setup();
    const w: Want = { source: S, system: 'chr1', start: 0, end: 1_000_000, level: 1000 };
    expect(store.keysOf({ ...w, items: true })[0]).not.toBe(store.keysOf(w)[0]);
  });
});

describe('the whole-assembly level (Juicebox’s “All”)', () => {
  it('is cut from the blocks’ planes in one want, a symmetric pair mirrored', async () => {
    const n1 = 64;
    const n2 = 40;
    const cov: Coverage = {
      systems: [
        { name: 'chr1', length: n1 },
        { name: 'chr2', length: n2 },
        { name: 'all', length: n1 + n2 },
      ],
      shapes: [
        {
          system: 'chr1',
          shape: 'grid',
          y: ['chr1', 'chr2'],
          symmetric: true,
          aggregate: 'sum',
          whole: { system: 'all', sequences: ['chr1', 'chr2'] },
        },
        { system: 'chr2', shape: 'grid', y: ['chr2'], symmetric: true, aggregate: 'sum' },
        { system: 'all', shape: 'grid', y: ['all'], symmetric: true, aggregate: 'sum' },
      ],
      levels: [1],
    };
    const dense = (nx: number, ny: number, diagonal: boolean): GridDense => ({
      form: 'dense',
      x: 0,
      y: 0,
      nx,
      ny,
      value: contactCells(nx, ny, diagonal),
    });
    const pair = (x: string, nx: number, y: string, ny: number, diagonal: boolean) => ({
      request: {
        system: x,
        start: 0,
        end: nx,
        level: 1,
        shapes: ['grid' as const],
        y: { system: y, start: 0, end: ny },
      },
      answer: { ok: true as const, tables: { grid: dense(nx, ny, diagonal) } },
    });
    const src = withSummaries(
      openSeq({
        coverage: cov,
        pairs: [
          pair('chr1', n1, 'chr1', n1, true),
          pair('chr1', n1, 'chr2', n2, false),
          pair('chr2', n2, 'chr2', n2, true),
        ],
      }),
      undefined,
      undefined,
      Now,
    );
    const store = new DataStore({ clock: new FakeClock() });
    store.attach(S, src);
    const w: Want = {
      source: S,
      system: 'all',
      start: 0,
      end: n1 + n2,
      level: 8,
      y: { system: 'all', start: 0, end: n1 + n2 },
    };
    expect(store.tilesOf(w).length).toBe(1);
    store.want([w]);
    await store.whenSettled();
    const e = store.get(store.keysOf(w)[0]);
    expect(e?.status).toBe('loaded');
    const g = e?.status === 'loaded' ? (e.tables.grid as GridDense) : undefined;
    const total = (v: Float32Array) => v.reduce((a, x) => a + (Number.isNaN(x) ? 0 : x), 0);
    // every count once, chr1 × chr2 twice (both triangles)
    const sum =
      total(contactCells(n1, n1, true)) + 2 * total(contactCells(n1, n2, false)) + total(contactCells(n2, n2, true));
    expect(total(g!.value)).toBe(sum);
    expect(g!.nx).toBe(Math.ceil((n1 + n2) / 8));
  });
});
