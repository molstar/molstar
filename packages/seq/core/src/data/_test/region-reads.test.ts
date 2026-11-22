/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { itWithin } from '@molstar/test-harness/timing';
import { describe, expect, it } from 'vitest';
import { dataKey, foldTile } from '../keys.js';
import type { GridDense, OpenedSource, PairsBins, PairsList, RegionRequest, Tables } from '../protocol.js';
import { decodeSeq, encodeSeq } from '../seq-format.js';
import { openSeq } from '../seq-runtime.js';
import { coarsenBins, withSummaries } from '../summaries.js';
import { binPairs } from '../summaries-pairs.js';
import type { SliceClock } from '../summaries-region.js';
import { transposeGrid } from '../summaries-grid.js';
import { answerBytes, checkAnswer, checkTables, joinTables, sliceTable } from '../tables.js';
import { checkRegionReads } from './conformance.js';
import {
  Chr1,
  genotypes,
  leveledPairsSource,
  matrixSource,
  oneLevelGrid,
  pairsOf,
  pairsSource,
  rowsSource,
  sidedPairsSource,
  transcripts,
  twoLevelGrid,
} from './region-fixtures.js';

// The region shapes (_spec/README.md "Data": Protocol, Tables, Summaries): their keys, checks, cuts and joins, the
// `seq` format's round trips, the pairs level and the grid levels, and every source's region reads against its whole
// read clipped.

const signal = new AbortController().signal;
const Now: SliceClock = { now: () => 0, yield: () => Promise.resolve() };

/** A bin pair's key, for comparing bins as sets. */
function binKeys(b: PairsBins): Map<string, number> {
  const out = new Map<string, number>();
  for (let i = 0; i < b.count.length; i++) {
    const k = `${b.aSystem?.get(i)}:${b.aBin[i]}|${b.bSystem?.get(i)}:${b.bBin[i]}|${b.strand?.[i]}|${b.category?.get(i)}`;
    out.set(k, (out.get(k) ?? 0) + b.count[i]);
  }
  return out;
}

describe('keys', () => {
  it('keep a 2D key, a sided system and an items tile apart from the 1D keys', () => {
    const one = dataKey('s', 'chr1', 16, 3);
    expect(dataKey('s', 'chr1', 16, 3, { y: { system: 'chr2', tile: 0 } })).not.toBe(one);
    expect(dataKey('s', 'chr1', 16, 3, { side: 'a' })).toBe('s|a:chr1|16|3');
    // an unsided system named like a sided part never meets it: `:` is escaped
    expect(dataKey('s', 'a:chr1', 16, 3)).not.toBe(dataKey('s', 'chr1', 16, 3, { side: 'a' }));
    expect(dataKey('s', 'chr1', 16, 3, { items: true })).toBe('s|chr1|i16|3');
  });

  it('fold a symmetric tile below the diagonal to the one served', () => {
    const order = (s: string) => (s === 'chr1' ? 0 : 1);
    expect(foldTile({ system: 'chr1', tile: 3 }, { system: 'chr1', tile: 1 }, order)).toEqual({
      x: { system: 'chr1', tile: 1 },
      y: { system: 'chr1', tile: 3 },
      transposed: true,
    });
    expect(foldTile({ system: 'chr2', tile: 0 }, { system: 'chr1', tile: 5 }, order).transposed).toBe(true);
    expect(foldTile({ system: 'chr1', tile: 1 }, { system: 'chr1', tile: 1 }, order).transposed).toBe(false);
  });
});

describe('tables of the region shapes', () => {
  it('check their forms and columns', () => {
    const p = pairsOf(10);
    expect(checkTables({ pairs: p })).toBeUndefined();
    expect(checkTables({ pairs: binPairs(p, 100) })).toBeUndefined();
    expect(checkTables({ pairs: { ...p, aEnd: new Int32Array(3) } })).toBe(
      'the pairs table: columns of unequal lengths',
    );
    // a ranges list is no pairs table, though both are `list`s
    expect(checkTables({ pairs: transcripts() } as unknown as Tables)).toMatch(/no form of its shape/);
    const dense: GridDense = { form: 'dense', x: 0, y: 0, nx: 2, ny: 3, value: new Float32Array(6) };
    expect(checkTables({ grid: dense })).toBeUndefined();
    expect(checkTables({ grid: { ...dense, nx: 3 } })).toMatch(/nx × ny/);
    expect(checkTables({ matrix: genotypes().matrix })).toBeUndefined();
    expect(checkTables({ rows: (rowsSource().pairs[0].answer as { tables: Tables }).tables.rows })).toBeUndefined();
  });

  it('weigh a dense tile at 4 bytes a cell and a sparse one at 12', () => {
    const dense: GridDense = { form: 'dense', x: 0, y: 0, nx: 256, ny: 256, value: new Float32Array(256 * 256) };
    expect(answerBytes({ grid: dense }, undefined)).toBe(256 * 1024);
    const sparse = {
      form: 'sparse' as const,
      x: 0,
      y: 0,
      xBin: new Int32Array(10),
      yBin: new Int32Array(10),
      value: new Float32Array(10),
      fill: 0,
    };
    expect(answerBytes({ grid: sparse }, undefined)).toBe(120);
  });

  it('refuse a grid answer outside its x and y regions, and bins asked as items', () => {
    const r: RegionRequest = {
      system: 'chr1',
      start: 0,
      end: 64,
      level: 4,
      shapes: ['grid'],
      y: { system: 'chr1', start: 0, end: 32 },
    };
    const tile = (ny: number): GridDense => ({
      form: 'dense',
      x: 0,
      y: 0,
      nx: 16,
      ny,
      value: new Float32Array(16 * ny),
    });
    expect(checkAnswer({ ok: true, tables: { grid: tile(8) } }, r)).toBeUndefined();
    expect(checkAnswer({ ok: true, tables: { grid: tile(9) } }, r)).toMatch(/outside the region/);
    const p: RegionRequest = { system: 'chr1', start: 0, end: 1000, level: 16, shapes: ['pairs'], items: true };
    expect(checkAnswer({ ok: true, tables: { pairs: binPairs(pairsOf(5), 16) } }, p)).toMatch(/bins asked as items/);
  });

  it('cut pairs by either end and join the cuts into each pair once (EC-62, EC-208)', () => {
    const p = pairsOf(200);
    const cut = { start: 1000, end: 3000, system: 'chr1' };
    const whole = sliceTable(p, 1, cut) as PairsList;
    // every pair kept has an end on chr1 meeting the region
    for (let i = 0; i < whole.aStart.length; i++) {
      const a = whole.aSystem!.get(i) === 'chr1' && whole.aStart[i] < 3000 && whole.aEnd[i] > 1000;
      const b = whole.bSystem!.get(i) === 'chr1' && whole.bStart[i] < 3000 && whole.bEnd[i] > 1000;
      expect(a || b).toBe(true);
    }
    const parts = [1000, 1500, 2200, 3000];
    const joined = joinTables(
      parts.slice(1).map((to, i) => sliceTable(p, 1, { ...cut, from: parts[i], to })),
      1,
    ) as PairsList;
    expect([...joined.id!.codes].map((c) => joined.id!.dict[c]).sort()).toEqual(
      [...whole.id!.codes].map((c) => whole.id!.dict[c]).sort(),
    );
  });

  it('emit a span row for an item with no piece in the cut, so a region read and the whole clipped agree (EC-229)', () => {
    const t = transcripts();
    // inside t1's intron: no exon, its extent over the region
    const inIntron = sliceTable(t, 1, { start: 1000, end: 2000 });
    expect([...inIntron.start]).toEqual([1000]);
    expect([...inIntron.end]).toEqual([1000]);
    expect(inIntron.itemStart?.[0]).toBe(100);
    expect(inIntron.itemEnd?.[0]).toBe(3100);
    // a region over an exon of t1: no span row for it; t2 by its piece
    const over = sliceTable(t, 1, { start: 150, end: 450 });
    expect([...over.start]).toEqual([100, 400]);
    // joined cuts hold the span row once
    const joined = joinTables(
      [
        sliceTable(t, 1, { start: 1000, end: 2000, from: 1000, to: 1500 }),
        sliceTable(t, 1, { start: 1000, end: 2000, from: 1500, to: 2000 }),
      ],
      1,
    );
    expect(joined).toEqual(inIntron);
  });
});

describe('the seq format', () => {
  it('round-trips pairs, bins, sided systems, grids, a matrix and keyed rows', () => {
    for (const data of [
      pairsSource(30),
      leveledPairsSource(30),
      sidedPairsSource(9),
      oneLevelGrid(20, 12),
      matrixSource(),
      rowsSource(),
    ]) {
      const back = decodeSeq(encodeSeq(data));
      expect(back.issues).toEqual([]);
      expect(back.data).toEqual(data);
    }
  });

  it('refuses a system listed twice on one side, and takes one listed once a side', () => {
    const sided = sidedPairsSource(3);
    expect(decodeSeq(encodeSeq(sided)).issues).toEqual([]);
    const twice = {
      ...sided,
      coverage: { ...sided.coverage, systems: [sided.coverage.systems[0], sided.coverage.systems[0]] },
    };
    expect(decodeSeq(encodeSeq(twice)).issues[0]?.message).toMatch(/listed twice/);
  });
});

describe('region reads equal the whole read clipped, at every level and shape (conformance)', () => {
  it('pairs with their own levels, found from either end (EC-62, EC-208)', async () => {
    await checkRegionReads(openSeq(leveledPairsSource(300)), 64);
  });

  it('pairs at one level with the summaries’ levels', async () => {
    await checkRegionReads(withSummaries(openSeq(pairsSource(300)), undefined, undefined, Now), 64);
  });

  it('sided pairs, each side read on its own', async () => {
    await checkRegionReads(openSeq(sidedPairsSource()), 64);
  });

  it('a grid with two levels, and its one-level twin with the summaries’ (EC-223)', async () => {
    await checkRegionReads(openSeq(twoLevelGrid(96, 64)), 8);
    await checkRegionReads(withSummaries(openSeq(oneLevelGrid(96, 64)), undefined, undefined, Now), 8);
  });

  it('a matrix and keyed rows (EC-220, EC-224)', async () => {
    await checkRegionReads(openSeq(matrixSource()), 256);
    await checkRegionReads(openSeq(rowsSource()), 8);
  });
});

describe('the pairs level (EC-208)', () => {
  it('counts sum to the items at every level, and per bin pair to the items inside it', () => {
    const p = pairsOf(2000);
    let below = binPairs(p, 16, { category: true });
    expect(below.count.reduce((a, b) => a + b, 0)).toBe(2000);
    for (const level of [256, 4096]) {
      const coarse = coarsenBins(below, level / (level === 256 ? 16 : 256));
      expect(coarse.count.reduce((a, b) => a + b, 0)).toBe(2000);
      // each bin pair's count is the items inside it, binned straight from the items
      expect(binKeys(coarse)).toEqual(binKeys(binPairs(p, level, { category: true })));
      below = coarse;
    }
  });

  it('keeps straight and reversed bin pairs apart, and sums a field (EC-206)', () => {
    const p = pairsOf(500);
    const b = binPairs(p, 1_000_000, { field: 'score' });
    // one bin a system at this width: at most two systems squared by two strands
    expect(b.count.length).toBeLessThanOrEqual(8);
    expect(new Set(b.strand)).toEqual(new Set([1, -1]));
    const total = (p.fields![0].column as Float32Array).reduce((a, v) => a + v, 0);
    expect(b.sum!.reduce((a, v) => a + v, 0)).toBeCloseTo(total, 0);
  });

  it('binned from held pairs equals binned by the store’s summaries from region reads, at every kept level', async () => {
    const N = 400;
    const s = withSummaries(openSeq(pairsSource(N)), undefined, undefined, Now);
    const level = s.coverage.levels[s.coverage.levels.length - 1];
    // the pyramid built by the first read: every kept level after it equals the held pairs binned there
    await s.read({ system: 'chr1', start: 0, end: Chr1, level, shapes: ['pairs'] }, signal);
    const kept = s.coverage.levels.filter((l) => l > 1);
    expect(kept.length).toBeGreaterThan(0);
    for (const l of kept) {
      const r = await s.read({ system: 'chr1', start: 0, end: Chr1, level: l, shapes: ['pairs'] }, signal);
      const held = sliceTable(binPairs(pairsOf(N), l, { category: true, systems: true }), l, {
        start: 0,
        end: Chr1,
        system: 'chr1',
      });
      expect(binKeys((r.ok ? r.tables.pairs : undefined) as PairsBins)).toEqual(binKeys(held as PairsBins));
    }
    const a = await s.read({ system: 'chr1', start: 0, end: Chr1, level, shapes: ['pairs'] }, signal);
    const fromStore = a.ok ? a.tables.pairs : undefined;
    expect(fromStore?.form).toBe('bins');
    const held = sliceTable(binPairs(pairsOf(N), level, { category: true, systems: true }), level, {
      start: 0,
      end: Chr1,
      system: 'chr1',
    });
    expect(binKeys(fromStore as PairsBins)).toEqual(binKeys(held as PairsBins));
    // the shape's numeric field summed into the bins (Decided 2), as binned from the held pairs
    const total = (c: Float32Array | undefined) => (c ?? new Float32Array(0)).reduce((a, v) => a + v, 0);
    const heldSum = sliceTable(binPairs(pairsOf(N), level, { category: true, systems: true, field: 'score' }), level, {
      start: 0,
      end: Chr1,
      system: 'chr1',
    }) as PairsBins;
    expect((fromStore as PairsBins).sum).toBeDefined();
    expect(total((fromStore as PairsBins).sum)).toBeCloseTo(total(heldSum.sum), 2);
    // a level whose bins do not halve the rows below leaves the coverage once built (the reader takes the next kept)
    expect(s.coverage.levels).not.toContain(16);
    expect(s.coverage.levels).toContain(level);
    // asked meanwhile, it is answered as items
    const fine = await s.read({ system: 'chr1', start: 0, end: Chr1, level: 16, shapes: ['pairs'] }, signal);
    expect(fine.ok && fine.tables.pairs?.form).toBe('list');
  });

  it('past the summaries’ bytes, refuses a coarse tile before reading unless its read cost fits', async () => {
    const reads: RegionRequest[] = [];
    const counted = (cost?: number): OpenedSource => {
      const inner = openSeq(pairsSource(400));
      return {
        coverage: inner.coverage,
        read: (r, sig) => {
          reads.push(r);
          return inner.read(r, sig);
        },
        ...(cost !== undefined && { readCost: () => cost }),
        dispose: () => inner.dispose(),
      };
    };
    const tiny = { factor: 16, maxBytes: 100, sliceMs: 4 };
    const tile: RegionRequest = { system: 'chr1', start: 0, end: 4096, level: 256, shapes: ['pairs'] };
    const blind = withSummaries(counted(), tiny, undefined, Now);
    const a = await blind.read(tile, signal);
    expect(a.ok).toBe(false);
    // only the whole reads the pyramid tried: no tile read to be refused after
    expect(reads.every((r) => r.end - r.start > 4096)).toBe(true);
    reads.length = 0;
    const cheap = withSummaries(counted(64), tiny, undefined, Now);
    const b = await cheap.read(tile, signal);
    expect(b.ok).toBe(true);
    expect(reads.some((r) => r.level === 1 && r.end === 4096)).toBe(true);
  });

  it('answers an items tile at a coarse level as the level-1 items over it', async () => {
    const s = withSummaries(openSeq(pairsSource(400)), undefined, undefined, Now);
    const r: RegionRequest = { system: 'chr1', start: 0, end: 4096, level: 256, shapes: ['pairs'], items: true };
    const a = await s.read(r, signal);
    expect(checkAnswer(a, r)).toBeUndefined();
    expect(a.ok && a.tables.pairs?.form).toBe('list');
  });
});

describe('grid levels and symmetric tiles (EC-223)', () => {
  it('a symmetric tile transposed equals the tile asked across the diagonal, every tile pair at every level', async () => {
    const s = openSeq(twoLevelGrid(64, 40));
    let pairs = 0;
    for (const level of [1, 4])
      for (let x = 0; x < 64; x += 16 * level)
        for (let y = 0; y < 64; y += 16 * level) {
          const w = 16 * level;
          const up: RegionRequest = {
            system: 'chr1',
            start: x,
            end: Math.min(64, x + w),
            level,
            shapes: ['grid'],
            y: { system: 'chr1', start: y, end: Math.min(64, y + w) },
          };
          const down: RegionRequest = {
            ...up,
            start: y,
            end: Math.min(64, y + w),
            y: { system: 'chr1', start: x, end: Math.min(64, x + w) },
          };
          const a = await s.read(up, signal);
          const b = await s.read(down, signal);
          expect(a.ok && b.ok).toBe(true);
          if (a.ok && b.ok) expect(transposeGrid(a.tables.grid!)).toEqual(b.tables.grid);
          pairs++;
        }
    // 4 × 4 tile pairs at level 1 (the diagonal and both sides of it), one at level 4
    expect(pairs).toBe(17);
  });

  it('sums counts at factor 2 a side, every level summing to the whole', async () => {
    const s = withSummaries(openSeq(oneLevelGrid(64, 40)), undefined, undefined, Now);
    expect(s.coverage.levels.slice(0, 4)).toEqual([1, 2, 4, 8]);
    const sum = (g: GridDense) => g.value.reduce((a, v) => a + (Number.isNaN(v) ? 0 : v), 0);
    const y = { system: 'chr1', start: 0, end: 64 };
    const one = await s.read({ system: 'chr1', start: 0, end: 64, level: 1, shapes: ['grid'], y }, signal);
    const total = one.ok ? sum(one.tables.grid as GridDense) : -1;
    for (const level of [2, 4, 8]) {
      const a = await s.read({ system: 'chr1', start: 0, end: 64, level, shapes: ['grid'], y }, signal);
      expect(a.ok && sum(a.tables.grid as GridDense)).toBe(total);
    }
  });
});

describe('the pairs pyramid’s cost (performance)', () => {
  /**
   * A one-level source of `n` pairs answering each system whole at once (its lists made before: the bound is the
   * pyramid's, not the source's), and a read of chr1 at the coarsest level (which builds the pyramid).
   */
  function pyramid(n: number, clock: SliceClock = Now) {
    // no ids or types: 26 bytes a pair, so 600,000 are within the summaries' 16 MB
    const data = pairsSource(n, true);
    const lists = new Map(data.pairs.map((p) => [p.request.system, p.answer]));
    const inner: OpenedSource = {
      coverage: data.coverage,
      read: (r) => Promise.resolve(lists.get(r.system)!),
      dispose() {},
    };
    const s = withSummaries(inner, undefined, undefined, clock);
    const level = s.coverage.levels[s.coverage.levels.length - 1];
    return () => s.read({ system: 'chr1', start: 0, end: Chr1, level, shapes: ['pairs'] }, signal);
  }

  itWithin('of 12,000 pairs is built under 10 ms', async (within) => {
    // the code warmed first
    for (let i = 0; i < 3; i++) await pyramid(12_000)();
    const read = pyramid(12_000);
    const a = await within(10, read);
    expect(a.ok && a.tables.pairs?.form).toBe('bins');
  });

  itWithin(
    'of 600,000 pairs (the summaries’ bytes) is built in slices, none a long task (merge and bins alike)',
    async (within) => {
      await pyramid(20_000)();
      // a clock of real time whose yields record each slice's length (the work between two yields)
      let last = performance.now();
      let longest = 0;
      const clock: SliceClock = {
        now: () => performance.now(),
        yield: () => {
          const t = performance.now();
          longest = Math.max(longest, t - last);
          return Promise.resolve().then(() => {
            last = performance.now();
          });
        },
      };
      const read = pyramid(600_000, clock);
      // from the read on (the fixture's pairs are made before)
      last = performance.now();
      const a = await read();
      longest = Math.max(longest, performance.now() - last);
      expect(a.ok && a.tables.pairs?.form).toBe('bins');
      // the longest slice, replayed under the bound (a frame: a system's whole read checked is one task, 15 ms at
      // 300,000 rows; the merge and the bins `sliceMs` 4 and a step's grain): what the heavy lane times
      within(20, () => {
        const t0 = performance.now();
        while (performance.now() - t0 < longest);
      });
    },
    60_000,
  );
});
