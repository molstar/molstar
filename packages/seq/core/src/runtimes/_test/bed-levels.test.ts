/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import { readTabixIndex, regionBytes } from '@molstar/io/common/tabix';
import { describe, expect, it } from 'vitest';
import {
  DefaultOpenLimits,
  type OpenContext,
  type PairsBins,
  type PairsList,
  type RangesList,
  type RegionAnswer,
  type Shape,
  type SourceEntry,
} from '../../data/protocol.js';
import { sourceSuite } from '../../data/_test/source-suite.js';
import { binPairs } from '../../data/summaries-pairs.js';
import { binValues } from '../../data/summaries.js';
import { checkAnswer } from '../../data/tables.js';
import { BedRuntime, bedRuntime } from '../bed.js';
import { densityOf } from '../internal/levels.js';
import {
  generateBed6,
  generateBed12,
  generateBedGraph,
  generateBedpe,
  indexed,
  memoryAssets,
  plainRows,
} from './bed-fixtures.js';

// BED and bedGraph by region past level 1 (_spec/README.md "Data", BED family; Genome runtimes' internals, Levels): the
// levels the index's bytes allow, a tile's read cost the index's; at every level a region read equals the whole read
// clipped: items tiles row for row, a BED's density counted from the tile's items (BED12's blocks one item), a
// bedGraph's bins equal to the whole read's values binned (a later interval over an earlier).

const signal = new AbortController().signal;
const names = ['chr1', 'chr2'];

const context = (files: Record<string, Uint8Array | string>): OpenContext => ({
  assets: memoryAssets(files),
  signal,
  limits: DefaultOpenLimits,
});
const entry = (format: string, more: Partial<SourceEntry> = {}): SourceEntry => ({
  address: '/data/s',
  format,
  uri: 'f',
  ...more,
});

describe('BED and bedGraph by region, at every level', () => {
  for (const [format, kind, gen] of [
    ['bed', 'BED12', (seed: number) => generateBed12(300 + seed * 37, names, seed, 3_000_000)],
    ['bed', 'BED6 with items of no length', (seed: number) => generateBed6(300 + seed * 37, names, seed, 3_000_000)],
    ['bedgraph', 'bedGraph', (seed: number) => generateBedGraph(400 + seed * 41, names, seed, 3_000_000)],
  ] as const)
    it(`${kind}: items tiles and coarse answers equal the whole read clipped (random tiles)`, async () => {
      const rng = new PCG(5);
      const shape: Shape = format === 'bed' ? 'ranges' : 'values';
      for (let seed = 1; seed <= 3; seed++) {
        const g = gen(seed);
        const { file, index } = indexed(g, 2000);
        // the whole read's level 1 over a 4 Mb tile: a position each, past the default rows (the test's oracle only)
        const whole = await BedRuntime.open(entry(format), {
          ...context({ f: g.text }),
          limits: { ...DefaultOpenLimits, maxRows: 8_000_000 },
        });
        const region = await BedRuntime.open(entry(format, { index: 'i' }), context({ f: file, i: index }));
        expect(region.coverage.levels).toEqual([1, 16, 256, 4096]);
        const ix = await readTabixIndex(index);
        for (const level of region.coverage.levels)
          for (let n = 0; n < 5; n++) {
            const system = names[Math.floor(rng.float() * names.length)];
            const k = Math.floor((rng.float() * 3_000_000) / (level * 1024));
            const t = { system, start: k * level * 1024, end: (k + 1) * level * 1024, level, shapes: [shape] };
            expect(region.readCost!(t)).toBe(regionBytes(ix, system, Math.max(0, t.start - 1), t.end));
            const coarse = await region.read(t, signal);
            expect(checkAnswer(coarse, t)).toBeUndefined();
            if (!coarse.ok) throw new Error(`${level}: ${coarse.issue.message}`);
            if (format === 'bed') {
              const w = await whole.read({ ...t, level: 1 }, signal);
              const items = await region.read({ ...t, items: true }, signal);
              expect(checkAnswer(items, { ...t, items: true })).toBeUndefined();
              expect(plainRows(items), `${level} ${system} ${k}`).toEqual(plainRows(w));
              if (!w.ok) throw new Error('not read');
              if (level === 1) expect(plainRows(coarse)).toEqual(plainRows(w));
              else expect(coarse.tables.ranges).toEqual(densityOf(w.tables.ranges as RangesList, t));
            } else if (level === 1) {
              expect(plainRows(coarse)).toEqual(plainRows(await whole.read(t, signal)));
            } else {
              // the whole read's values of the tile (a position each), binned as the store's summaries bin them
              const w = await whole.read({ ...t, level: 1 }, signal);
              if (!w.ok || w.tables.values?.form !== 'positions') throw new Error('not read');
              const bins = binValues(w.tables.values, level, t.end, t.start);
              expect(coarse.tables.values).toEqual(bins);
            }
          }
        whole.dispose();
        region.dispose();
      }
    });

  it('a BED12 item is counted once in each bin its extent touches, its introns too', async () => {
    const text = 'chr1\t1000\t40000\tg\t0\t+\t1000\t40000\t0\t2\t100,100,\t0,38900,\n';
    const s = await BedRuntime.open(entry('bed'), context({ f: text }));
    const w = await s.read({ system: 'chr1', start: 0, end: 65_536, level: 1, shapes: ['ranges'] }, signal);
    if (!w.ok) throw new Error('not read');
    const d = densityOf(w.tables.ranges as RangesList, { start: 0, end: 65_536, level: 16_384 });
    expect([Array.from(d.bin), Array.from(d.count)]).toEqual([
      [0, 1, 2],
      [1, 1, 1],
    ]);
  });

  it('`maxLevel` caps the levels; a level not listed is not served and has no read cost', async () => {
    const g = generateBed6(200, names, 3, 3_000_000);
    const { file, index } = indexed(g, 2000);
    const s = await bedRuntime({ maxLevel: 16 }).open(entry('bed', { index: 'i' }), context({ f: file, i: index }));
    expect(s.coverage.levels).toEqual([1, 16]);
    const t = { system: 'chr1', start: 0, end: 256 * 1024, level: 256, shapes: ['ranges' as const] };
    expect(s.readCost!(t)).toBeUndefined();
    expect(await s.read(t, signal)).toMatchObject({ ok: false, kind: 'not-served' });
  });
});

describe('BEDPE through the store’s protocol (EC-208)', () => {
  it('lists `pairs` on each system and passes the source suite at every level', async () => {
    const text = generateBedpe(600, names, 4, 2_000_000);
    const s = await BedRuntime.open(entry('bedpe', { uri: undefined, value: text as never }), context({}));
    expect(s.coverage.shapes).toEqual(names.map((system) => ({ system, shape: 'pairs' })));
    expect(s.coverage.levels).toEqual([1, 16, 256, 4096]);
    await sourceSuite(BedRuntime, entry('bedpe', { uri: undefined, value: text as never }), context({}));
  });

  it('a region’s pairs from either end, the items at any level; bin pairs tile by tile, each bin pair whole', async () => {
    const text = generateBedpe(800, names, 6, 2_000_000);
    const s = await BedRuntime.open(entry('bedpe', { uri: undefined, value: text as never }), context({}));
    const all = (await s.read(
      { system: 'chr1', start: 0, end: 4_194_304, level: 4096, shapes: ['pairs'], items: true },
      signal,
    )) as RegionAnswer & { ok: true };
    const list = all.tables.pairs as PairsList;
    for (const level of [16, 256, 4096]) {
      // the bins of every tile over chr1, deduped by their whole key as a reader does: the items' bins on chr1
      const seen = new Map<string, number>();
      for (let k = 0; k * level * 1024 < 2_100_000; k++) {
        const t = {
          system: 'chr1',
          start: k * level * 1024,
          end: (k + 1) * level * 1024,
          level,
          shapes: ['pairs' as const],
        };
        const a = await s.read(t, signal);
        expect(checkAnswer(a, t)).toBeUndefined();
        const b = (a as RegionAnswer & { ok: true }).tables.pairs as PairsBins;
        for (let i = 0; i < b.count.length; i++) {
          const key = `${b.aSystem!.get(i)}|${b.bSystem!.get(i)}|${b.aBin[i]}|${b.bBin[i]}|${b.strand?.[i] ?? 0}`;
          const was = seen.get(key);
          if (was !== undefined) expect(b.count[i]).toBe(was);
          seen.set(key, b.count[i]);
        }
      }
      const whole = binPairs(list, level);
      const want = new Map<string, number>();
      for (let i = 0; i < whole.count.length; i++)
        want.set(
          `${whole.aSystem!.get(i)}|${whole.bSystem!.get(i)}|${whole.aBin[i]}|${whole.bBin[i]}|${whole.strand?.[i] ?? 0}`,
          whole.count[i],
        );
      expect(seen).toEqual(want);
    }
  });
});
