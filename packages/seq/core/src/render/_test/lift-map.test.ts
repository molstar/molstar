/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A map read by region (_spec/README.md "Resolve", Maps read by region; EC-212): its runs a sided pairs source's rows,
// read by either end (`RegionRequest.side`), one read for the ribbons and the liftover.

import { openSeq, SeqRuntime } from '@molstar/seqx-core/data';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { beforeAll, describe, expect, it } from 'vitest';
import { loadBridges } from '../bridge-load.js';
import { settle } from './helpers.js';
import { LiftMap, liftSource, liftSpec, Runs } from './lift-fixtures.js';

beforeAll(async () => {
  await loadGenome();
  await loadBridges();
});

const Size = { width: 900, height: 400 };
const chrOf = (system: string) => /:name (chr\d)\)/.exec(system)?.[1] ?? system;

describe('a map read by region (`map.source`)', () => {
  it('its runs read by `to` equal the whole read clipped to the window (EC-212)', async () => {
    const src = openSeq(liftSource());
    const signal = new AbortController().signal;
    const read = async (side: 'a' | 'b', system: string, start: number, end: number) => {
      const a = await src.read({ system, start, end, level: 1, shapes: ['pairs'], side, items: true }, signal);
      if (!a.ok || a.tables.pairs?.form !== 'list') throw new Error('not read');
      return Array.from(a.tables.pairs.row!).sort((x, y) => x - y);
    };
    for (const [start, end] of [
      [0, 80_000],
      [15_000, 45_000],
      [61_000, 62_000],
    ])
      for (const chr of ['chr1', 'chr2']) {
        const whole = Runs.flatMap((r, i) => (r.b === chr && r.to < end && r.to + r.length > start ? [i] : []));
        expect(await read('b', chr, start, end)).toEqual(whole);
      }
    // by `from`: the duplicated stretch's two runs
    expect(await read('a', 'chr1', 72_000, 78_000)).toEqual([3, 4]);
  });

  it('is a pairs annotation over both assemblies, its params read', async () => {
    const s = await settle(liftSpec(), [SeqRuntime]);
    const m = s.model.annotations.get(LiftMap)!;
    expect(m.shape).toBe('pairs');
    expect(m.bySystem?.get('/data/new/(sequence :name chr1)')).toMatchObject({ system: 'chr1', side: 'b' });
    expect(m.bySystem?.get('/data/old/(sequence :name chr1)')).toMatchObject({ system: 'chr1', side: 'a' });
    expect(s.model.issues.filter((i) => i.code === 'seq.resolve.param.unread')).toEqual([]);
    expect(s.model.maps.has(LiftMap)).toBe(false);
  });

  it('is drawn as a bridge’s ribbons, a run each; unaligned stretches gaps (EC-212)', async () => {
    const s = await settle(liftSpec(), [SeqRuntime]);
    const f = await s.frame({ size: Size });
    const ribbons = f.layout.bridges!.list[0].ribbons();
    const got = ribbons
      .filter((r) => r.bin === undefined)
      .map((r) => [chrOf(r.aSystem), r.aStart, r.aEnd, chrOf(r.bSystem), r.bStart, r.bEnd, r.reversed])
      .sort();
    const want = Runs.map((r) => [r.a, r.from, r.from + r.length, r.b, r.to, r.to + r.length, r.strand < 0]).sort();
    expect(got).toEqual(want);
  });
});
