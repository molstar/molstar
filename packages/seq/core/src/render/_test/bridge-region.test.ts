/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A bridge over pairs read by region (_spec/README.md "Render", Bridges): each side read at its own level, each pair
// drawn once by ownership; an alignment's blocks under a pixel apart one ribbon; pending tiles the hatch (EC-169).

import { delayedSource, type FormatRuntime, openSeq, SeqRuntime } from '@molstar/seqx-core/data';
import { describe, expect, it, vi } from 'vitest';
import { hitView, TracksLayout } from '../index.js';
import { RibbonsType } from '../pairs-ribbons.js';
import type { RibbonsLay } from '../pairs-ribbons.js';
import { groupsOf } from '../pairs-ribbons-lay.js';
import {
  type Alignment,
  AlnLength,
  alignmentList,
  alignmentSource,
  randomAlignments,
  regionBridgeSpec,
} from './bridge-region-fixtures.js';
import { manualClock } from './bridge-sides-fixtures.js';
import { settle } from './helpers.js';

const Size = { width: 720, height: 500 };
const layOf = (f: { layout: { bridges?: { list: readonly { lays: readonly { lay: unknown }[] }[] } } }) =>
  f.layout.bridges?.list[0]?.lays[0]?.lay as RibbonsLay;
const windows = (top: [number, number], bottom: [number, number]) => ({
  size: Size,
  viewports: {
    '/views/v/top': { start: top[0], end: top[1] },
    '/views/v/bottom': { start: bottom[0], end: bottom[1] },
  },
});

/** How many drawn ribbons cover each block: a ribbon of its alignment over it, or a bin pair holding both its ends. */
function coverage(lay: RibbonsLay, alignments: readonly Alignment[]): number[] {
  const list = alignmentList(alignments);
  const n = list.aStart.length;
  const out = new Array<number>(n).fill(0);
  for (let i = 0; i < lay.n; i++) {
    const prep = lay.reads[lay.src[i]];
    const d = prep.data!;
    const k = lay.pair[i];
    const h = lay.hull.subarray(4 * i, 4 * i + 4);
    const binned = prep.level !== undefined;
    const group = groupsOf(d)?.[k] ?? -1;
    for (let j = 0; j < n; j++) {
      const a = list.aStart[j];
      const b = list.bStart[j];
      if (!(a >= h[0] && a < h[1] && b >= h[2] && b < h[3])) continue;
      if (binned ? d.bStrand[k] !== list.bStrand![j] : group !== list.row![j]) continue;
      out[j]++;
    }
  }
  return out;
}

describe('a bridge over pairs read by region', () => {
  it('reads each side at its own level and draws each pair once, by the finer side’s ownership (EC-209)', async () => {
    const alignments = randomAlignments(500, 3);
    const s = await settle(regionBridgeSpec(alignments), [SeqRuntime], { pairs: { maxArcs: 1024, binPx: 40 } });
    const f = await s.frame(windows([400_000, 450_000], [0, AlnLength]));
    const lay = layOf(f);
    expect(lay.reads.length).toBe(2);
    // the zoomed side on items, the whole side on bins
    expect(lay.reads[0].level).toBeUndefined();
    expect(lay.reads[1].level).toBe(16_384);
    expect(new Set(lay.src.subarray(0, lay.n))).toEqual(new Set([0, 1]));
    // the whole of B in view: every block drawn exactly once, by its items or its bin
    const c = coverage(lay, alignments);
    expect(c.filter((x) => x !== 1)).toEqual([]);
  });

  it('draws nothing twice when both sides are on items, the window of `a` its own (EC-209)', async () => {
    const alignments = randomAlignments(200, 5);
    const s = await settle(regionBridgeSpec(alignments), [SeqRuntime]);
    const f = await s.frame(windows([0, 600_000], [300_000, AlnLength]));
    const lay = layOf(f);
    expect(lay.reads.map((r) => r.level)).toEqual([undefined, undefined]);
    const list = alignmentList(alignments);
    const c = coverage(lay, alignments);
    for (let j = 0; j < c.length; j++) {
      const inA = list.aStart[j] < 600_000;
      const inB = list.bStart[j] >= 300_000;
      expect(c[j]).toBe(inA || inB ? 1 : 0);
    }
  });
});

describe('an alignment’s blocks in a bridge', () => {
  // three blocks 100 apart, then one 50 kb on: at the fit a px is ~1.4 kb
  const one: Alignment[] = [
    {
      a: 100_000,
      b: 200_000,
      blocks: [
        [500, 100, 100],
        [500, 100, 100],
        [500, 50_000, 50_000],
        [500, 0, 0],
      ],
    },
  ];

  it('merges consecutive blocks under a device px apart into one ribbon, their hull, its hit the run', async () => {
    const s = await settle(regionBridgeSpec(one), [SeqRuntime]);
    const f = await s.frame(windows([0, AlnLength], [0, AlnLength]));
    const lay = layOf(f);
    expect(lay.n).toBe(2);
    const m = Array.from(lay.merged.subarray(0, 2)).indexOf(1);
    expect(m).toBeGreaterThanOrEqual(0);
    expect(Array.from(lay.hull.subarray(4 * m, 4 * m + 4))).toEqual([100_000, 101_700, 200_000, 201_700]);
    // its hit is the drawn run: both hulls, on A and on B
    const r = lay.rib.subarray(6 * m, 6 * m + 6);
    const x = (r[0] + r[1] + r[3] + r[4]) / 4;
    const h = hitView(f.layout, x, (r[2] + r[5]) / 2);
    expect(h?.bridge).toBe(0);
    expect(h?.extent?.parts.map((p) => Array.from(p.ranges))).toEqual([
      [100_000, 101_700],
      [200_000, 201_700],
    ]);
  });

  it('parts the ribbon where a gap reaches a px: zoomed in, each block its own', async () => {
    const s = await settle(regionBridgeSpec(one), [SeqRuntime]);
    const lay = layOf(await s.frame(windows([99_000, 103_000], [199_000, 203_000])));
    expect(lay.n).toBe(3);
    expect(Array.from(lay.merged.subarray(0, 3))).toEqual([0, 0, 0]);
  });

  it('EC-169: a side’s tiles pending draw the hatch along its edge, never an empty band', async () => {
    const s = await settle(regionBridgeSpec(randomAlignments(50, 7)), [SeqRuntime]);
    const input = windows([0, AlnLength], [0, AlnLength]);
    const lay = layOf(s.renderer.frame(input));
    expect(lay.n).toBe(0);
    const ys = new Set<number>();
    for (let i = 0; i + 3 < lay.hatch.length; i += 4) ys.add(lay.hatch[i + 2]);
    expect(ys).toEqual(new Set([lay.edgeA, lay.edgeB]));
    const settled = layOf(await s.frame(input));
    expect(settled.n).toBeGreaterThan(0);
    expect(settled.hatch).toEqual([]);
  });
});

/** Every block anchored in either window drawn exactly once, none twice. */
function expectDrawnOnce(lay: RibbonsLay, alignments: readonly Alignment[], top: number[], bottom: number[]) {
  const list = alignmentList(alignments);
  const c = coverage(lay, alignments);
  const wrong: string[] = [];
  for (let j = 0; j < c.length; j++) {
    const shown =
      (list.aStart[j] >= top[0] && list.aStart[j] < top[1]) ||
      (list.bStart[j] >= bottom[0] && list.bStart[j] < bottom[1]);
    if (c[j] > 1 || (shown && c[j] !== 1)) wrong.push(`${list.aStart[j]}|${list.bStart[j]} ×${c[j]}`);
  }
  expect(wrong.slice(0, 5)).toEqual([]);
}

describe('the near side’s `R`, read whole or not at all (EC-209)', () => {
  const alignments = randomAlignments(500, 3);
  const config = { pairs: { maxArcs: 1024, binPx: 40 } };
  const top = [400_000, 450_000];
  const bottom = [0, AlnLength];

  it('past `maxTilesPerWant`: `R` never asked, the window shrunk to the far grid owned, no hole, no double, no flicker', async () => {
    // the window's items tiles (50) fit, `R`'s (64) do not
    const s = await settle(regionBridgeSpec(alignments), [SeqRuntime], config, {}, { maxTilesPerWant: 60 });
    const at = windows(top as [number, number], bottom as [number, number]);
    const first = layOf(await s.frame(at));
    expect(first.reads.map((r) => r.level)).toEqual([undefined, 16_384]);
    expect([...first.ownR.values()].filter((r) => r.length)).toEqual([[409_600, 442_368]]);
    expectDrawnOnce(first, alignments, top, bottom);
    const n = first.n;
    for (let i = 0; i < 10; i++) {
      const lay = layOf(await s.frame(at));
      expect([lay.n, lay.readSide, s.store.layerRefused('/views/v/br/aln')]).toEqual([n, [0, 1], undefined]);
    }
  });

  it('the near side landed, its widening not: no hole, no double; then `R` whole', async () => {
    const clock = manualClock();
    const slow: FormatRuntime = {
      name: 'slow-widening',
      formats: ['slow-widening'],
      // the window's items tiles at once, the rest (the widening, the margins) on the clock
      open: async () =>
        delayedSource(openSeq(alignmentSource(alignments)), {
          clock,
          delayMs: (r) => (r.items && (r.start < 399_360 || r.end > 450_560) ? 1 : 0),
        }),
    };
    const s = await settle(regionBridgeSpec(alignments, { format: 'slow-widening' }), [slow], config);
    const at = windows(top as [number, number], bottom as [number, number]);
    for (let i = 0; i < 6; i++) {
      s.store.want(s.renderer.frame(at).wants);
      await new Promise((r) => setTimeout(r, 0));
    }
    const landed = layOf(s.renderer.frame(at));
    expect(landed.reads.map((r) => r.level)).toEqual([undefined, 16_384]);
    expect([...landed.ownR.values()].filter((r) => r.length)).toEqual([[409_600, 442_368]]);
    expectDrawnOnce(landed, alignments, top, bottom);
    clock.flush();
    const whole = layOf(await s.frame(at));
    expect([...whole.ownR.values()].filter((r) => r.length)).toEqual([[393_216, 458_752]]);
    expectDrawnOnce(whole, alignments, top, bottom);
  });
});

describe('the cull, a side at a time', () => {
  it('with the default `maxArcs`, a side at the fit past it never pushes the zoomed side’s items out', async () => {
    // thousands of bin pairs at the fit; in the zoomed window, besides the random blocks, 50 one-position blocks (the
    // least aligned length: culled first when the two sides are culled together)
    const short: Alignment[] = Array.from({ length: 50 }, (_, i) => ({
      a: 401_000 + 300 * i,
      b: 7 * i * 1000,
      blocks: [[1, 10, 10]],
    }));
    const alignments = [...randomAlignments(12_000, 12), ...short];
    const s = await settle(regionBridgeSpec(alignments), [SeqRuntime], { pairs: { binPx: 4 } });
    const top = [400_000, 420_000];
    const lay = layOf(await s.frame(windows(top as [number, number], [0, AlnLength])));
    expect(lay.reads[0].level).toBeUndefined();
    expect(lay.densityN).toBeGreaterThan(0);
    // every block anchored in the zoomed window drawn, by its item
    const list = alignmentList(alignments);
    const c = coverage(lay, alignments);
    for (let j = 0; j < c.length; j++) if (list.aStart[j] >= top[0] && list.aStart[j] < top[1]) expect(c[j]).toBe(1);
  });
});

describe('a bridge over a collapsed axis', () => {
  it('EC-33: a side of 150,000 single-position pieces owns a few merged ranges, each pair drawn once (EC-209)', async () => {
    const alignments = randomAlignments(500, 4);
    const s = await settle(regionBridgeSpec(alignments, { pieces: 150_000 }), [SeqRuntime], {
      pairs: { maxArcs: 4096 },
    });
    // the top's whole axis (A's first 300 kb, every other position), the bottom zoomed in: the top the far side
    const lay = layOf(
      await s.frame({ size: Size, viewports: { '/views/v/bottom': { start: 500_000, end: 520_000 } } }),
    );
    expect(lay.reads.length).toBe(2);
    const ranges = [...lay.ownR.values(), ...lay.ownFar.values()];
    expect(Math.max(...ranges.map((r) => r.length))).toBeLessThanOrEqual(4);
    expect(coverage(lay, alignments).filter((x) => x > 1)).toEqual([]);
  });
});

describe('a bridge’s pending hatch while a side reads items', () => {
  it('EC-169: a side on items panned to tiles not loaded hatches its own edge alone, then its items again', async () => {
    const s = await settle(regionBridgeSpec(randomAlignments(300, 9)), [SeqRuntime]);
    const zoomed = windows([400_000, 450_000], [400_000, 450_000]);
    const before = layOf(await s.frame(zoomed));
    expect(before.reads.map((r) => r.level)).toEqual([undefined, undefined]);
    expect(before.hatch).toEqual([]);
    // the top panned far: its tiles pending, the bottom's loaded
    const panned = windows([800_000, 850_000], [400_000, 450_000]);
    const lay = layOf(s.renderer.frame(panned));
    const ys = new Set<number>();
    for (let i = 0; i + 3 < lay.hatch.length; i += 4) ys.add(lay.hatch[i + 2]);
    expect(ys).toEqual(new Set([lay.edgeA]));
    // the bottom still draws what it reads, each once
    expect(lay.readSide).toEqual([1]);
    expect(lay.n).toBeGreaterThan(0);
    const after = layOf(await s.frame(panned));
    expect(after.hatch).toEqual([]);
    expect(after.reads.map((r) => r.level)).toEqual([undefined, undefined]);
  });
});

describe('a bridge over a region source, kept', () => {
  it('a pan of a third section lays out that section alone; a pan of a side lays out the bridge again', async () => {
    const s = await settle(regionBridgeSpec(randomAlignments(200, 7), { third: true }), [SeqRuntime]);
    const at = (third: number, top = 0) => ({
      size: Size,
      viewports: {
        '/views/v/top': { start: top, end: AlnLength / 4 + top },
        '/views/v/bottom': { start: 0, end: AlnLength / 4 },
        '/views/v/third': { start: third, end: AlnLength / 2 + third },
      },
    });
    await s.frame(at(0));
    await s.frame(at(0));
    expect(layOf(await s.frame(at(0))).n).toBeGreaterThan(0);
    const place = vi.spyOn(TracksLayout, 'place');
    const ribbons = vi.spyOn(RibbonsType, 'layout');
    try {
      s.renderer.frame(at(1000));
      expect(place.mock.calls.map((c) => c[0].address)).toEqual(['/views/v/third']);
      expect(ribbons).not.toHaveBeenCalled();
      s.renderer.frame(at(1000, 500));
      expect(ribbons).toHaveBeenCalledTimes(1);
    } finally {
      place.mockRestore();
      ribbons.mockRestore();
    }
  });
});
