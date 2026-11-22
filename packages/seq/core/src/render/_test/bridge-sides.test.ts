/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A bridge between two genomes (_spec/README.md "Render", Bridges: reads): each side reads its own blocks at its own
// level (a sided source on the end it shows, `chr1` on both genomes), each pair drawn exactly once by ownership over
// the sides' systems, the near side's `R` read at its own level (EC-209).

import { SeqRuntime } from '@molstar/seqx-core/data';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { beforeAll, describe, expect, it } from 'vitest';
import { loadBridges } from '../bridge-load.js';
import { regionState } from '../pairs-region.js';
import type { RibbonsLay } from '../pairs-ribbons.js';
import type { BridgeRibbon } from '../tracks-types.js';
import {
  BottomSection,
  ChrLength,
  type SideBlock,
  sideBlocks,
  sidesBridgeSpec,
  TopSection,
} from './bridge-sides-fixtures.js';
import { type Settled, settle } from './helpers.js';

const Size = { width: 720, height: 500 };

beforeAll(async () => {
  await loadGenome();
  await loadBridges();
});

type Window = readonly [number, number];
const Fit: Window = [0, 2 * ChrLength];

const input = (top: Window, bottom: Window) => ({
  size: Size,
  viewports: {
    [TopSection]: { start: top[0], end: top[1] },
    [BottomSection]: { start: bottom[0], end: bottom[1] },
  },
});

/** A position of a chromosome in a window of a genome section (`chr1` then `chr2`, end to end). */
const inWindow = (w: Window, chr: string, p: number) => {
  const u = (chr === 'chr2' ? ChrLength : 0) + p;
  return u >= w[0] && u < w[1];
};

const chrOf = (system: string) => /:name (chr\d)\)/.exec(system)?.[1] ?? system;

async function frameOf(s: Settled, top: Window, bottom: Window) {
  const f = await s.frame(input(top, bottom));
  const bridge = f.layout.bridges!.list[0];
  return { f, ribbons: bridge.ribbons(), lay: bridge.lays[0].lay as RibbonsLay };
}

/** How many drawn ribbons cover each block: its item, or a bin pair of its strand holding both its starts. */
function coverage(ribbons: readonly BridgeRibbon[], blocks: readonly SideBlock[]): number[] {
  const out = new Array<number>(blocks.length).fill(0);
  for (const r of ribbons) {
    const a = chrOf(r.aSystem);
    const b = chrOf(r.bSystem);
    for (let j = 0; j < blocks.length; j++) {
      const x = blocks[j];
      if (x.a !== a || x.b !== b || r.reversed !== x.strand < 0) continue;
      const hit =
        r.bin === undefined
          ? r.aStart === x.as && r.bStart === x.bs
          : x.as >= r.aStart && x.as < r.aEnd && x.bs >= r.bStart && x.bs < r.bEnd;
      if (hit) out[j]++;
    }
  }
  return out;
}

/** Every block anchored in either window drawn exactly once, and none twice. */
function expectDrawnOnce(ribbons: readonly BridgeRibbon[], blocks: readonly SideBlock[], top: Window, bottom: Window) {
  const c = coverage(ribbons, blocks);
  const wrong: string[] = [`${top} ${bottom}`];
  for (let j = 0; j < blocks.length; j++) {
    const x = blocks[j];
    const shown = inWindow(top, x.a, x.as) || inWindow(bottom, x.b, x.bs);
    if (c[j] > 1 || (shown && c[j] !== 1)) wrong.push(`${x.a}:${x.as} ${x.b}:${x.bs} ×${c[j]}`);
  }
  expect(wrong.length > 1 ? wrong.slice(0, 6) : []).toEqual([]);
}

describe('a bridge between two genomes, each side reading its own blocks', () => {
  // about 1,200 blocks: past maxArcs at the fit (bins), items zoomed in, never culled
  const blocks = sideBlocks(400, 3);
  const config = { pairs: { maxArcs: 1000, binPx: 40, itemsHysteresis: 1 } };
  const Zoomed: Window = [1_300_000, 1_340_000];
  const Middle: Window = [600_000, 1_400_000];
  const windows: Window[] = [Fit, Middle, Zoomed, [100_000, 140_000]];

  it('a sided source read on the end each side shows: chr1 of both genomes apart (EC-209)', async () => {
    const s = await settle(sidesBridgeSpec(blocks), [SeqRuntime], { pairs: { maxArcs: 4096 } });
    const { ribbons, lay } = await frameOf(s, Fit, Fit);
    expect(lay.reads.length).toBe(2);
    expect(ribbons.length).toBe(blocks.length);
    expectDrawnOnce(ribbons, blocks, Fit, Fit);
    // each end on its own genome's sequence
    expect(new Set(ribbons.map((r) => r.aSystem.slice(0, 8)))).toEqual(new Set(['/data/P/']));
    expect(new Set(ribbons.map((r) => r.bSystem.slice(0, 8)))).toEqual(new Set(['/data/G/']));
  });

  it('each pair drawn exactly once at every pair of windows and levels, items and bins, overlapping or not (EC-209)', async () => {
    const s = await settle(sidesBridgeSpec(blocks), [SeqRuntime], config);
    const levels = new Set<string>();
    for (const top of windows)
      for (const bottom of windows) {
        const { ribbons, lay } = await frameOf(s, top, bottom);
        expect(lay.reads.length).toBe(2);
        levels.add(lay.reads.map((r) => r.level ?? 'items').join('/'));
        expect(lay.densityN, JSON.stringify([top, bottom, lay.n, lay.reads.map((r) => [r.level, r.data!.n])])).toBe(0);
        expectDrawnOnce(ribbons, blocks, top, bottom);
      }
    // both sides on bins, both on items, and one of each, both ways
    expect(levels).toContain('16384/16384');
    expect(levels).toContain('items/items');
    expect(levels).toContain('items/16384');
    expect(levels).toContain('16384/items');
  });

  it('a zoomed-in side owns its window widened to the far bins, read at its own level (EC-209)', async () => {
    const s = await settle(sidesBridgeSpec(blocks), [SeqRuntime], config);
    const { lay } = await frameOf(s, Fit, Zoomed);
    // the near side (the bottom, on items) reads its `R`: chr2 from 288 to 344 kb (the 16 kb grid)
    expect(lay.reads[1].level).toBeUndefined();
    expect(lay.reads[0].level).toBe(16_384);
    const r = [...lay.ownR].filter(([, x]) => x.length);
    expect(r.map(([k, x]) => [chrOf(k), ...x])).toEqual([['chr2', 294_912, 344_064]]);
  });

  it('a region read per side equals the whole read clipped to that side’s window', async () => {
    const s = await settle(sidesBridgeSpec(blocks), [SeqRuntime], { pairs: { maxArcs: 4096 } });
    const top: Window = [200_000, 1_500_000];
    const bottom: Window = [900_000, 1_100_000];
    const { lay } = await frameOf(s, top, bottom);
    for (const [i, w, end] of [
      [0, top, 'a'],
      [1, bottom, 'b'],
    ] as const) {
      const prep = regionState(lay.sides[i]).prep!;
      const d = prep.data!;
      const got = new Set<string>();
      for (let k = 0; k < d.n; k++) {
        const sys = chrOf((i ? prep.endUnits![d.bUnit![k]] : prep.endUnits![d.aUnit![k]]) ?? '');
        const p = i ? d.bStart[k] : d.aStart[k];
        if (inWindow(w, sys, p)) got.add(`${chrOf(prep.endUnits![d.aUnit![k]])}:${d.aStart[k]}|${d.bStart[k]}`);
      }
      const want = new Set(
        blocks
          .filter((x) => inWindow(w, end === 'a' ? x.a : x.b, end === 'a' ? x.as : x.bs))
          .map((x) => `${x.a}:${x.as}|${x.bs}`),
      );
      expect(got).toEqual(want);
    }
  });
});
