/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The dot plot's grids at scale (seqspec:format#layer-types `dot_plot`): a 1 × 1 grid's one-bin level, a small grid on
// a long block (tiles over the grid, not the block), many block pairs on a progressive backend (the coarse levels at
// once, then a budget of samples a frame), stacked lanes sharing one frame count and budget (they settle), only the
// plot's entry drawn again while tiles build, per-tile eviction past the soft cap, the pyramid reused when a patch
// changes params alone and checked against the cells by brute force, and the hovered cell's ink from the track's first
// dot plot with a grid there.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import { afterEach, describe, expect, it } from 'vitest';
import { DotPlotLimits, type DotPlotState, DotPlotVisits, dotPlotSample } from '../src/draw/dot-plot.js';
import { BackendCaps, Dirty, Pass } from '../src/draw/index.js';
import { Op, type RecordingBackend } from '../src/testing/index.js';
import { recording } from './composer-fixtures.js';
import { directed } from './conformance-fixtures.js';
import { dotPlotDoc } from './dot-plot-fixtures.js';
import { lit, passRects } from './overlay-fixtures.js';

const limits = { ...DotPlotLimits };
afterEach(() => Object.assign(DotPlotLimits, limits));

type Lit = Awaited<ReturnType<typeof lit>>;
function prepared(t: Lit, layer = 'l-plot'): DotPlotState {
  const m = (
    t.ctx.graph.assembled.value as unknown as {
      view: {
        resolved: { layers: { ref: Int32Array } };
        set: { strings: { get(i: number): string } };
        prepared: { state: unknown }[];
      };
    }
  ).view;
  const k = Array.from(m.resolved.layers.ref).findIndex((x) => m.set.strings.get(x) === layer);
  return m.prepared[k].state as DotPlotState;
}
function images(r: RecordingBackend) {
  let n = 0;
  for (let k = 0; k < r.opCount; k++) if (r.ops[4 * k] === Op.image) n++;
  return n;
}
/** One sequence of n residues and a grid on it. */
function oneSeqGrid(n: number, grid: Record<string, unknown>, stacked = false) {
  const b = createSeqViewSpec({ title: 'one' });
  const a = b.data.sequence({ ref: 'a', value: 'M'.repeat(n) });
  const g = b.data.annotation({ ref: 'g', label: 'G', geometry: 'grid', value: grid });
  b.view({ ref: 'v', label: 'V' })
    .section({ ref: 'v.s', arrangement: stacked ? 'stacked' : 'side' })
    .block({ sequence: a })
    .track({ ref: 't', label: 'P' })
    .dot_plot({ ref: 'l-plot', data: g });
  return b.getDocument() as SeqViewSpecDoc;
}
/** n chains of L residues, a grid per pair (`self`: on the diagonal only), side by side or stacked. */
function chains(n: number, L: number, o: { self?: boolean; stacked?: boolean } = {}) {
  const b = createSeqViewSpec({ title: 'chains' });
  const cs = Array.from({ length: n }, (_, k) => b.data.sequence({ ref: `c${k}`, value: 'M'.repeat(L) }));
  const cols = { target: [] as string[], b_target: [] as string[], value: [] as Float32Array[] };
  for (let p = 0; p < n; p++)
    for (let q = 0; q < n; q++) {
      if (o.self && p !== q) continue;
      cols.target.push(`c${p}`);
      cols.b_target.push(`c${q}`);
      cols.value.push(new Float32Array(L * L).map((_, x) => (x % 31) / 3));
    }
  const g = b.data.annotation({ ref: 'g', label: 'G', geometry: 'grid', value: cols });
  const sec = b.view({ ref: 'v', label: 'V' }).section({ ref: 'v.s', arrangement: o.stacked ? 'stacked' : 'side' });
  for (const c of cs) sec.block({ sequence: c });
  sec.track({ ref: 't', label: 'P' }).dot_plot({ ref: 'l-plot', data: g });
  return b.getDocument() as SeqViewSpecDoc;
}
/** Frames on a progressive backend until none asks for another (at most `max`): their count and the tiles built. */
function settle(t: Lit, be: RecordingBackend, max = 200, whole = true) {
  const before = DotPlotVisits.tiles;
  let frames = 0;
  do {
    if (whole) t.composer.redraw();
    t.frame(be, Dirty.all);
    frames++;
  } while (t.composer.again && frames < max);
  return { frames, built: DotPlotVisits.tiles - before };
}

describe('grids at scale', () => {
  it('a 1 × 1 grid has its one-bin level and draws zoomed out', async () => {
    const t = await lit(
      oneSeqGrid(2400, {
        target: ['a'],
        start: [5],
        stop: [6],
        b_target: ['a'],
        b_start: [7],
        b_stop: [8],
        value: [[3]],
      }),
      { width: 900, height: 600, dpr: 2 },
    );
    const s = prepared(t);
    expect([s.grids[0].nx, s.grids[0].ny, s.grids[0].levels?.length]).toEqual([1, 1, 1]);
    const be = t.frame(recording(), Dirty.all);
    expect(images(be)).toBe(1);
    expect(t.ctx.issues.current.value.filter((i) => i.code.startsWith('seqv.render'))).toEqual([]);
  });

  it('a small grid on a long block keeps tiles over the grid, not the block (2 Mb, zoomed in on it)', async () => {
    const value = [new Float32Array(100 * 100).fill(5)];
    const doc = oneSeqGrid(2_000_000, {
      target: ['a'],
      start: [1000],
      stop: [1100],
      b_target: ['a'],
      b_start: [1000],
      b_stop: [1100],
      value,
    });
    const t = await lit(doc, { width: 900, height: 600, dpr: 2 });
    t.vp?.setView(4, 1000 * 4, 'api');
    const before = process.memoryUsage().arrayBuffers;
    const tiles = DotPlotVisits.tiles;
    const be = t.frame(recording(BackendCaps.Progressive), Dirty.all);
    expect(DotPlotVisits.tiles - tiles).toBeGreaterThan(0);
    expect(images(be)).toBeGreaterThan(0);
    // Tile ids per level over the block would be some 61M entries at level 0 (244 MB).
    expect(process.memoryUsage().arrayBuffers - before).toBeLessThan(8e6);
  });

  it('many block pairs: the first frame builds the coarse levels only, then a budget of samples a frame', async () => {
    const t = await lit(chains(20, 135), { width: 1400, height: 900, dpr: 2 });
    expect(prepared(t).grids.length).toBe(400);
    const be = recording(BackendCaps.Progressive);
    const s0 = DotPlotVisits.samples;
    t.frame(be, Dirty.all);
    // Each block about 17 samples a side (level 3): 400 forced tiles of 17 × 17.
    expect(DotPlotVisits.samples - s0).toBeLessThanOrEqual(400 * 17 * 17 + DotPlotLimits.samplesPerFrame);
    expect(t.composer.again).toBe(true);
    let worst = 0;
    let frames = 1;
    while (t.composer.again && frames < 200) {
      const s = DotPlotVisits.samples;
      t.frame(be, Dirty.all);
      worst = Math.max(worst, DotPlotVisits.samples - s);
      frames++;
    }
    expect(t.composer.again).toBe(false);
    // A tile past what is left: at most one tile's samples over the budget.
    expect(worst).toBeLessThanOrEqual(DotPlotLimits.samplesPerFrame + 256 * 256);
  });

  /** `n` stacked lanes of `L` residues settled: the frames and tiles it took, and the steady frames after. */
  async function stacked(n: number, L: number) {
    const t = await lit(chains(n, L, { self: true, stacked: true }), { width: 1400, height: 9000, dpr: 2 });
    const be = recording(BackendCaps.Progressive);
    const { frames, built } = settle(t, be, 100);
    expect(t.composer.again).toBe(false);
    expect(frames).toBeLessThan(100);
    const after = settle(t, be, 5);
    expect([after.frames, after.built]).toEqual([1, 0]);
    return built;
  }

  it('stacked lanes share one frame count and budget: they settle, and steady frames build nothing (6 lanes)', async () => {
    expect(await stacked(6, 500)).toBeLessThan(300);
  });

  it('stacked lanes share one frame count and budget: they settle, and steady frames build nothing', {
    tags: ['slow'],
  }, async () => {
    // Before: each lane counted a frame and evicted the others' levels (2,262 tiles in 100 frames, never settled).
    expect(await stacked(12, 1000)).toBeLessThan(300);
  });

  it("while tiles build, only the plot's entry is drawn again, not the whole section", async () => {
    const t = await lit(dotPlotDoc('MKTAYIAKQR'.repeat(240), true), { width: 900, height: 600, dpr: 2 });
    DotPlotLimits.samplesPerFrame = 256 * 256;
    const be = recording(BackendCaps.Progressive);
    t.frame(be, Dirty.all);
    expect(t.composer.again).toBe(true);
    const full = t.composer.stats.fullDraws;
    const { frames } = settle(t, be, 100, false);
    expect(frames).toBeGreaterThan(1);
    expect(t.composer.again).toBe(false);
    expect(t.composer.stats.fullDraws).toBe(full);
  });

  it('past the cap, tiles not drawn in the last two frames go, one by one; a view in use keeps its own', async () => {
    const t = await lit(dotPlotDoc('MKTAYIAKQR'.repeat(240), true), { width: 900, height: 600, dpr: 2 });
    DotPlotLimits.tiles = 8;
    const be = recording();
    const vp = t.vp;
    if (!vp) throw new Error('no viewport');
    // The fit: 2,400 residues at level 1, some 25 tiles, over the cap of 8.
    const fit = vp.live.current.cell;
    const tiles = settle(t, be, 3).built;
    expect(tiles).toBeGreaterThan(8);
    // In use: no rebuilds in steady frames, though over the cap.
    expect(settle(t, be, 3).built).toBe(0);
    // Zoomed in elsewhere for three frames, then back: the fit's tiles past the cap went and are built again.
    vp.setView(8, 8 * 1200, 'api');
    for (let k = 0; k < 3; k++) settle(t, be, 3);
    vp.setView(fit, 0, 'api');
    expect(settle(t, be, 3).built).toBeGreaterThanOrEqual(tiles - 8);
  });

  it("a patch to a plot's color reuses its cells and pyramid (they key on the data, not the params)", async () => {
    const t = await lit(dotPlotDoc('MKTAYIAKQR'.repeat(20)), { width: 900, height: 600 });
    const s0 = prepared(t, 'l-dist');
    await t.ctx.applyPatch({
      format: 'seq_view_spec_patch',
      ops: [{ op: 'set', ref: 'l-dist', params: { color: { theme: 'uniform', params: { color: '#336699' } } } }],
    } as never);
    const s1 = prepared(t, 'l-dist');
    expect(s1).not.toBe(s0);
    expect([s1.grids[0].values === s0.grids[0].values, s1.grids[0].levels === s0.grids[0].levels]).toEqual([
      true,
      true,
    ]);
  });
});

describe('the pyramid', () => {
  for (const dir of [1, -1] as const)
    for (const ranged of [false, true])
      it(`is the mean of each sample's cells at every level, by brute force (direction ${dir}${ranged ? ', ranged' : ''})`, async () => {
        const n = 77;
        const rnd = (k: number) => ((k * 2654435761) % 1000) / 10;
        // a[13, 70) × a[5, 64): odd offsets and sizes, a fifth missing.
        const [s, t1, u, w] = [13, 70, 5, 64];
        const value = Array.from({ length: (t1 - s) * (w - u) }, (_, k) => (k % 5 === 2 ? null : rnd(k)));
        const grid = {
          target: ['a'],
          start: [s],
          stop: [t1],
          b_target: ['a'],
          b_start: [u],
          b_stop: [w],
          value: [value],
        };
        const t = await lit(directed(oneSeqGrid(n, grid), dir, ranged), { width: 900, height: 600 });
        const g = prepared(t).grids[0];
        const ax = t.ctx.layout.axis(t.section)?.axis;
        if (!ax) throw new Error('no axis');
        const N = ax.shown[0];
        let bad = 0;
        for (let L = 0; L <= 8; L++) {
          const step = 2 ** L;
          for (let sx = 0; sx * step < N; sx++)
            for (let sy = 0; sy * step < N; sy++) {
              let [sum, c] = [0, 0];
              for (let kx = sx * step; kx < Math.min(N, (sx + 1) * step); kx++)
                for (let ky = sy * step; ky < Math.min(N, (sy + 1) * step); ky++) {
                  const i = ax.indexOf(0, kx as never);
                  const j = ax.indexOf(0, ky as never);
                  if (i < s || i >= t1 || j < u || j >= w) continue;
                  const x = value[(i - s) * (w - u) + (j - u)];
                  if (x !== null) [sum, c] = [sum + Math.fround(x), c + 1];
                }
              const want = c > 0 ? sum / c : Number.NaN;
              const got = dotPlotSample(ax, g, L, sx, sy);
              if (Number.isNaN(want) !== Number.isNaN(got) || (!Number.isNaN(want) && Math.abs(got - want) > 1e-3))
                bad++;
            }
        }
        expect(bad).toBe(0);
      });
});

describe('the hovered cell', () => {
  it("is outlined once, in its grid's ink, whether the grid's plot or another dot plot comes first on the track", async () => {
    const outline = async (wordsFirst: boolean) => {
      const doc = dotPlotDoc();
      if (wordsFirst) {
        const track = (doc.root.children ?? [])
          .flatMap(function walk(n): (typeof n)[] {
            return [n, ...(n.children ?? []).flatMap(walk)];
          })
          .find((n) => n.ref === 't-plot');
        track?.children?.reverse();
      }
      const t = await lit(doc, { width: 900, height: 600 });
      t.frame(recording(), Dirty.all);
      const l = t.ctx.layout.get(t.section)?.current;
      const v = t.ctx.layout.axis(t.section);
      const a = t.ctx.graph.assembled.value;
      if (!l || !v || !a) throw new Error('no layout');
      const track = a.view?.resolved.tracks.ref.findIndex((r) => a.set.strings.get(r) === 't-plot') ?? -1;
      const e = Array.from({ length: l.entryCount }, (_, k) => k).find((k) => l.entryTrack[k] === track) ?? -1;
      const top = l.scrolledTop + l.entryTop[e] - v.scrollY;
      const xOf = (k: number) => v.x(0, k as never) + v.cell / 2;
      t.hover(xOf(5), top + ((xOf(35) - v.left) / (v.right - v.left)) * l.entryHeight[e]);
      return passRects(t.overlay(), Pass.overlay).filter((r) => r[6] > 0);
    };
    const gridFirst = await outline(false);
    expect(gridFirst.length).toBeGreaterThan(0);
    expect(await outline(true)).toEqual(gridFirst);
  });
});
