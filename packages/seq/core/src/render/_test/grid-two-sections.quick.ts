/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A grid over two sections of blocks (EC-209): a warm frame of the plot (a strip a block pair, its cells and its dots)
// and its crosshair hover allocate nothing; a categorical grid's (EC-35) draw with its labels and its hover neither. By hand, never in CI: `pnpm -C packages/seq/core test:perf`.

import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import { drawView, type FrameInput, GridPlacement, type Hit, hitView } from '../index.js';
import { kinship } from './grid-categorical-fixtures.js';
import { twoSectionsSpec } from './grid-two-sections-fixtures.js';
import { settle } from './helpers.js';
import { ScaledControlWhy, scaledControl } from './scaled-control.js';
import { reusedHover } from './scenarios.js';

/** A categorical pan frame's bound, bytes: its layout's (a positions grid of as many measured about 50 kB). */
const CategoricalPanBound = 40_000;

const Input: FrameInput = {
  size: { width: 1440, height: 1600 },
  dpr: 2,
  settled: true,
  viewports: { '/views/v/a': { start: 3, end: 300 }, '/views/v/b': { start: 10, end: 260 } },
};

describe('a grid over two sections of blocks draws and hovers with nothing allocated once warm', () => {
  it('a warm frame, a warm draw and the crosshair hover', async ({ task }) => {
    const s = await settle(
      twoSectionsSpec({
        xs: [120, 90, 100],
        ys: [80, 110, 70],
        reversed: [1],
        pairs: {
          x: 1,
          y: 0,
          rows: [
            [2, 40, 4, 42, 1],
            [50, 80, 10, 40, -1],
            [85, 86, 60, 61, 1],
          ],
        },
      }),
    );
    const f = await s.frame(Input);
    const frame = quickCase(task, 'grid two sections frame', () => void s.renderer.frame(Input), {});
    console.info(`grid two sections frame: ${frame.r.bytesPerStep.toFixed(0)} bytes a frame (${frame.r.note})`);
    expect(frame.r.bytesPerStep, frame.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
    const d = quickCase(task, 'grid two sections draw', () => drawView(f.layout, f.lists), {});
    console.info(`grid two sections draw: ${d.r.bytesPerStep.toFixed(0)} bytes a frame (${d.r.note})`);
    expect(d.r.bytesPerStep, d.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
    const plot = f.layout.sections[2].placement;
    const hits: Hit[] = [];
    for (let k = 0; k < 40; k++) {
      const x = plot.box.x + 200 + ((k * 37) % 900);
      const y = plot.box.y + 10 + ((k * 53) % Math.max(1, plot.box.height - 20));
      const h = hitView(f.layout, x, y);
      if (h) hits.push({ ...h });
    }
    expect(hits.length).toBeGreaterThan(10);
    const hover = reusedHover(hits);
    const at = { ...Input, emphasis: hover.emphasis };
    const h = quickCase(
      task,
      'grid two sections hover',
      (i) => {
        hover.set(i);
        s.renderer.frame(at);
      },
      {},
    );
    console.info(`grid two sections hover: ${h.r.bytesPerStep.toFixed(0)} bytes a frame (${h.r.note})`);
    expect(h.r.bytesPerStep, h.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });

  it('a categorical grid (EC-35): a warm draw with its members’ labels, and its crosshair hover', async ({ task }) => {
    const s = await settle(kinship(400, 3, { field: 'score' }, { field: 'score', direction: 'descending' }).spec);
    const input: FrameInput = { size: { width: 1200, height: 1400 }, dpr: 2, settled: true };
    const f = await s.frame(input);
    const d = quickCase(task, 'categorical grid draw', () => drawView(f.layout, f.lists), {});
    console.info(`categorical grid draw: ${d.r.bytesPerStep.toFixed(0)} bytes a frame (${d.r.note})`);
    expect(d.r.bytesPerStep, d.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
    const plot = f.layout.sections[0].placement as GridPlacement;
    const hits: Hit[] = [];
    for (let k = 0; k < 40; k++) {
      const h = hitView(f.layout, plot.strips.axes[0].x0 + 5 + ((k * 37) % 800), plot.plotTop + 5 + ((k * 53) % 800));
      if (h) hits.push({ ...h });
    }
    expect(hits.length).toBeGreaterThan(10);
    const hover = reusedHover(hits);
    const at = { ...input, emphasis: hover.emphasis };
    const h = quickCase(
      task,
      'categorical grid hover',
      (i) => {
        hover.set(i);
        s.renderer.frame(at);
      },
      {},
    );
    console.info(`categorical grid hover: ${h.r.bytesPerStep.toFixed(0)} bytes a frame (${h.r.note})`);
    expect(h.r.bytesPerStep, h.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });

  it('a categorical grid of 2,000 members (EC-35): a pan frame a bounded size, below a positions grid’s', async ({
    task,
  }) => {
    // a pan relays the section (a layout a bounded size a frame, as any grid's: a 2,000-position grid's pan measured
    // about 50 kB); the categorical axis adds nothing: its members and labels are kept, measured 35 kB
    const cat = await settle(kinship(2000, 5, { field: 'score' }, { field: 'score', direction: 'descending' }).spec);
    const at = (k: number): FrameInput => ({
      size: { width: 1200, height: 1400 },
      dpr: 2,
      settled: true,
      viewports: { '/views/v/plot': { start: 100 + (k % 2), end: 1100 + (k % 2), y: { start: 0, end: 2000 } } },
    });
    for (const k of [0, 1]) await cat.frame(at(k));
    const pan = (i: number) => void cat.renderer.frame(at(i));
    // the 88-byte control read 56 over the case in one run (a window of 800 bytes a step more): scaledControl's 4 KB
    const c = quickCase(task, 'categorical grid pan', pan, { control: false, why: ScaledControlWhy });
    console.info(`categorical grid pan: ${c.r.bytesPerStep.toFixed(0)} bytes a frame (${c.r.note})`);
    expect(c.r.bytesPerStep, c.r.note).toBeLessThan(CategoricalPanBound);
    scaledControl('categorical grid pan', pan, c.r);
  });
});
