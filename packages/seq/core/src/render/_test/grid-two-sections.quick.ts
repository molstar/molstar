/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A grid over two sections of blocks (EC-209): a warm frame of the plot (a strip a block pair, its cells and its dots)
// and its crosshair hover allocate nothing. By hand, never in CI: `pnpm -C packages/seq/core test:perf`.

import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import { drawView, type FrameInput, type Hit, hitView } from '../index.js';
import { twoSectionsSpec } from './grid-two-sections-fixtures.js';
import { settle } from './helpers.js';
import { reusedHover } from './scenarios.js';

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
});
