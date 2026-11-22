/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Region reads' allocation targets (_spec/README.md "Render", The grid layer: by region; The matrix layer; Rings
// layout: images as runs): a warm draw of a region grid's store tiles, of a matrix's items and of its bins, and of a
// ring of 40 matrix states as runs, allocates nothing; a hit on a region grid's cell or a matrix cell allocates no hit
// object. By hand, never in CI: `pnpm -C packages/seq/core test:perf`.

import { SeqRuntime } from '@molstar/seqx-core/data';
import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import { drawView, type FrameInput, GridType, MatrixType } from '../index.js';
import type { GridLay } from '../grid.js';
import type { MatrixLay } from '../matrix.js';
import { settle } from './helpers.js';
import { gridRegionSpec, Many, sourced } from './region-render-fixtures.js';

/** The point a layer's hit reads. */
const At = new Float64Array(2);

describe('region reads draw with nothing allocated once warm', () => {
  it('a region grid: a warm draw of its store tiles, a hit on a cell', async ({ task }) => {
    const s = await settle(gridRegionSpec(), [SeqRuntime]);
    const input: FrameInput = { size: { width: 1200, height: 2400 }, dpr: 2, viewport: { start: 512, end: 1536 } };
    const f = await s.frame(input);
    const lay = f.layout.sections[1].tracks[0].layers[0].lay as GridLay;
    expect(lay.images.length).toBeGreaterThanOrEqual(16);
    const d = quickCase(task, 'region grid draw', () => drawView(f.layout, f.lists), {});
    console.info(`region grid draw: ${d.r.bytesPerStep.toFixed(0)} bytes a frame (${d.r.note})`);
    expect(d.r.bytesPerStep, d.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
    const b = lay.box;
    const h = quickCase(
      task,
      'region grid hit',
      (i) => {
        At[0] = b.x + ((i * 37) % 100) * (b.width / 100);
        At[1] = b.y + ((i * 53) % 100) * (b.height / 100);
        GridType.hit(lay, At);
      },
      {},
    );
    console.info(`region grid hit: ${h.r.bytesPerStep.toFixed(0)} bytes a hit (${h.r.note})`);
    expect(h.r.bytesPerStep, h.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });

  for (const [name, viewport] of [
    ['items', { start: 0, end: 2000 }],
    ['bins', undefined],
  ] as const)
    it(`a matrix's ${name}: a warm draw, a hit on a cell`, async ({ task }) => {
      const s = await settle(sourced(), [SeqRuntime]);
      const input: FrameInput = { size: { width: 800, height: 200 }, ...(viewport ? { viewport } : {}) };
      const f = await s.frame(input);
      const lay = f.layout.sections[0].tracks[0].layers[0].lay as MatrixLay;
      expect(lay.n + lay.images.length).toBeGreaterThan(0);
      const d = quickCase(task, `matrix ${name} draw`, () => drawView(f.layout, f.lists), {});
      console.info(`matrix ${name} draw: ${d.r.bytesPerStep.toFixed(0)} bytes a frame (${d.r.note})`);
      expect(d.r.bytesPerStep, d.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
      const b = lay.box;
      const h = quickCase(
        task,
        `matrix ${name} hit`,
        (i) => {
          At[0] = b.x + ((i * 37) % 100) * (b.width / 100);
          At[1] = b.y + ((i * 53) % 100) * (b.height / 100);
          MatrixType.hit(lay, At);
        },
        {},
      );
      console.info(`matrix ${name} hit: ${h.r.bytesPerStep.toFixed(0)} bytes a hit (${h.r.note})`);
      expect(h.r.bytesPerStep, h.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
    });

  // open (RG's report): a warm draw of the ring's runs allocates about 1.5 kB a frame (40 rows, 41,000 sectors), the
  // cause not found yet; skipped until it is, its bound not raised
  it('a ring of 40 matrix states: a warm draw of its runs', async ({ task }) => {
    const s = await settle(sourced(Many, 'rings'), [SeqRuntime]);
    const f = await s.frame({ size: { width: 640, height: 900 } });
    drawView(f.layout, f.lists);
    const d = quickCase(task, 'matrix ring draw', () => drawView(f.layout, f.lists), {});
    console.info(`matrix ring draw: ${d.r.bytesPerStep.toFixed(0)} bytes a frame (${d.r.note})`);
    expect(d.r.bytesPerStep, d.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });
});
