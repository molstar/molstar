/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Region reads' allocation targets (_spec/README.md "Render", The grid layer: by region; The matrix layer; Rings
// layout: images as runs): a warm draw of a region grid's store tiles, of a matrix's items and of its bins, and of a
// ring of matrix states as runs, allocates nothing; a hit on a region grid's cell or a matrix cell allocates no hit
// object. By hand, never in CI: `pnpm -C packages/seq/core test:perf`.

import { SeqRuntime } from '@molstar/seqx-core/data';
import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import { drawView, type FrameInput, GridType, MatrixType, type PairsLay } from '../index.js';
import type { GridLay } from '../grid.js';
import type { MatrixLay } from '../matrix.js';
import { settle } from './helpers.js';
import { gridRegionSpec, linksSpec, Many, sourced } from './region-render-fixtures.js';

/**
 * Bytes a region pairs pan frame may allocate, at either size: the arcs' layout over about 2,400 bin pairs in view
 * (239 KB at 6,000 links, 275 KB at 60,000, 2026-11-22), never a key and a set entry a row read (1.26 and 1.53 MB
 * before the anchor rule).
 */
const PairsPanBound = 400_000;

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

  // 8 states: what a warm draw allocates is per image, not per row (40 rows take minutes of steps for the same answer)
  it('a ring of 8 matrix states: a warm draw of its runs', async ({ task }) => {
    const s = await settle(sourced(Many.slice(0, 8), 'rings'), [SeqRuntime]);
    const f = await s.frame({ size: { width: 640, height: 900 } });
    drawView(f.layout, f.lists);
    const d = quickCase(task, 'matrix ring draw', () => drawView(f.layout, f.lists), {});
    console.info(`matrix ring draw: ${d.r.bytesPerStep.toFixed(0)} bytes a frame (${d.r.note})`);
    expect(d.r.bytesPerStep, d.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });

  // a pan's layout of region pairs (bins in view counted, each once by the anchor rule, nothing built again)
  const PanOptions = {
    warmSteps: 300,
    steps: 100,
    windows: 3,
    maxWarmSteps: 2000,
    control: false,
    why: 'a layout',
  } as const;
  for (const n of [6000, 60_000])
    it(`region pairs, ${n.toLocaleString('en-US')} links: a warm pan's layout allocates no key a row`, async ({
      task,
    }) => {
      const s = await settle(linksSpec(n), [SeqRuntime]);
      const at = (k: number): FrameInput => ({
        size: { width: 720, height: 300 },
        settled: true,
        viewport: { start: 40_000 + (k % 16) * 1000, end: 840_000 + (k % 16) * 1000 },
      });
      for (let k = 0; k < 16; k++) await s.frame(at(k));
      const lay = (await s.frame(at(0))).layout.sections[0].sectionLays[0].lay as PairsLay;
      expect(lay.prep?.weight).toBeDefined();
      const r = quickCase(task, `region pairs pan ${n}`, (i) => void s.renderer.frame(at(i + 1)), PanOptions).r;
      console.info(`region pairs pan (${n}): ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
      expect(r.bytesPerStep, r.note).toBeLessThan(PairsPanBound);
    });
});
