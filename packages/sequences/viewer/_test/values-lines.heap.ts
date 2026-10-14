/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Values colored by their values through the conformance kit's heap guard (seqspec:format#layer-types;
// seqspec:viewer#the-render-backend): a frame of `linesDoc`'s line in four bands (a discrete palette, its color changing
// halfway between two positions), its area under a continuous palette and its plain line allocates nothing after
// warm-up, at each of the kit's zooms (one below a device px per position: columns' means and bands). The positive
// control goes through the kit's `heap` hook: at each of those frames the same warm-up and measurement of the frame
// allocating an object each (`leaking`) must catch it.
//
// A file of its own: in draws.heap.ts, after "values on 40,000 positions" the line's frame at 0.62 px (tiny) measured
// about 52 bytes a step, and HEAD's uncolored line the same (68), a draw V8 left in its middle tier under the vitest
// module runner (draws.heap.ts's header); alone it is flat.

import { expect, it, vi } from 'vitest';
import type { LayerRenderer } from '../src/index.js';
import {
  controlCaught,
  describeScaling,
  FRAME_N,
  FRAME_WARM,
  FRAME_WARM_MS,
  heapScaling,
  isFlat,
  leaking,
} from './heap.js';

/** The guard's measurement: warmed by count and by time, then `FRAME_N` frames and 4× measured. */
const measure = (frame: (i: number) => void) => {
  const t0 = performance.now();
  for (let i = 0; i < FRAME_WARM || performance.now() - t0 < FRAME_WARM_MS; i++) frame(i);
  return heapScaling(FRAME_N, frame);
};
const guard = (frame: () => void) => {
  const r = measure(frame);
  if (!isFlat(r)) throw new Error(describeScaling(r, FRAME_N));
};
/** The positive control as the kit's `heap` hook: a frame's check fails unless its leaking frame is caught. */
const control = (frame: () => void) => {
  if (!controlCaught(() => measure(leaking(frame)))) throw new Error('the leaking frame was not caught');
};

async function run(heap = guard) {
  vi.resetModules();
  const { KernelRenderers } = await import('../src/index.js');
  const { rendererConformance } = await import('../src/testing/index.js');
  const { opened } = await import('./layout-fixtures.js');
  const { linesDoc } = await import('./prepare-fixtures.js');
  const all = KernelRenderers as LayerRenderer[];
  const entry = all.find((r) => r.name === 'values') as LayerRenderer;
  const checks = await rendererConformance(entry, {
    open: async ({ entry: has }) =>
      (await opened(linesDoc(20_000), { renderers: has ? all : all.filter((r) => r !== entry) })).ctx,
    heap,
    items: false,
  });
  return checks.filter((c) => c.check === 'heap' && c.geometry === 'linear');
}

it('a frame of values colored by their values allocates nothing after warm-up, at every zoom', async () => {
  const heap = await run();
  expect(heap.length).toBeGreaterThan(0);
  expect(heap.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
}, 240_000);

it('a values frame allocating an object each is caught, at every zoom (positive control)', async () => {
  const heap = await run(control);
  expect(heap.length).toBeGreaterThan(0);
  expect(heap.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
}, 240_000);
