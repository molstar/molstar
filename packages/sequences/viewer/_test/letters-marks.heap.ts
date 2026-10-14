/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A frame of a text-only letters row's marks (dots, then bars), another polymer's (`of`) and the plain row's under an
// empty `rules` theme, allocates nothing after warm-up (seqspec:viewer#the-render-backend); each case loads the viewer
// afresh and warms its frame until it stops allocating (see draws.heap.ts, `warmFrame`). The positive control: the
// design row's dots allocating an object a frame (`leaking`), through the same warm-up and measurement, are caught.

import { describe, expect, it, vi } from 'vitest';
import { controlCaught, describeScaling, FRAME_N, heapScaling, isFlat, leaking, warmFrame } from './heap.js';

/** The guard's measurement: warmed until the frame stops allocating (`warmFrame`), then `FRAME_N` frames and 4×. */
const measure = (step: (i: number) => void) => {
  let i = 0;
  const frame = () => step(i++);
  warmFrame(frame);
  return heapScaling(FRAME_N, frame);
};
const guard = (frame: () => void) => {
  const r = measure(frame);
  if (!isFlat(r)) throw new Error(describeScaling(r, FRAME_N));
};

async function frames(label: string, cell: number, dpr: number, heap = guard) {
  vi.resetModules();
  const { KernelRenderers } = await import('../src/index.js');
  const { LayerFrame, RecordingBackend } = await import('../src/testing/index.js');
  const { drawn, layerOf } = await import('./draw-fixtures.js');
  const { marksDoc } = await import('./letters-marks-fixtures.js');
  const t = await drawn(marksDoc(false), { width: 600, height: 200, dpr });
  const { k, section } = layerOf(t.ctx, label);
  t.ctx.viewports.get(section)?.setView(cell, 0, 'api');
  await t.ctx.graph.update();
  const entry = KernelRenderers.find((r) => r.name === 'letters');
  if (!entry) throw new Error('no letters renderer');
  const f = new LayerFrame(t.ctx, entry as never, section, k);
  const direct = new RecordingBackend();
  const n = f.draw(direct);
  const marked = direct.glyphItems.n + direct.rectItems.n;
  heap(() => f.draw(f.direct));
  t.ctx.dispose();
  return { n, marked };
}

describe("a text-only letters row's marks allocate nothing", () => {
  // The plain row (`blockMarks`): its 48 bytes a frame under a count warm-up were the draw's climb to TurboFan.
  for (const label of ['Design', 'Plain'])
    for (const [cell, dpr, what] of [
      [5, 2, 'dots'],
      [5, 1, 'bars'],
      [1.2, 1, 'bars, a pixel a position'],
    ] as const)
      it(`${label}: ${what}`, async () => {
        const r = await frames(label, cell, dpr);
        // Not vacuous: the frame draws the marks.
        expect(r.marked).toBeGreaterThan(0);
      }, 240_000);

  it('Design: dots, an object allocated a frame (the positive control)', async () => {
    let caught = false;
    await frames('Design', 5, 2, (frame) => {
      caught = controlCaught(() => measure(leaking(frame)));
    });
    expect(caught).toBe(true);
  }, 240_000);
});
