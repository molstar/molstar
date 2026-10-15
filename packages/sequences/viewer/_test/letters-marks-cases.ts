/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A text-only letters row's marks cases, for their heap guard and its quick twin (letters-marks.heap.ts,
// letters-marks.quick.ts): the design row's and the plain row's letters layer drawn directly at a cell width and a
// device pixel ratio (dots, bars, a pixel a position), a fresh viewer for each (`vi.resetModules`).

import { vi } from 'vitest';

/** The rows: the design row's marks and the plain row's (`blockMarks`). */
export const LABELS = ['Design', 'Plain'] as const;

/** The marks: a cell width, a device pixel ratio, and what they draw. */
export const MARKS = [
  [5, 2, 'dots'],
  [5, 1, 'bars'],
  [1.2, 1, 'bars, a pixel a position'],
] as const;

/** The row's letters layer at `cell` px a position and `dpr`: its frame (a direct draw), and what the first drew. */
export async function marksFrame(label: string, cell: number, dpr: number) {
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
  return { frame: () => f.draw(f.direct), n, marked, dispose: () => t.ctx.dispose() };
}
