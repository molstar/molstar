/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A drag-select's work a move and a frame, a drag kind each (the viewer's _spec, "SelectionModel and HighlightModel",
// Drags; `drag-kinds.ts`): a rect on the 4,000 × 4,000 dot plot, a range on its ruler, a box on the PAE grid, a sweep
// across 12 chromosomes on the 3,000-contig ring and on its row. Each reports the work a pointer move does (ms, bytes,
// selection writes), a move and its frame, the layouts, and the SVG twin's export (what a story beside the canvas draws
// again): the profile, no bound (`drags.quick.ts` holds them; a file of its own, so its runs leave the bounds' fork as
// it was). By hand, never in CI: `pnpm -s -C packages/seq/storybook test:perf`.

import { loadGenome } from '@molstar/seqx-core/resolve';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { measureBytesPerStep } from '@molstar/test-harness/quick';
import { afterEach, beforeAll, describe, it, vi } from 'vitest';
import { DragKinds, openDragView, type Point } from './drag-kinds.js';

afterEach(() => {
  vi.useRealTimers();
});

beforeAll(async () => {
  await loadGenome();
});

/** The mean ms of `n` calls of `f`. */
function msOf(f: (i: number) => void, n: number): number {
  const t0 = performance.now();
  for (let i = 0; i < n; i++) f(i);
  return (performance.now() - t0) / n;
}

/**
 * A drag along `points` (its first the press): what a move does alone (several a frame, as a pointer sends them), and
 * a move with its frame; frames slower than 20 ms are timed only (their bytes would take minutes to read).
 */
async function profile(name: string, v: SeqViewer, points: readonly Point[]): Promise<void> {
  let writes = 0;
  const off = v.selection.value.listen(() => writes++);
  const n = points.length;
  v.input.down(...points[0]);
  v.input.move(...points[1]);
  v.draw();
  const layouts = v.counts.layouts;
  const w0 = writes;
  const moveMs = msOf((i) => v.input.move(...points[i]), n);
  const perMove = (writes - w0) / n;
  const frame = (i: number) => {
    v.input.move(...points[i % n]);
    v.draw();
  };
  const frameMs = msOf(frame, n);
  const moveBytes = measureBytesPerStep((i) => v.input.move(...points[i % n])).bytesPerStep;
  const frameBytes = frameMs < 20 ? `${measureBytesPerStep(frame).bytesPerStep.toFixed(0)} B` : '(not read)';
  const laid = v.counts.layouts - layouts;
  v.input.up();
  off();
  const t = performance.now();
  const svg = v.renderSvg({ idPrefix: 'd-' });
  await vi.advanceTimersByTimeAsync(400);
  await svg;
  const svgMs = performance.now() - t;
  console.info(
    `${name}: a move ${moveMs.toFixed(3)} ms, ${moveBytes.toFixed(0)} B, ${perMove} selection writes; a move and its frame ${frameMs.toFixed(2)} ms, ${frameBytes}; ${laid} layouts; the SVG twin's export ${svgMs.toFixed(0)} ms`,
  );
}

describe('a drag-select’s work a move and a frame', () => {
  for (const kind of DragKinds)
    it(kind.name, async () => {
      const v = await openDragView(kind.view, true);
      await profile(kind.name, v, kind.path(v, 64));
      v.retire();
    });
});
