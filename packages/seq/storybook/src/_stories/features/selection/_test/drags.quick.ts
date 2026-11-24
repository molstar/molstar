/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A drag-select's frame, bounded, a drag kind each (the viewer's _spec, "SelectionModel and HighlightModel", Drags;
// `drag-kinds.ts`; the profile is `drag-profile.quick.ts`): a drag's warm frame (the pointer held, its moves and frames
// going on: the selection's end unchanged) allocates nothing, the move, its frame and the hit once a frame included; a
// frame that moves the end allocates the new selection value and what derives from it, bounded. The wall-clock bounds
// are `drag-time.test.ts`'s (heavy lane). By hand, never in CI: `pnpm -s -C packages/seq/storybook test:perf`.

import { loadGenome } from '@molstar/seqx-core/resolve';
import { quickCase } from '@molstar/test-harness/quick';
import { afterEach, beforeAll, describe, expect, inject, it, vi } from 'vitest';
import { DragKinds, openDragView } from './drag-kinds.js';

afterEach(() => {
  vi.useRealTimers();
});

beforeAll(async () => {
  await loadGenome();
});

/**
 * The bytes a frame moving the drag's end may allocate: the new value and its projections (measured 3.5 to 11.3 kB, the
 * sweeps' twelve parts the most).
 */
const MovingFrameBytes = 16 * 1024;

/**
 * What a frame of the genome's views allocates whatever runs (48 B: two sections' totals of 3·10⁹ units returned boxed
 * by the rails' idle check, and the window's start): the frame's, not the drag's (the frame alone is measured beside).
 */
const GenomeFrameBytes = 48;

describe('a drag-select’s frame, bounded', () => {
  for (const kind of DragKinds)
    it(kind.name, async ({ task }) => {
      const v = await openDragView(kind.view);
      const path = kind.path(v, 64);
      v.input.down(...path[0]);
      v.input.move(...path[32]);
      v.draw();
      // the end unchanged: the pointer held where it is, its moves and frames going on (at these scales a pixel is
      // more than a cell, 0.25 px a cell on the dot plot, 100 kb on the genome)
      const [x, y] = path[32];
      const e = {};
      const warm = quickCase(task, `drag frame warm: ${kind.name}`, () => {
        v.input.move(x, y, e);
        v.draw();
      }).r;
      // the frame alone, no move (the emphasis of the selection drawn again: what any frame over it costs)
      const still = quickCase(task, `drag frame still: ${kind.name}`, () => void v.draw()).r;
      // across the path: each frame a new end
      const moving = quickCase(
        task,
        `drag frame moving: ${kind.name}`,
        (i) => {
          const p = path[16 + (i % 32)];
          v.input.move(p[0], p[1]);
          v.draw();
        },
        { control: false },
      ).r;
      v.input.up();
      v.retire();
      console.info(
        `${kind.name}: a warm frame ${warm.bytesPerStep.toFixed(0)} B (${warm.note}); the frame alone ${still.bytesPerStep.toFixed(0)} B; a frame moving the end ${moving.bytesPerStep.toFixed(0)} B`,
      );
      const residue = kind.view.dataset === 'genome' ? GenomeFrameBytes : 0;
      expect(warm.bytesPerStep - residue, warm.note).toBeLessThanOrEqual(inject('perfThreshold'));
      expect(still.bytesPerStep - residue, still.note).toBeLessThanOrEqual(inject('perfThreshold'));
      expect(moving.bytesPerStep, moving.note).toBeLessThanOrEqual(MovingFrameBytes);
    });
});
