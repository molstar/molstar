/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A drag-select's frame within its time, a drag kind each (the viewer's _spec, "SelectionModel and HighlightModel",
// Drags; `drag-kinds.ts`): a pointer's moves and the frame that applies them, the end moving every frame, at most 4 ms
// a frame on average (a quarter of a 60 Hz frame: the canvas's own drawing is not in it, headless). The wall-clock half
// of `drags.quick.ts`, whose allocation bounds stay in the perf tier. Heavy lane only: a view of 16 M values or 3,000
// blocks opens in seconds.

// @module-tag heavy

import { loadGenome } from '@molstar/seqx-core/resolve';
import { itWithin } from '@molstar/test-harness/timing';
import { afterEach, beforeAll, describe, expect, vi } from 'vitest';
import { DragKinds, openDragView } from './drag-kinds.js';

beforeAll(async () => {
  await loadGenome();
});

afterEach(() => {
  vi.useRealTimers();
});

/** A frame's mean bound, ms. */
const FrameMs = 4;
const Frames = 200;
const MovesAFrame = 3;

describe('a drag-select’s frame within its time', () => {
  for (const kind of DragKinds)
    itWithin(
      `${kind.name}: at most ${FrameMs} ms a frame`,
      async (within) => {
        const v = await openDragView(kind.view, true);
        const path = kind.path(v, 64);
        v.input.down(...path[0]);
        v.input.move(...path[1]);
        v.draw();
        // warm: one pass along the path
        for (let i = 0; i < 64; i++) {
          v.input.move(...path[i]);
          v.draw();
        }
        const ms = within(FrameMs * Frames, () => {
          const t0 = performance.now();
          for (let f = 0; f < Frames; f++) {
            for (let m = 0; m < MovesAFrame; m++) v.input.move(...path[(f * MovesAFrame + m) % 64]);
            v.draw();
          }
          return (performance.now() - t0) / Frames;
        });
        v.input.up();
        v.retire();
        expect(ms).toBeGreaterThanOrEqual(0);
      },
      120_000,
    );
});
