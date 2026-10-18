/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The case of the recording backend's guard (seqspec:viewer#the-render-backend; recording.heap.ts) and its quick twin
// (recording.quick.ts): a recorded frame of 4,000 rects and 1,000 text runs with every other batch kind, an image and a
// pan's copy, the backend reset each step. The heap guard and the quick twin share it.

import { RecordingBackend } from '../../src/testing/index.js';
import { drawFrame, Frame, PIXELS } from '../support/frame.js';

/** The quick case's counts: a frame takes about a millisecond, so windows of 250 steps (`DRAW`, wrapped-lanes-cases.ts). */
export const FRAMES = { warmSteps: 1_000, steps: 250, windows: 3 };

/** The step: one frame recorded; `rects` is the last frame's rect count (4,001). */
export function recorded(): { step: (i: number) => void; rects: () => number } {
  const r = new RecordingBackend();
  const f = new Frame(4000, 1000);
  const img = r.images.fromPixels(2, 2, PIXELS);
  return {
    step: (i) => {
      r.reset();
      drawFrame(r, f, i, img);
    },
    rects: () => r.rectItems.n,
  };
}
