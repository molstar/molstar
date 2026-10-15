/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The overlay's progressive-loading cases, for the heap guards and quick twins of streams-progress.heap.ts and
// streams-progress-landed.heap.ts (streams-progress.quick.ts, streams-progress-landed.quick.ts): a live composer over
// the streamed genome (genes, variants, letters, a signal) at 100 kb, DPR 2, with every region pending (skeletons
// shimmering, the clock advancing a frame each step, the surface's next frame decided from `animating`), every region
// failed (error bands, hatches and hints) or every region landed after a second of pending frames, the clock held so
// the fade never ends (`overlayFrame`, stream-progress-fixtures.ts). A step is an overlay frame. Each call loads the
// viewer afresh (`vi.resetModules`).

import { vi } from 'vitest';
import type { ProgressCase } from './stream-progress-fixtures.js';

/** The overlay frame over the genome in state `c`, and the fixture's `rectsOf` and `SkinToken`. */
export async function progressCase(c: ProgressCase) {
  vi.resetModules();
  const { overlayFrame, rectsOf, SkinToken } = await import('./stream-progress-fixtures.js');
  return { ...(await overlayFrame(c)), rectsOf, SkinToken };
}
