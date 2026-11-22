/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Region reads' frame times (_spec/README.md "Render", The grid layer: by region; Rings layout: images as runs): a
// region grid's pan frame over 4 × 4 store tiles, and a ring of 40 matrix states drawn as runs, each under 8 ms. Each
// bound is the heavy lane's (`itWithin`, best of three); the quick lane runs them unbounded.

import { SeqRuntime } from '@molstar/seqx-core/data';
import { itWithin } from '@molstar/test-harness/timing';
import { describe, expect } from 'vitest';
import type { GridLay } from '../grid.js';
import type { FrameInput } from '../index.js';
import { settle } from './helpers.js';
import { gridRegionSpec, Many, sourced } from './region-render-fixtures.js';

const Frames = 20;

/** The median ms of `frame(k)` over `Frames` frames. */
function median(frame: (k: number) => void): number {
  const t: number[] = [];
  for (let k = 0; k < Frames; k++) {
    const t0 = performance.now();
    frame(k);
    t.push(performance.now() - t0);
  }
  t.sort((a, b) => a - b);
  return t[t.length >> 1];
}

describe('region reads in a frame', () => {
  itWithin('a region grid pans over 4 × 4 store tiles (EC-223)', async (within) => {
    const s = await settle(gridRegionSpec(), [SeqRuntime]);
    const size = { width: 1200, height: 2400 };
    const at = (k: number): FrameInput => ({ size, dpr: 2, viewport: { start: 512 + (k % 8), end: 1536 + (k % 8) } });
    // every tile of the pan read and held
    for (let k = 0; k < 8; k++) await s.frame(at(k));
    const lay = (await s.frame(at(0))).layout.sections[1].tracks[0].layers[0].lay as GridLay;
    expect(lay.binLevel).toBe(1);
    expect(lay.images.length).toBeGreaterThanOrEqual(16);
    const pan = () => median((k) => s.renderer.frame(at(k + 1)));
    pan();
    const ms = within(8 * Frames, pan);
    console.info(`region grid pan, ${lay.images.length} tiles: median ${ms.toFixed(3)} ms`);
  });

  itWithin('a ring of 40 matrix states drawn as runs (EC-108)', async (within) => {
    const s = await settle(sourced(Many, 'rings'), [SeqRuntime]);
    // a pixel wider each frame: each one laid out and drawn again
    const at = (k: number): FrameInput => ({ size: { width: 640 + (k & 1), height: 900 } });
    await s.frame(at(0));
    await s.frame(at(1));
    const ring = () =>
      median((k) => {
        if (!s.renderer.frame(at(k)).rebuilt) throw new Error('a frame at another width lays the ring out again');
      });
    ring();
    const ms = within(8 * Frames, ring);
    console.info(`ring of 40 matrix states: median ${ms.toFixed(3)} ms`);
  });
});
