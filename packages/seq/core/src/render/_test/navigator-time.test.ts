/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The navigators' targets (performance.md, navigators): a window move (the picture's key held: the windows and marks
// drawn again over it) under 0.5 ms; the picture laid out again: a strip of p53's tracks under 4 ms, a chromosome's
// ideogram of about 60 bands under 1 ms, the 10,000-row overview from its summary grid under 8 ms. Each bound the
// heavy lane's (`itWithin`, the median of a run, best of three); the quick lane runs them unbounded.

import { overviewNavigator, stripNavigator } from '@molstar/seqx-core/lang';
import { itWithin } from '@molstar/test-harness/timing';
import { describe, expect } from 'vitest';
import type { NavigatorFrameInput, ViewLayout } from '../index.js';
import { p53TracksSpec } from './fixtures.js';
import { bandsSpec, settleNavigator, withNavigator } from './navigator-fixtures.js';
import { rowsSpec, Synthetic } from './rows-fixtures.js';

const Frames = 40;
/** The warm-up's frames: enough to compile the paths, a fifth of a timed run. */
const WarmFrames = 8;

/** The median of `frames` (`Frames`) navigator frames, each laid out again or not as `input(k)` says. */
function median(frame: (k: number) => void, frames = Frames): number {
  const t: number[] = [];
  for (let k = 0; k < frames; k++) {
    const t0 = performance.now();
    frame(k);
    t.push(performance.now() - t0);
  }
  t.sort((a, b) => a - b);
  return t[t.length >> 1];
}

describe('the navigators’ time', () => {
  itWithin('a window move: the picture kept, its windows and marks drawn again', async (within) => {
    const s = await settleNavigator(
      withNavigator(p53TracksSpec(), stripNavigator('/views/main/protein', { content: 'tracks' })),
    );
    const views: ViewLayout[] = [];
    for (let k = 0; k < 8; k++) {
      const r = await s.frames(
        { size: { width: 1200, height: 300 }, viewport: { start: 20 * k, end: 20 * k + 100 } },
        { size: { width: 1200, height: 48 } },
      );
      views.push(r.view.layout);
    }
    const input: NavigatorFrameInput = { size: { width: 1200, height: 48 }, moves: { '/views/main/protein': 1 } };
    const run = (frames = Frames) =>
      median((k) => {
        const f = s.nav.frame({ ...input, view: views[k % views.length] });
        if (f.rebuilt) throw new Error('a window move lays the picture out again');
      }, frames);
    run(WarmFrames);
    const m = within(0.5 * Frames, () => run());
    console.info(`navigator, a window move: median ${m.toFixed(3)} ms`);
  });

  itWithin('a strip of p53’s tracks, laid out again', async (within) => {
    const s = await settleNavigator(
      withNavigator(p53TracksSpec(), stripNavigator('/views/main/protein', { content: 'tracks' })),
    );
    const { view } = await s.frames({ size: { width: 1200, height: 300 } }, { size: { width: 1200, height: 48 } });
    const run = (frames = Frames) =>
      median((k) => s.nav.frame({ size: { width: 1200 - (k & 1), height: 48 }, view: view.layout }), frames);
    expect(s.nav.frame({ size: { width: 1199, height: 48 }, view: view.layout }).rebuilt).toBe(true);
    run(WarmFrames);
    const m = within(4 * Frames, () => run());
    console.info(`navigator, p53's tracks strip: median ${m.toFixed(3)} ms`);
  });

  itWithin('a chromosome’s ideogram of 60 bands, laid out again', async (within) => {
    const s = await settleNavigator(bandsSpec());
    const { view, nav } = await s.frames({ size: { width: 1200, height: 100 } }, { size: { width: 1200, height: 24 } });
    expect(nav.issues.filter((i) => i.severity !== 'info')).toEqual([]);
    expect(s.nav.frame({ size: { width: 1199, height: 24 }, view: view.layout }).rebuilt).toBe(true);
    const run = (frames = Frames) =>
      median((k) => s.nav.frame({ size: { width: 1200 - (k & 1), height: 24 }, view: view.layout }), frames);
    run(WarmFrames);
    const m = within(1 * Frames, () => run());
    console.info(`navigator, an ideogram of 60 bands: median ${m.toFixed(3)} ms`);
  });

  itWithin('the 10,000-row overview from its summary grid, laid out again', async (within) => {
    const s = await settleNavigator(withNavigator(rowsSpec(10_000, true), overviewNavigator('/views/v/s')), 0, [
      Synthetic,
    ]);
    const { view } = await s.frames({ size: { width: 1200, height: 600 } }, { size: { width: 96, height: 600 } });
    expect(s.nav.frame({ size: { width: 95, height: 600 }, view: view.layout }).rebuilt).toBe(true);
    expect(s.nav.pictureSection!.tracks.length).toBe(1);
    const run = (frames = Frames) =>
      median((k) => s.nav.frame({ size: { width: 96 - (k & 1), height: 600 }, view: view.layout }), frames);
    run(WarmFrames);
    const m = within(8 * Frames, () => run());
    console.info(`navigator, the 10,000-row overview: median ${m.toFixed(3)} ms`);
  });
});
