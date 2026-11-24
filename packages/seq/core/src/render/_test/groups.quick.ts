/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A section with track groups (_spec/README.md "Render", Track groups): a warm frame and its hover redraw allocate
// nothing (a fold changes the layout once, never a warm frame's work); a pan frame (a layout and a draw) allocates at
// most a bound. By hand, never in CI: `pnpm -C packages/seq/core test:perf`.

import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import { groupsSpec } from '../../view/_test/groups-spec.js';
import { drawView, hitView } from '../index.js';
import { settle } from './helpers.js';

const Section = '/views/v/p53';
const Size = { width: 800, height: 600 };
const LayoutOptions = { warmSteps: 500, steps: 200, windows: 3, maxWarmSteps: 4000 } as const;
/**
 * A pan frame's bytes: 62.7 KB measured (2026-11-24; 67.8 KB before the folding and the headers were kept across
 * layouts), as the seven tracks' layout without groups ('p53 small layers': 65 KB), with half again as headroom.
 */
const PanBound = 94_000;

describe('a section with track groups', () => {
  it('Features open, Variants and Signals folded: a warm frame and a hover redraw allocate nothing; a pan is bounded', async ({
    task,
  }) => {
    const s = await settle(groupsSpec());
    const groups = ['variants', 'signals'];
    const emphasis: { hover?: ReturnType<typeof hitView> } = {};
    const input = { size: Size, settled: true, viewports: { [Section]: { start: 0, end: 40, groups } }, emphasis };
    const f = await s.frame(input);
    const sec = f.layout.sections[0];
    expect(sec.groups!.map((g) => g.folded)).toEqual([false, true, true]);
    // a hover on a domain, then on the folded Variants row's summary, then none
    const t = sec.tracks[1].box;
    const v = sec.tracks[3].box;
    const hits = [hitView(f.layout, t.x + 40, t.y + t.height / 2), hitView(f.layout, v.x + 40, v.y + v.height / 2)];
    const warm = quickCase(
      task,
      'groups warm frame',
      (i) => {
        emphasis.hover = i % 3 === 2 ? undefined : hits[i % 3];
        s.renderer.frame(input);
      },
      {},
    );
    console.info(`groups warm frame: ${warm.r.bytesPerStep.toFixed(0)} bytes a frame (${warm.r.note})`);
    expect(warm.r.bytesPerStep, warm.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
    const { layout, lists } = s.renderer.frame(input);
    const d = quickCase(task, 'groups draw', () => drawView(layout, lists), {});
    console.info(`groups draw: ${d.r.bytesPerStep.toFixed(0)} bytes a frame (${d.r.note})`);
    expect(d.r.bytesPerStep, d.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
    // a pan: two windows in turn, each a layout and a draw
    const pans = [0, 6].map((start) => ({
      size: Size,
      settled: true,
      viewports: { [Section]: { start, end: start + 40, groups } },
      emphasis: {},
    }));
    const pan = quickCase(task, 'groups pan frame', (i) => void s.renderer.frame(pans[i & 1]), LayoutOptions);
    console.info(`groups pan frame: ${pan.r.bytesPerStep.toFixed(0)} bytes a frame (${pan.r.note})`);
    expect(pan.r.bytesPerStep, pan.r.note).toBeLessThanOrEqual(PanBound);
  });
});
