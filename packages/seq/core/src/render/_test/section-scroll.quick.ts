/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A tracks section past its cap (_spec/README.md "Render", Section scroll): a warm frame of the scrolled section (no
// layout: its draw and its emphasis redraw) allocates nothing. By hand, never in CI: `pnpm -C packages/seq/core
// test:perf`.

import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import { drawView, hitView } from '../index.js';
import { settle } from './helpers.js';
import { geneTracks } from './section-scroll-fixtures.js';

describe('a scrolled section allocates nothing warm', () => {
  it('30 gene tracks scrolled under the ruler: a warm frame and a hover redraw', async ({ task }) => {
    const s = await settle(geneTracks());
    const size = { width: 800, height: 2000 };
    const f0 = await s.frame({ size, settled: true });
    const a = f0.layout.sections[0].axis!;
    const emphasis: { hover?: ReturnType<typeof hitView> } = {};
    const input = { size, settled: true, viewport: { start: a.u0, end: a.u1, top: 200 }, emphasis };
    const f = await s.frame(input);
    expect(f.layout.sections[0].placement.top).toBe(200);
    const t = f.layout.sections[0].tracks[8].box;
    const hits = [hitView(f.layout, t.x + 40, t.y + 8), undefined];
    const { r } = quickCase(
      task,
      'section scrolled warm frame',
      (i) => {
        emphasis.hover = hits[i & 1];
        s.renderer.frame(input);
      },
      {},
    );
    console.info(`section scrolled warm frame: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
    const { layout, lists } = s.renderer.frame(input);
    const d = quickCase(task, 'section scrolled draw', () => drawView(layout, lists), {});
    console.info(`section scrolled draw: ${d.r.bytesPerStep.toFixed(0)} bytes a frame (${d.r.note})`);
    expect(d.r.bytesPerStep, d.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });
});
