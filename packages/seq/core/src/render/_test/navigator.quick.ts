/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The navigators' allocation targets (performance.md, Memory): a window move redraws the navigator over its kept
// picture with nothing allocated once warm (the windows carried from the view's layouts, the marks); a rail drawn with
// its window lit; a window move on a zoom strip. By hand, never in CI: `pnpm -C packages/seq/core test:perf`.

import { stripNavigator } from '@molstar/seqx-core/lang';
import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import { type FrameInput, hitView, type ViewLayout } from '../index.js';
import { type NavigatorFrameInput } from '../navigator.js';
import { p53TracksSpec } from './fixtures.js';
import { settle } from './helpers.js';
import { bandsSpec, settleNavigator, withNavigator } from './navigator-fixtures.js';

describe('a navigator redraws with nothing allocated once warm', () => {
  it('a window move over the strip of p53’s tracks', async ({ task }) => {
    const s = await settleNavigator(
      withNavigator(p53TracksSpec(), stripNavigator('/views/main/protein', { content: 'tracks' })),
    );
    // each view layout of its own renderer: a window move is the view's frame changing under the navigator
    const views: ViewLayout[] = [];
    for (let k = 0; k < 4; k++) {
      const v = await settle(p53TracksSpec());
      views.push(
        (await v.frame({ size: { width: 1200, height: 300 }, viewport: { start: 30 * k, end: 30 * k + 90 } })).layout,
      );
    }
    await s.frames({ size: { width: 1200, height: 300 } }, { size: { width: 1200, height: 48 } });
    const inputs: NavigatorFrameInput[] = views.map((view) => ({
      size: { width: 1200, height: 48 },
      view,
      moves: { '/views/main/protein': 1 },
      lit: { section: '/views/main/protein', part: 'body' as const },
    }));
    const { r } = quickCase(task, 'navigator window move', (i) => s.nav.frame(inputs[i & 3]), {});
    console.info(`navigator window move: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });

  it('a window move on a zoom strip (a region a sliver of its chromosome)', async ({ task }) => {
    const s = await settleNavigator(bandsSpec());
    const views: ViewLayout[] = [];
    for (let k = 0; k < 4; k++) {
      const v = await settle(bandsSpec());
      views.push(
        (await v.frame({ size: { width: 1200, height: 300 }, viewport: { start: 100 * k, end: 100 * k + 300 } }))
          .layout,
      );
    }
    await s.frames({ size: { width: 1200, height: 300 } }, { size: { width: 1200, height: 24 } });
    expect(s.nav.zoomHeight).toBeGreaterThan(0);
    const inputs: NavigatorFrameInput[] = views.map((view) => ({
      size: { width: 1200, height: 24 },
      view,
      moves: { '/views/v/region': 1 },
      lit: { section: '/views/v/region', part: 'body' as const },
      pointerX: 300,
      pointerY: 5,
    }));
    const { r } = quickCase(task, 'zoom strip window move', (i) => s.nav.frame(inputs[i & 3]), {});
    console.info(`zoom strip window move: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });

  it('a rail hovered: its window lit and grown', async ({ task }) => {
    const s = await settle(p53TracksSpec());
    const rails = { '/views/main/protein': true };
    const input: FrameInput = {
      size: { width: 1200, height: 300 },
      viewport: { start: 50, end: 150 },
      rails,
      railsActive: rails,
    };
    const f = await s.frame(input);
    const rail = f.layout.sections[0].rail!;
    const hits = [0.1, 0.12, 0.3, 0.7].map((t) => ({ ...hitView(f.layout, rail.x + t * rail.width, rail.top + 1)! }));
    const emphasis: { hover?: (typeof hits)[number] } = {};
    const at = { ...input, emphasis };
    const { r } = quickCase(
      task,
      'rail hover',
      (i) => {
        emphasis.hover = hits[i & 3];
        s.renderer.frame(at);
      },
      {},
    );
    console.info(`rail hover: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });
});
