/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { loadGenome } from '@molstar/seqx-core/resolve';
import { itWithin } from '@molstar/test-harness/timing';
import { beforeAll, describe, expect } from 'vitest';
import type { ViewLayout } from '../index.js';
import type { NavigatorFrameInput } from '../navigator.js';
import { genomeSpec, toyGenome } from './genome-fixtures.js';
import { settle } from './helpers.js';
import { settleNavigator, withNavigator } from './navigator-fixtures.js';

// A drag along 3,000 contigs: the strip over the touched blocks laid out again for most frames, each within 2 ms (the
// wall-clock half of navigator-blocks.quick.ts's drag case, whose allocation bound stays in the perf tier). Quick lane:
// the frames run unbounded; heavy lane: the bound.

const Section = '/views/v/genome';
const Size = { width: 1200, height: 300 };
const NavSize = { width: 1200, height: 28 };

beforeAll(async () => {
  await loadGenome();
});

describe('a drag across many small blocks', () => {
  itWithin('3,000 contigs, a window of 1 Mb a few contigs on each frame: at most 2 ms a frame', async (within) => {
    const spec = () =>
      withNavigator(genomeSpec(toyGenome(25, 3000, 7)), { of: Section, axes: 'positions', content: 'axis' });
    const s = await settleNavigator(spec());
    const blocks = s.view.sections[0].blocks!;
    const first = blocks.findIndex((b) => b.axis.system.includes('chrUn'));
    const views: ViewLayout[] = [];
    for (let k = 0; k < 32; k++) {
      const v = await settle(spec());
      const start = blocks[first + k * 20].offset!;
      views.push((await v.frame({ size: Size, viewport: { start, end: start + 1e6 } })).layout);
    }
    await s.frames({ size: Size }, { size: NavSize });
    const inputs: NavigatorFrameInput[] = views.map((view) => ({
      size: NavSize,
      view,
      moves: { [Section]: 1 },
      lit: { section: Section, part: 'body' as const },
    }));
    for (const i of inputs) s.nav.frame(i);
    expect(s.nav.zoomHeight).toBeGreaterThan(0);
    const n = 500;
    const ms = within(2 * n, () => {
      const t0 = performance.now();
      for (let i = 0; i < n; i++) s.nav.frame(inputs[i & 31]);
      return (performance.now() - t0) / n;
    });
    expect(ms).toBeGreaterThanOrEqual(0);
  }, 60_000);
});
