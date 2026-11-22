/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A navigator over a genome's blocks (performance.md, Memory and Time; _spec/README.md "Render", Navigators, A navigator
// over blocks): a window move over 25 chromosomes, the window carried by identity in units onto the picture's row,
// redraws with nothing allocated once warm and in at most 0.5 ms a frame, on its zoom strip (always over blocks, its
// viewport the touched blocks' hull: a move across blocks that changes the hull lays the strip out again, not timed here).
// By hand, never in CI: `pnpm -C packages/seq/core test:perf`.

import { quickCase } from '@molstar/test-harness/quick';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { beforeAll, describe, expect, inject, it } from 'vitest';
import type { ViewLayout } from '../index.js';
import type { NavigatorFrameInput } from '../navigator.js';
import { genomeSpec, toyGenome } from './genome-fixtures.js';
import { settle } from './helpers.js';
import { settleNavigator, withNavigator } from './navigator-fixtures.js';

const Section = '/views/v/genome';
const Size = { width: 1200, height: 300 };
const NavSize = { width: 1200, height: 28 };

beforeAll(async () => {
  await loadGenome();
});

/** The genome with its navigator, and four view layouts of windows `w` units wide moving across chr3 into chr4. */
async function navigatorOver(w: number, turned = false) {
  // turned: chr3 reversed and chr4 at twice its length (the window moves across a reversed block into a scaled one)
  const options = turned ? { reversed: ['chr3'], scaled: { chr4: 2 } } : {};
  const spec = () =>
    withNavigator(genomeSpec(toyGenome(25, 0), options), { of: Section, axes: 'positions', content: 'axis' });
  const s = await settleNavigator(spec());
  const gap = s.view.sections[0].blocks![3].offset!;
  const views: ViewLayout[] = [];
  for (let k = 0; k < 4; k++) {
    const v = await settle(spec());
    const start = gap - w / 2 + (k - 2) * (w / 8);
    views.push((await v.frame({ size: Size, viewport: { start, end: start + w } })).layout);
  }
  await s.frames({ size: Size }, { size: NavSize });
  const inputs: NavigatorFrameInput[] = views.map((view) => ({
    size: NavSize,
    view,
    moves: { [Section]: 1 },
    lit: { section: Section, part: 'body' as const },
  }));
  for (const i of inputs) s.nav.frame(i);
  return { s, inputs };
}

/** The best of three runs' mean ms a frame over `n` frames. */
function msPerFrame(frame: (i: number) => void, n = 2000): number {
  let best = Number.POSITIVE_INFINITY;
  for (let t = 0; t < 3; t++) {
    const t0 = performance.now();
    for (let i = 0; i < n; i++) frame(i);
    best = Math.min(best, (performance.now() - t0) / n);
  }
  return best;
}

describe('a navigator over 25 chromosomes moves its window with nothing allocated, in at most 0.5 ms', () => {
  it('a window of 100 Mb across a block gap (the strip over both blocks, whose hull holds)', async ({ task }) => {
    const { s, inputs } = await navigatorOver(1e8);
    // over blocks the chromosome is always the control row: the same two rows at any width
    expect(s.nav.zoomHeight).toBeGreaterThan(0);
    const { r } = quickCase(task, 'navigator over blocks window move', (i) => s.nav.frame(inputs[i & 3]), {});
    const ms = msPerFrame((i) => s.nav.frame(inputs[i & 3]));
    console.info(`navigator over blocks: ${r.bytesPerStep.toFixed(0)} bytes, ${ms.toFixed(3)} ms a frame (${r.note})`);
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
    expect(ms).toBeLessThanOrEqual(0.5);
  });

  it('a window of 2 Mb: the zoom strip of the chromosomes it touches', async ({ task }) => {
    const { s, inputs } = await navigatorOver(2e6);
    expect(s.nav.zoomHeight).toBeGreaterThan(0);
    const { r } = quickCase(
      task,
      'navigator over blocks zoom strip window move',
      (i) => s.nav.frame(inputs[i & 3]),
      {},
    );
    const ms = msPerFrame((i) => s.nav.frame(inputs[i & 3]));
    console.info(
      `navigator over blocks, zoom strip: ${r.bytesPerStep.toFixed(0)} bytes, ${ms.toFixed(3)} ms (${r.note})`,
    );
    // 16 B a frame were the mark's `wraps` read through a getter returning the genome's total (boxed): inline now
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
    expect(ms).toBeLessThanOrEqual(0.5);
  });

  it('over a reversed and a scaled block (chr3 reversed, chr4 ×2), with and without the zoom strip', async ({
    task,
  }) => {
    for (const w of [1e8, 2e6]) {
      const { s, inputs } = await navigatorOver(w, true);
      const { r } = quickCase(task, `navigator over turned blocks ${w}`, (i) => s.nav.frame(inputs[i & 3]), {});
      const ms = msPerFrame((i) => s.nav.frame(inputs[i & 3]));
      console.info(
        `navigator over turned blocks (${w}): ${r.bytesPerStep.toFixed(0)} bytes, ${ms.toFixed(3)} ms (${r.note})`,
      );
      expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
      expect(ms).toBeLessThanOrEqual(0.5);
    }
  });
});

describe('a drag across many small blocks: the strip’s viewport follows the hull of the blocks the window touches', () => {
  it('a window of 1 Mb over 3,000 contigs, a different hull most frames: the strip laid out again, bounded', async ({
    task,
  }) => {
    const spec = () =>
      withNavigator(genomeSpec(toyGenome(25, 3000, 7)), { of: Section, axes: 'positions', content: 'axis' });
    const s = await settleNavigator(spec());
    const blocks = s.view.sections[0].blocks!;
    const first = blocks.findIndex((b) => b.axis.system.includes('chrUn'));
    const views: ViewLayout[] = [];
    // 32 windows, each a few contigs on from the last (a drag along the run): every frame's hull is another
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
    const { r } = quickCase(
      task,
      'navigator over blocks drag across small blocks',
      (i) => s.nav.frame(inputs[i & 31]),
      {},
    );
    const ms = msPerFrame((i) => s.nav.frame(inputs[i & 31]), 500);
    console.info(
      `navigator drag across small blocks: ${r.bytesPerStep.toFixed(0)} bytes, ${ms.toFixed(3)} ms (${r.note})`,
    );
    // a new hull is a new layout of the strip's section (~110 KB for twenty contigs: the layout's own, measured): a
    // bound, so it never grows with the genome's 3,000 blocks (the context is not laid out again)
    // its 2 ms bound is a wall-clock one: in the heavy lane, navigator-drag-time.test.ts
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(256 * 1024);
  });
});
