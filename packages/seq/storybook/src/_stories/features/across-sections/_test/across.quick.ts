/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A bridge on the pointer's path, over a region source (the stories' peach and grape, read by region): a pan carried
// through a lock allocates what both sides panned by hand do (the carrying kept a resolved view, never found again a move), and a
// hover change over the ribbons builds the hit ribbon's card, never the list of every ribbon. By hand, never in CI:
// `npx vitest run -c packages/seq/storybook/vitest.perf.config.ts`.

import { AcrossViews, GrapeSection, PeachSection } from '@molstar/seqx-examples/across/views';
import { describerFor } from '@molstar/seqx-viewer/describe';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { measureBytesPerStep, quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';

const fit = AcrossViews.find((v) => v.id === 'across-fit')!;

async function open(lock: boolean): Promise<SeqViewer> {
  const base = DefaultSeqViewerSpec();
  const v = createSeqViewer({
    ...base,
    ui: undefined,
    registry: [...(base.registry ?? []), ...fit.registry!({ delayMs: 0 })],
  });
  v.resize(1200);
  await v.load(fit.build());
  await v.settle();
  if (lock) await v.execute('bridge.lock', { bridge: '/views/main/peach-grape' });
  await v.settle();
  v.view.zoom(PeachSection, 4);
  v.view.zoom(GrapeSection, 4);
  await v.settle();
  v.draw();
  return v;
}

/** A pan of peach by a px, back and forth (its frame not drawn: the move alone). */
const panner = (v: SeqViewer) => (i: number) => void v.view.panPx(PeachSection, i & 1 ? 1 : -1);
/** Both sides panned by hand, each alone: what a carried pan must cost at most. */
const both = (v: SeqViewer) => (i: number) => {
  v.view.panPx(PeachSection, i & 1 ? 1 : -1);
  v.view.panPx(GrapeSection, i & 1 ? 1 : -1);
};

describe('a bridge on the pointer’s path', () => {
  it('a pan through a lock allocates what both sides panned by hand do', async ({ task }) => {
    const free = await open(false);
    const locked = await open(true);
    expect(locked.view.resolved.value!.bridges![0].lock).toBe('px');
    const reads = [both(free), panner(locked), both(free), panner(locked)].map(
      (s) => measureBytesPerStep(s).bytesPerStep,
    );
    const { r } = quickCase(task, 'bridge pan through a lock', panner(locked));
    free.retire();
    locked.retire();
    const plain = Math.min(reads[0], reads[2]);
    const carried = Math.min(reads[1], reads[3]);
    console.info(
      `pan: ${carried.toFixed(0)} bytes a step locked, ${plain.toFixed(0)} both by hand; the case ${r.note}`,
    );
    // the other side's move is its own pan's bytes: nothing for the carrying
    expect(carried - plain).toBeLessThanOrEqual(inject('perfThreshold'));
  });

  it('a hover change over the ribbons builds the hit ribbon’s card, never the list', async ({ task }) => {
    const v = await open(false);
    const bl = v.frame!.layout.bridges!.list[0];
    const [r0, r1] = bl.ribbons().filter((r) => !r.stub && r.ax1 - r.ax0 > 2 && r.ax0 > bl.box.x);
    expect(r1).toBeDefined();
    const y = bl.box.y + 2;
    const describer = describerFor({ ...v, layout: () => v.frame?.layout });
    const owner = {};
    const step = (i: number) => {
      const r = i & 1 ? r1 : r0;
      v.hover.set(v.hit((r.ax0 + r.ax1) / 2, y), owner);
      describer.describe(v.hover.current);
    };
    step(0);
    expect(describer.describe(v.hover.current)?.title).toContain('↔');
    const { r } = quickCase(task, 'bridge hover a ribbon and its card', step);
    v.retire();
    // a card's words and one ribbon: a few hundred bytes, never every ribbon's object (hundreds of them)
    expect(r.bytesPerStep, r.note).toBeLessThan(4096);
  });
});
