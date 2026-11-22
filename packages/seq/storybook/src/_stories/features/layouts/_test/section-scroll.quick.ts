/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Layouts story's section past its cap, composed as the story composes it: a wheel step over the section (its
// vertical delta scrolling the tracks below the ruler: the input's chain, the section's top set) allocates nothing, and
// a warm frame of the scrolled section neither. By hand, never in CI: `npx vitest run -c
// packages/seq/storybook/vitest.perf.config.ts`.

import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import { Caps, capConfig, capSpec, Genes } from '../compose.js';

async function open(): Promise<SeqViewer> {
  const cap = Caps[0];
  const base = DefaultSeqViewerSpec();
  const v = createSeqViewer({ ...base, ui: undefined, config: [...(base.config ?? []), ...capConfig(cap)] });
  v.resize(1200);
  await v.load(capSpec(cap));
  await v.settle();
  v.draw();
  return v;
}

describe('a section past its cap: the wheel and a warm frame', () => {
  it('a wheel step down and up over a gene (no frame): 0 bytes', async ({ task }) => {
    const v = await open();
    const s = v.frame!.layout.sections.find((x) => x.section.address === Genes)!;
    const g = s.tracks[s.tracks.length - 1].box;
    const x = g.x + 20;
    const y = s.placement.box.y + s.placement.box.height - 10;
    // the section a little scrolled, so each step moves it (down, then back up)
    v.view.setTop(Genes, 100);
    v.draw();
    const down = { deltaX: 0, deltaY: 1, deltaMode: 0 };
    const up = { deltaX: 0, deltaY: -1, deltaMode: 0 };
    const step = (i: number) => void v.input.wheel(x, y, (i & 1) === 0 ? down : up);
    for (let i = 0; i < 4; i++) step(i);
    expect(v.view.get(Genes)?.top).toBe(100);
    const { r } = quickCase(task, 'section wheel step', step, {});
    console.info(`section wheel step: ${r.bytesPerStep.toFixed(0)} bytes a step (${r.note})`);
    v.retire();
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });

  it('a warm frame of the scrolled section: 0 bytes', async ({ task }) => {
    const v = await open();
    v.view.setTop(Genes, 150);
    v.draw();
    const { r } = quickCase(task, 'section scrolled viewer frame', () => v.draw(), {});
    console.info(`section scrolled viewer frame: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
    v.retire();
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });
});
