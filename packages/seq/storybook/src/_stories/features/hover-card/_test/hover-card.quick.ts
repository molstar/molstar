/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The hover card costs the pointer nothing: with the card listed (the default spec), a warm pointer move over the
// letters and the hover's frame drawn again allocate what they do without it (the card is made on rest, never per move;
// one timer waits for the rest however many moves come). The frame's own bytes are the emphasis case's
// (`selection/_test/emphasis.quick.ts`), measured here on a viewer without the card and taken off. By hand, never in CI:
// `npx vitest run -c packages/seq/storybook/vitest.perf.config.ts`.

import { P53Views } from '@molstar/seqx-examples';
import { HoverCardProvider } from '@molstar/seqx-viewer/behaviors';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { measureBytesPerStep, quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import { cellPoint } from '../../selection/compose.js';

const tracks = P53Views.find((v) => v.id === 'tracks')!;

async function open(card: boolean): Promise<SeqViewer> {
  const base = DefaultSeqViewerSpec();
  const v = createSeqViewer({
    ...base,
    ui: undefined,
    behaviors: card ? base.behaviors : base.behaviors.filter((b) => b.provider !== HoverCardProvider),
  });
  v.resize(1440);
  await v.load(tracks.build());
  await v.settle();
  v.draw();
  return v;
}

/** A pointer move within residue 120 (the hover stays) and the hover's frame drawn again. */
function mover(v: SeqViewer): (i: number) => void {
  const [x, y] = cellPoint(v, 0, 0, 120);
  v.input.move(x, y);
  v.draw();
  return (i) => {
    v.input.move(x + (i & 1) * 0.5, y);
    v.draw();
  };
}

describe('the hover card on the pointer’s path', () => {
  it('p53 tracks: a pointer move within a residue and the hover frame, with the card and without', async ({ task }) => {
    // the same view and step, without the card and with it, each measured twice in turn (a second viewer in one
    // process reads a few tens of bytes apart from the first whichever it is: the lower of each is compared)
    const without = await open(false);
    const v = await open(true);
    expect(v.behaviors.has('hoverCard')).toBe(true);
    const plainStep = mover(without);
    const step = mover(v);
    const reads = [plainStep, step, plainStep, step].map((s) => measureBytesPerStep(s).bytesPerStep);
    const plain = { bytesPerStep: Math.min(reads[0], reads[2]) };
    const card = { bytesPerStep: Math.min(reads[1], reads[3]), note: reads.map((n) => n.toFixed(0)).join(' / ') };
    const { r } = quickCase(task, 'viewer pointer move with the hover card', step);
    without.retire();
    v.retire();
    console.info(
      `pointer move: ${card.bytesPerStep.toFixed(0)} bytes a step with the card, ${plain.bytesPerStep.toFixed(0)} without; the case ${r.bytesPerStep.toFixed(0)} (${r.note})`,
    );
    expect(card.bytesPerStep - plain.bytesPerStep, card.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });
});
