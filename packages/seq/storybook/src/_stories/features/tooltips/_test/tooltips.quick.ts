/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The spec's tooltips cost the pointer nothing (the viewer's _spec "Hover card", Tooltips; core "Render", Tooltips): a
// warm pointer move and the hover's frame, with a tooltip table of 1,000,000 ranges over the sequence, allocate what
// they do without it (the card's providers run at rest, never per move). The card's show with the tooltips is recorded
// cold (the first: the entry's interval index built, the parse made) and warm (the same place again: the index and the
// parse cached). By hand, never in CI: `npx vitest run -c packages/seq/storybook/vitest.perf.config.ts`.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { describerFor } from '@molstar/seqx-viewer/describe';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { measureBytesPerStep, quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';

const Length = 2_000_000;
const Rows = 1_000_000;

/** A long protein-like sequence with a domains track, and (`tips`) a tooltip table of a million short ranges on it. */
function spec(tips: boolean): Node {
  const b = createSeqSpec({ title: 'tooltips' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(Length / 4) });
  const d = b.data.ranges('d', { target: s, rows: { start: [100, 5000], end: [900, 9000], label: ['a', 'b'] } });
  if (tips) {
    const start = new Int32Array(Rows);
    const end = new Int32Array(Rows);
    for (let i = 0; i < Rows; i++) {
      start[i] = (i * 2) % Length;
      end[i] = start[i] + 3;
    }
    const tooltip = Array.from({ length: Rows }, (_, i) =>
      i % 3 === 0 ? 'a **note**' : i % 3 === 1 ? 'another' : 'third',
    );
    b.views.tooltip('notes', { target: s, rows: { start, end, tooltip }, label: 'Note' });
  }
  b.views.view('v', (v) =>
    v.section('s', { of: s, range: [0, 600] }, (sec) => {
      sec.track('letters').layer('l', { type: 'letters' });
      sec.track('d').layer('d', { type: 'ranges', data: d });
    }),
  );
  return b.build();
}

async function open(tips: boolean): Promise<SeqViewer> {
  const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
  v.resize(1440);
  await v.load(spec(tips));
  await v.settle();
  v.draw();
  return v;
}

/** A pointer move within one position of the letters (the hover stays) and the hover's frame drawn again. */
function mover(v: SeqViewer): (i: number) => void {
  const s = v.frame!.layout.sections[0];
  const t = s.tracks[0].box;
  const x = (s.axis!.x(300) + s.axis!.x(301)) / 2;
  const y = t.y + t.height / 2;
  v.input.move(x, y);
  v.draw();
  return (i) => {
    v.input.move(x + (i & 1) * 0.5, y);
    v.draw();
  };
}

describe('tooltips on the pointer’s path', () => {
  it('a warm pointer move with a million-range tooltip table, and a cold and a warm card show', async ({ task }) => {
    const without = await open(false);
    const v = await open(true);
    const plainStep = mover(without);
    const step = mover(v);
    const reads = [plainStep, step, plainStep, step].map((s) => measureBytesPerStep(s).bytesPerStep);
    const plain = Math.min(reads[0], reads[2]);
    const tips = Math.min(reads[1], reads[3]);
    const { r } = quickCase(task, 'viewer pointer move with a million-range tooltip table', step);
    // the card at the resting place: cold (index built, parse made), then warm (the same place again)
    const d = describerFor({ ...v, layout: () => v.frame?.layout });
    const s = v.frame!.layout.sections[0];
    const place = { section: s.section.address, layer: s.tracks[0].layers[0].key, system: '/data/s', position: 300 };
    const before = process.memoryUsage().heapUsed;
    const cold = d.describe(place);
    const coldBytes = process.memoryUsage().heapUsed - before;
    expect(cold?.rows.some((r) => r.label === 'Note')).toBe(true);
    const warm = measureBytesPerStep(() => void d.describe(place)).bytesPerStep;
    const plainCard = measureBytesPerStep(
      () =>
        void describerFor({ ...without, layout: () => without.frame?.layout }).describe({
          ...place,
          layer: without.frame!.layout.sections[0].tracks[0].layers[0].key,
        }),
    ).bytesPerStep;
    without.retire();
    v.retire();
    console.info(
      `pointer move: ${tips.toFixed(0)} bytes a step with a million-range table, ${plain.toFixed(0)} without; the case ${r.bytesPerStep.toFixed(0)} (${r.note}). Card show: cold about ${Math.max(0, coldBytes).toFixed(0)} bytes (the index built), warm ${warm.toFixed(0)} a show with the tooltips, ${plainCard.toFixed(0)} without`,
    );
    expect(tips - plain).toBeLessThanOrEqual(inject('perfThreshold'));
  });
});
