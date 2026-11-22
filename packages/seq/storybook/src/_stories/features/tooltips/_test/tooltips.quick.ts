/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The spec's tooltips cost the pointer nothing and a warm card about nothing (the viewer's _spec "Hover card",
// Tooltips; core "Render", Tooltips): with a tooltip table of 1,000,000 ranges from a BED source read by region, a warm
// pointer move and the hover's frame allocate what they do without it (the card's providers run at rest, never per
// move), and a warm card show held at one place allocates within a few hundred bytes of the card without the tooltip
// (its rows and texts made once a row, the lookup's buffers reused). The card's first show is recorded (the index
// built, the texts made). By hand, never in CI: `npx vitest run -c packages/seq/storybook/vitest.perf.config.ts`.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { BedRuntime } from '@molstar/seqx-core/runtimes/bed';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { type Place, describerFor } from '@molstar/seqx-viewer/describe';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { measureBytesPerStep, quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';

const Length = 2_000_000;
const Rows = 1_000_000;
/**
 * A warm card with the tooltip's rows against the card without (measured 416 B): the live region's longer sentence (the
 * rows' phrases joined) and the lookup's few small temporaries; its rows, texts, keys and buffers are kept.
 */
const WarmCardSlack = 512;

/** A BED of a million short ranges over `chr1`, each named (its id: what the tooltip reads), Markdown in it escaped. */
function bed(): string {
  const lines: string[] = [];
  for (let i = 0; i < Rows; i++) {
    const s = (i * 2) % Length;
    lines.push(`chr1\t${s}\t${s + 3}\tn${i}_**note**`);
  }
  return lines.join('\n');
}

/** A long DNA sequence with a domains track, and (`tips`) a tooltip over a BED source of a million ranges on it. */
function spec(tips: boolean, text: string): Node {
  const b = createSeqSpec({ title: 'tooltips' });
  const s = b.data.sequence('s', { alphabet: 'dna', name: 'chr1', letters: 'ACGT'.repeat(Length / 4) });
  const d = b.data.ranges('d', { target: s, rows: { start: [100, 5000], end: [900, 9000], label: ['a', 'b'] } });
  if (tips) {
    const notes = b.data.source('notes', { format: 'bed', value: text as never });
    b.views.tooltip('notes', { target: s, source: notes, field: 'id', label: 'Note' });
  }
  b.views.view('v', (v) =>
    v.section('s', { of: s, range: [0, 600] }, (sec) => {
      sec.track('letters').layer('l', { type: 'letters' });
      sec.track('d').layer('d', { type: 'ranges', data: d });
    }),
  );
  return b.build();
}

async function open(tips: boolean, text: string): Promise<SeqViewer> {
  const base = DefaultSeqViewerSpec();
  const v = createSeqViewer({
    ...base,
    ui: undefined,
    registry: [...(base.registry ?? []), { kind: 'format-runtime', runtime: BedRuntime }],
  });
  v.resize(1440);
  await v.load(spec(tips, text));
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

/** The card at the letters' position 300, held there: its describer and place. */
function held(v: SeqViewer): { show: () => void; card: () => ReturnType<ReturnType<typeof describerFor>['describe']> } {
  const d = describerFor({ ...v, layout: () => v.frame?.layout });
  const s = v.frame!.layout.sections[0];
  const place: Place = {
    section: s.section.address,
    layer: s.tracks[0].layers[0].key,
    system: '/data/s',
    position: 300,
  };
  return { show: () => void d.describe(place), card: () => d.describe(place) };
}

describe('tooltips on the pointer’s path', () => {
  it('a million ranges from a source read by region: a warm pointer move, a warm card held', async ({ task }) => {
    const text = bed();
    const without = await open(false, text);
    const v = await open(true, text);
    const plainStep = mover(without);
    const step = mover(v);
    const reads = [plainStep, step, plainStep, step].map((s) => measureBytesPerStep(s).bytesPerStep);
    const plain = Math.min(reads[0], reads[2]);
    const tips = Math.min(reads[1], reads[3]);
    const { r } = quickCase(task, 'viewer pointer move with a million-range tooltip table', step);
    // the card held at the place: the first show (the index built, the texts made), then warm, against the plain card
    const card = held(v);
    const plainCard = held(without);
    const first = card.card();
    expect(first?.rows.filter((r) => r.label === 'Note').map((r) => r.value)).toEqual([
      'n149_**note**@298-301',
      'n150_**note**@300-303',
    ]);
    const before = process.memoryUsage().heapUsed;
    card.show();
    const coldBytes = process.memoryUsage().heapUsed - before;
    const shows = [plainCard.show, card.show, plainCard.show, card.show].map(
      (s) => measureBytesPerStep(s).bytesPerStep,
    );
    const plainShow = Math.min(shows[0], shows[2]);
    const tipShow = Math.min(shows[1], shows[3]);
    const { r: c } = quickCase(task, 'viewer warm card show held, with a million-range tooltip table', card.show);
    without.retire();
    v.retire();
    console.info(
      `pointer move: ${tips.toFixed(0)} bytes a step with a million-range table, ${plain.toFixed(0)} without; the case ${r.bytesPerStep.toFixed(0)} (${r.note}). Card held: warm ${tipShow.toFixed(0)} a show with the tooltip, ${plainShow.toFixed(0)} without (the case ${c.bytesPerStep.toFixed(0)}); a second show about ${Math.max(0, coldBytes).toFixed(0)}`,
    );
    expect(tips - plain).toBeLessThanOrEqual(inject('perfThreshold'));
    expect(tipShow - plainShow).toBeLessThanOrEqual(WarmCardSlack);
  });
});
