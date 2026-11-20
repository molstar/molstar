/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The headers behavior costs the pointer nothing: a warm pointer moving down the label column of 200 tracks (each step
// onto the next track's label, so the actions of another header are drawn) and its frame allocate what they do without
// the behavior (its slots are an emphasis redraw, its geometry compared by the layout's identity). By hand, never in
// CI: `npx vitest run -c packages/seq/storybook/vitest.perf.config.ts`.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { HeadersProvider } from '@molstar/seqx-viewer/behaviors/headers';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { measureBytesPerStep, quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';

const Tracks = 200;

function spec() {
  const b = createSeqSpec({ title: '200 tracks' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQETFSDLWKLLPENN'.repeat(10) });
  b.views.view('main', (v) => {
    const s = v.section('protein', { of: p });
    for (let i = 0; i < Tracks; i++) s.track(`t${i}`, { label: `Track ${i}` }).layer('l', { type: 'letters' });
  });
  return b.build();
}

async function open(headers: boolean): Promise<SeqViewer> {
  const base = DefaultSeqViewerSpec();
  const v = createSeqViewer({
    ...base,
    ui: undefined,
    behaviors: headers ? base.behaviors : base.behaviors.filter((b) => b.provider !== HeadersProvider),
  });
  v.resize(1440);
  await v.load(spec());
  await v.settle();
  v.draw();
  return v;
}

/** A pointer onto the next track's label each step, and its frame drawn. */
function mover(v: SeqViewer): (i: number) => void {
  const tracks = v.frame!.layout.sections[0].tracks;
  const ys = new Float64Array(tracks.map((t) => t.well.y + t.well.height / 2));
  const x = v.context.renderConfig.inset + 4;
  return (i) => {
    v.input.move(x, ys[i % ys.length]);
    v.draw();
  };
}

describe('the headers on the pointer’s path', () => {
  it('200 tracks: a pointer down the label column and the actions’ frame, with the behavior and without', async ({
    task,
  }) => {
    const without = await open(false);
    const v = await open(true);
    expect(v.behaviors.has('headers')).toBe(true);
    const plainStep = mover(without);
    const step = mover(v);
    const layouts = v.counts.layouts;
    const reads = [plainStep, step, plainStep, step].map((s) => measureBytesPerStep(s).bytesPerStep);
    expect(v.counts.layouts).toBe(layouts);
    const plain = Math.min(reads[0], reads[2]);
    const headers = { bytesPerStep: Math.min(reads[1], reads[3]), note: reads.map((n) => n.toFixed(0)).join(' / ') };
    const { r } = quickCase(task, 'viewer pointer down 200 track headers', step);
    without.retire();
    v.retire();
    console.info(
      `header pointer move: ${headers.bytesPerStep.toFixed(0)} bytes a step with the behavior, ${plain.toFixed(0)} without; the case ${r.bytesPerStep.toFixed(0)} (${r.note})`,
    );
    expect(headers.bytesPerStep - plain, headers.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });
});
