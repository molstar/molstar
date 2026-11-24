/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Features › Coordinates › Protein Over Its Coding DNA: a hover moving over the domain's connector (its hairline drawn
// in the hover ink) and the frame drawn again allocate what the same hover on an item without one does (the exons').
// By hand, never in CI: `npx vitest run -c packages/seq/storybook/vitest.perf.config.ts`.

import { CoordinatesViews } from '@molstar/seqx-examples';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { measureBytesPerStep, quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';

const view = CoordinatesViews.find((v) => v.id === 'coding-dna')!;
const Owner = {};

async function open(): Promise<SeqViewer> {
  const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
  v.resize(1440);
  await v.load(view.build());
  await v.settle();
  v.draw();
  return v;
}

/** The hover moving along a gap or an item (a new version each move) and the frame drawn again. */
function mover(v: SeqViewer, layer: string, key: string, positions: [number, number]): (i: number) => void {
  const hit = { layer, system: '/data/gene', element: { owner: `/data/${layer.split('/').pop()}`, key }, position: 0 };
  const at = (i: number) => {
    hit.position = positions[i & 1];
    v.hover.set(hit, Owner);
    v.draw();
  };
  at(0);
  return at;
}

describe('the connector on the hover’s path', () => {
  it('a hover moving over the domain’s connector and its frame allocate what an item’s hover does', async ({
    task,
  }) => {
    const v = await open();
    // the intron's stretch holds the connector (46–76): two positions inside it; the exons' own marks for the control
    const connector = mover(v, '/views/main/gene/domain/domain', 'D', [55, 56]);
    const item = mover(v, '/views/main/gene/exons/exons', 'E1', [20, 21]);
    const reads = [item, connector, item, connector].map((s) => measureBytesPerStep(s).bytesPerStep);
    const plain = Math.min(reads[0], reads[2]);
    const joined = Math.min(reads[1], reads[3]);
    const { r } = quickCase(task, 'viewer hover over a connector', connector);
    v.retire();
    console.info(
      `hover over a connector: ${joined.toFixed(0)} bytes a step, ${plain.toFixed(0)} over an item; the case ${r.bytesPerStep.toFixed(0)} (${r.note})`,
    );
    expect(joined - plain, reads.map((n) => n.toFixed(0)).join(' / ')).toBeLessThanOrEqual(inject('perfThreshold'));
  });
});
