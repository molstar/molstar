/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `linkedPairs` behavior costs the hover nothing a frame past its first: on Complex › Contacts with the cyclin box
// hovered (its contacts lit), a hover moving within the box (a new hover version a move) and the frame drawn again
// allocate what they do without the behavior (its query cached by the element's key, the channel's value kept). By
// hand, never in CI: `npx vitest run -c packages/seq/storybook/vitest.perf.config.ts`.

import { ComplexViews } from '@molstar/seqx-examples';
import { nodeFixtureAssets } from '@molstar/seqx-examples/node-assets';
import { LinkedPairsChannel, LinkedPairsProvider } from '@molstar/seqx-viewer/behaviors';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { measureBytesPerStep, quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';

const contacts = ComplexViews.find((v) => v.id === 'contacts')!;

async function open(linked: boolean): Promise<SeqViewer> {
  const base = DefaultSeqViewerSpec();
  const v = createSeqViewer({
    ...base,
    ui: undefined,
    ports: { ...base.ports, assets: () => nodeFixtureAssets() },
    behaviors: linked ? base.behaviors : base.behaviors.filter((b) => b.provider !== LinkedPairsProvider),
  });
  v.resize(1440);
  await v.load(contacts.build());
  await v.settle();
  v.draw();
  return v;
}

const Owner = {};
/** A hit on the cyclin box at a position, as one record reused (the layer's hit is). */
const hit = {
  layer: '/views/main/complex/features/cyclin',
  system: '/data/cyclin',
  element: { owner: '/data/cyclin-features', key: 'box-n' },
  position: 40,
};

/** The hover moving within the cyclin box (a new version each move) and the frame drawn again. */
function mover(v: SeqViewer): (i: number) => void {
  v.hover.set(hit, Owner);
  v.draw();
  return (i) => {
    hit.position = 40 + (i & 1);
    v.hover.set(hit, Owner);
    v.draw();
  };
}

describe('linkedPairs on the hover’s path', () => {
  it('Complex › Contacts: a hover within the cyclin box and its frame, with the behavior and without', async ({
    task,
  }) => {
    const without = await open(false);
    const v = await open(true);
    const plainStep = mover(without);
    const step = mover(v);
    expect(v.highlight.get(LinkedPairsChannel)?.elements.keys.length).toBeGreaterThan(0);
    const reads = [plainStep, step, plainStep, step].map((s) => measureBytesPerStep(s).bytesPerStep);
    const plain = Math.min(reads[0], reads[2]);
    const linked = Math.min(reads[1], reads[3]);
    const { r } = quickCase(task, 'viewer hover with linked pairs', step);
    without.retire();
    v.retire();
    console.info(
      `hover within an element: ${linked.toFixed(0)} bytes a step with linked pairs, ${plain.toFixed(0)} without; the case ${r.bytesPerStep.toFixed(0)} (${r.note})`,
    );
    expect(linked - plain, reads.map((n) => n.toFixed(0)).join(' / ')).toBeLessThanOrEqual(inject('perfThreshold'));
  });
});
