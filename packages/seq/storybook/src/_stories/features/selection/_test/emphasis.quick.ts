/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The viewer's warm emphasis frame allocates nothing: with a selection projected onto the gene, a highlight channel and
// a hover projected from the protein, a frame drawn again (the viewer's emphasis records updated, the projections read
// from their caches, the renderer's emphasis redraw) makes no object. A change of what is selected or hovered projects
// once (bounded, not here). By hand, never in CI: `npx vitest run -c packages/seq/storybook/vitest.perf.config.ts`.

import { P53Views } from '@molstar/seqx-examples';
import { nodeFixtureAssets } from '@molstar/seqx-examples/node-assets';
import { GeneFile } from '@molstar/seqx-examples/p53/data';
import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import { cellPoint, clickAt, composeViewer } from '../compose.js';

const onGene = P53Views.find((v) => v.id === 'protein-on-gene')!;

describe('the viewer’s warm emphasis frame', () => {
  it('p53 on its gene: selection, highlight and hover, projected, drawn again', async ({ task }) => {
    const gene = await nodeFixtureAssets().reader.readText(GeneFile, new AbortController().signal);
    const v = composeViewer(onGene, {
      ports: {
        assets: () => ({
          reader: { readText: async () => gene, readBytes: async () => new TextEncoder().encode(gene) },
        }),
      },
    });
    v.resize(1440);
    await v.load(onGene.build());
    await v.settle();
    v.view.setRange('/views/main/protein', 150, 250);
    await v.settle();
    clickAt(v, cellPoint(v, 0, 0, 170));
    v.addresses.highlightAddress('/data/p53/(range :from 190 :to 200)');
    v.input.move(...cellPoint(v, 0, 0, 180));
    v.draw();
    const { r } = quickCase(task, 'viewer emphasis frame', () => void v.draw(), {});
    console.info(`viewer emphasis frame: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
    v.retire();
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });
});
