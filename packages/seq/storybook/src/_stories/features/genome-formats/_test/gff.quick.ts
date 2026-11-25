/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// TP53's gene models from Ensembl's GFF3 (Features/Genome formats/GFF3), composed as the story's view: a warm draw of
// the gene line (boxes, hairlines, chevrons) and a warm hover on a transcript allocate nothing; a pan frame (a layout
// and a draw) a bounded size. By hand, never in CI: `npx vitest run -c packages/seq/storybook/vitest.perf.config.ts`.

import { BuiltInLayerTypes, BuiltInLayouts } from '@molstar/seqx-core/render/built-ins';
import { DataStore, SeqRuntime } from '@molstar/seqx-core/data';
import {
  bundledTextMetrics,
  drawView,
  type FrameInput,
  hitView,
  LightTokens,
  openModel,
  Renderer,
} from '@molstar/seqx-core/render';
import { GenomeFormatsViews } from '@molstar/seqx-examples';
import { nodeFixtureAssets } from '@molstar/seqx-examples/node-assets';
import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';

const view = GenomeFormatsViews.find((v) => v.id === 'gff-tp53')!;
const Size = { width: 1440, height: 640 };
const Steps = 64;

const LayoutOptions = {
  warmSteps: 300,
  steps: 128,
  windows: 3,
  maxWarmSteps: 2000,
} as const;

async function settled() {
  const runtimes = [
    SeqRuntime,
    ...(view.registry?.({ delayMs: 0 }) ?? []).flatMap((e) => (e.kind === 'format-runtime' ? [e.runtime] : [])),
  ];
  const store = new DataStore();
  const reader = nodeFixtureAssets().reader;
  const opened = await openModel(view.build(), store, { runtimes, assets: reader });
  const renderer = new Renderer({
    types: BuiltInLayerTypes,
    layouts: BuiltInLayouts,
    model: opened.model,
    view: opened.view!,
    data: store,
    tokens: LightTokens,
    text: bundledTextMetrics(),
  });
  // a pan over the section's 30 kb (its own positions from 0): 20 kb wide, out and back
  const inputs: FrameInput[] = Array.from({ length: Steps }, (_, i) => {
    const t = i < Steps / 2 ? i : Steps - i;
    return { size: Size, dpr: 2, settled: true, viewport: { start: 2000 + t * 100, end: 22_000 + t * 100 } };
  });
  for (const input of [inputs[0], inputs[Steps / 2]])
    for (let k = 0; k < 4; k++) {
      store.want(renderer.frame(input).wants);
      await store.whenSettled();
    }
  return { renderer, inputs };
}

describe('TP53’s gene models: warm draws and hovers allocate nothing, a pan a bounded size', () => {
  it('pan (layout and draw)', async ({ task }) => {
    const { renderer, inputs } = await settled();
    const { r } = quickCase(task, 'gff tp53 pan', (i) => void renderer.frame(inputs[i % Steps]), LayoutOptions);
    console.info(`gff tp53 pan: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
    // 279 KB (2026-11-21, M-series, Node 26) with half again as headroom: three tracks of about 50 transcripts and their
    // 800 pieces read from the tiles (the items' reading, shared by every items layer), their lanes, boxes, hairlines
    // and chevrons; the same data drawn as blocks 395 KB
    expect(r.bytesPerStep, r.note).toBeLessThan(420_000);
  });

  it('draw: 0 bytes warm (B6)', async ({ task }) => {
    const { renderer, inputs } = await settled();
    const { layout, lists } = renderer.frame(inputs[0]);
    const { r } = quickCase(task, 'gff tp53 draw', () => drawView(layout, lists), {});
    console.info(`gff tp53 draw: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });

  it('hover on a transcript and off: 0 bytes warm (B6)', async ({ task }) => {
    const { renderer, inputs } = await settled();
    const f = renderer.frame(inputs[0]);
    const track = f.layout.sections[0].tracks.find((t) => t.track.label === 'Transcripts')!;
    const lay = track.layers[0].lay as {
      m: number;
      mx0: Float32Array;
      mx1: Float32Array;
      mlane: Int32Array;
      laneY: number;
      laneHeight: number;
    };
    const hits = [0, 1, 2, 3].map((k) => {
      const j = Math.floor((k * lay.m) / 4);
      return hitView(
        f.layout,
        (lay.mx0[j] + lay.mx1[j]) / 2,
        lay.laneY + lay.mlane[j] * lay.laneHeight + lay.laneHeight / 2,
      );
    });
    expect(hits.every((h) => h?.element)).toBe(true);
    const emphasis: { hover?: (typeof hits)[number] } = {};
    const at = { ...inputs[0], emphasis };
    const { r } = quickCase(
      task,
      'gff tp53 hover',
      (i) => {
        emphasis.hover = hits[i & 3];
        renderer.frame(at);
      },
      {},
    );
    console.info(`gff tp53 hover: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });
});
