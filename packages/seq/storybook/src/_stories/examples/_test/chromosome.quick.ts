/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// p53's chromosome: a warm pan's frame (a layout and a draw) allocates at most a bound that follows the screen, not the
// chromosome's 83 Mb, at the fit and at a density between two of its levels (512 positions a device pixel, between
// 256 and 1,024). By hand, never in CI: `npx vitest run -c packages/seq/storybook/vitest.perf.config.ts`.

import { DataStore, SeqRuntime } from '@molstar/seqx-core/data';
import { bundledTextMetrics, type FrameInput, LightTokens, openModel, Renderer } from '@molstar/seqx-core/render';
import { P53Views } from '@molstar/seqx-examples';
import { chromosomeRuntime, DefaultChromosome } from '@molstar/seqx-examples/p53/chromosome-data';
import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, it } from 'vitest';

const view = P53Views.find((v) => v.id === 'chromosome')!;
const Length = DefaultChromosome.length;
const Size = { width: 1440, height: 400 };
const Steps = 64;

// a layout allocates kilobytes by design (the wants, the tiles' keys, the genes' rows read): the control's 88 bytes are
// lost in that
const LayoutOptions = {
  warmSteps: 300,
  steps: 128,
  windows: 3,
  maxWarmSteps: 2000,
  control: false,
  why: 'a layout allocates kilobytes: they hide the control',
} as const;

/**
 * Bytes a pan frame may allocate, measured (2026-11-19, M-series, Node 26) with half again as headroom: at the fit
 * 36 KB (the genes a density, a count a bin: 348 KB when they came as 1,177 rows), between levels 37 KB, the signal's
 * tiles and their keys.
 */
const LayoutBound: Record<string, number> = { fit: 54_000, 'between two levels': 56_000 };

/** The view settled over a pan of `span` positions from `start`: the inputs of its walk (a pan out and back). */
async function pan(start: number, span: number): Promise<{ renderer: Renderer; inputs: FrameInput[] }> {
  const runtimes = [SeqRuntime, chromosomeRuntime()];
  const store = new DataStore();
  const opened = await openModel(view.build(), store, { runtimes });
  const renderer = new Renderer({
    model: opened.model,
    view: opened.view!,
    data: store,
    tokens: LightTokens,
    text: bundledTextMetrics(),
  });
  const step = span / 2000;
  const inputs: FrameInput[] = Array.from({ length: Steps }, (_, i) => {
    const t = i < Steps / 2 ? i : Steps - i;
    return { size: Size, dpr: 2, settled: true, viewport: { start: start + t * step, end: start + span + t * step } };
  });
  for (const input of [inputs[0], inputs[Steps / 2]])
    for (let k = 0; k < 4; k++) {
      store.want(renderer.frame(input).wants);
      await store.whenSettled();
    }
  return { renderer, inputs };
}

describe('p53’s chromosome: a warm pan allocates a bounded size', () => {
  const cases = [
    { name: 'fit', start: 0, span: Length - Length / 30 },
    { name: 'between two levels', start: 20_000_000, span: 512 * Size.width * 2 },
  ];
  for (const c of cases)
    it(`${c.name}: pan (layout and draw)`, async ({ task }) => {
      const { renderer, inputs } = await pan(c.start, c.span);
      const { r } = quickCase(
        task,
        `chromosome ${c.name} pan`,
        (i) => void renderer.frame(inputs[i % Steps]),
        LayoutOptions,
      );
      console.info(`chromosome ${c.name} pan: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
      expect(r.bytesPerStep, r.note).toBeLessThan(LayoutBound[c.name]);
    });
});
