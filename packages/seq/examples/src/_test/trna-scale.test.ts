/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BuiltInLayerTypes, BuiltInLayouts } from '@molstar/seqx-core/render/built-ins';
import type { Node } from '@molstar/core/repr/tree';
import { DataStore } from '@molstar/seqx-core/data';
import {
  attachSources,
  bundledTextMetrics,
  type FrameInput,
  type GridLay,
  LightTokens,
  Renderer,
} from '@molstar/seqx-core/render';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { resolveViews } from '@molstar/seqx-core/view';
import { itWithin } from '@molstar/test-harness/timing';
import { describe, expect } from 'vitest';
import { ScaleSize } from '../trna/data.js';
import { scale } from '../trna/trna.js';

// The tRNA dot plot's scale variant (examples.md, tRNA: the dot plot's scale variant; performance.md, Scale): a generated
// 4,000 × 4,000 grid drawn at the frame target, a pan frame under 8 ms once its tiles are built, as the core's 2,000 × 2,000
// bound. The quick lane runs the same at 300 × 300, unbounded; the bound is the heavy lane's (`itWithin`, best of three).

const Size = { width: 1440, height: 1600 };
const Leader = '/views/main/rna';

const at = (start: number, end: number): FrameInput => ({
  size: Size,
  dpr: 2,
  settled: true,
  viewports: { [Leader]: { start, end } },
});

async function settled(spec: Node) {
  const model = resolveSpec(spec);
  const store = new DataStore();
  await attachSources(model, store);
  const renderer = new Renderer({
    types: BuiltInLayerTypes,
    layouts: BuiltInLayouts,
    model,
    view: resolveViews(model)[0],
    data: store,
    tokens: LightTokens,
    text: bundledTextMetrics(),
  });
  return {
    renderer,
    async frame(input: FrameInput) {
      for (let i = 0; i < 4; i++) {
        const f = renderer.frame(input);
        store.want(f.wants);
        await store.whenSettled();
        if (!renderer.frame(input).rebuilt) break;
      }
      return renderer.frame(input);
    },
  };
}

function panFrame(n: number) {
  itWithin('a pan frame once its tiles are built', async (within) => {
    const s = await settled(scale(n));
    const inputs = [0, 1, 2, 3].map((k) => at(n * 0.1 + 3 * k, n * 0.7 + 3 * k));
    for (const i of inputs) await s.frame(i);
    const lay = (await s.frame(inputs[0])).layout.sections[1].lays[0][0].lay as GridLay;
    expect(lay.images.length).toBeGreaterThan(0);
    let k = 0;
    const pan = () => {
      const f = s.renderer.frame(inputs[++k % inputs.length]);
      if (!f.rebuilt)
        throw new Error('a pan frame must lay out again, and did not: the bound would time a cached frame');
    };
    for (let i = 0; i < 8; i++) pan();
    within(8, pan);
  });
}

describe('the dot plot at 300 × 300', () => panFrame(300));

// the 4,000 × 4,000 plot is 16M values to build: the heavy lane's, the quick lane keeps the small one
describe('the dot plot at 4,000 × 4,000', { tags: ['heavy'] }, () => panFrame(ScaleSize));
