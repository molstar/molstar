/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BuiltInLayerTypes, BuiltInLayouts } from '@molstar/seqx-core/render/built-ins';
import { readFileSync } from 'node:fs';
import { DataStore } from '@molstar/seqx-core/data';
import {
  attachSources,
  bundledTextMetrics,
  type FrameInput,
  LightTokens,
  Renderer,
  type RowsLay,
} from '@molstar/seqx-core/render';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { resolveViews } from '@molstar/seqx-core/view';
import { generatedRuntime } from '@molstar/seqx-examples/globins/data';
import { scale } from '@molstar/seqx-examples/globins/family';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import { itWithin } from '@molstar/test-harness/timing';
import { describe, expect } from 'vitest';

// The Globins' scale view against the frame targets (performance.md: open a 10,000 × 300 alignment under 300 ms to the
// first frame; a pan frame under 8 ms): the seed's profile sampled to 10,000 rows by a balanced tree, a first frame
// and pan frames at the fit (image tiles) and at letters. The bounds are the heavy lane's at 10,000 rows; the quick
// lane runs the same frames unbounded over 1,000.

const Rows = ['heavy', 'all'].includes((import.meta as { env?: { MODE?: string } }).env?.MODE ?? '') ? 10_000 : 1_000;
const Section = '/views/main/rows';
const Track = `${Section}/rows`;
const Size = { width: 1440, height: 800 };
const at = (bandTop: number, rowScale: number, start = 0, end = 300): FrameInput => ({
  size: Size,
  dpr: 2,
  settled: true,
  viewports: { [Section]: { start, end, bands: { [Track]: { bandTop, rowScale } } } },
});

const assets = {
  readText: async (uri: string) => readFileSync(`${PackageFolder}${uri}`, 'utf8'),
  readBytes: async (uri: string) => new Uint8Array(readFileSync(`${PackageFolder}${uri}`)),
};

async function settled() {
  const model = resolveSpec(scale({ rows: Rows }), { runtimes: [generatedRuntime()] });
  const view = resolveViews(model)[0];
  const store = new DataStore();
  await attachSources(model, store, { assets });
  const renderer = new Renderer({
    types: BuiltInLayerTypes,
    layouts: BuiltInLayouts,
    model,
    view,
    data: store,
    tokens: LightTokens,
    text: bundledTextMetrics(),
  });
  const frame = async (input: FrameInput) => {
    let f = renderer.frame(input);
    for (let i = 0; i < 8; i++) {
      store.want(f.wants);
      await store.whenSettled();
      f = renderer.frame(input);
      if (!f.rebuilt) break;
    }
    return f;
  };
  return { model, view, store, renderer, frame };
}

const rowsLay = (f: Awaited<ReturnType<Awaited<ReturnType<typeof settled>>['frame']>>) => {
  for (const s of f.layout.sections)
    for (const b of s.lays) for (const l of b ?? []) if (l.type.name === 'rows') return l.lay as RowsLay;
  throw new Error('no rows');
};

/** The median of `n` frames, each its own input (a pan). */
function median(renderer: Renderer, inputs: readonly FrameInput[], n = 30): number {
  const t: number[] = [];
  for (let k = 0; k < n; k++) {
    const t0 = performance.now();
    renderer.frame(inputs[k % inputs.length]);
    t.push(performance.now() - t0);
  }
  t.sort((a, b) => a - b);
  return t[t.length >> 1];
}

describe('the Globins scale view’s time', () => {
  itWithin('a first frame of the rows sampled from the seed, by tree', async (within) => {
    const s = await settled();
    const fit = at(0, 800 / Rows);
    await s.frame(fit);
    const first = () => {
      const r = new Renderer({
        types: BuiltInLayerTypes,
        layouts: BuiltInLayouts,
        model: s.model,
        view: s.view,
        data: s.store,
        tokens: LightTokens,
        text: bundledTextMetrics(),
      });
      const t0 = performance.now();
      r.frame(fit);
      return performance.now() - t0;
    };
    expect(within(300, first)).toBeGreaterThan(0);
  });

  itWithin('a pan frame at the fit (image tiles) and at letters', async (within) => {
    const s = await settled();
    const fit = [0, 1, 2, 3].map((k) => at(k, 0.06));
    const letters = [0, 20, 40, 60].map((k) => at(Rows / 2 + k / 20, 20, 100, 160));
    for (const i of [...fit, ...letters]) await s.frame(i);
    expect(rowsLay(await s.frame(letters[0])).letters).toBe(true);
    for (const inputs of [fit, letters]) {
      median(s.renderer, inputs);
      within(8 * 30, () => median(s.renderer, inputs));
    }
  });
});
