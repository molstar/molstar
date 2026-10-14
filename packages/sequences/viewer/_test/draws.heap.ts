/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// RN1's draws through the conformance kit's heap guard (seqspec:viewer#extension-points;
// seqspec:viewer#the-render-backend): a frame of `letters`, `color`, `values`, `ruler` and `emphasis` allocates nothing
// after warm-up, at each of the kit's zooms and below the row's minimum, linear and wrapped, on documents that reach
// every path: every type on two chains, case 14's residues over their codons, 40,000 positions zoomed out (classes,
// values' pyramids), dense values and a row field under a continuous palette (the mean), many emphasized positions
// (density).
//
// Each case loads the viewer afresh (`vi.resetModules`): in one process, after other documents deoptimized it, V8 can
// leave a draw in its middle tier (Maglev) for good, where the vitest module runner's export getters and Math builtins
// box doubles (about 110 to 145 bytes a frame for the ruler and emphasis; RN1). The guard measures the optimized frame;
// the same draws under `--no-turbofan`, outside vitest, are measured flat as well (RN1).
// The draws are large functions, called once a frame: each guard warms its frame until it stops allocating
// (`warmFrame`, heap.ts says why). The positive control goes through the kit's `heap` hook on letters over every type:
// at each frame the same warm-up and measurement of the frame allocating an object each (`leaking`) must catch it.

import { describe, expect, it, vi } from 'vitest';
import type { LayerRenderer } from '../src/draw/api.js';
import { controlCaught, describeScaling, FRAME_N, heapScaling, isFlat, leaking, warmFrame } from './heap.js';

const N = FRAME_N;
/** The guards' measurement: warmed until the frame stops allocating (`warmFrame`), then `N` frames and `4N` measured. */
const measure = (step: (i: number) => void) => {
  let i = 0;
  const frame = () => step(i++);
  warmFrame(frame);
  return heapScaling(N, frame);
};
const guard = (frame: () => void) => {
  const r = measure(frame);
  if (!isFlat(r)) throw new Error(describeScaling(r, N));
};
/** The positive control as the kit's `heap` hook: a frame's check fails unless its leaking frame is caught. */
const control = (frame: () => void) => {
  if (!controlCaught(() => measure(leaking(frame)))) throw new Error('the leaking frame was not caught');
};

const DOCS = ['all types', 'case 14', '40,000 positions', 'means'] as const;
/** Each renderer on the documents that reach its paths. */
const CASES: Record<string, readonly (typeof DOCS)[number][]> = {
  letters: ['all types', 'case 14', '40,000 positions'],
  color: ['all types', 'means'],
  values: ['all types', '40,000 positions'],
  ruler: ['all types', 'case 14', 'means'],
  emphasis: ['all types', 'means'],
};

async function run(name: string, doc: (typeof DOCS)[number], heap = guard) {
  vi.resetModules();
  const { KernelRenderers } = await import('../src/index.js');
  const { rendererConformance } = await import('../src/testing/index.js');
  const { asText } = await import('./conformance-fixtures.js');
  const { meanDoc } = await import('./draw-fixtures.js');
  const { opened } = await import('./layout-fixtures.js');
  const { allTypes, gene, scaleDoc } = await import('./prepare-fixtures.js');
  const make = {
    'all types': () => allTypes(),
    'case 14': gene,
    '40,000 positions': () => scaleDoc({ length: 20_000 }),
    means: () => meanDoc(10_000),
  }[doc];
  const all = KernelRenderers as LayerRenderer[];
  const entry = all.find((r) => r.name === name) as LayerRenderer;
  return rendererConformance(entry, {
    open: async ({ geometry, entry: has }) =>
      (
        await opened(geometry === 'wrapped' ? asText(make()) : make(), {
          renderers: has ? all : all.filter((r) => r !== entry),
        })
      ).ctx,
    heap,
    items: false,
  });
}

describe("RN1's frames allocate nothing after warm-up", () => {
  for (const [name, docs] of Object.entries(CASES))
    for (const doc of docs)
      it(`${name} on ${doc}`, async () => {
        const heap = (await run(name, doc)).filter((c) => c.check === 'heap');
        expect(heap.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
      }, 240_000);
});

describe('the guard catches a draw that allocates', () => {
  it('letters on all types, an object allocated a frame (the positive control)', async () => {
    const heap = (await run('letters', 'all types', control)).filter((c) => c.check === 'heap');
    expect(heap.length).toBeGreaterThan(0);
    expect(heap.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
  }, 240_000);
});
