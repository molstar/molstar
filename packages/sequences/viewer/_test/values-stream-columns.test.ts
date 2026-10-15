/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// VT-06's draw side through the `values` renderer (LC-03; seqspec:sources#the-viewers-side): a streamed values track
// on a cyclic block shown from the middle, forward and reversed, at less than a device px a position, asks its stream
// each device column's runs, never the hull of a column's two ends (across the origin, the whole block but the
// column).

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { expect, it } from 'vitest';
import type { StreamView } from '../src/draw/api.js';
import { KernelRenderers, type LayerRenderer } from '../src/index.js';
import { LayerFrame } from '../src/testing/index.js';
import { opened } from './layout-fixtures.js';

const L = 4000;
function doc(direction: 1 | -1) {
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: 'a', alphabet: 'dna', cyclic: true, value: 'ACGT'.repeat(L / 4) });
  const value = Array.from({ length: L }, (_, i) => i % 7);
  const vals = b.data.annotation({ ref: 'vals', geometry: 'values', target: a, value: { value } } as never);
  const v = b.view({ ref: 'v', label: 'V' }).section({ ref: 'v.main', direction });
  v.block({ sequence: a, origin: L / 2 });
  v.track({ ref: 't', label: 'Values' }).values({ data: vals });
  return b.getDocument();
}

for (const direction of [1, -1] as const)
  it(`${direction < 0 ? 'reversed' : 'forward'}: each column asks its runs, the origin's two apart`, async () => {
    const { ctx } = await opened(doc(direction), { renderers: [...KernelRenderers] });
    const view = ctx.graph.assembled.value?.view;
    if (!view) throw new Error('no view');
    const k = view.prepared.findIndex((p) => p.renderer === 'values');
    // Streamed: the layer reads its stream's regions (a stream over the block's system, every bin loaded).
    (view.prepared[k].state as { value: unknown[] }).value[0] = null;
    const asked: [number, number][] = [];
    const stream = {
      source: { systems: [0] },
      version: 0,
      budget: 0,
      values: (_b: number, i0: number, i1: number, out: Float64Array) => {
        asked.push([i0, i1]);
        out[0] = out[1] = 1;
        out[2] = out[3] = 1;
        return true;
      },
    } as unknown as StreamView;
    ctx.streams.get = () => stream;
    const vp = ctx.viewports.get('v.main');
    vp?.fit();
    const f = new LayerFrame(ctx, KernelRenderers.find((r) => r.name === 'values') as LayerRenderer, 'v.main', k);
    f.draw(f.direct);
    // A column holds 4 positions at the fit (0.25 px a position).
    expect(asked.length).toBeGreaterThan(500);
    const widest = Math.max(...asked.map(([i0, i1]) => i1 - i0));
    expect(widest).toBeLessThanOrEqual(8);
    // The column across the origin asked its two runs, at either end of the system.
    expect(asked.some(([i0]) => i0 === 0) && asked.some(([, i1]) => i1 === L)).toBe(true);
    ctx.dispose();
  });
