/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The cases of RN2's hot paths, for their heap guard and its quick twin (rn2-hot-paths.heap.ts,
// rn2-hot-paths.quick.ts): the shared density (points stacked and majority, raw and from the count pyramid; ranges'
// coverage; a block counted through its items beside one through its pyramid) through the test entries `ranges` and
// `points`, `pairs`' draw and the matrix's tiles through the conformance kit (`kitRun`, the frame handed to the kit's
// `heap` hook at each of its zooms and geometries), the pairs' hit over a grid of points (`pairsHit`) and the overlay's
// lit pairs through the channels' runs and the end index (`litPairsFrame`). Each call loads the viewer afresh
// (`vi.resetModules`; draws.heap.ts says why: Maglev and the module runner).

import { vi } from 'vitest';
import type { LayerRenderer } from '../src/index.js';

export type Case = { entry: 'points' | 'ranges' | 'pairs' | 'matrix'; doc: string; mixed?: boolean };
export const CASES: Case[] = [
  { entry: 'points', doc: 'points, raw' },
  { entry: 'points', doc: 'points, pyramid, majority' },
  { entry: 'points', doc: 'points, pyramid, majority', mixed: true },
  { entry: 'ranges', doc: 'ranges, raw' },
  { entry: 'ranges', doc: 'ranges, pyramid' },
  { entry: 'pairs', doc: 'pairs' },
  { entry: 'matrix', doc: 'matrix' },
];

/** A case's name, as the guards and the quick twin give it. */
export const caseName = (c: Case) => `${c.entry} on ${c.doc}${c.mixed ? ', a block through its items' : ''}`;

/**
 * The case through the conformance kit, `heap` getting each frame the kit checks and where it is (the geometry and the
 * frame's place in the kit's order there: `linear #3`). Omitted: the kit's other checks only.
 */
export async function kitRun(c: Case, heap?: (frame: () => void, where: string) => void | Promise<void>) {
  vi.resetModules();
  const { KernelRenderers } = await import('../src/index.js');
  const { rendererConformance } = await import('../src/testing/index.js');
  const { opened } = await import('./layout-fixtures.js');
  const fx = await import('./rn2-fixtures.js');
  const { asText } = await import('./conformance-fixtures.js');
  const { OrderedSet } = await import('@molstar/seq-core/model');
  const make = {
    'points, raw': () => fx.densityDoc(),
    'points, pyramid, majority': () => fx.densityDoc({ n: 200_000, blocks: 2, items: 20_000, density: 'majority' }),
    'ranges, raw': () => fx.densityDoc({ type: 'ranges', items: 6000, len: 300 }),
    'ranges, pyramid': () => fx.densityDoc({ n: 400_000, type: 'ranges', items: 20_000, len: 3000 }),
    pairs: () => fx.pairsDoc({ count: 2000 }),
    matrix: () => fx.wideMatrix(),
  }[c.doc] as () => ReturnType<typeof fx.densityDoc>;
  const small =
    c.entry === 'pairs'
      ? () => fx.pairsDoc({ n: 2000, count: 300 })
      : () =>
          fx.densityDoc({
            n: 2000,
            items: 300,
            blocks: c.mixed ? 2 : 1,
            type: c.entry === 'ranges' ? 'ranges' : 'points',
          });
  const all =
    c.entry === 'points' || c.entry === 'ranges' ? fx.withDensity(c.entry) : ([...KernelRenderers] as LayerRenderer[]);
  const entry = all.find((r) => r.name === c.entry) as LayerRenderer;
  fx.forced.positions = c.mixed ? [null, OrderedSet.ofRange(0, 200_000)] : null;
  let at = 'linear';
  let k = 0;
  try {
    return await rendererConformance(entry, {
      open: async ({ geometry, entry: has }) => {
        at = geometry;
        k = 0;
        // Wrapped reads items (no pyramid there): a small document, so the text layout's lines stay few.
        return (
          await opened(geometry === 'wrapped' ? asText(c.entry === 'matrix' ? make() : small()) : make(), {
            renderers: has ? all : all.filter((r) => r !== entry),
          })
        ).ctx;
      },
      heap: heap && ((frame) => heap(frame, `${at} #${++k}`)),
      items: false,
    });
  } finally {
    fx.forced.positions = null;
  }
}

/** A grid of 64 × 8 fractional points beside one long pair, hit-tested each step. */
export async function pairsHit() {
  vi.resetModules();
  const { KernelRenderers, Hit } = await import('../src/index.js');
  const { opened } = await import('./layout-fixtures.js');
  const { layerOf, entryOf } = await import('./draw-fixtures.js');
  const { pairsDoc } = await import('./rn2-fixtures.js');
  const t = await opened(pairsDoc(), { renderers: [...KernelRenderers] as LayerRenderer[] });
  const { section } = layerOf(t.ctx, 'Contacts');
  t.ctx.viewports.get(section)?.setView(4, 4 * 15_000, 'api');
  await t.ctx.graph.update();
  const e = entryOf(t.ctx, 'Contacts');
  const h = Hit.create();
  const xs = Float64Array.from({ length: 64 }, (_, i) => e.left + 3.3 + i * ((e.right - e.left - 8) / 64));
  const ys = Float64Array.from({ length: 8 }, (_, i) => e.top + 0.7 + i * (e.height / 8));
  return () => {
    for (let i = 0; i < xs.length; i++)
      for (let j = 0; j < ys.length; j++) t.ctx.layout.hitTest(section, xs[i], ys[j], h);
  };
}

/** The overlay's lit pairs: 40 lit runs through the channels' runs and the end index, each step. */
export async function litPairsFrame() {
  vi.resetModules();
  const { KernelRenderers } = await import('../src/index.js');
  const { litPairs, MarkBit } = await import('../src/draw/api.js');
  const { MarkSet } = await import('../src/marks.js');
  const { opened } = await import('./layout-fixtures.js');
  const { layerOf } = await import('./draw-fixtures.js');
  const { pairsDoc } = await import('./rn2-fixtures.js');
  const t = await opened(pairsDoc(), { renderers: [...KernelRenderers] as LayerRenderer[] });
  const { k, section } = layerOf(t.ctx, 'Contacts');
  const st = t.ctx.graph.assembled.value?.view?.prepared[k].state as never;
  const v = t.ctx.layout.axis(section);
  if (!v) throw new Error('no axis');
  const set = new MarkSet();
  set.reset(1);
  set.begin(0);
  for (let i = 0; i < 40; i++) set.add(i * 400, i * 400 + 30);
  set.close(0, 20_000);
  const marks = {
    bitsAt: () => 0,
    any: () => false,
    runs: () => set,
    primary: () => -1,
    item: () => -1,
    itemAnnotation: () => -1,
    element: () => false,
  };
  const out = new Int32Array(5000);
  return () => {
    litPairs(v, st, marks, MarkBit.selection | MarkBit.hover, out);
  };
}
