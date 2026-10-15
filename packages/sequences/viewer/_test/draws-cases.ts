/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// RN1's draw cases through the conformance kit, for their heap guard and its quick twin (draws.heap.ts,
// draws.quick.ts): `letters`, `color`, `values`, `ruler` and `emphasis`, each on the documents that reach its paths
// (every type on two chains, case 14's residues over their codons, 40,000 positions zoomed out, dense values and a row
// field under a continuous palette, the mean), linear and wrapped; the kit hands each frame it draws (a layer at a
// zoom, `tiny` too) to `heap` with its geometry (kit-cases.ts `KitHeap`), or with `fit`, only each geometry's
// first full draw (its first layer at the fit). Each run loads the viewer afresh (`vi.resetModules`).

import { vi } from 'vitest';
import type { AllocOptions } from '../../core/_test/model/alloc.js';
import type { LayerRenderer } from '../src/draw/api.js';
import type { KitHeap } from './kit-cases.js';

/**
 * The quick twin's counts: windows of 250 frames. Under the quick flags emphasis over the means, at the fit, allocated
 * 144 bytes a frame from about 1,000 frames to 3,500, then none (2026-10-14, Node 24, M4: RN1's middle tier); the
 * warm-up runs on until it reads quiet (alloc.ts).
 */
export const DRAW: AllocOptions = { warmSteps: 1_000, steps: 250, windows: 3 };

export const DOCS = ['all types', 'case 14', '40,000 positions', 'means'] as const;
export type DrawDoc = (typeof DOCS)[number];
/** Each renderer on the documents that reach its paths. */
export const CASES: Record<string, readonly DrawDoc[]> = {
  letters: ['all types', 'case 14', '40,000 positions'],
  color: ['all types', 'means'],
  values: ['all types', '40,000 positions'],
  ruler: ['all types', 'case 14', 'means'],
  emphasis: ['all types', 'means'],
};

/** `name`'s draws on `doc` through the kit; `heap` gets every frame, or with `fit` each geometry's first full draw. */
export async function drawsKit(name: string, doc: DrawDoc, heap: KitHeap, fit = false) {
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
  // A layer's frames come zoom by zoom, `tiny` first where the entry has it: the first full draw is the second then.
  const first = entry.tiny && entry.draw ? 1 : 0;
  let at = '';
  let k = 0;
  return rendererConformance(entry, {
    open: async ({ geometry, entry: has }) => {
      at = geometry;
      k = 0;
      return (
        await opened(geometry === 'wrapped' ? asText(make()) : make(), {
          renderers: has ? all : all.filter((r) => r !== entry),
        })
      ).ctx;
    },
    heap: (frame) => {
      if (!fit || k++ === first) heap(frame, at);
    },
    items: false,
  });
}
