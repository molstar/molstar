/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Values colored by their values through the conformance kit (seqspec:format#layer-types;
// seqspec:viewer#the-render-backend), for the heap guard and its quick twin (values-lines.heap.ts,
// values-lines.quick.ts): `linesDoc`'s line in four bands (a discrete palette, its color changing halfway between two
// positions), its area under a continuous palette and its plain line over 20,000 positions, with the frame handed to the
// kit's `heap` hook at each of its zooms (one below a device px per position: columns' means and bands). The viewer is
// loaded afresh per call (`vi.resetModules`).

import { vi } from 'vitest';
import type { LayerRenderer } from '../src/index.js';

/**
 * The kit over the values entry, `heap` getting each linear frame the kit checks and where it is (its place in the
 * kit's order: `linear #3`). Returns the heap checks of the linear geometry.
 */
export async function valuesLines(heap?: (frame: () => void, where: string) => void | Promise<void>) {
  vi.resetModules();
  const { KernelRenderers } = await import('../src/index.js');
  const { rendererConformance } = await import('../src/testing/index.js');
  const { opened } = await import('./layout-fixtures.js');
  const { linesDoc } = await import('./prepare-fixtures.js');
  const all = KernelRenderers as LayerRenderer[];
  const entry = all.find((r) => r.name === 'values') as LayerRenderer;
  let at = 'linear';
  let k = 0;
  const checks = await rendererConformance(entry, {
    open: async ({ geometry, entry: has }) => {
      at = geometry;
      k = 0;
      return (await opened(linesDoc(20_000), { renderers: has ? all : all.filter((r) => r !== entry) })).ctx;
    },
    // The linear frames only: the other geometries' results are not judged (nor measured by the twin).
    heap: heap && ((frame) => (at === 'linear' ? heap(frame, `linear #${++k}`) : undefined)),
    items: false,
  });
  return checks.filter((c) => c.check === 'heap' && c.geometry === 'linear');
}
