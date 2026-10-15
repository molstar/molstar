/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// AX1's reversed and ranged pieces' cases, for their heap guard and its quick twin (pieces.heap.ts, pieces.quick.ts):
// the conformance kit run on every built-in renderer's kind over documents that reach their paths with every section
// `direction: -1` and every block of a sequence a `range`, with the kit's `heap` hook given by the caller (a deep
// guard's measurement, or a quick case); and the overlay's hover frames on reversed pieces. Each call loads the viewer
// afresh (`vi.resetModules`).

import { vi } from 'vitest';
import type { AllocOptions } from '../../core/_test/model/alloc.js';
import type { LayerRenderer } from '../src/index.js';

/**
 * The quick twin's counts for each frame the kit draws (there are a dozen a case, each measured with its control): a
 * frame of a few tens of microseconds, windows of 500 steps. A warm-up allowed to end at 1,500 steps stopped on a
 * plateau of 144 bytes a step in one frame (2026-10-15), so it runs at least 4,000.
 */
export const KIT: AllocOptions = { warmSteps: 1_000, steps: 250, windows: 3 };

// A reversed piece goes through the same code for every renderer of a kind (the axis's run conversions), so one
// renderer is kept for each of the three: letters, values and color (values), ranges (ranges, points, variants,
// secondary); the renderers with a code path of their own stay. The forward cases of every renderer are in
// draws.heap.ts and rn2-draws.heap.ts. Dropped on 2026-10-12 (the tier's run time): the `DEEP_CASES`.
export const CASES: [string, string][] = [
  ['letters', 'all types'],
  ['values', '40,000 positions'],
  ['pairs', 'all types'],
  ['ranges', 'all types'],
  ['rows', 'family'],
  ['logo', 'family'],
  ['matrix', 'matrix'],
  ['seqv:map', 'miniature'],
];

/** The cases the trim dropped: `HEAP_DEEP` runs them (the quick twin does not). */
export const DEEP_CASES: [string, string][] = [
  ['letters', '40,000 positions'],
  ['color', 'means'],
  ['ruler', 'all types'],
  ['emphasis', 'all types'],
  ['ranges', 'density'],
  ['points', 'density'],
  ['variants', 'all types'],
  ['secondary', 'all types'],
];

export async function run(name: string, doc: string, heap: (frame: () => void) => void | Promise<void>) {
  vi.resetModules();
  const { KernelRenderers } = await import('../src/index.js');
  const { rendererConformance } = await import('../src/testing/index.js');
  const { asText, directed, matrixDoc, miniatureDoc } = await import('./conformance-fixtures.js');
  const { meanDoc } = await import('./draw-fixtures.js');
  const { opened } = await import('./layout-fixtures.js');
  const { allTypes, scaleDoc } = await import('./prepare-fixtures.js');
  const { densityDoc } = await import('./rn2-fixtures.js');
  const { case3 } = await import('../../core/_test/format/cases.js');
  const make = {
    'all types': () => allTypes(),
    '40,000 positions': () => scaleDoc({ length: 20_000 }),
    means: () => meanDoc(10_000),
    density: () => densityDoc({ type: name === 'ranges' ? 'ranges' : 'points', items: 6000, len: 300 }),
    family: case3,
    matrix: matrixDoc,
    miniature: () => miniatureDoc(),
  }[doc] as () => ReturnType<typeof allTypes>;
  const all = KernelRenderers as LayerRenderer[];
  const entry = all.find((r) => r.name === name) as LayerRenderer;
  return rendererConformance(entry, {
    open: async ({ geometry, entry: has }) => {
      const d = directed(make(), -1, true);
      return (
        await opened(geometry === 'wrapped' ? asText(d) : d, {
          renderers: has ? [...all] : all.filter((r) => r !== entry),
        })
      ).ctx;
    },
    heap,
    items: false,
  });
}

/** The overlay's hover frames on reversed pieces (overlay-fixtures.ts `hoverFrames`), over the scale or the rows. */
export async function reversedHover(kind: 'scale' | 'rows') {
  vi.resetModules();
  const { hoverFrames } = await import('./overlay-fixtures.js');
  return hoverFrames(kind, -1);
}
