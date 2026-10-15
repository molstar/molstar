/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The conformance kit's heap cases (XR1b, CG1b, RN3), for their heap guard and its quick twin (conformance.heap.ts,
// conformance.quick.ts): `org.example:stars` (a star per item through `mark`, colored by a field provider's value) at
// each of the kit's zooms and on rings; the rows (single rows at every zoom the kit draws, then the dense band) and
// `seqv:map` (its image built once). The kit hands each frame it draws (a layer at a zoom, `tiny` too) to `heap`, with
// the geometry it opened (kit-cases.ts `KitHeap`).

import type { LayerRenderer } from '../src/draw/api.js';
import { KernelRenderers } from '../src/index.js';
import { DefaultPresets } from '../src/presets.js';
import { MapRenderer } from '../src/renderers/map.js';
import { RowsRenderer } from '../src/renderers/rows.js';
import { rendererConformance, StarsFormat, StarsRenderer } from '../src/testing/index.js';
import { miniatureDoc, ringed, starsDoc } from './conformance-fixtures.js';
import { familyWithTree } from './families.js';
import type { KitHeap } from './kit-cases.js';
import { opened } from './layout-fixtures.js';

/** The stars through the kit, linear and on rings; `heap` gets every frame. */
export function starsKit(heap: KitHeap) {
  let at = '';
  const entry = StarsRenderer as LayerRenderer;
  return rendererConformance(entry, {
    open: async ({ geometry, entry: has }) => {
      at = geometry;
      return (
        await opened(geometry === 'circular' ? ringed(starsDoc()) : starsDoc(), {
          renderers: [...KernelRenderers, ...(has ? [entry] : [])],
          registry: [{ svs: StarsFormat }, { presets: DefaultPresets.filter((p) => p.name === 'circular') }],
        })
      ).ctx;
    },
    heap: (frame) => heap(frame, at),
  });
}

/** RN3's entries through the kit: the rows on a family with its tree, and `seqv:map` on the miniature. */
export const RN3 = [
  [RowsRenderer, () => familyWithTree(120, 12, { columns: true, sites: true })],
  [MapRenderer, miniatureDoc],
] as const;

/** One of RN3's entries through the kit; `heap` gets every frame. */
export function rn3Kit([entry, doc]: (typeof RN3)[number], heap: KitHeap) {
  let at = '';
  return rendererConformance(entry as LayerRenderer, {
    open: async ({ geometry, entry: has }) => {
      at = geometry;
      return (
        await opened(doc(), { renderers: has ? [...KernelRenderers] : KernelRenderers.filter((r) => r !== entry) })
      ).ctx;
    },
    heap: (frame) => heap(frame, at),
  });
}
