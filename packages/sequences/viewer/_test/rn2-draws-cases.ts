/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The cases of RN2's feature, pair-overlay and logo draws, for their heap guard and its quick twin (rn2-draws.heap.ts,
// rn2-draws.quick.ts): an entry (`ranges`, `points`, `variants`, `secondary`, `logo`) through the conformance kit on a
// document that reaches one of its paths (`kitRun`, the frame handed to the kit's `heap` hook at each of its zooms and
// geometries), and the `pairs` overlay in both of its modes over many lit ends (`pairsOverlay`, one frame). Each call
// loads the viewer afresh (`vi.resetModules`; draws.heap.ts says why: Maglev and the module runner).

import { vi } from 'vitest';
import type { LayerRenderer } from '../src/index.js';

export type Case = { entry: 'ranges' | 'points' | 'variants' | 'secondary' | 'logo'; doc: string };
export const CASES: Case[] = [
  { entry: 'ranges', doc: 'features' },
  { entry: 'ranges', doc: 'boundaries over letters' },
  { entry: 'ranges', doc: 'transcripts' },
  { entry: 'ranges', doc: 'density, raw' },
  { entry: 'ranges', doc: 'density, pyramid' },
  { entry: 'points', doc: 'features' },
  { entry: 'points', doc: 'density, raw' },
  { entry: 'points', doc: 'density, pyramid' },
  { entry: 'variants', doc: 'features' },
  { entry: 'variants', doc: 'scale' },
  { entry: 'secondary', doc: 'features' },
  { entry: 'logo', doc: 'logo' },
  { entry: 'logo', doc: 'wide logo' },
];

/** The pairs overlay's modes: how the arcs are drawn, and whether an item is hovered in the `hover` mode. */
export const PAIRS_MODES = ['always', 'hover', 'hover, an item hovered'] as const;

/**
 * The case through the conformance kit, `heap` getting each frame the kit checks and where it is (the geometry and the
 * frame's place in the kit's order there: `linear #3`). Omitted: the kit's other checks only.
 */
export async function kitRun(c: Case, heap?: (frame: () => void, where: string) => void | Promise<void>) {
  vi.resetModules();
  const { KernelRenderers } = await import('../src/index.js');
  const { rendererConformance } = await import('../src/testing/index.js');
  const { asText } = await import('./conformance-fixtures.js');
  const { opened } = await import('./layout-fixtures.js');
  const { scaleDoc, transcripts } = await import('./prepare-fixtures.js');
  const fx = await import('./rn2-fixtures.js');
  const mine = await import('./rn2b-fixtures.js');
  const { pinsDoc } = await import('./boundary-pins-fixtures.js');
  const make = {
    features: () => mine.featuresDoc({ n: 300 }),
    'boundaries over letters': () => pinsDoc(),
    transcripts: () => transcripts(),
    'density, raw': () => fx.densityDoc({ type: c.entry === 'ranges' ? 'ranges' : 'points', items: 6000, len: 300 }),
    'density, pyramid': () =>
      fx.densityDoc({ n: 400_000, type: c.entry === 'ranges' ? 'ranges' : 'points', items: 20_000, len: 3000 }),
    scale: () => scaleDoc({ length: 20_000 }),
    logo: () => mine.logoDoc(),
    'wide logo': () => mine.logoDoc(12, 8000),
  }[c.doc] as () => ReturnType<typeof mine.featuresDoc>;
  const all = [...KernelRenderers] as LayerRenderer[];
  const entry = all.find((r) => r.name === c.entry) as LayerRenderer;
  let at = 'linear';
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
    heap: heap && ((frame) => heap(frame, `${at} #${++k}`)),
    items: false,
  });
}

/** The pairs overlay of `mode`: a frame drawing it on the recording backend, and what it drew (not vacuous). */
export async function pairsOverlay(mode: (typeof PAIRS_MODES)[number]) {
  const arcs = mode.startsWith('hover') ? 'hover' : 'always';
  vi.resetModules();
  const { KernelRenderers } = await import('../src/index.js');
  const { opened } = await import('./layout-fixtures.js');
  const { layerOf } = await import('./draw-fixtures.js');
  const { RecordingBackend } = await import('../src/testing/index.js');
  const { pairsDoc } = await import('./rn2-fixtures.js');
  const { lit, MarkBit, overlayFrame } = await import('./rn2b-fixtures.js');
  const t = await opened(pairsDoc({ arcs }), { renderers: [...KernelRenderers] as LayerRenderer[] });
  const { section } = layerOf(t.ctx, 'Contacts');
  t.ctx.viewports.get(section)?.setView(0.05, 0, 'api');
  await t.ctx.graph.update();
  const runs: [number, number][] = [];
  for (let i = 0; i < 40; i++) runs.push([i * 400, i * 400 + 30]);
  const marks = lit(MarkBit.hover, runs, 1, 20_000) as ReturnType<typeof lit> & {
    hoverLayer: number;
    hoverAnnotation: number;
    hoverItem: number;
  };
  if (mode !== arcs) {
    // The hovered pair alone for the hover (`marks.item` of its annotation), through its index by item.
    const { k } = layerOf(t.ctx, 'Contacts');
    const s = t.ctx.graph.assembled.value?.view?.prepared[k].state as { annotation: Int32Array };
    marks.hoverLayer = k;
    marks.hoverAnnotation = s.annotation[0];
    marks.hoverItem = 7;
  }
  const frame = overlayFrame(t.ctx, 'Contacts', marks);
  const rec = new RecordingBackend();
  frame.draw(rec);
  return { frame: () => frame.draw(rec), paths: () => rec.pathItems.n };
}
