/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The fingerprint's heap cases, for its heap guard and its quick twin (fingerprint.heap.ts, fingerprint.quick.ts):
// `fingerprintKit`, the conformance kit on a scene (case 6, 20 poses × 8 kinds, a short track, mixed data, the absent
// property, 20,000 positions, case 6 and the wide one reversed), each frame it draws handed to `heap` with its geometry
// (../kit-cases.ts `KitHeap`); `hitSweep`, a hit a step over every row and stripe of a scene, each row and each kind
// hit somewhere first. Each call loads the viewer afresh (`vi.resetModules`): a frame built in this module would box
// from the second case on (047, RN3).

import { expect, vi } from 'vitest';
import type { KitHeap } from '../kit-cases.js';

export const SCENES = [
  'case 6',
  'many',
  'compact',
  'mixed',
  'absent',
  'wide',
  'case 6 reversed',
  'wide reversed',
] as const;
export type Scene = (typeof SCENES)[number];
export const HIT_SCENES = ['case 6', 'compact', 'mixed'] as const;

async function load() {
  vi.resetModules();
  const fx = await import('./fixtures.js');
  const kit = await import('../../src/testing/index.js');
  const viewer = await import('../../src/index.js');
  const draws = await import('../conformance-fixtures.js');
  const layout = await import('../layout-fixtures.js');
  const cases = await import('../../../core/_test/format/cases.js');
  const src = await import('../../src/interactions/viewer/index.js');
  return { fx, kit, viewer, draws, layout, cases, src };
}

/** The kit on `name`'s scene; `heap` gets every frame it draws. */
export async function fingerprintKit(name: Scene, heap: KitHeap) {
  const { fx, kit, viewer, draws, layout, cases, src } = await load();
  const scene = {
    'case 6': () => ({ doc: cases.case6(), section: 'poses-view.main', view: 'poses-view' }),
    many: () => ({ ...fx.manyScene(), view: 'v' }),
    compact: () => ({ ...fx.manyScene(20, 8, 1), view: 'v' }),
    mixed: () => ({ ...fx.mixedScene(), view: 'v' }),
    absent: () => ({ ...fx.absentScene(), view: 'v' }),
    wide: () => ({ ...fx.wideScene(), view: 'v' }),
    // AX1: every piece reversed.
    'case 6 reversed': () => ({
      doc: draws.directed(cases.case6(), -1),
      section: 'poses-view.main',
      view: 'poses-view',
    }),
    'wide reversed': () => ({ ...fx.wideScene(), doc: draws.directed(fx.wideScene().doc, -1), view: 'v' }),
  }[name]();
  const entry = src.FingerprintRenderer as never;
  let at = '';
  return kit.rendererConformance(entry, {
    open: async ({ geometry, entry: has }) => {
      at = geometry;
      return (
        await layout.opened(geometry === 'wrapped' ? draws.asText(scene.doc) : scene.doc, {
          renderers: [...viewer.KernelRenderers],
          view: scene.view,
          registry: has
            ? [src.InteractionsViewerRegistry, fx.ShortFingerprint]
            : [{ svs: src.InteractionsViewerRegistry.svs }],
        })
      ).ctx;
    },
    heap: (frame) => heap(frame, at),
    items: false,
  });
}

/**
 * A hit a step (`step(i)`) over every row and stripe of `name`'s scene, points every 3 px across and 0.5 px down the
 * track; before it, every row (state) and every kind (part) is checked hit somewhere, and on mixed data the block
 * without the property too.
 */
export async function hitSweep(name: (typeof HIT_SCENES)[number]) {
  const { fx, layout, viewer, cases, src } = await load();
  const scene = {
    'case 6': () => ({
      doc: cases.case6(),
      section: 'poses-view.main',
      track: 'poses-view.fingerprint',
      view: 'poses-view',
    }),
    compact: () => ({ ...fx.manyScene(20, 8, 1), view: 'v' }),
    mixed: () => ({ ...fx.mixedScene(), view: 'v' }),
  }[name]();
  const t = await layout.opened(scene.doc, {
    renderers: [...viewer.KernelRenderers],
    view: scene.view,
    registry: [src.InteractionsViewerRegistry, fx.ShortFingerprint],
  });
  for (const ref of t.ctx.viewports.sections.value) t.ctx.viewports.get(ref)?.setSize(1000, 600, 2);
  await t.ctx.graph.update();
  const l = t.layout(scene.section);
  const e = Array.from(l.entryTrack.subarray(0, l.entryCount)).indexOf(t.trackIndex(scene.track));
  const top = e < l.pinnedCount ? l.entryTop[e] : l.scrolledTop + l.entryTop[e];
  const points: number[] = [];
  for (let x = l.dataLeft; x < l.dataRight; x += 3)
    for (let y = top - 2; y < top + l.entryHeight[e] + 2; y += 0.5) points.push(x, y);
  const hit = viewer.Hit.create();
  // Every row (state) and every kind (part) is hit somewhere.
  const states = new Set<number>();
  const parts = new Set<number>();
  const blocks = new Set<number>();
  for (let k = 0; k < points.length; k += 2) {
    t.ctx.layout.hitTest(scene.section, points[k], points[k + 1], hit);
    states.add(hit.state);
    parts.add(hit.part);
    if (hit.state >= 0) blocks.add(hit.block);
  }
  // Mixed: the block without the property (the third) is hit as bands too.
  if (name === 'mixed') expect([...blocks].sort()).toEqual([0, 1, 2]);
  const s = t.ctx.graph.assembled.value?.view?.prepared.find((p) => p.renderer === 'interactions:fingerprint')
    ?.state as {
    matrix: { states: number };
    kinds: Uint8Array;
  };
  for (let q = 0; q < s.matrix.states; q++) expect(states.has(q)).toBe(true);
  for (const k of s.kinds) if (name !== 'mixed') expect(parts.has(k)).toBe(true);
  expect(states.has(-1)).toBe(true);
  const n = points.length / 2;
  const si = t.ctx.layout.sectionIndex(scene.section);
  // The pointer-rate call (`hitTestAt`, the point in the record), as the dispatcher's: `hitTest`'s fractional x and y
  // are boxed doubles where V8 does not inline it into the caller (8 bytes a hit on the half-pixel rows, 2026-10-14).
  const step = (i: number) => {
    const k = (i * 7) % n;
    hit.x = points[2 * k];
    hit.y = points[2 * k + 1];
    t.ctx.layout.hitTestAt(si, hit);
  };
  return { step };
}
