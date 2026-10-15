/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The cells of several kinds' case, for its heap guard and its quick twin (pairs-kinds.heap.ts, pairs-kinds.quick.ts):
// pairs between hubs every 12 residues, one with 3 kinds, the next with 11 (more than fit: a "+n" marker), and a frame
// of a pan, a hover on one of the first hubs' glyphs and the overlay; a fresh viewer (`vi.resetModules`).

import { vi } from 'vitest';

async function setup() {
  vi.resetModules();
  const { lit } = await import('./overlay-fixtures.js');
  const { recording, COPY_BACK } = await import('./composer-fixtures.js');
  const { layerOf } = await import('./draw-fixtures.js');
  const { createSeqViewSpec } = await import('@molstar/seq-core/builder');
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: 'a', value: 'MKVLATGHEKWRSPQIDNCY'.repeat(100) });
  // Hubs every 12 residues: one with 3 kinds, the next with 11 (more than fit: a marker), partners on the next cells.
  const start: number[] = [];
  const bStart: number[] = [];
  const category: string[] = [];
  for (let hub = 0; hub < 150; hub++) {
    const kinds = hub % 2 === 0 ? 3 : 11;
    for (let k = 0; k < kinds; k++) {
      start.push(hub * 12);
      bStart.push(hub * 12 + 1 + (k % 8));
      category.push(`k${k}`);
    }
  }
  const pairs = b.data.annotation({
    ref: 'pp',
    geometry: 'pairs',
    target: a,
    value: { start, b_target: start.map(() => 'a'), b_start: bStart, category },
  } as never);
  const v = b.view({ ref: 'v', label: 'V', height: 'box' }).block({ sequence: a });
  v.track({ label: 'Letters' }).letters();
  v.track({ label: 'Contacts' }).pairs({
    data: pairs,
    arcs: 'hover',
    glyph: { theme: 'field_map', params: { map: { k0: 'diamond', k1: 'circle', k2: 'square' } } },
  } as never);
  const t = await lit(b.getDocument(), { width: 1000, height: 300, dpr: 2 });
  const be = recording(COPY_BACK);
  t.vp?.setView(4, 0, 'api');
  t.frame(be, 7);
  // The hover label's strings are one allocation per hovered element by design (IN2's guard): the overlay alone.
  t.ctx.behaviors.get('hoverLabel')?.remove();
  const { k } = layerOf(t.ctx, 'Contacts');
  const s = t.ctx.graph.assembled.value?.view?.prepared[k].state as { frame: Float64Array; cellKinds: Int32Array };
  return { t, be, s };
}

/** The case: a frame (it gets a running index) of a pan, a hover on a glyph and its overlay; and whether a cell holds kinds. */
export async function kindsFrame() {
  const { t, be, s } = await setup();
  const into = (await import('./composer-fixtures.js')).recording();
  // The glyphs of the first hubs, from a draw of the layer.
  const { record, layerOf } = await import('./draw-fixtures.js');
  const glyphs = Array.from(record(t.ctx, 'Contacts').rec.glyphItems.x.subarray(0, 24), (x) => x / 2);
  const y = s.frame[0];
  const frame = (i: number) => {
    t.vp?.panBy((i & 1) === 0 ? 5 : -5, 0, 'pointer');
    t.frame(be, 7);
    t.hover(glyphs[i % glyphs.length], y);
    t.overlay(into);
  };
  return { frame, multi: s.cellKinds.some((n) => n > 1), inSection: layerOf(t.ctx, 'Contacts').section === t.section };
}
