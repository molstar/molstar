/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The case of arcs whose partners lie beyond the window, for its heap guard and its quick twin
// (arcs-off-window.heap.ts, arcs-off-window.quick.ts): case 7 zoomed into barstar, every contact's partner off-window
// in barnase; a step pans the view (`Arcs.frame` recomputed), draws the frame, hovers an end and draws the overlay.
// Each call loads the viewer afresh (`vi.resetModules`).

import { vi } from 'vitest';

export async function arcsOffWindow() {
  vi.resetModules();
  const { lit } = await import('./overlay-fixtures.js');
  const { recording, COPY_BACK } = await import('./composer-fixtures.js');
  const { layerOf, record } = await import('./draw-fixtures.js');
  const { case7 } = await import('../../core/_test/format/cases.js');
  const t = await lit(case7(), { width: 1000, height: 400, dpr: 2 });
  const be = recording(COPY_BACK);
  t.vp?.zoomTo(1, 34, 86);
  await t.ctx.graph.update();
  t.frame(be, 7);
  // The hover label's strings are one allocation per hovered element by design (IN2's guard): the overlay alone.
  t.ctx.behaviors.get('hoverLabel')?.remove();
  const { k, section } = layerOf(t.ctx, 'Contacts');
  if (section !== t.section) throw new Error('the contacts are not in the first section');
  const s = t.ctx.graph.assembled.value?.view?.prepared[k].state as { frame: Float64Array };
  const glyphs = Array.from(record(t.ctx, 'Contacts').rec.glyphItems.x.subarray(0, 12), (x) => x / 2);
  const into = recording();
  const y = s.frame[0];
  const frame = (i: number) => {
    t.vp?.panBy((i & 1) === 0 ? 5 : -5, 0, 'pointer');
    t.frame(be, 7);
    t.hover(glyphs[i % glyphs.length], y);
    t.overlay(into);
  };
  return { frame, glyphs };
}
