/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The wrapped lanes' cases, for their heap guard and its quick twin (wrapped-lanes.heap.ts, wrapped-lanes.quick.ts): on
// the text preset, a section of letters, underline lanes, a `values` strip in four bands with missing values, a `color`
// strip alone and a `color` layer beside the letters; `draw` scrolls the rows down and back a step each, the lines
// redrawn whole; `hover` moves the pointer over the strips and the underlines through the real dispatcher, each move a
// hit on a strip's layer or an item, the overlay outlining the hovered strip cell. Each call loads the viewer afresh
// (`vi.resetModules`).

import { vi } from 'vitest';
import type { AllocOptions } from '../../core/_test/model/alloc.js';

/**
 * The draw case's quick counts: a frame redraws every line (about 0.7 ms on an M4), so windows of 250 steps; exact
 * accounting reads a byte a step in 250 bytes, and tier-up is done by 1,000 steps under the quick flags.
 */
export const DRAW: AllocOptions = { warmSteps: 1_000, steps: 250, windows: 3 };

export async function wrappedLanes(mode: 'draw' | 'hover') {
  vi.resetModules();
  const { driven } = await import('./input-fixtures.js');
  const { recording, COPY_BACK } = await import('./composer-fixtures.js');
  const { BANDS, BAND_SCORES } = await import('./prepare-fixtures.js');
  const { createSeqViewSpec } = await import('@molstar/seq-core/builder');
  const { DefaultPresets, Hit, HoverKind, Inputs, KernelRenderers, PointerType } = await import('../src/index.js');
  const { runX } = await import('../src/runs.js');
  const { Dirty, SectionComposer } = await import('../src/draw/index.js');
  const S = 'v.main';
  const n = 2000;
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: 'a', value: 'ACDEFGHIKL'.repeat(n / 10) });
  const f = b.data.annotation({
    ref: 'f',
    geometry: 'ranges',
    target: a,
    value: Array.from({ length: 60 }, (_, i) => ({ start: i * 30, stop: i * 30 + 45, label: `f${i}` })),
  });
  const value = Array.from({ length: n }, (_, i) => (i % 23 === 7 ? null : BAND_SCORES[Math.floor(i / 7) % 4]));
  const d = b.data.annotation({ ref: 'd', geometry: 'values', target: a, value: { value } } as never);
  const palette = { theme: 'field_palette', params: { field: 'value', palette: BANDS } };
  const v = b.view({ ref: 'v', label: 'V', layout: 'text', height: 'box' });
  v.block({ sequence: a });
  v.track({ ref: 'v.letters' }).letters();
  v.track({ ref: 'v.features' }).ranges({ data: f, display: 'underline' } as never);
  v.track({ ref: 'v.plddt' }).values({ data: d, domain: [0, 100], color: palette } as never);
  v.track({ ref: 'v.shade' }).color({ data: d, color: palette } as never);
  const behind = v.track({ ref: 'v.behind' });
  behind.letters();
  behind.color({ data: d, color: palette } as never);
  const text = DefaultPresets.find((p) => p.name === 'wrapped');
  const t = await driven(b.getDocument(), {
    size: { [S]: [500, 1200] },
    registry: [{ presets: [{ ...(text as NonNullable<typeof text>), name: 'lab:text' }] }],
    renderers: KernelRenderers as never,
    labels: false,
  });
  t.ctx.ui.setStyle(S, 'lab:text');
  await t.ctx.graph.update();
  // The frames alone: the hover label's strings are one allocation per hovered element by design (IN2's guard).
  t.ctx.behaviors.get('hoverLabel')?.remove();
  const l = t.layout(S);
  const ax = t.ctx.layout.axis(S);
  const run0 = ax?.runs[0];
  if (!ax || !run0 || l.lines.count < 6) throw new Error('no lines');
  const composer = new SectionComposer(t.ctx, S);
  // The draw's frames redraw the tracks whole (no copy to reuse); the hover's, the overlay over a copied back buffer.
  const be = mode === 'draw' ? recording() : recording(COPY_BACK);
  composer.frame(be, Dirty.all);
  // Points over lines 1 to 4: on the values strip and the color strip (each its layer's hit), on the features'
  // underline (an item's), on the letters beside their color.
  const bottom = (track: string, line: number) => {
    const tr = t.trackIndex(track);
    for (let e = 0; e < l.entryCount; e++)
      if (l.entryTrack[e] === tr && l.entryRow[e] === line) return l.scrolledTop + l.entryTop[e] + l.entryHeight[e];
    throw new Error(`no ${track} band on line ${line}`);
  };
  const pts = new Float64Array(64);
  for (let i = 0; i < 32; i++) {
    const line = 1 + (i % 4);
    const k = l.lines.start[line] + ((i * 7) % (l.lines.stop[line] - l.lines.start[line]));
    pts[2 * i] = l.dataLeft + runX(run0, l.lines.start[line], k, ax.gap) + 3.25;
    const where = i % 4;
    pts[2 * i + 1] =
      where === 0
        ? bottom('v.plddt', line) - 1.5
        : where === 1
          ? bottom('v.shade', line) - 1.5
          : where === 2
            ? bottom('v.features', line) - 3.75
            : bottom('v.behind', line) - 10.5;
  }
  // The strips' points hit their layers (the hover path the guard covers), the underline's an item.
  const h = Hit.create();
  let strips = 0;
  let items = 0;
  for (let i = 0; i < 32; i++) {
    if (!t.ctx.layout.hitTest(S, pts[2 * i], pts[2 * i + 1], h)) continue;
    if (h.kind === HoverKind.item) items++;
    else if (i % 4 < 2 && h.kind === HoverKind.seq && h.track === t.trackIndex(i % 4 ? 'v.shade' : 'v.plddt')) strips++;
  }
  if (strips < 12 || items < 2) throw new Error(`the points hit ${strips} strips and ${items} items`);
  const p = Inputs.pointer();
  p.section = t.ctx.layout.sectionIndex(S);
  p.type = PointerType.move;
  p.button = -1;
  p.buttons = 0;
  const vp = t.vp(S);
  let hovered = 0;
  const step = (i: number) => {
    be.reset();
    if (mode === 'draw') {
      // Scrolled down and back a step each (a frame with nothing moved draws nothing): the lines redrawn whole.
      vp.panBy(0, i & 1 ? -11 : 11, 'pointer');
      composer.frame(be, Dirty.tracks);
      return;
    }
    p.x = pts[2 * (i & 31)];
    p.y = pts[2 * (i & 31) + 1];
    t.ctx.input.pointer(p);
    if (t.ctx.hover.live.current.index >= 0) hovered++;
    composer.frame(be, Dirty.overlay);
  };
  const state = () => ({ hovered, rects: be.rectItems.n });
  return { step, state };
}
