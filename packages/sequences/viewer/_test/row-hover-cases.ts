/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The data hover's row band (the `hoverRow` behavior; seqspec:interaction#hover, seqspec:viewer#performance) as a step,
// for its heap guard and its quick twin (row-hover.heap.ts, row-hover.quick.ts): the hover document (hover-fixtures.ts)
// lit at DPR 2, the pointer sweeping down its tracks (letters, ranges, points, values, a matrix's bands through its
// renderer's `rowHover`), the element changing on most moves, an overlay frame per move on the recording backend. Each
// call loads the viewer afresh (`vi.resetModules`).

import { vi } from 'vitest';

/**
 * The hover document lit at DPR 2, and a frame: the pointer on a track's row, an overlay frame, the band counted
 * (`hovers`: the element changed; `drawn`: the frame drew the band).
 */
export async function rowHover() {
  vi.resetModules();
  const { Dirty, SkinToken } = await import('../src/draw/index.js');
  const { Inputs, PointerType } = await import('../src/index.js');
  const { COPY_BACK, recording } = await import('./composer-fixtures.js');
  const { hoverDoc } = await import('./hover-fixtures.js');
  const { lit } = await import('./overlay-fixtures.js');
  const t = await lit(hoverDoc(100), { width: 1000, height: 600, dpr: 2 });
  // The overlay alone: the hover label's strings are one allocation per hovered element by design (IN2's guard).
  t.ctx.behaviors.get('hoverLabel')?.remove();
  t.vp?.setView(12, 0, 'api');
  await t.ctx.graph.update();
  const be = recording(COPY_BACK);
  t.frame(be, Dirty.all);
  const tracks = ['t.letters', 't.ranges', 't.points', 't.values', 't.matrix'];
  const ys: number[] = [];
  for (const ref of tracks) ys.push(t.yOf(t.section, ref));
  // A matrix's three bands, each its own row band.
  const m = ys[ys.length - 1];
  ys.push(m - 12, m + 12);
  // Eight positions against seven rows (coprime): the sweep meets every pair.
  const xs = [3, 4, 5, 6, 7, 8, 9, 10].map((i) => t.xOf(t.section, 0, i));
  const p = Inputs.pointer();
  p.section = t.ctx.layout.sectionIndex(t.section);
  p.type = PointerType.move;
  const HOVER = SkinToken.hover;
  let drawn = 0;
  const frame = (i: number) => {
    p.x = xs[i % xs.length];
    p.y = ys[i % ys.length];
    t.ctx.input.pointer(p);
    t.frame(be, Dirty.overlay);
    const r = be.rectItems;
    let band = 0;
    for (let j = 0; j < r.n && band === 0; j++) if (r.x[j] === 0 && r.color[j] === HOVER) band++;
    drawn += band;
  };
  const channels = (t.composer as unknown as { channels: { rebuilds: { hover: number } } }).channels;
  return { frame, hovers: () => channels.rebuilds.hover, drawn: () => drawn };
}
