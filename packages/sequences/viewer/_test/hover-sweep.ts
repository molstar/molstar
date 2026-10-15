/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Hover sweeps for the hover emphasis' heap guards (HV1; seqspec:interaction#hover), imported afresh per case: the
// pointer moving over one track of the hover document at DPR 2, an overlay frame per move on the recording backend,
// alternating between the track's items (or a matrix's cells) so the hovered element changes on most moves; and the
// positive control, the frames whose overlay drew an emphasis.

import { Dirty } from '../src/draw/index.js';
import { Inputs, PointerType } from '../src/index.js';
import { COPY_BACK, recording } from './composer-fixtures.js';
import { hoverDoc } from './hover-fixtures.js';
import { lit } from './overlay-fixtures.js';

/** Per track: the zoom (CSS px a position; 0: the fit), the positions swept and the CSS px below the entry's top. */
const SWEEPS = {
  't.ranges': { cell: 0, at: [12, 42, 18, 65, 35, 75], dy: [6, 6, 15, 6, 15, 6] },
  't.points': { cell: 0, at: [12, 50, 70], dy: [11] },
  't.variants': { cell: 0, at: [12, 50], dy: [8] },
  't.matrix': { cell: 12, at: [3, 4, 5, 6, 7, 8, 9], dy: [4, 20, 36] },
} as const;

export async function hoverSweep(track: keyof typeof SWEEPS) {
  const s = SWEEPS[track];
  const t = await lit(hoverDoc(100), { width: 1000, height: 600, dpr: 2 });
  // The overlay alone: the hover label's strings are one allocation per hovered element by design (IN2's guard).
  t.ctx.behaviors.get('hoverLabel')?.remove();
  if (s.cell > 0) t.vp?.setView(s.cell, 0, 'api');
  await t.ctx.graph.update();
  const be = recording(COPY_BACK);
  t.frame(be, Dirty.all);
  const l = t.layout(t.section);
  const k = t.trackIndex(track);
  let e = 0;
  while (l.entryTrack[e] !== k) e++;
  const top = (e < l.pinnedCount ? 0 : l.scrolledTop) + l.entryTop[e];
  const xs = s.at.map((i) => t.xOf(t.section, 0, i));
  const p = Inputs.pointer();
  p.section = t.ctx.layout.sectionIndex(t.section);
  p.type = PointerType.move;
  let drawn = 0;
  const frame = (i: number) => {
    p.x = xs[i % xs.length];
    p.y = top + s.dy[i % s.dy.length];
    t.ctx.input.pointer(p);
    t.frame(be, Dirty.overlay);
    let lit = be.glyphItems.n + be.pathItems.n;
    for (let j = 0; j < be.rectItems.n && lit === 0; j++) if (be.rectItems.stroke[j] > 0) lit++;
    if (lit > 0) drawn++;
  };
  const channels = (t.composer as unknown as { channels: { rebuilds: { hover: number } } }).channels;
  return { frame, hovers: () => channels.rebuilds.hover, drawn: () => drawn };
}
