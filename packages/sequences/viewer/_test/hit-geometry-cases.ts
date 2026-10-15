/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The hover sweeps over sized shapes, for their heap guard and its quick twin (hit-geometry.heap.ts,
// hit-geometry.quick.ts): `sizedShapes`, the pointer over dots, bands, boundaries' glyphs (the end boundary's too, the
// pointer around its edge), pair ends and wide arcs (hit-geometry-fixtures.ts' document); `ownRow`, over an own unit's
// row in `UI.OwnRowHover` "residue" mode (split-units-fixtures.ts' document). Each move a hit, the hover rebuilt and an
// overlay frame on the recording backend. Each call loads the viewer afresh (`vi.resetModules`).

import { vi } from 'vitest';
import type { AllocOptions } from '../../core/_test/model/alloc.js';

/**
 * The sized shapes' quick counts: the sweep crosses six kinds of hit, each its own functions climbing to TurboFan, so
 * even under the quick flags it allocates until about 32,000 moves (40 bytes a move over the first 6,000, 13 at
 * 18,000, under 0.1 from 36,000; 2026-10-14, Node 24, M4): a higher cap, the warm-up running on until it reads quiet.
 */
export const SHAPES: AllocOptions = { maxWarmSteps: 60_000 };

export async function sizedShapes() {
  vi.resetModules();
  const { lit } = await import('./overlay-fixtures.js');
  const { recording, COPY_BACK } = await import('./composer-fixtures.js');
  const { hitGeometryDoc, S } = await import('./hit-geometry-fixtures.js');
  const { Dirty } = await import('../src/draw/index.js');
  const { Inputs, PointerType } = await import('../src/index.js');
  const t = await lit(hitGeometryDoc(), { width: 900, height: 600, dpr: 2 });
  t.ctx.behaviors.get('hoverLabel')?.remove();
  // The whole sequence in view (fit), so both ends' boundaries are on screen.
  t.vp?.setView(1, 0, 'api');
  const be = recording(COPY_BACK);
  t.frame(be, Dirty.all);
  const p = Inputs.pointer();
  p.section = t.ctx.layout.sectionIndex(S);
  p.type = PointerType.move;
  const ys = ['t-dots', 't-bars', 't-ticks', 't-pairs', 't-arcs'].map((r) => t.yOf(S, r));
  const ends = ['t-ticks', 't-tick-points', 't-tick-variants'].map((r) => t.yOf(S, r));
  const x0 = t.xOf(S, 0, 8);
  // The edge after the last position: the end boundary's glyph, the pointer on either side of it.
  const vp = t.ctx.layout.axis(S);
  const xEnd = vp ? vp.x(0, (vp.shownOf(0, 399) + 1) as never) : 0;
  const frame = (i: number) => {
    if (i % 4 === 3) {
      p.x = xEnd + ((i % 9) - 5) * 1.5;
      p.y = ends[i % 3] + ((i % 5) - 2) * 2;
    } else {
      p.x = x0 + (i % 60) * 2.3;
      p.y = ys[i % 5] + ((i % 7) - 3) * 2;
    }
    t.ctx.input.pointer(p);
    t.frame(be, Dirty.overlay);
  };
  const channels = (t.composer as unknown as { channels: { rebuilds: { hover: number } } }).channels;
  return { frame, hovers: () => channels.rebuilds.hover };
}

export async function ownRow() {
  vi.resetModules();
  const { lit } = await import('./overlay-fixtures.js');
  const { recording, COPY_BACK } = await import('./composer-fixtures.js');
  const { splitUnitsDoc } = await import('./split-units-fixtures.js');
  const { Dirty } = await import('../src/draw/index.js');
  const { Inputs, PointerType, SeqViewerConfig } = await import('../src/index.js');
  const t = await lit(splitUnitsDoc({ highlights: true }), { width: 1000, height: 400, dpr: 2 });
  t.ctx.behaviors.get('hoverLabel')?.remove();
  t.ctx.config.set(SeqViewerConfig.UI.OwnRowHover, 'residue');
  const be = recording(COPY_BACK);
  t.frame(be, Dirty.all);
  const p = Inputs.pointer();
  p.section = t.ctx.layout.sectionIndex(t.section);
  p.type = PointerType.move;
  const ys = [0, 1, 2].map((r) => t.yOf(t.section, 't-chains', r));
  const frame = (i: number) => {
    p.x = t.xOf(t.section, 0, 8 + (i % 44));
    p.y = ys[i % 3];
    t.ctx.input.pointer(p);
    t.frame(be, Dirty.overlay);
  };
  return { frame };
}
