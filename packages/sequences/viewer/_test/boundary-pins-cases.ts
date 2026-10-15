/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The pinned boundaries' hover sweep, for its heap guard and its quick twin (boundary-pins.heap.ts,
// boundary-pins.quick.ts): the pointer moving across a design row's pinned boundaries and the letters beside them
// (boundary-pins-fixtures.ts' document), each move a hit and an overlay frame outlining the hovered pin; `boundary`
// counts the frames that hovered a boundary. Each call loads the viewer afresh (`vi.resetModules`).

import { vi } from 'vitest';

export async function boundaryPins() {
  vi.resetModules();
  const { lit } = await import('./overlay-fixtures.js');
  const { recording, COPY_BACK } = await import('./composer-fixtures.js');
  const { pinsDoc, S } = await import('./boundary-pins-fixtures.js');
  const { Dirty } = await import('../src/draw/index.js');
  const { HoverKind, Inputs, PointerType } = await import('../src/index.js');
  const ITEM = HoverKind.item;
  const t = await lit(pinsDoc(), { width: 900, height: 400, dpr: 2 });
  t.ctx.behaviors.get('hoverLabel')?.remove();
  t.vp?.setView(16, 0, 'api');
  const be = recording(COPY_BACK);
  t.frame(be, Dirty.all);
  const p = Inputs.pointer();
  p.section = t.ctx.layout.sectionIndex(S);
  p.type = PointerType.move;
  const a = t.ctx.layout.axis(S);
  const at = a ? a.x(0, a.shownOf(0, 20)) : 0;
  // On the middle boundary and a few px off it (its pin), then the residues beside it; on both rows.
  const xs = [at, at - 2, at + 3, at - 8, at + 8, at + 1];
  const ys = [t.yOf(S, 't-of'), t.yOf(S, 't-own')];
  let boundary = 0;
  const frame = (i: number) => {
    p.x = xs[i % xs.length];
    p.y = ys[(i / xs.length) & 1];
    t.ctx.input.pointer(p);
    t.frame(be, Dirty.overlay);
    if (t.ctx.hover.live.current.kind === ITEM) boundary++;
  };
  return { frame, boundary: () => boundary };
}
