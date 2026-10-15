/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The carets' hover sweep, for its heap guard and its quick twin (insert-carets.heap.ts, insert-carets.quick.ts): the
// pointer over the design's inserts and the residues beside them, the wild type's deleted 29–31 and a point past its
// block (insert-carets-fixtures.ts' document), an overlay frame of both sections a move, with the wild type's 30–32
// selected (a caret every frame). `share` is the fraction of frames that drew two carets or more (the selection's and a
// hover's). Each call loads the viewer afresh (`vi.resetModules`).

import { vi } from 'vitest';

export async function insertCarets() {
  vi.resetModules();
  const { lit } = await import('./overlay-fixtures.js');
  const { recording, COPY_BACK } = await import('./composer-fixtures.js');
  const { caretsDoc } = await import('./insert-carets-fixtures.js');
  const { Dirty, SectionComposer } = await import('../src/draw/index.js');
  const { Inputs, PointerType } = await import('../src/index.js');
  const { OrderedSet } = await import('@molstar/seq-core/model');
  const t = await lit(caretsDoc(), { width: 2400, height: 400, dpr: 2 });
  t.ctx.behaviors.get('hoverLabel')?.remove();
  // WT 30–32, deleted in the design, selected: a caret at the design's junction (its 33), drawn every frame.
  const set = t.set;
  const unit = set.units.key.findIndex((k) => set.strings.get(k) === 'wt/unit');
  t.ctx.selection.set(
    { kind: 'seq-loci', set, elements: [{ unit, indices: OrderedSet.ofRange(29, 32) }] } as never,
    'api',
  );
  const be = recording(COPY_BACK);
  t.frame(be, Dirty.all);
  // The design's section too (the fixture's composer draws the first, the wild type's).
  const dc = new SectionComposer(t.ctx, 'v.d');
  const bd = recording(COPY_BACK);
  dc.frame(bd, Dirty.all);
  const p = Inputs.pointer();
  p.type = PointerType.move;
  const D = t.ctx.layout.sectionIndex('v.d');
  const W = t.ctx.layout.sectionIndex('v.wt');
  const wt = t.ctx.layout.axis('v.wt');
  const end = wt ? wt.x(0, (wt.shownOf(0, 39) + 1) as never) : 0;
  if (end + 5 >= t.layout('v.wt').dataRight) throw new Error("no room past the wild type's block");
  // The design's inserted residues (20–23), its residues beside them and past the junction; the wild type's 29–31,
  // and past its end.
  const points: [number, number, number][] = [
    ...[20, 21, 22, 23, 18, 33, 22, 34].map((i): [number, number, number] => [
      D,
      t.xOf('v.d', 0, i),
      t.yOf('v.d', 't-d'),
    ]),
    ...[29, 30, 31].map((i): [number, number, number] => [W, t.xOf('v.wt', 0, i), t.yOf('v.wt', 't-wt')]),
    [W, end + 5, t.yOf('v.wt', 't-wt')],
  ];
  let carets = 0;
  let frames = 0;
  const frame = (i: number) => {
    const q = points[i % points.length];
    p.section = q[0];
    p.x = q[1];
    p.y = q[2];
    t.ctx.input.pointer(p);
    t.frame(be, Dirty.overlay);
    bd.reset();
    dc.frame(bd, Dirty.overlay);
    frames++;
    let n = 0;
    for (let j = 0; j < be.rectItems.n; j++) if (be.rectItems.w[j] === 6) n++;
    for (let j = 0; j < bd.rectItems.n; j++) if (bd.rectItems.w[j] === 6) n++;
    // The selection's caret and a hover's.
    if (n >= 2) carets++;
  };
  return { frame, share: () => carets / Math.max(1, frames) };
}
