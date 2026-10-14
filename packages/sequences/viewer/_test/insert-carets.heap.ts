/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Carets allocate nothing per hover frame after warm-up (seqspec:interaction#channels, seqspec:viewer#performance):
// the pointer sweeping the design's residues, its inserts among them, each move a hover whose projection onto the wild
// type's section finds the insert's next pair and lights a caret; then the wild type's 29–31, which the design deletes
// (the forward leg's insert: a caret at the design's junction), and a point just past the wild type's block (the end
// reach, its boundary marks asked); an overlay frame each, with a selected insert drawn as a caret every frame.
// Positive controls: most measured frames drew a hover's caret beside the selection's, and the same frames allocating an
// object each (`leaking`), through the same warm-up and measurement, are caught.

import { describe, expect, it, vi } from 'vitest';
import {
  controlCaught,
  describeScaling,
  FRAME_N,
  FRAME_SETTLE_MS,
  FRAME_WARM,
  FRAME_WARM_MS,
  heapScaling,
  isFlat,
  leaking,
  settleFrame,
} from './heap.js';

async function setup() {
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
  let measuring = false;
  const frame = (i: number) => {
    const q = points[i % points.length];
    p.section = q[0];
    p.x = q[1];
    p.y = q[2];
    t.ctx.input.pointer(p);
    t.frame(be, Dirty.overlay);
    bd.reset();
    dc.frame(bd, Dirty.overlay);
    if (!measuring) return;
    frames++;
    let n = 0;
    for (let j = 0; j < be.rectItems.n; j++) if (be.rectItems.w[j] === 6) n++;
    for (let j = 0; j < bd.rectItems.n; j++) if (bd.rectItems.w[j] === 6) n++;
    // The selection's caret and a hover's.
    if (n >= 2) carets++;
  };
  return {
    frame,
    measure: () => {
      measuring = true;
    },
    share: () => carets / Math.max(1, frames),
  };
}

/** The guard's measurement: warmed by count and by time, then `heapScaling` over `n` frames and `4n`. */
function run(n: number, frame: (i: number) => void, measure = () => {}) {
  const t0 = performance.now();
  for (let i = 0; i < FRAME_WARM || performance.now() - t0 < FRAME_WARM_MS; i++) frame(i);
  // Settled after the timed warm-up, as every timed warm-up is (lib/seqstar/AGENTS.md, Tests; LC-20).
  let settled = 0;
  settleFrame(() => frame(settled++), FRAME_SETTLE_MS);
  measure();
  let i = 0;
  return heapScaling(n, () => frame(i++));
}

describe('carets over a hover sweep', () => {
  it('allocate nothing; most frames draw a caret; an object allocated per frame is caught', async () => {
    const { frame, measure, share } = await setup();
    const r = run(FRAME_N, frame, measure);
    expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
    // 8 of the 12 points light a caret (the inserts, the deleted stretch).
    expect(share()).toBeGreaterThan(0.5);
    expect(controlCaught(() => run(FRAME_N, leaking(frame)))).toBe(true);
  }, 240_000);
});
