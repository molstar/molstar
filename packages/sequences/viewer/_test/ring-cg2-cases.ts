/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// CG2's ring case, for its heap guard and its quick twin (ring-cg2.heap.ts, ring-cg2.quick.ts), on a viewer loaded
// afresh (`vi.resetModules`): a circular section of `n` positions with items' labels outside the ring (leaders and texts
// from the kept placement), arrowheads by a direction field and pairs as chords, and the step, the pointer moved round
// the rings (a polar hit and the hover) and a whole frame on the recording backend.

import { vi } from 'vitest';
import type { AllocOptions } from '../../core/_test/model/alloc.js';

/**
 * The quick twin's counts: a frame takes about 1.4 ms (M4), so windows of 500 steps; under the quick flags the frames
 * still climb for several thousand steps (hundreds of bytes a step falling to a few, 2026-10-14, Node 24), so the
 * warm-up is longer than a draw's.
 */
export const CG2: AllocOptions = {
  warmSteps: 2_000,
  steps: 500,
  windows: 3,
  maxWarmSteps: 16_000,
};

/** The ring composed over `n` positions, and a step: a pointer move round the rings and a whole frame. */
export async function ringCg2(n: number) {
  vi.resetModules();
  const { composed, recording } = await import('./composer-fixtures.js');
  const { createSeqViewSpec } = await import('@molstar/seq-core/builder');
  const { DefaultPresets } = await import('../src/presets.js');
  const { Inputs, PointerType } = await import('../src/index.js');
  const b = createSeqViewSpec();
  const s = b.data.sequence({ ref: 's', value: 'ACGTACGTTA'.repeat(n / 10), alphabet: 'dna', cyclic: true });
  const at = (x: number) => Math.floor(x * n);
  const value = {
    start: [at(0.8), at(1 / 6), at(0.55)],
    stop: [at(1 / 8), at(1 / 2), at(0.7)],
    label: ['across', 'forward', 'reverse'],
    strand: ['+', '+', '-'],
  };
  const f = b.data.annotation({ ref: 'f', geometry: 'ranges', target: s, value } as never);
  const cuts = { boundary: [at(0.3), at(0.31), at(0.32), at(0.9)], label: ['A', 'B', 'C', 'D'] };
  const c = b.data.annotation({ ref: 'c', geometry: 'ranges', target: s, value: cuts } as never);
  const pv = { start: [at(0.1), at(0.4)], b_target: ['s', 's'], b_start: [at(0.6), at(0.45)] };
  const pp = b.data.annotation({ ref: 'pp', geometry: 'pairs', target: s, value: pv } as never);
  const v = b.view({ ref: 'v', label: 'V', style: 'circular' }).block({ sequence: s });
  v.track({ label: '' }).ruler();
  v.track({ label: 'F', height: 2 }).ranges({ data: f, labels: true, direction_field: 'strand' } as never);
  v.track({ label: 'C' }).ranges({ data: c, labels: true } as never);
  v.track({ label: 'P' }).pairs({ data: pp, display: 'both' } as never);
  const registry = [{ presets: DefaultPresets.filter((p) => p.name === 'circular') }];
  const t = await composed(b.getDocument(), { width: 600, height: 600, dpr: 2, registry });
  const l = t.ctx.layout.get(t.section)?.current;
  if (!l) throw new Error('no layout');
  const be = recording();
  // The frames alone: the hover label's strings are one allocation per hovered element by design (IN2's guard).
  t.ctx.behaviors.get('hoverLabel')?.remove();
  const p = Inputs.pointer();
  p.section = t.ctx.layout.sectionIndex(t.section);
  p.type = PointerType.move;
  const step = (i: number) => {
    // Round the rings, in and out across the entries: a polar hit and the hover each move.
    const a = (i * 0.37) % (2 * Math.PI);
    const r = l.ringR * (0.5 + 0.5 * ((i * 0.13) % 1));
    p.x = l.ringCx + r * Math.cos(a);
    p.y = l.ringCy + r * Math.sin(a);
    t.ctx.input.pointer(p);
    t.composer.redraw();
    t.frame(be, 7);
  };
  return { step, be };
}
