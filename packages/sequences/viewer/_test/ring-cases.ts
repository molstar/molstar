/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ring's case (CG1a, CG1b), for its heap guard and its quick twin (ring.heap.ts, ring.quick.ts), on a viewer loaded
// afresh (`vi.resetModules`): a circular section of `n` positions (a ruler, letters, labelled ranges, values) with a
// selection and a highlight across the origin and the keyboard cursor's ring, and the step, the pointer moved round the
// rings (a polar hit and the hover) and a whole frame (main and overlay passes) on the recording backend.

import { vi } from 'vitest';
import type { AllocOptions } from '../../core/_test/model/alloc.js';

/**
 * The quick twin's counts: a frame takes 1 to 3 ms (M4), so windows of 500 steps; under the quick flags 240 positions'
 * frames still climb for several thousand steps (hundreds of bytes a step falling to a few, 2026-10-14, Node 24), so
 * the warm-up is longer than a draw's. 120,000 positions (3 ms a frame) are quiet by 3,500 steps: `RING_LARGE`.
 */
export const RING: AllocOptions = {
  warmSteps: 2_000,
  steps: 500,
  windows: 3,
  maxWarmSteps: 16_000,
};
export const RING_LARGE: AllocOptions = { warmSteps: 2_000, steps: 250, windows: 3 };

/** The ring composed over `n` positions with its marks and cursor, and a step: a pointer move and a whole frame. */
export async function ring(n: number) {
  vi.resetModules();
  const { composed, recording } = await import('./composer-fixtures.js');
  const { createSeqViewSpec } = await import('@molstar/seq-core/builder');
  const { DefaultPresets } = await import('../src/presets.js');
  const { Inputs, PointerType } = await import('../src/index.js');
  const b = createSeqViewSpec();
  const s = b.data.sequence({ ref: 's', value: 'ACGTACGTTA'.repeat(n / 10), alphabet: 'dna', cyclic: true });
  const value = { start: [Math.floor(0.8 * n), Math.floor(n / 6)], stop: [Math.floor(n / 8), Math.floor(n / 2)] };
  const f = b.data.annotation({ ref: 'f', geometry: 'ranges', target: s, value });
  const vals = b.data.annotation({ ref: 'vals', geometry: 'values', target: s, value: { value: [0.5, 1, 0.2] } });
  const v = b.view({ ref: 'v', label: 'V', style: 'circular' }).block({ sequence: s });
  v.track({ label: '' }).ruler();
  v.track({ label: 'S' }).letters();
  v.track({ label: 'F' }).ranges({ data: f, labels: true });
  v.track({ label: 'V' }).values({ data: vals });
  const registry = [{ presets: DefaultPresets.filter((p) => p.name === 'circular') }];
  const t = await composed(b.getDocument(), { width: 600, height: 600, dpr: 2, registry });
  const l = t.ctx.layout.get(t.section)?.current;
  if (!l) throw new Error('no layout');
  const be = recording();
  // The frames alone: the hover label's strings are one allocation per hovered element by design (IN2's guard).
  t.ctx.behaviors.get('hoverLabel')?.remove();
  // Marks across the origin, and the keyboard cursor on the letters (set once: a cursor move publishes an object).
  t.ctx.selection.set({ systems: ['seqstar:s/unit'], runs: Int32Array.of(0, n - n / 12, n, 0, 0, n / 24) }, 'api');
  t.ctx.highlight.set('host', { systems: ['seqstar:s/unit'], runs: Int32Array.of(0, n / 3, n / 2) });
  const m = t.ctx.graph.assembled.value?.view;
  const track = m ? m.set.strings.get(m.resolved.tracks.ref[l.entryTrack[1]]) : '';
  const row = l.entryRow[1];
  t.ctx.ui.setKeyCursor({ section: t.section, track, header: false, block: 0, index: 7, row: -1, displayRow: row });
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
