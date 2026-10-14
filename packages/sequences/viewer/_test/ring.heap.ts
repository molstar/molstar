/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Rings allocate nothing per frame after warm-up (CG1a; seqspec:viewer#the-render-backend): whole frames of a circular
// section (main and overlay passes through the polar map: every type's batches mapped to sectors, paths, glyphs and
// text) with the pointer moving round the rings (polar hits and the hover, the hover label's strings left out as
// overlay-fixtures.ts does), over 240 positions and over 120,000 (reductions per angular column). CG1b: with a
// selection and a highlight across the origin (marks as sectors), the keyboard cursor's ring (its rect mapped each
// frame), the hover as a radial line or the hovered item's sector and outline, and the ruler's minor ticks. The
// positive control: the same steps allocating an object each (`leaking`), through the same warm-up, settle and
// measurement, are caught.
//
// Each case loads the viewer afresh (`vi.resetModules`), as composer.heap.ts does.

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

/** The guards' measurement: the shared warm-up, settled (`settleFrame`), then `FRAME_N` steps and 4× measured. */
function measure(step: (i: number) => void) {
  // The shared warm-up (2,000 frames and 3 s, restored after the trim): at 1 s, 240 positions failed most runs (about 45
  // to 150 bytes a step; flat with HEAP_DEEP), the pointer still meeting new items, positions and rings.
  const t0 = performance.now();
  let i = 0;
  for (; i < FRAME_WARM || performance.now() - t0 < FRAME_WARM_MS; i++) step(i);
  // Then until a frame stops allocating: over 120,000 positions 3 s were about 4,700 frames, and the climb to TurboFan
  // ran about 6,000 more at 20 to 210 bytes a frame (V8 13.6, the pinned Node, 2026-10-13: 48.8 bytes a step measured).
  settleFrame(() => step(i++), FRAME_SETTLE_MS);
  return heapScaling(FRAME_N, () => step(i++));
}

/** The ring composed over `n` positions with its marks and cursor, and a step: a pointer move and a whole frame. */
async function setup(n: number) {
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

describe('ring frames allocate nothing after warm-up', () => {
  it.each([240, 120_000])(
    'whole frames with polar hovers round the rings, %i positions',
    async (n) => {
      const { step, be } = await setup(n);
      const r = measure(step);
      expect(be.sectorItems.n).toBeGreaterThan(0);
      expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
    },
    120_000,
  );

  it('the same steps allocating an object each are caught (positive control)', async () => {
    const { step } = await setup(240);
    expect(controlCaught(() => measure(leaking(step)))).toBe(true);
  }, 120_000);
});
