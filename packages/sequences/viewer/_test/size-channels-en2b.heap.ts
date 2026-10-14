/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// EN2b's channels in frames allocate nothing after warm-up: pairs at a height theme's apex with their labels, variants
// and secondary structure at a thickness, through pans.
// The positive control runs the same frames allocating an object each (`leaking`) through the same warm-up and
// measurement, which must catch it. The case loads the viewer afresh (`vi.resetModules`), as the composer's guards do.

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
  const { composed, recording, COPY_BACK } = await import('./composer-fixtures.js');
  const { createSeqViewSpec } = await import('@molstar/seq-core/builder');
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: 'a', value: 'MKVLATGHEKWRSPQIDNCY'.repeat(100) });
  const n = 300;
  const scale = (field: string, range: [number, number]) => ({ theme: 'field_scale', params: { field, range } });
  const start = Array.from({ length: n }, (_, i) => i * 6);
  const pairs = b.data.annotation({
    ref: 'pp',
    geometry: 'pairs',
    target: a,
    value: {
      start,
      b_target: start.map(() => 'a'),
      b_start: start.map((s) => s + 4 + (s % 5)),
      category: start.map((_, i) => `k${i % 4}`),
      count: start.map((_, i) => 1 + (i % 7)),
    },
  } as never);
  const items = b.data.annotation({
    ref: 'f',
    geometry: 'ranges',
    target: a,
    value: {
      start,
      stop: start.map((s) => s + 4),
      alt: start.map((_, i) => 'ACGT'[i % 4]),
      kind: start.map((_, i) => ['helix', 'strand', 'turn'][i % 3]),
      conf: start.map((_, i) => 0.2 + (i % 5) / 5),
      category: start.map((_, i) => `c${i % 3}`),
    },
  } as never);
  const v = b.view({ ref: 'v', label: 'V', height: 'box' }).block({ sequence: a });
  v.track({ label: 'Pairs', height: 6 }).pairs({
    data: pairs,
    display: 'arcs',
    labels: true,
    width: scale('count', [1, 3]),
    height: scale('count', [0.5, 1]),
  } as never);
  v.track({ label: 'Variants' }).variants({ data: items, thickness: scale('conf', [0.4, 1]) } as never);
  v.track({ label: 'Secondary' }).secondary({ data: items, thickness: scale('conf', [0.3, 1]) } as never);
  const t = await composed(b.getDocument(), { width: 1000, height: 400, dpr: 2 });
  const be = recording(COPY_BACK);
  t.vp?.setView(4, 0, 'api');
  t.frame(be, 7);
  // The labels are drawn (digits under the arcs' apexes), so the guard covers them.
  const digits = Array.from(be.textItems.text.slice(0, be.textItems.n)).filter((x) => /^[1-7]$/.test(x));
  if (digits.length === 0) throw new Error('no pair labels drawn');
  return { t, be };
}

/** The guard's measurement: warmed by count and by time, then `heapScaling` over `n` frames and `4n`. */
function run(n: number, frame: (i: number) => void) {
  const t0 = performance.now();
  for (let i = 0; i < FRAME_WARM || performance.now() - t0 < FRAME_WARM_MS; i++) frame(i);
  // Settled after the timed warm-up, as every timed warm-up is (lib/seqstar/AGENTS.md, Tests; LC-20).
  let settled = 0;
  settleFrame(() => frame(settled++), FRAME_SETTLE_MS);
  let i = 0;
  return heapScaling(n, () => frame(i++));
}

describe("frames of EN2b's size channels", () => {
  it('pans allocate nothing; an object allocated per frame is caught', async () => {
    const { t, be } = await setup();
    const frame = (i: number) => {
      t.vp?.panBy((i & 1) === 0 ? 5 : -5, 0, 'pointer');
      t.frame(be, 7);
    };
    const r = run(FRAME_N, frame);
    expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
    expect(controlCaught(() => run(FRAME_N, leaking(frame)))).toBe(true);
  }, 240_000);
});
