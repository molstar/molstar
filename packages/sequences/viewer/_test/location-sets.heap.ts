/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// LS1's frames allocate nothing after warm-up: a section's `highlights` (a `band` set of 120 items in the main pass,
// a `highlight` set on the overlay) and a stepped bookmark's accent, through pans that copy and draw strips and the
// overlay, and a band set on a wrapped section (painted as a highlight there, LS2). The positive control: the pans
// allocating an object a frame (`leaking`), warmed and measured alike, are caught; a channel set on every frame (a
// rebuild per frame) is not flat either. Each case loads the viewer afresh (`vi.resetModules`), as the composer's
// guards do.

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
  const starts = Array.from({ length: 120 }, (_, k) => k * 16);
  const sites = b.data.annotation({
    ref: 'sites',
    geometry: 'ranges',
    target: a,
    value: { start: starts, stop: starts.map((s) => s + 5), id: starts.map((_, k) => `s${k}`) },
  });
  const highlights = [{ data: 'sites' }, { data: 'sites', items: ['s3', 's40', 's77'], style: 'highlight' }];
  const v = b.view({ ref: 'v', label: 'V', highlights } as never);
  v.block({ sequence: a });
  v.track({ label: 'L' }).letters();
  v.track({ label: 'R' }).ranges({ data: sites });
  const t = await composed(b.getDocument(), { width: 1000, height: 300, dpr: 2 });
  const be = recording(COPY_BACK);
  t.vp?.setView(4, 200, 'api');
  t.frame(be, 7);
  return { t, be };
}

/** The file's warm-up and measurement of a frame, as its guards and its control run them. */
function run(n: number, frame: (i: number) => void) {
  const t0 = performance.now();
  for (let i = 0; i < FRAME_WARM || performance.now() - t0 < FRAME_WARM_MS; i++) frame(i);
  // Settled after the timed warm-up, as every timed warm-up is (lib/seqstar/AGENTS.md, Tests; LC-20).
  let settled = 0;
  settleFrame(() => frame(settled++), FRAME_SETTLE_MS);
  let i = 0;
  return heapScaling(n, () => frame(i++));
}

describe("location sets' frames", () => {
  it("pans with a band set, a highlight set and a bookmark's accent allocate nothing", async () => {
    const { t, be } = await setup();
    const set = t.ctx.graph.assembled.value?.set;
    t.ctx.actions.run('bookmark.add', { section: t.section, loci: t.ctx.refs.toLoci('sites', { items: ['s5'] }) });
    await t.ctx.graph.update();
    expect(t.ctx.graph.assembled.value?.set).not.toBe(set);
    expect(t.ctx.actions.run('bookmark.step', { section: t.section, dir: 1 })).toBe(true);
    const r = run(FRAME_N, (i) => {
      t.vp?.panBy((i & 1) === 0 ? 7 : -7, 0, 'pointer');
      t.frame(be, 7);
    });
    expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
  }, 120_000);

  it('a band set in a wrapped section, painted as a highlight on the overlay (LS2), allocates nothing per frame', async () => {
    vi.resetModules();
    const { lit, twoBlocks, seqLoci } = await import('./overlay-fixtures.js');
    const { recording, COPY_BACK } = await import('./composer-fixtures.js');
    const t = await lit(twoBlocks({ wrapped: true }), { width: 600, height: 300, dpr: 2 });
    const be = recording(COPY_BACK);
    t.ctx.highlight.set('site', seqLoci(t.set, 'a', 20, 60), { style: 'band' });
    t.frame(be, 7);
    const r = run(FRAME_N, (i) => {
      t.vp?.panBy(0, (i & 1) === 0 ? 9 : -9, 'pointer');
      t.frame(be, 7);
    });
    expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
  }, 120_000);

  it('the same pans allocating an object a frame are caught (positive control)', async () => {
    const { t, be } = await setup();
    const pan = (i: number) => {
      t.vp?.panBy((i & 1) === 0 ? 7 : -7, 0, 'pointer');
      t.frame(be, 7);
    };
    expect(controlCaught(() => run(FRAME_N, leaking(pan)))).toBe(true);
  }, 120_000);

  it('a channel set on every frame (a rebuild per frame) is not flat', async () => {
    const { t, be } = await setup();
    const loci = t.ctx.refs.toLoci('sites', { items: ['s9'] });
    const r = run(2_000, (i) => {
      t.ctx.highlight.set('host', loci, { style: 'band' });
      t.vp?.panBy((i & 1) === 0 ? 7 : -7, 0, 'pointer');
      t.frame(be, 7);
    });
    expect(isFlat(r), describeScaling(r, 2_000)).toBe(false);
  }, 120_000);
});
