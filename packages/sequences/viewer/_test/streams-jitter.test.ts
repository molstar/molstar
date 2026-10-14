/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A grace on cancellation (seqspec:sources#the-viewers-side; VT-02): a request no view needs is aborted only once no
// pass has needed it for `Stream.CancelGraceMs`, or at once when its slot is wanted, the oldest first. A pan jittering
// across a region's edge (at the view's, with no margin; at the margin's, with the defaults) asks it once and it lands,
// where it was asked and aborted every other frame, never landing while a read outlasted a frame.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { describe, expect, it } from 'vitest';
import { GRCH38 } from '../../core/_test/format/stream-provider.js';
import { RegionState, SeqViewerConfig } from '../src/index.js';
import type { Streams } from '../src/streams.js';
import { layer, opened, preparedRenderers } from './layout-fixtures.js';
import { GENOME, Margin, openGenome } from './stream-fixtures.js';

const flush = () => new Promise((ok) => setTimeout(ok, 0));
const SPAN = 4096;
/** A data region of chr2. */
const K = 24_420;
const { CancelGraceMs, MaxInFlight, Prefetch } = SeqViewerConfig.Stream;

/**
 * 40 frames of 16 ms, the pan alternating ±9 px (90 positions) with the view's right edge `before` positions short of
 * region K (or still), every read one frame long: K's requests and aborts, and its state after.
 */
async function jitter(o: { before: number; margin?: boolean; grace?: number; still?: boolean }) {
  const t = await openGenome({ chromosomes: GRCH38.slice(0, 3), gated: true });
  if (!o.margin) {
    t.ctx.config.set(Margin, 0);
    t.ctx.config.set(Prefetch, false);
  }
  if (o.grace !== undefined) t.ctx.config.set(CancelGraceMs, o.grace);
  t.ctx.config.set(MaxInFlight, 64);
  const vp = t.vp(GENOME);
  vp.zoomTo(1, K * SPAN - o.before - 10_000, K * SPAN - o.before);
  if (o.still) vp.panBy(9, 0, 'pointer');
  t.layout(GENOME);
  const p = t.streams.phylop;
  const of = (list: readonly { level?: number; start?: number }[]) =>
    list.filter((r) => r.level === -1 && r.start === K * SPAN).length;
  for (let f = 0; f < 40; f++) {
    const before = p.held;
    t.scheduler.advance(16);
    if (!o.still) vp.panBy(f % 2 === 0 ? 9 : -9, 0, 'pointer');
    t.layout(GENOME);
    t.manager.pass();
    // A read one frame long: what was asked before this frame's pass is answered now.
    p.release(before);
    await flush();
  }
  return { asked: of(p.requests), aborted: of(p.aborted), state: t.manager.state(t.layer, 1, K * SPAN + 10), t };
}

describe('a grace on cancellation (VT-02)', () => {
  it("a region at the view's edge, no margin: asked once, never aborted, and it lands", async () => {
    const r = await jitter({ before: 50 });
    expect([r.asked, r.aborted, r.state]).toEqual([1, 0, RegionState.loaded]);
  });

  it("a region at the margin's edge, the defaults: asked once, never aborted", async () => {
    const r = await jitter({ before: 50 + 10_000, margin: true });
    expect([r.asked, r.aborted]).toEqual([1, 0]);
  });

  it('`Stream.CancelGraceMs` 0 aborts at once, as before: the region churns and never lands', async () => {
    const r = await jitter({ before: 50, grace: 0 });
    expect(r.asked).toBeGreaterThan(10);
    expect(r.aborted).toBeGreaterThan(10);
    expect(r.state).toBe(RegionState.pending);
  });

  it('a still view asks once; a region left for longer than the grace is aborted', async () => {
    const r = await jitter({ before: 50, still: true });
    expect([r.asked, r.aborted, r.state]).toEqual([1, 0, RegionState.loaded]);
    const t = await openGenome({ chromosomes: GRCH38.slice(0, 3), gated: true });
    t.ctx.config.set(Margin, 0);
    t.ctx.config.set(Prefetch, false);
    t.ctx.config.set(MaxInFlight, 64);
    const vp = t.vp(GENOME);
    vp.zoomTo(1, K * SPAN, K * SPAN + 10_000);
    t.layout(GENOME);
    t.manager.pass();
    const asked = t.streams.phylop.held;
    expect(asked).toBeGreaterThan(0);
    // Far away: the old view's requests stay within the grace, and go after it.
    vp.zoomTo(1, 2 * K * SPAN, 2 * K * SPAN + 10_000);
    t.layout(GENOME);
    t.manager.pass();
    expect(t.streams.phylop.aborted.length).toBe(0);
    t.scheduler.advance(300);
    t.manager.pass();
    // The old view's data regions (a coarse one the new view needs too stays).
    const old = (list: readonly { level?: number; start?: number }[]) =>
      list.filter((r) => r.level === -1 && (r.start ?? 0) < K * SPAN + 20_000).length;
    expect(old(t.streams.phylop.aborted)).toBe(old(t.streams.phylop.requests));
    expect(old(t.streams.phylop.aborted)).toBeGreaterThan(0);
  });
});

describe("the grace's order and a local job (LB-16, LB-17)", () => {
  it('when the bound forces aborts, the oldest stale region goes first', async () => {
    const t = await openGenome({ chromosomes: GRCH38.slice(0, 3), gated: true });
    t.ctx.config.set(Margin, 0);
    t.ctx.config.set(Prefetch, false);
    t.ctx.config.set(CancelGraceMs, 100_000);
    t.ctx.config.set(MaxInFlight, 64);
    const vp = t.vp(GENOME);
    const p = t.streams.phylop;
    const views = [10_000_000, 60_000_000, 110_000_000, 160_000_000];
    const sets: number[][] = [];
    for (let v = 0; v < 3; v++) {
      const before = p.requests.length;
      t.scheduler.advance(16);
      vp.zoomTo(1, views[v], views[v] + 10_000);
      t.layout(GENOME);
      t.manager.pass();
      sets.push(p.requests.slice(before).map((r) => r.start));
    }
    expect(p.aborted.length).toBe(0);
    // Full: the fourth view's regions take the slots of the first view's, the oldest stale, before any other's.
    t.ctx.config.set(MaxInFlight, t.manager.stats.current.inFlight);
    t.scheduler.advance(16);
    vp.zoomTo(1, views[3], views[3] + 10_000);
    t.layout(GENOME);
    t.manager.pass();
    const order = p.aborted.map((r) => sets.findIndex((s) => s.includes(r.start)));
    expect(order.length).toBeGreaterThan(0);
    expect(order.every((v) => v === 0)).toBe(true);
    t.ctx.dispose();
  });

  it("a graced local job pauses: the new view's tile is computed first, as with no grace", async () => {
    /** A held alignment large enough that one tile is more than a frame's work. */
    const doc = () => {
      const AA = 'ACDEFGHIKLMNPQRSTVWY';
      const lines = ['# STOCKHOLM 1.0'];
      for (let k = 0; k < 1000; k++) {
        let r = '';
        for (let c = 0; c < 8000; c++) r += (c * 7 + k) % 9 === 0 ? '-' : AA[(c * 7 + k * 3) % 20];
        lines.push(`R${k}/1-8000 ${r}`);
      }
      const b = createSeqViewSpec({ title: 'large' });
      const value = [...lines, '//', ''].join('\n');
      const aln = b.data.alignment({ ref: 'big', source: { value, format: 'stockholm', max_bytes: 1 << 27 } });
      const v = b.view({ ref: 'large', label: 'Large', style: 'bare', height: 'box' });
      v.block({ alignment: aln });
      v.track({ label: 'Rows' }).rows();
      return b.getDocument();
    };
    const passes: number[] = [];
    // One document, opened twice: its 8,000,000 letters written once.
    const large = doc();
    for (const grace of [250, 0]) {
      const t = await opened(structuredClone(large), {
        renderers: preparedRenderers(),
        size: { 'large.main': [600, 400] },
      });
      t.ctx.config.set(CancelGraceMs, grace);
      const s = t.ctx.streams as Streams;
      const { k, m } = layer(t.ctx, 'Rows');
      const ref = m.set.strings.get(m.resolved.layers.ref[k]);
      const vp = t.vp('large.main');
      // The whole view's coarse job starts (more cells than a frame's budget), then the view jumps to the far end.
      vp.zoomTo(0, 0, 8000);
      vp.setRowScale(0.01, 0);
      t.layout('large.main');
      s.pass();
      vp.zoomTo(0, 4000, 8000);
      vp.setRowScale(0.02, 0);
      t.layout('large.main');
      let n = 0;
      while (s.state(ref, 0, 6000, 10) !== RegionState.loaded && n < 40) {
        s.pass();
        n++;
      }
      passes.push(n);
      t.ctx.dispose();
    }
    expect(passes[0]).toBe(passes[1]);
  }, 120_000);
});
