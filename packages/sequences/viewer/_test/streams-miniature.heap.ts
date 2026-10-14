/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A miniature over streamed 1-D layers allocates nothing per pan once its regions are held (GM-04;
// seqspec:sources#the-viewers-side): over the streamed genome with the `strip` navigator, a pan's pass plans the map's
// regions every pass (each 1-D layer's level over its block's whole range) and the frame draws the tracks, the overlay
// and the navigator's map, whose layers read the map's level. Not vacuous: the map is drawn from its regions (no
// placeholder band in it), and nothing is asked while measured. The positive control: the same frame allocating an
// object a step (`leaking`), warmed, settled and measured alike, is caught. Each case loads the viewer afresh
// (`vi.resetModules`).

import { describe, expect, it, vi } from 'vitest';
import {
  controlCaught,
  describeScaling,
  FRAME_SETTLE_MS,
  FRAME_WARM,
  FRAME_WARM_MS,
  heapScaling,
  isFlat,
  leaking,
  settleFrame,
} from './heap.js';

/** 1,000 steps, as the other stream guards: they catch about 22 bytes a step. */
const N = 1_000;
/** The file's warm-up, settle and measurement, as its guard and its control run them. */
const measure = (frame: () => void) => {
  const t0 = performance.now();
  for (let i = 0; i < FRAME_WARM || performance.now() - t0 < FRAME_WARM_MS; i++) frame();
  // Settled after the timed warm-up (lib/seqstar/AGENTS.md, Tests).
  settleFrame(frame, FRAME_SETTLE_MS);
  return heapScaling(N, frame);
};

/** The streamed genome (two chromosomes, letters, values and items) under the `strip` miniature at 200 kb, settled. */
async function miniatureFrame() {
  vi.resetModules();
  const { hg38Streams, streamedDoc } = await import('../../core/_test/format/stream-provider.js');
  const { recording } = await import('./composer-fixtures.js');
  const { SkinToken } = await import('../src/draw/backend.js');
  const { Dirty, SectionComposer } = await import('../src/draw/composer.js');
  const { RecordingBackend, fixtureMeasurer } = await import('../src/testing/index.js');
  const viewer = await import('../src/index.js');
  const s = hg38Streams({
    letters: true,
    chromosomes: [
      ['chr1', 248_956_422],
      ['chr2', 242_193_529],
    ],
  });
  const presets = viewer.DefaultPresets.map((p) =>
    p.geometry === 'linear' && p.name !== 'strip' ? { ...p, navigator: 'strip' } : p,
  );
  const scheduler = viewer.manualScheduler();
  const ctx = viewer.createSeqViewer(
    { ...viewer.DefaultSeqViewerSpec(), registry: [{ ...viewer.DefaultRegistry, presets }] },
    {
      scheduler,
      loaders: { load: (r) => Promise.reject(new Error(`no ${r.uri}`)), stream: s.stream },
      text: fixtureMeasurer(),
    },
  );
  const size = () => {
    for (const ref of ctx.viewports.sections.value) ctx.viewports.get(ref)?.setSize(1000, 300, 1);
  };
  const done = ctx.setDocument(streamedDoc({ items: true }));
  size();
  await done;
  await ctx.graph.update();
  const section = ctx.viewports.sections.value[0];
  const vp = ctx.viewports.get(section);
  if (!vp) throw new Error('no viewport');
  vp.zoomTo(1, 100_000_000, 100_200_000);
  const manager = ctx.streams as import('../src/streams.js').Streams;
  const composer = new SectionComposer(ctx, section);
  const be = recording();
  const flush = () => new Promise((ok) => setTimeout(ok, 0));
  // Both pan positions' regions asked and landed before measuring (a region asked at one position only stays in
  // flight within `Stream.CancelGraceMs`).
  for (let k = 0; k < 30; k++) {
    for (const dx of [9, -9]) {
      vp.panBy(dx, 0, 'pointer');
      manager.pass();
      await flush();
    }
    scheduler.advance(16);
    await ctx.graph.update();
    be.reset();
    composer.frame(be, Dirty.all);
  }
  await manager.settled();
  be.reset();
  composer.frame(be, Dirty.all);
  /** The pending bands in the map's images. */
  const pending = () => {
    let maps = 0;
    let n = 0;
    for (const c of be.images.content)
      if (c instanceof RecordingBackend) {
        maps++;
        for (let i = 0; i < c.rectItems.n; i++) if (c.rectItems.color[i] === SkinToken.pending) n++;
      }
    return { maps, n };
  };
  let i = 0;
  return {
    pending,
    requests: () => manager.stats.current.requests,
    frame: () => {
      // A pan (small ints: a double passed to a call is boxed), its pass, and the frame with the navigator.
      vp.panBy((i++ & 1) === 0 ? 9 : -9, 0, 'pointer');
      manager.pass();
      be.reset();
      composer.frame(be, Dirty.all);
    },
    dispose: () => ctx.dispose(),
  };
}

describe('a miniature over streamed layers allocates nothing per pan once its regions are held (GM-04)', () => {
  it("a pass and a frame with the strip's map, 200 kb", async () => {
    const g = await miniatureFrame();
    try {
      const p = g.pending();
      expect(p.maps).toBeGreaterThan(0);
      expect(p.n).toBe(0);
      const asked = g.requests();
      const r = measure(g.frame);
      expect(isFlat(r), describeScaling(r, N)).toBe(true);
      expect(g.requests()).toBe(asked);
    } finally {
      g.dispose();
    }
  }, 180_000);

  it('the same frame allocating an object a step is caught (positive control)', async () => {
    const g = await miniatureFrame();
    try {
      expect(controlCaught(() => measure(leaking(g.frame) as () => void))).toBe(true);
    } finally {
      g.dispose();
    }
  }, 180_000);
});
