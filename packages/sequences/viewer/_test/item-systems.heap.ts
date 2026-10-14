/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST7's item coordinate systems' heap guard (the memory rule): once a transcript's coordinate system is held, a pan's
// pass and frame over its own residues (P1's domains through a document mapping onto it, the genome's variants and
// phyloP through its map) and over the genome (P1's domains along its exons) allocate nothing. A frame allocating an
// object a step fails the same guard (the positive control).

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

const N = 1_000;
const WARM = 6 * FRAME_WARM;
/** The file's warm-up, settle and measurement, as its guards and its control run them. */
const measure = (frame: () => void) => {
  const t0 = performance.now();
  for (let i = 0; i < WARM || performance.now() - t0 < FRAME_WARM_MS; i++) frame();
  // Settled after the timed warm-up (GS-02): a climb to TurboFan can outlast it.
  settleFrame(frame, FRAME_SETTLE_MS);
  return heapScaling(N, frame);
};
const guard = (frame: () => void) => {
  const r = measure(frame);
  return isFlat(r) ? '' : describeScaling(r, N);
};

/** The transcript story's document (DPR 2) on `section` at block 0's [start, stop), every region held. */
async function transcriptFrame(section: string, block: number, start: number, stop: number) {
  const { composed, recording } = await import('./composer-fixtures.js');
  const { Dirty } = await import('../src/draw/composer.js');
  const f = await import('../../core/_test/format/stream-provider.js');
  const three = f.GRCH38.slice(0, 3);
  const keys = three.map((c) => c[0]);
  const s = f.hg38Streams({ chromosomes: three, gated: true });
  const genes = f.scripted(f.itemsManifest(three, f.GENE_LEVELS, f.VOCABULARY.genes), (r) =>
    f.transcriptItems(r, keys),
  );
  genes.gated = true;
  const stream = (src: { stream: string | null }) => (src.stream === 'genes' ? genes.provider : s.stream(src));
  const loaders = { load: () => Promise.reject(new Error('none')), stream };
  const t = await composed(f.transcriptDoc({ item: true }), { width: 1000, height: 400, dpr: 2, loaders, section });
  const m = t.ctx.streams as import('../src/streams.js').Streams;
  const settle = async () => {
    for (let k = 0; k < 8; k++) {
      for (const sec of t.ctx.viewports.sections.value) t.layout(sec);
      m.pass();
      for (const p of [s.genome, s.phylop, genes, s.variants]) p.release();
      await new Promise((ok) => setTimeout(ok, 0));
      await t.ctx.graph.update();
    }
  };
  await settle();
  const vp = t.ctx.viewports.get(section);
  if (!vp) throw new Error('no viewport');
  vp.zoomTo(block, start, stop);
  const be = recording();
  for (let k = 0; k < 30; k++) {
    for (const dx of [9, -9]) vp.panBy(dx, 0, 'pointer');
    await settle();
    t.frame(be, Dirty.all);
  }
  /** On T1's own block: the variants' features and phyloP's values read through its map (not vacuous). */
  const throughMap = () => {
    const out = new Int32Array(7);
    let features = 0;
    const st = m.get('tx.variants-2');
    for (let k = 0; st?.items; k++) {
      const r = st.items(0, 450, 650, k, out);
      if (r === null) break;
      features += r.features?.count ?? 0;
    }
    const v = new Float64Array(4);
    m.get('tx.values')?.values?.(0, 500, 501, v);
    return [features > 0, v[3] > 0];
  };
  let i = 0;
  return {
    throughMap,
    backend: be,
    dispose: () => t.ctx.dispose(),
    frame: () => {
      vp.panBy((i++ & 1) === 0 ? 9 : -9, 0, 'pointer');
      m.pass();
      t.frame(be, Dirty.all);
    },
  };
}

describe("an item coordinate system's tracks allocate nothing per frame once held", () => {
  for (const [label, section, block, start, stop] of [
    ["the transcript's own residues, 450–650", 'tx.item', 0, 450, 650],
    ['the genome along its exons', 'tx.genome', 1, 100_009_000, 100_095_000],
  ] as const)
    it(`a pass and a frame: ${label}`, async () => {
      vi.resetModules();
      const g = await transcriptFrame(section, block, start, stop);
      g.frame();
      // Not vacuous: features drawn (rects); on T1's residues, the genome's data read through its map.
      expect(g.backend.rectItems.n).toBeGreaterThan(2);
      if (section === 'tx.item') expect(g.throughMap()).toEqual([true, true]);
      expect(guard(g.frame)).toBe('');
      g.dispose();
    }, 120_000);

  it('the same frame allocating an object a step is caught (positive control)', async () => {
    vi.resetModules();
    const g = await transcriptFrame('tx.item', 0, 450, 650);
    expect(controlCaught(() => measure(leaking(g.frame) as () => void))).toBe(true);
    g.dispose();
  }, 120_000);
});
