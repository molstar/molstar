/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST3a-2's heap guards' frame (streams-items.heap.ts, streams-items-glyphs.heap.ts): over the streamed genome with
// genes, variants and letters and without the values track (ST1's values draw allocates about 23 bytes a frame at 100 kb
// and DPR 2, under its own guard's bound, which would hide what these draws add). Two files: a case run after another in
// the same isolate allocated where it alone did not (code shared across the fresh module instances, deoptimized).

/**
 * The genome composed (DPR 2, 1,000 × 300) at [start, stop) of chr2 (0: the whole view), held or every region pending;
 * its variants drawn as `variants` (stacked) or `points`.
 */
export async function itemsFrame(
  start: number,
  stop: number,
  pending = false,
  variantsAs: 'variants' | 'points' = 'variants',
) {
  const { composed, recording } = await import('./composer-fixtures.js');
  const { Dirty } = await import('../src/draw/composer.js');
  const { GRCH38, hg38Streams } = await import('../../core/_test/format/stream-provider.js');
  const { GENOME, streamedDoc } = await import('./stream-fixtures.js');
  const s = hg38Streams({ chromosomes: GRCH38.slice(0, 3), letters: true, gated: true });
  const loaders = { load: () => Promise.reject(new Error('none')), stream: s.stream };
  const t = await composed(streamedDoc({ items: true, sequence: true, values: false, variantsAs }), {
    width: 1000,
    height: 300,
    dpr: 2,
    loaders,
  });
  const m = t.ctx.streams as import('../src/streams.js').Streams;
  const vp = t.vp;
  if (!vp) throw new Error('no viewport');
  if (stop > start) vp.zoomTo(1, start, stop);
  else vp.fit();
  const be = recording();
  const asked = () => s.genes.requests.length + s.variants.requests.length + s.genome.requests.length;
  // Pans of ±9 px with every request answered, until a round asks for nothing more (ST6: a prefetch's lead is as many
  // regions as the view's, which the in-flight bound lets through a few a pass); 8 rounds where all stay pending.
  for (let k = 0, quiet = 0; k < (pending ? 8 : 60) && (pending || quiet < 2); k++) {
    const before = asked();
    for (const dx of [9, -9]) {
      vp.panBy(dx, 0, 'pointer');
      t.layout(GENOME);
      m.pass();
      if (!pending) for (const p of [s.genome, s.phylop, s.genes, s.variants]) p.release();
      await new Promise((ok) => setTimeout(ok, 0));
    }
    t.frame(be, Dirty.all);
    quiet = asked() === before ? quiet + 1 : 0;
  }
  let i = 0;
  return {
    backend: be,
    asked,
    dispose: () => t.ctx.dispose(),
    /** A pan, its pass (none while every region is pending: the manager's churn is not the draws'), and the frame. */
    frame: () => {
      vp.panBy((i++ & 1) === 0 ? 9 : -9, 0, 'pointer');
      if (!pending) m.pass();
      t.frame(be, Dirty.all);
    },
  };
}
