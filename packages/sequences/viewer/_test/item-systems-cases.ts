/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The item coordinate systems' cases, for their heap guard and its quick twin (item-systems.heap.ts,
// item-systems.quick.ts): the transcript story's document with every region held, a pan's pass and frame over the
// transcript's own residues and over the genome along its exons.

import { vi } from 'vitest';

/** The cases: a label, the section, and block 0 or 1's [start, stop). */
export const SYSTEMS = [
  ["the transcript's own residues, 450–650", 'tx.item', 0, 450, 650],
  ['the genome along its exons', 'tx.genome', 1, 100_009_000, 100_095_000],
] as const;

/** The transcript story's document (DPR 2) on `section` at block 0's [start, stop), every region held. */
export async function transcriptFrame(section: string, block: number, start: number, stop: number) {
  vi.resetModules();
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
