/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The folded genes' case, for its heap guard and its quick twin (streams-items-folded.heap.ts,
// streams-items-folded.quick.ts): a pan's pass and composer frame over the demo's GENCODE slices at TP53 (about 30
// transcripts folded into two lanes), every region held. Each call loads the viewer afresh (`vi.resetModules`).

import { vi } from 'vitest';

export async function foldedFrame() {
  vi.resetModules();
  const { composed, recording } = await import('./composer-fixtures.js');
  const { Dirty } = await import('../src/draw/composer.js');
  const { demoSources } = await import('./jbrowse-sources.js');
  const { JBROWSE_SECTION, jbrowseDoc } = await import('./jbrowse-story-fixtures.js');
  const src = demoSources('nav');
  const t = await composed(jbrowseDoc('chr17', 7_661_779, 7_687_538), {
    width: 1000,
    height: 300,
    dpr: 2,
    loaders: src.loaders,
  });
  const m = t.ctx.streams as import('../src/streams.js').Streams;
  const vp = t.vp;
  if (!vp) throw new Error('no viewport');
  t.layout(JBROWSE_SECTION);
  let done = false;
  void m.settled().then(() => {
    done = true;
  });
  for (let k = 0; k < 2000 && !done; k++) {
    t.scheduler.advance(16);
    await new Promise((ok) => setTimeout(ok, 0));
  }
  const be = recording();
  let i = 0;
  const frame = () => {
    vp.panBy((i++ & 1) === 0 ? 9 : -9, 0, 'pointer');
    m.pass();
    t.frame(be, Dirty.all);
  };
  // Both pan positions' regions asked before measuring: a region asked at one position and not at the other stays in
  // flight within `Stream.CancelGraceMs` (VT-02), never asked and cancelled every other frame, so nothing is asked
  // while measured.
  for (let k = 0; k < 8; k++) {
    frame();
    await new Promise((ok) => setTimeout(ok, 0));
  }
  return { frame, backend: be, dispose: () => t.ctx.dispose(), requests: () => m.stats.current.requests };
}
