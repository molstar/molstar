/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The prefetch case, for its heap guard and its quick twin (streams-prefetch.heap.ts, streams-prefetch.quick.ts): a pass
// with prefetch on, once what it asks is held: the lead past the margin as a pan alternates sides, the next level as a
// zoom alternates in and out, and the snapshots' focus regions (two hints, held). A step is the motion and its pass. Each
// call loads the viewer afresh (`vi.resetModules`).

import { vi } from 'vitest';

export const MOTIONS = ['pan', 'zoom'] as const;
export type Motion = (typeof MOTIONS)[number];

/** The case's name, as both tiers list it. */
export const motionName = (motion: Motion) =>
  `a pass as the view ${motion === 'pan' ? 'pans to and fro' : 'zooms in and out'}, hints held`;

/** The pass over the `motion`, every region it asks held; `pass` is the step. */
export async function prefetchCase(motion: Motion) {
  vi.resetModules();
  return prefetchFrame(motion);
}

async function prefetchFrame(motion: 'pan' | 'zoom') {
  const { createSeqViewSpec } = await import('@molstar/seq-core/builder');
  const { GRCH38, hg38Streams } = await import('../../core/_test/format/stream-provider.js');
  const { composed } = await import('./composer-fixtures.js');
  const { GENOME } = await import('./stream-fixtures.js');
  const AT = 100_000_000;
  const b = createSeqViewSpec({ title: 'hg38' });
  const genome = b.data.source({ ref: 'genome', stream: 'genome', data_kind: 'letters' });
  const hg38 = b.data.complex({ ref: 'hg38', source_ref: genome, label: 'GRCh38' });
  b.data.source({ ref: 'phylop_src', stream: 'phylop', data_kind: 'values', coordinate_system: 'hg38' });
  const phylop = b.data.annotation({ ref: 'phylop', geometry: 'values', source_ref: 'phylop_src' });
  const v = b.view({ ref: 'browser', label: 'hg38' });
  v.block({ complex: hg38 });
  v.focus({ ref: 'browser.focus', target: 'hg38/chr2', start: AT, stop: AT + 10_000 });
  v.track({ label: 'Ruler' }).ruler();
  v.track({ label: 'phyloP' }).values({ data: phylop });
  for (const [ref, at] of [
    ['a', AT],
    ['b', 150_000_000],
    ['c', 20_000_000],
  ] as const)
    b.snapshot({
      ref,
      duration_ms: 1000,
      transition_ms: 0,
      ops: [{ op: 'set', ref: 'browser.focus', params: { start: at, stop: at + 10_000 } }],
    });
  const s = hg38Streams({ chromosomes: GRCH38.slice(0, 3) });
  const t = await composed(b.getDocument(), {
    width: 1000,
    height: 200,
    dpr: 2,
    loaders: { load: () => Promise.reject(new Error('none')), stream: s.stream },
  });
  const m = t.ctx.streams as import('../src/streams.js').Streams;
  const vp = t.vp;
  if (!vp) throw new Error('no viewport');
  await t.ctx.snapshots.apply('a');
  vp.zoomTo(1, AT, AT + 10_000);
  const flush = () => new Promise((ok) => setTimeout(ok, 0));
  // The cell at the first view and at a view a third as wide (a zoom in), set through the viewport's own verb.
  const wide = vp.live.current.cell;
  const narrow = wide * 3;
  const x = vp.live.current.scrollX;
  for (let k = 0; k < 8; k++) {
    for (const [dx, cell] of [
      [9, wide],
      [-9, narrow],
    ] as const) {
      if (motion === 'pan') vp.panBy(dx, 0, 'pointer');
      else vp.setView(cell, x, 'api');
      t.layout(GENOME);
      m.pass();
      await flush();
    }
  }
  let i = 0;
  return {
    manager: m,
    asked: () => s.phylop.requests.length,
    prefetched: () => s.phylop.requests.filter((r) => r.priority === 'prefetch').length,
    pass: () => {
      // A step of the motion (small ints: a double passed to a call is boxed), and its pass.
      if (motion === 'pan') vp.panBy((i++ & 1) === 0 ? 9 : -9, 0, 'pointer');
      else vp.setView((i++ & 1) === 0 ? wide : narrow, x, 'api');
      m.pass();
    },
    dispose: () => t.ctx.dispose(),
  };
}
