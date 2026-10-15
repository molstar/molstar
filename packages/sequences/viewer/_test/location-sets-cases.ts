/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// LS1's location-sets cases, for their heap guard and its quick twin (location-sets.heap.ts, location-sets.quick.ts):
// a section's `highlights` (a `band` set of 120 items in the main pass, a `highlight` set on the overlay) and a stepped
// bookmark's accent through pans that copy and draw strips and the overlay (`pans`, `bandPans`), a band set on a
// wrapped section painted as a highlight there (`wrappedBand`, LS2), and a channel set on every frame (`rebuilding`,
// not flat by design: the guard's negative case). Each call loads the viewer afresh (`vi.resetModules`) and returns
// the case's step (it gets a running index).

import { expect, vi } from 'vitest';

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

/** Pans over the section with its highlights, a step a pan and a whole frame. */
export async function pans() {
  const { t, be } = await setup();
  return (i: number) => {
    t.vp?.panBy((i & 1) === 0 ? 7 : -7, 0, 'pointer');
    t.frame(be, 7);
  };
}

/** The same pans with a bookmark added and stepped to (its accent in the frames). */
export async function bandPans() {
  const { t, be } = await setup();
  const set = t.ctx.graph.assembled.value?.set;
  t.ctx.actions.run('bookmark.add', { section: t.section, loci: t.ctx.refs.toLoci('sites', { items: ['s5'] }) });
  await t.ctx.graph.update();
  expect(t.ctx.graph.assembled.value?.set).not.toBe(set);
  expect(t.ctx.actions.run('bookmark.step', { section: t.section, dir: 1 })).toBe(true);
  return (i: number) => {
    t.vp?.panBy((i & 1) === 0 ? 7 : -7, 0, 'pointer');
    t.frame(be, 7);
  };
}

/** A band set in a wrapped section, painted as a highlight on the overlay (LS2), through vertical pans. */
export async function wrappedBand() {
  vi.resetModules();
  const { lit, twoBlocks, seqLoci } = await import('./overlay-fixtures.js');
  const { recording, COPY_BACK } = await import('./composer-fixtures.js');
  const t = await lit(twoBlocks({ wrapped: true }), { width: 600, height: 300, dpr: 2 });
  const be = recording(COPY_BACK);
  t.ctx.highlight.set('site', seqLoci(t.set, 'a', 20, 60), { style: 'band' });
  t.frame(be, 7);
  return (i: number) => {
    t.vp?.panBy(0, (i & 1) === 0 ? 9 : -9, 'pointer');
    t.frame(be, 7);
  };
}

/** A channel set on every frame (a rebuild per frame) with the same pans: allocates, by design. */
export async function rebuilding() {
  const { t, be } = await setup();
  const loci = t.ctx.refs.toLoci('sites', { items: ['s9'] });
  return (i: number) => {
    t.ctx.highlight.set('host', loci, { style: 'band' });
    t.vp?.panBy((i & 1) === 0 ? 7 : -7, 0, 'pointer');
    t.frame(be, 7);
  };
}
