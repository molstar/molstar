/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// VL-14: a section that left the view draws nothing and asks for no more frames, whatever it was animating, so a
// surface kept mounted for it stops polling; nor does one that turned wrapped (GS-10's scenario).

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { expect, it } from 'vitest';
import { AxisGeometry } from '../src/axis.js';
import { composed, recording } from './composer-fixtures.js';

it('a frame of a section out of the current view clears `animating` and `again`', async () => {
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: 'a', value: 'MKVLATGHEKWRSPQIDNCY'.repeat(5) });
  const v1 = b.view({ ref: 'v1', label: 'V1' });
  v1.section({ ref: 's1' }).block({ sequence: a });
  v1.track({ label: 'L' }).letters();
  b.view({ ref: 'v2', label: 'V2' }).block({ sequence: a }).track({ label: 'R' }).ruler();
  const t = await composed(b.getDocument(), { width: 600, height: 200, dpr: 1 });
  expect(t.frame(recording(), 7).opCount).toBeGreaterThan(0);
  t.ctx.views.setCurrent('v2');
  await t.ctx.graph.update();
  expect(t.ctx.layout.axis(t.section)?.axis ?? null).toBeNull();
  // As if a region was fading when the view switched.
  t.composer.animating = true;
  t.composer.again = true;
  t.composer.redraw();
  expect(t.frame(recording(), 7).opCount).toBe(0);
  expect([t.composer.animating, t.composer.again]).toEqual([false, false]);
});

it('a frame of an animating section turned wrapped clears `animating` (LC-40: only linear regions fade)', async () => {
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: 'a', value: 'MKVLATGHEKWRSPQIDNCY'.repeat(5) });
  const v1 = b.view({ ref: 'v1', label: 'V1' });
  v1.section({ ref: 's1' }).block({ sequence: a });
  v1.track({ label: 'L' }).letters();
  const t = await composed(b.getDocument(), { width: 600, height: 200, dpr: 1 });
  expect(t.frame(recording(), 7).opCount).toBeGreaterThan(0);
  await t.ctx.applyPatch({
    format: 'seq_view_spec_patch',
    ops: [{ op: 'set', ref: 'v1', params: { layout: 'text' } }],
  });
  await t.ctx.graph.update();
  expect(t.layout(t.section).geometry).toBe(AxisGeometry.wrapped);
  // As if a region was fading when the section turned wrapped.
  t.composer.animating = true;
  t.composer.redraw();
  t.frame(recording(), 7);
  expect(t.composer.animating).toBe(false);
});
