/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A header's hover is hover-rate (seqspec:interaction#hover, seqspec:viewer#performance): the hover alternating between
// two track headers and a group's header, each landing described by the `headerHover` behavior (its items counted, its
// runs written) and projected by the channels, allocates nothing after warm-up. Not vacuous: the described count
// alternates as the headers' items do, and the channels rebuilt the hover every step. The positive control: the same
// steps allocating an object each (`leaking`), warmed and measured alike, are caught. `toLoci` (item loci, allocating)
// is read only by a host, never here. The viewer is opened afresh per case (header-hover-cases.ts, shared with the
// quick twin, header-hover.quick.ts); the step is warmed and settled by count (`frameScaling`), then `N` and `4N`
// steps are measured.

import { expect, it } from 'vitest';
import { controlCaught, describeScaling, frameScaling, isFlat, leaking } from './heap.js';
import { headers, N } from './header-hover-cases.js';

it('headers alternating, each described and projected: flat; the count alternates, the channels rebuild', async () => {
  const h = await headers();
  const before = h.ch.rebuilds.hover;
  const r = frameScaling(h.step, N);
  expect(isFlat(r), describeScaling(r, N)).toBe(true);
  expect(h.counts).toEqual([40, 10, 5]);
  expect(h.ch.rebuilds.hover - before).toBeGreaterThanOrEqual(5 * N);
  expect(h.sink).toBeGreaterThan(0);
  h.t.ctx.dispose();
}, 300_000);

it('the same steps allocating an object each are caught (positive control)', async () => {
  const h = await headers();
  expect(controlCaught(() => frameScaling(leaking(h.step), N))).toBe(true);
  h.t.ctx.dispose();
}, 300_000);
