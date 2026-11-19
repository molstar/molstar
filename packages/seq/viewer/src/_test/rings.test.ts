/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A ring and its linear detail through the ViewModel (_spec/README.md "ViewModel", Viewports: a cyclic viewport wraps):
// the ring's window round the origin, never clamped at the total; the detail following it over the origin in two
// pieces; a pan or a zoom asked of the detail moving the ring's window (`follow.followerMoves: lead`, EC-114); a ring
// following a detail; a rotated linear view panned across its origin.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { SpecModel } from '../models/spec.js';
import { ViewModel } from '../models/view.js';

const Ring = '/views/v/ring';
const Detail = '/views/v/detail';

function ringSpec(ringFollows = false) {
  const b = createSeqSpec({ title: 'ring' });
  const c = b.data.sequence('c', { alphabet: 'dna', letters: 'ACGT'.repeat(150), cyclic: true });
  b.views.view('v', (v) => {
    if (ringFollows) {
      const d = v.section('detail', { of: c });
      d.track('l').layer('l', { type: 'letters' });
      v.section('ring', { of: c, layout: 'rings', follow: d }).track('l').layer('l', { type: 'letters' });
      return;
    }
    const ring = v.section('ring', { of: c, layout: 'rings' });
    ring.track('l').layer('l', { type: 'letters' });
    v.section('detail', { of: c, follow: ring }).track('l').layer('l', { type: 'letters' });
  });
  return b.build();
}

function setup(ringFollows = false) {
  const view = new ViewModel(new SpecModel(ringSpec(ringFollows)));
  view.setWidth(Ring, 400);
  view.setWidth(Detail, 400);
  return view;
}

describe('a ring and its linear detail in the ViewModel', () => {
  it('the ring’s window pans round the origin; the detail follows over it in two pieces', () => {
    const view = setup();
    view.setRange(Ring, 100, 160);
    expect(view.get(Detail)).toMatchObject({ start: 100, end: 160 });
    view.pan(Ring, 480);
    // past the total: the window wraps, never clamped
    expect(view.get(Ring)).toMatchObject({ start: 580, end: 640 });
    expect(view.get(Detail)).toMatchObject({ start: 580, end: 640, pieces: [580, 600, 0, 40] });
    view.pan(Ring, 30);
    expect(view.get(Ring)).toMatchObject({ start: 10, end: 70 });
    expect(view.get(Detail)?.pieces).toBeUndefined();
  });

  it('a pan or a zoom asked of the detail moves the ring’s window (the follower leads)', () => {
    const view = setup();
    view.setRange(Ring, 200, 260);
    view.pan(Detail, 20);
    expect(view.get(Ring)).toMatchObject({ start: 220, end: 280 });
    view.zoom(Detail, 2, 250);
    const r = view.get(Ring)!;
    expect(r.end - r.start).toBeCloseTo(30, 6);
    expect(view.get(Detail)).toMatchObject({ start: r.start, end: r.end });
  });

  it('a ring following a detail: the detail leads, the ring shows its whole axis with the window the detail’s range', () => {
    const view = setup(true);
    view.setRange(Detail, 590, 630);
    expect(view.get(Detail)).toMatchObject({ start: 590, end: 630 });
    // the ring follows: its viewport is the leader's range carried over the origin
    expect(view.get(Ring)).toMatchObject({ start: 590, end: 630, pieces: [590, 600, 0, 30] });
  });

  it('fit shows the whole ring; the detail with it', () => {
    const view = setup();
    view.setRange(Ring, 580, 640);
    view.fit(Ring);
    expect(view.get(Ring)).toMatchObject({ start: 0, end: 600 });
    expect(view.get(Detail)).toMatchObject({ start: 0, end: 600 });
  });
});
