/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for the axis and the viewport at pointer rate (VC3; seqspec:viewer#performance), on `heapScaling` (a step
// is one move): per-move axis lookups (`indexOf`, `shownOf`, `breakAt`, the view's `x`) and both run converters on a
// forward ranged unit, a filtered, a reversed and a cyclic piece; viewport pans and anchored zooms with a live
// subscriber; programmatic (`api`) moves, each committed at once, with no `committed` subscriber (a subscriber's
// records are its own, and show). The steps are axis-cases.ts's, shared with the quick twin (axis.quick.ts). The
// positive control runs the axis lookups' step allocating an object each (`leaking`) through the same measurement,
// which must catch it.

import { describe, expect, it } from 'vitest';
import { api, lookups, N, panning, runs, subscribed, sunk } from './axis-cases.js';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from './heap.js';

describe('the axis allocates nothing per move', () => {
  it('indexOf, shownOf, breakAt and x over forward, filtered, reversed and cyclic pieces', async () => {
    const step = await lookups();
    const r = heapScaling(N, step);
    expect(sunk()).not.toBe(0.5);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(controlCaught(() => heapScaling(N, leaking(step)))).toBe(true);
  });

  it('runToShown and shownToRuns', async () => {
    const r = heapScaling(N, await runs());
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });
});

describe('the viewport allocates nothing per move', () => {
  it('pans and anchored zooms at pointer rate, with a live subscriber', async () => {
    const { step, v, moves } = await panning();
    const r = heapScaling(N, step);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(moves()).toBeGreaterThan(N);
    expect(v.live.current.interacting).toBe(true); // one settle timer for the whole gesture
  });
});

describe('an api move allocates nothing while nobody reads `committed`', () => {
  it('pans and views set with origin api: committed each, built only when read', async () => {
    const { step, v, before } = await api();
    const r = heapScaling(N, step);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    // Read, the last move's record; still at once for api (no interaction).
    expect(v.committed.value).not.toBe(before);
    expect([v.committed.value.origin, v.committed.value.scrollX, v.live.current.interacting]).toEqual([
      'api',
      v.live.current.scrollX,
      false,
    ]);
  });
  it("with a subscriber, each move's record is delivered and shows", async () => {
    const { step, delivered } = await subscribed();
    // Short runs: the records fill the young generation within a few thousand moves.
    const M = 1_000;
    const r = heapScaling(M, step);
    expect(delivered()).toBeGreaterThan(M);
    expect(isFlat(r), describeScaling(r, M)).toBe(false);
  });
});
