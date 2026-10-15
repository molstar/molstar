/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of axis.heap.ts (`test:perf`): the axis lookups and run converters, the viewport's pans and zooms and
// the api moves, the guard's steps (axis-cases.ts), each with its positive control. Not here: the api move with a
// `committed` subscriber, which allocates its records by design (the heap guard checks that it does).

import { describe, expect, it } from 'vitest';
import { api, lookups, panning, runs, sunk } from './axis-cases.js';
import { quickCase } from './heap.js';

describe('the axis and the viewport, quick', () => {
  it('indexOf, shownOf, breakAt and x', async ({ task }) => {
    quickCase(task, 'axis: indexOf, shownOf, breakAt and x', await lookups());
    expect(sunk()).not.toBe(0.5);
  });

  it('runToShown and shownToRuns', async ({ task }) => {
    quickCase(task, 'axis: runToShown and shownToRuns', await runs());
  });

  it('pans and anchored zooms with a live subscriber', async ({ task }) => {
    const { step, moves } = await panning();
    quickCase(task, 'axis: pans and anchored zooms', step);
    expect(moves()).toBeGreaterThan(1_000);
  });

  it('api pans and views, nobody reading committed', async ({ task }) => {
    const { step, v, before } = await api();
    quickCase(task, 'axis: api moves, committed unread', step);
    expect(v.committed.value).not.toBe(before);
  });
});
