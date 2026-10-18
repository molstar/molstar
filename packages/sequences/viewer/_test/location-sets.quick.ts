/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of location-sets.heap.ts (`test:perf`): pans with a band set, a highlight set and a bookmark's accent,
// and a band set in a wrapped section (the guard's cases, location-sets-cases.ts), each with its control. The guard's
// last case, a channel set on every frame, allocates by design (a rebuild per frame): it asserts that it is not flat,
// which a bytes-a-step baseline cannot, so it has no twin.

import { describe, it } from 'vitest';
import { quickCase } from './heap.js';
import { bandPans, wrappedBand } from './location-sets-cases.js';

/** A frame of several draws (about 0.4 ms) still read 1 to 3 bytes a step in some windows at 5,000 steps. */
const LONG = { warmSteps: 2_000, steps: 1_000, windows: 3, maxWarmSteps: 60_000 };

describe("location sets' frames, quick", () => {
  it("pans with a band set, a highlight set and a bookmark's accent", async ({ task }) => {
    quickCase(task, 'location-sets: band, highlight and bookmark', await bandPans(), LONG);
  });

  it('a band set in a wrapped section, painted as a highlight on the overlay (LS2)', async ({ task }) => {
    quickCase(task, 'location-sets: wrapped band', await wrappedBand());
  });
});
