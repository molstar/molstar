/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of fingerprint.heap.ts (`test:perf`): the same cases (fingerprint-cases.ts), each on a viewer loaded
// afresh with its positive control: every frame the kit draws on each scene a quick case (`KIT` counts, ids by
// ../kit-cases.ts `frameIds`), and a hit a step over every row and stripe of three scenes (`HIT` counts).

import { describe, expect, it } from 'vitest';
import type { AllocOptions } from '../../../core/_test/model/alloc.js';
import { quickCase } from '../heap.js';
import { frameIds, KIT } from '../kit-cases.js';
import { fingerprintKit, HIT_SCENES, hitSweep, SCENES } from './fingerprint-cases.js';

/** A hit is about a microsecond: windows of 10,000 hits. */
const HIT: AllocOptions = { warmSteps: 2_000, steps: 10_000, windows: 3, maxWarmSteps: 400_000 };

describe("the fingerprint's frames, quick", () => {
  it.for(SCENES)('%s', async (name, { task }) => {
    const id = frameIds(`fingerprint: ${name}`);
    const heap = (await fingerprintKit(name, (frame, at) => quickCase(task, id(at), frame, KIT))).filter(
      (c) => c.check === 'heap',
    );
    expect(heap.length).toBeGreaterThanOrEqual(6);
    expect(heap.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
  });
});

describe("the fingerprint's hit, quick", () => {
  it.for(HIT_SCENES)('%s', async (name, { task }) => {
    const { step } = await hitSweep(name);
    quickCase(task, `fingerprint: hit on ${name}`, step, HIT);
  });
});
