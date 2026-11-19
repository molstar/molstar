/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A story's frame under its target (_spec/README.md "Stories", Sampling): `sample` of a plan of 200 tracks' channels
// (and a view, an angle, an entering and leaving track) is a small part of a frame: under 0.05 ms a call, timed over
// 1,000 calls, the best of three in the heavy lane; the quick lane runs it unbounded.

import { bestOf, itWithin } from '@molstar/test-harness/timing';
import { describe, expect } from 'vitest';
import { sample } from '../index.js';
import { manyChannels } from './plans.js';

describe('a story’s frame', () => {
  itWithin('sample: a plan of 200 tracks’ channels, 1,000 frames under 50 ms', (within) => {
    const p = manyChannels(200);
    expect(p.channels.length).toBeGreaterThan(300);
    for (let i = 0; i < 100; i++) sample(p, i / 99);
    bestOf(within, 50, () => {
      for (let i = 0; i < 1000; i++) sample(p, (i % 997) / 996);
    });
  });
});
