/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A story's frame allocates nothing once warm (_spec/README.md "Stories", Sampling): `sample` of a plan of every
// channel kind (view channels field by field, numbers, an angle, colors, an entering and a leaving track, cuts) at a
// moving time, into the plan's own frame. By hand, never in CI: `pnpm -C packages/seq/core test:perf`.

import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import { sample } from '../index.js';
import { manyChannels } from './plans.js';

describe('a story’s frame', () => {
  it('sample: 0 bytes warm', ({ task }) => {
    const p = manyChannels(40);
    expect(new Set(p.channels.map((c) => c.kind)).size).toBeGreaterThanOrEqual(6);
    // the times made once, held boxed (an array of mixed elements): a double the step computed and passed would be
    // boxed by the call, the step's bytes, not `sample`'s
    const ts: unknown[] = Array.from({ length: 997 }, (_, k) => k / 996 + 1e-9);
    ts.push('mixed');
    const { r } = quickCase(task, 'story sample', (i) => sample(p, ts[i % 997] as number), {});
    console.info(`story sample: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });
});
