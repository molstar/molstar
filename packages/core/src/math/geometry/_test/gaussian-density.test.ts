/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { OrderedSet } from '@molstar/core/data/int/ordered-set';
import type { RuntimeContext } from '@molstar/core/task';
import { getBoundary } from '../boundary.js';
import { GaussianDensityCPU } from '../gaussian-density/cpu.js';

describe('GaussianDensityCPU', () => {
  it('yields within the pass for small radii (perf-audit-2 02/B2)', async () => {
    // 27,000 atoms of radius 0.3 at resolution 0.5: about 23.5k atoms of work per chunk (1e5 / r⁹ gave billions)
    const n = 30;
    const count = n * n * n;
    const x = new Float32Array(count),
      y = new Float32Array(count),
      z = new Float32Array(count);
    for (let i = 0; i < count; i++) {
      x[i] = (i % n) * 2;
      y[i] = (Math.floor(i / n) % n) * 2;
      z[i] = Math.floor(i / (n * n)) * 2;
    }
    const position = { indices: OrderedSet.ofBounds(0, count), x, y, z };
    const { box } = getBoundary(position);

    let updates = 0;
    const ctx: RuntimeContext = {
      shouldUpdate: true,
      isSynchronous: false,
      update: async () => {
        updates++;
      },
    };
    await GaussianDensityCPU(ctx, position, box, () => 0.3, { resolution: 0.5, radiusOffset: 0, smoothness: 1.5 });
    expect(updates).toBe(2);
  });
});
