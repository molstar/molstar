/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { arrayMax, arrayMean, arrayMin, arrayRms } from '@molstar/core/util/array';
import type { DxFile } from '@molstar/io/reader/dx/parser';
import { volumeFromDx } from '../dx.js';

describe('volume stats', () => {
  it('a DX volume reads its values once for min, max, mean and sigma (perf-audit-2 04/P18)', async () => {
    const n = 4 * 5 * 6;
    const raw = Float64Array.from({ length: n }, (_, i) => Math.sin(i) * 10 + 3);
    let reads = 0;
    const values = new Proxy(raw, {
      get(t, k) {
        if (typeof k === 'string' && k !== 'length' && /^\d+$/.test(k)) reads++;
        const v = Reflect.get(t, k, t);
        return typeof v === 'function' ? v.bind(t) : v;
      },
    }) as Float64Array;
    const dx: DxFile = {
      name: 't',
      header: { dim: Vec3.create(4, 5, 6), min: Vec3(), h: Vec3.create(1, 1, 1) },
      values,
    };
    const volume = await volumeFromDx(dx).run();
    expect(reads).toBe(n);
    expect(volume.grid.stats).toEqual({
      min: arrayMin(raw),
      max: arrayMax(raw),
      mean: arrayMean(raw),
      sigma: arrayRms(raw),
    });
  });
});
