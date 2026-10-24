/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { Box3D } from '@molstar/core/math/geometry';
import { OrderedSet } from '@molstar/core/data/int';
import { computeGaussianDensityTexture2d, computeGaussianDensityTexture3d } from '../gaussian-density.js';

// which GPU path a dispatch takes, without a context
const gpu = vi.hoisted(() => ({ calls: [] as string[] }));
vi.mock('../gaussian-density/gpu.js', () => ({
  GaussianDensityTexture2d: () => (gpu.calls.push('2d'), '2d-texture'),
  GaussianDensityTexture3d: () => (gpu.calls.push('3d'), '3d-texture'),
}));

describe('computeGaussianDensityTexture', () => {
  it('the 3d variant computes a 3d texture (perf-audit-2 08/B5)', async () => {
    const position = { x: [0], y: [0], z: [0], indices: OrderedSet.ofSingleton(0) };
    const args = [position, Box3D(), () => 1, { resolution: 1, radiusOffset: 0, smoothness: 1.5 }, {} as any] as const;
    await computeGaussianDensityTexture3d(...args).run();
    await computeGaussianDensityTexture2d(...args).run();
    expect(gpu.calls).toEqual(['3d', '2d']);
  });
});
