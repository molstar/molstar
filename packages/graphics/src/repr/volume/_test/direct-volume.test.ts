/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Mat4, Tensor } from '@molstar/core/math/linear-algebra';
import { CustomProperties } from '@molstar/model/model/custom-property';
import type { Volume } from '@molstar/model/model/volume';
import { getDirectVolumeParams } from '../direct-volume.js';

/** 32³ cells of noise around 0 (σ about 1, a fixed sequence) and a few cells of 20, the maximum */
function createVolume(): Volume {
  const dim: [number, number, number] = [32, 32, 32];
  const space = Tensor.Space(dim, [2, 1, 0], Float32Array);
  const data = new Float32Array(32 * 32 * 32);
  let seed = 1;
  const random = () => (seed = (seed * 16807) % 2147483647) / 2147483647;
  for (let i = 0; i < data.length; ++i) data[i] = (random() + random() + random() - 1.5) * 2;
  for (let i = 0; i < 8; ++i) data[i * 997] = 20;
  let min = Infinity;
  let max = -Infinity;
  let sum = 0;
  for (const v of data) {
    min = Math.min(min, v);
    max = Math.max(max, v);
    sum += v;
  }
  const mean = sum / data.length;
  const sigma = Math.sqrt(data.reduce((s, v) => s + (v - mean) ** 2, 0) / data.length);
  return {
    label: 'test',
    grid: {
      transform: { kind: 'matrix', matrix: Mat4.identity() },
      cells: Tensor.create(space, Tensor.Data1(data)),
      stats: { min, max, mean, sigma },
    },
    instances: [{ transform: Mat4.identity() }],
    sourceData: { kind: 'test', name: 'test', data: {} } as any,
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
    _localPropertyData: Object.create(null),
  };
}

describe('direct-volume projection presets', () => {
  const volume = createVolume();
  const { min, max } = volume.grid.stats;
  const presets = getDirectVolumeParams({} as any, volume).controlPoints.getPresets!();
  const projections = presets.filter(([, name]) => name.startsWith('Projection'));

  it('offers the linear, bright and solid curves', () => {
    expect(projections.map(([, name]) => name)).toEqual([
      'Projection, linear (from ~+2σ)',
      'Projection, bright (~+2σ to +15σ)',
      'Projection, solid (from ~+3σ)',
    ]);
  });

  it('start dark above the mean and stay bright up to the map maximum', () => {
    const mean = (0 - min) / (max - min);
    for (const [points] of projections) {
      expect(points[0][1]).toBe(0);
      expect(points[0][0]).toBeGreaterThan(mean);
      expect(points[points.length - 1][0]).toBeCloseTo(0.999, 3);
      expect(points[points.length - 1][1]).toBeGreaterThan(0.5);
    }
  });
});
