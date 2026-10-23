/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Mat4, Tensor, Vec3 } from '@molstar/core/math/linear-algebra';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CustomProperties } from '@molstar/model/model/custom-property';
import { PositionLocation } from '@molstar/model/model/position-location';
import type { Grid, Volume } from '@molstar/model/model/volume';
import { ColorTheme } from '../color.js';
import { VolumeValueColorTheme, VolumeValueColorThemeParams } from '../color/volume-value.js';

// values 0..90 along x (10 per voxel), mean 20 and sigma 30 as given, not as computed: they only need to be not 0 and 1
const N = 10;
function volume(): Volume {
  const space = Tensor.Space([N, 2, 2], [2, 1, 0], Float32Array);
  const data = space.create();
  for (let k = 0; k < 2; k++) for (let j = 0; j < 2; j++) for (let i = 0; i < N; i++) space.set(data, i, j, k, 10 * i);
  return {
    grid: {
      transform: { kind: 'matrix', matrix: Mat4.identity() },
      cells: Tensor.create(space, data),
      stats: { min: 0, max: 90, mean: 20, sigma: 30 },
    } satisfies Grid,
    instances: [{ transform: Mat4.identity() }],
    sourceData: { kind: 'test', name: 'test', data: {} } as any,
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
    _localPropertyData: Object.create(null),
  };
}

const props = (p: Partial<PD.Values<VolumeValueColorThemeParams>>) => ({
  ...PD.getDefaultValues(VolumeValueColorThemeParams),
  ...p,
});

/** Where the vertex at x lands in the palette, 0..1 */
function fraction(theme: ColorTheme<any>, x: number) {
  return (theme.color(PositionLocation(Vec3.create(x, 0, 0)), false) as number) / ColorTheme.PaletteScale;
}

describe('volume-value color theme isRelative (perf-audit-2 09/B7)', () => {
  const v = volume();
  const auto = { name: 'auto' as const, params: { symmetric: false } };

  it('places a vertex in the palette as without isRelative on the auto domain', () => {
    const abs = VolumeValueColorTheme({ volume: v }, props({ domain: auto }));
    const rel = VolumeValueColorTheme({ volume: v }, props({ domain: auto, isRelative: true }));
    for (const x of [0, 2.5, 5, 9]) {
      expect(fraction(rel, x)).toBeCloseTo(fraction(abs, x), 5);
    }
    expect(fraction(rel, 5)).toBeCloseTo(50 / 90, 5);
  });

  it('reads a custom domain in sigma units', () => {
    // [-1, 1] sigma = values -10..50; 35 (x = 3.5) is 0.75 of the way
    const rel = VolumeValueColorTheme(
      { volume: v },
      props({ domain: { name: 'custom', params: [-1, 1] }, isRelative: true }),
    );
    expect(fraction(rel, 3.5)).toBeCloseTo(0.75, 5);
  });

  it('gives direct volume rendering the domain as normalized values', () => {
    const direct = { volume: v, locationKinds: ['direct-location' as const] };
    const rel = VolumeValueColorTheme(direct, props({ domain: auto, isRelative: true }));
    expect(rel.palette!.domain![0]).toBeCloseTo(0, 5);
    expect(rel.palette!.domain![1]).toBeCloseTo(1, 5);
    const custom = VolumeValueColorTheme(
      direct,
      props({ domain: { name: 'custom', params: [-1, 1] }, isRelative: true }),
    );
    expect(custom.palette!.domain![0]).toBeCloseTo(-10 / 90, 5);
    expect(custom.palette!.domain![1]).toBeCloseTo(50 / 90, 5);
  });
});
