/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { Box3D } from '@molstar/core/math/geometry';
import { Mat4, Tensor } from '@molstar/core/math/linear-algebra';
import { getPalette } from '@molstar/core/util/color/palette';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CustomProperties } from '@molstar/model/model/custom-property';
import { type Grid, Volume } from '@molstar/model/model/volume';
import { VolumeSegmentColorTheme, VolumeSegmentColorThemeParams } from '../color/volume-segment.js';

/** A segmentation of `count` segments, ids 2, 4, 6, ... (not 0..count - 1) */
function segmentation(count: number): Volume {
  const space = Tensor.Space([1, 1, 1], [2, 1, 0], Float32Array);
  const v: Volume = {
    grid: {
      transform: { kind: 'matrix', matrix: Mat4.identity() },
      cells: Tensor.create(space, space.create()),
      stats: { min: 0, max: 1, mean: 0, sigma: 1 },
    } satisfies Grid,
    instances: [{ transform: Mat4.identity() }],
    sourceData: { kind: 'test', name: 'test', data: {} } as any,
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
    _localPropertyData: Object.create(null),
  };
  const segments = new Map<Volume.SegmentIndex, Set<number>>();
  const bounds: { [k: Volume.SegmentIndex]: Box3D } = {};
  for (let i = 0; i < count; i++) {
    const id = (2 * i + 2) as Volume.SegmentIndex;
    segments.set(id, new Set([i + 1]));
    bounds[id] = Box3D();
  }
  Volume.Segmentation.set(v, { segments, sets: new Map(), bounds, labels: {} });
  return v;
}

describe('volume-segment color theme', () => {
  // a segment's palette index is its place among the segments: looked up, not searched for (perf-audit-2 09/P18)
  it('colors each segment by its place, without a scan of the segments per location', () => {
    const count = 500;
    const v = segmentation(count);
    const props = PD.getDefaultValues(VolumeSegmentColorThemeParams);
    const theme = VolumeSegmentColorTheme({ volume: v }, props);
    const loc = Volume.Segment.Location(v);
    const indexOf = vi.spyOn(Array.prototype, 'indexOf');
    const colors: number[] = [];
    for (let i = 0; i < count; i++) {
      loc.segment = (2 * i + 2) as Volume.SegmentIndex;
      colors.push(theme.color(loc, false) as number);
    }
    // (other code may call indexOf: only a search of the segments counts)
    expect(indexOf.mock.contexts.filter((c) => Array.isArray(c) && c.length === count).length).toBe(0);
    indexOf.mockRestore();

    // the palette's color of each segment's place, as before
    const palette = getPalette(count, props);
    expect(colors).toEqual(Array.from({ length: count }, (_, i) => palette.color(i)));
    expect(new Set(colors).size).toBeGreaterThan(1);
  });
});
