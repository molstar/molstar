/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Interval, OrderedSet } from '@molstar/core/data/int';
import { Box3D } from '@molstar/core/math/geometry';
import { Mat4, Tensor } from '@molstar/core/math/linear-algebra';
import { CustomProperties } from '@molstar/model/model/custom-property';
import { type Grid, Volume } from '@molstar/model/model/volume';
import { eachVolumeLoci } from '../util.js';

/** 4 cells labelled 1, 2, 1, 0 (sets 1 and 2 are segments 10 and 20), drawn twice: 2 instances */
function segmentation(): Volume {
  const space = Tensor.Space([4, 1, 1], [2, 1, 0], Float32Array);
  const data = space.create();
  data.set([1, 2, 1, 0]);
  const v: Volume = {
    grid: {
      transform: { kind: 'matrix', matrix: Mat4.identity() },
      cells: Tensor.create(space, data),
      stats: { min: 0, max: 2, mean: 1, sigma: 1 },
    } satisfies Grid,
    instances: [{ transform: Mat4.identity() }, { transform: Mat4.fromTranslation(Mat4(), [10, 0, 0]) }],
    sourceData: { kind: 'test', name: 'test', data: {} } as any,
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
    _localPropertyData: Object.create(null),
  };
  const segments = new Map<Volume.SegmentIndex, Set<number>>([
    [10 as Volume.SegmentIndex, new Set([1])],
    [20 as Volume.SegmentIndex, new Set([2])],
  ]);
  const sets = new Map<number, Set<Volume.SegmentIndex>>([
    [1, new Set([10 as Volume.SegmentIndex])],
    [2, new Set([20 as Volume.SegmentIndex])],
  ]);
  Volume.Segmentation.set(v, { segments, sets, bounds: { 10: Box3D(), 20: Box3D() }, labels: {} });
  return v;
}

function cellsOf(loci: Volume.Segment.Loci, volume: Volume) {
  const cells: number[] = [];
  eachVolumeLoci(loci, volume, undefined, (i) => {
    for (let k = Interval.start(i); k < Interval.end(i); k++) cells.push(k);
    return true;
  });
  return cells.sort((a, b) => a - b);
}

// the cells of a segment loci over the whole grid (a slice or direct volume of the segmentation): by the instances the
// loci names, each element its own segments (perf-audit-2 07/N3)
describe('eachVolumeLoci of segment loci', () => {
  const volume = segmentation();
  const one = (segment: number, instance: number) => ({
    segments: OrderedSet.ofSingleton(segment as Volume.SegmentIndex),
    instances: OrderedSet.ofSingleton(instance as Volume.InstanceIndex),
  });

  it('marks the instance the loci names', () => {
    // segment 10 (label 1: cells 0 and 2) of instance 1: cells 4 + 0 and 4 + 2
    expect(cellsOf(Volume.Segment.Loci(volume, [one(10, 1)]), volume)).toEqual([4, 6]);
  });

  it('marks each element its own segments', () => {
    // segment 10 of instance 0, segment 20 (label 2: cell 1) of instance 1
    expect(cellsOf(Volume.Segment.Loci(volume, [one(10, 0), one(20, 1)]), volume)).toEqual([0, 2, 5]);
  });
});
