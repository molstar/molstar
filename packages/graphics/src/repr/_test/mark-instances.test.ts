/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Interval } from '@molstar/core/data/int/interval';
import { MarkerAction } from '@molstar/core/util/marker-action';
import { ValueCell } from '@molstar/core/util/value-cell';
import { EveryLoci, type Loci } from '@molstar/model/model/loci';
import { createMarkers, setMarkerInstanceCounts } from '@molstar/graphics/geo/geometry/marker-data';
import type { GraphicsRenderObject } from '@molstar/graphics/gl/render-object';
import { Visual } from '../visual.js';

const GroupCount = 10;
const InstanceCount = 100;

type Values = {
  tMarker: ValueCell<{ array: Uint8Array }>;
  markerInstanceCounts: ValueCell<Uint32Array>;
  markedInstanceCount: ValueCell<number>;
};

/** The values `Visual.mark` reads: markers per group and instance, or per instance */
function renderObject(granularity: 'groupInstance' | 'instance') {
  const instance = granularity === 'instance';
  const values = {
    ...createMarkers(instance ? InstanceCount : GroupCount * InstanceCount, granularity),
    uGroupCount: ValueCell.create(GroupCount),
    instanceCount: ValueCell.create(InstanceCount),
    uInstanceGranularity: ValueCell.create(instance),
  };
  return { values } as unknown as GraphicsRenderObject;
}

/** A loci of the entries `[start, end)` */
const range = (start: number, end: number) => ({ kind: 'test-range', start, end }) as unknown as Loci;
const lociApply: Visual.LociApply = (loci, apply) => {
  const { start, end } = loci as unknown as { start: number; end: number };
  return apply(Interval.ofBounds(start, end));
};

const valuesOf = (ro: GraphicsRenderObject) => ro.values as unknown as Values;

/** The counts and the marked instances as counted anew from the markers */
function expectCounted(ro: GraphicsRenderObject, groupCount: number) {
  const { tMarker, markerInstanceCounts, markedInstanceCount } = valuesOf(ro);
  const counts = new Uint32Array(InstanceCount);
  const marked = setMarkerInstanceCounts(tMarker.ref.value.array, InstanceCount, groupCount, counts);
  expect(Array.from(markerInstanceCounts.ref.value.subarray(0, InstanceCount))).toEqual(Array.from(counts));
  expect(markedInstanceCount.ref.value).toBe(marked);
}

describe('marked instances', () => {
  it('counts the marked groups per instance through marks of ranges', () => {
    const ro = renderObject('groupInstance');
    const { markerInstanceCounts, markedInstanceCount } = valuesOf(ro);

    // sized by the first mark
    Visual.mark(ro, range(32, 36), MarkerAction.Highlight, lociApply);
    expect(markerInstanceCounts.ref.value[3]).toBe(4);
    expect(markedInstanceCount.ref.value).toBe(1);

    // across instances: 5 (5 groups), 6 (all), 7 (3 groups)
    Visual.mark(ro, range(55, 73), MarkerAction.Select, lociApply);
    expect(Array.from(markerInstanceCounts.ref.value.subarray(5, 8))).toEqual([5, 10, 3]);
    expect(markedInstanceCount.ref.value).toBe(4);

    // a group both highlighted and selected counts once
    Visual.mark(ro, range(55, 56), MarkerAction.Highlight, lociApply);
    expect(markerInstanceCounts.ref.value[5]).toBe(5);

    Visual.mark(ro, range(30, 40), MarkerAction.RemoveHighlight, lociApply);
    expect(markerInstanceCounts.ref.value[3]).toBe(0);
    expect(markedInstanceCount.ref.value).toBe(3);
    expectCounted(ro, GroupCount);
  });

  it('counts all instances where everything is marked or cleared', () => {
    const ro = renderObject('groupInstance');
    const { markerInstanceCounts, markedInstanceCount } = valuesOf(ro);

    Visual.mark(ro, range(0, 5), MarkerAction.Highlight, lociApply);
    Visual.mark(ro, EveryLoci, MarkerAction.Select, lociApply);
    expect(markedInstanceCount.ref.value).toBe(InstanceCount);
    expect(markerInstanceCounts.ref.value.every((c) => c === GroupCount)).toBe(true);
    expectCounted(ro, GroupCount);

    // removes one of two flags: the markers read again
    Visual.mark(ro, EveryLoci, MarkerAction.RemoveHighlight, lociApply);
    expectCounted(ro, GroupCount);

    Visual.mark(ro, EveryLoci, MarkerAction.Deselect, lociApply);
    expect(markedInstanceCount.ref.value).toBe(0);
    expectCounted(ro, GroupCount);

    // all covered by a range: as everything
    Visual.mark(ro, range(0, GroupCount * InstanceCount), MarkerAction.Highlight, lociApply);
    expect(markedInstanceCount.ref.value).toBe(InstanceCount);
    expectCounted(ro, GroupCount);
  });

  it('counts a sequence of marks as the markers do', () => {
    const ro = renderObject('groupInstance');
    const actions = [MarkerAction.Highlight, MarkerAction.Select, MarkerAction.RemoveHighlight, MarkerAction.Deselect];
    let seed = 1;
    const next = (n: number) => {
      seed = (seed * 16807) % 2147483647;
      return seed % n;
    };
    for (let i = 0; i < 200; ++i) {
      const start = next(GroupCount * InstanceCount);
      const end = Math.min(GroupCount * InstanceCount, start + 1 + next(40));
      Visual.mark(ro, range(start, end), actions[next(actions.length)], lociApply);
    }
    expectCounted(ro, GroupCount);
  });

  it('counts per instance with instance markers', () => {
    const ro = renderObject('instance');
    const { markerInstanceCounts, markedInstanceCount } = valuesOf(ro);
    Visual.mark(ro, range(10, 13), MarkerAction.Highlight, lociApply);
    expect(Array.from(markerInstanceCounts.ref.value.subarray(9, 14))).toEqual([0, 1, 1, 1, 0]);
    expect(markedInstanceCount.ref.value).toBe(3);
    expectCounted(ro, 1);
  });

  it('sizes the counts from the markers where something is marked already', () => {
    const ro = renderObject('groupInstance');
    const values = valuesOf(ro);
    Visual.mark(ro, range(32, 36), MarkerAction.Highlight, lociApply);
    // as values made before the counts were kept
    ValueCell.update(values.markerInstanceCounts, new Uint32Array(0));
    ValueCell.update(values.markedInstanceCount, 0);
    Visual.mark(ro, range(50, 51), MarkerAction.Highlight, lociApply);
    expect(values.markedInstanceCount.ref.value).toBe(2);
    expectCounted(ro, GroupCount);
  });
});
