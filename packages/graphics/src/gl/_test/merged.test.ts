/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import {
  createRenderObject,
  createMergedRenderObject,
  createMergedRecordsRenderObject,
  canMergeRenderObjects,
} from '../render-object.js';
import { canMergeValues, type MergeableValues, type MergeRecord, pickMergeRecord } from '../renderable/merged.js';
import { Scene } from '../scene.js';
import { SegmentedArray } from '../segmented-array.js';
import { getGLContext, tryGetGLContext } from './gl.js';
import { setDebugMode } from '@molstar/core/util/debug';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { Points } from '@molstar/graphics/geo/geometry/points/points';
import { ValueCell } from '@molstar/core/util';
import { Color } from '@molstar/core/util/color/color';
import { createObjectColor, createTextureColor } from '../../geo/geometry/color-data.js';
import { createDataArray } from '../renderable/util.js';
import { NullLocation } from '@molstar/model/model/location';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';
import { Cylinders, CylinderStride } from '@molstar/graphics/geo/geometry/cylinders/cylinders';
import { LocationIterator } from '@molstar/graphics/geo/util/location-iterator';
import { getDrawnBoundingSphere } from '../renderable/util.js';
import { SegmentTexels } from '../renderable/schema.js';
import { Mat4 } from '@molstar/core/math/linear-algebra';
import { Sphere3D } from '@molstar/core/math/geometry';
import { getDataWords } from '../../geo/geometry/center-group.js';

/** `count` items of a position group's array: xyz and the group (an int32 as its bits) */
function items(array: Float32Array, count: number) {
  const words = getDataWords(array);
  return Array.from({ length: count * 4 }, (_, i) => (i % 4 === 3 ? words[i] : array[i]));
}

function createSpheresRenderObject(centers: Float32Array, groups: Float32Array, count: number) {
  const spheres = Spheres.create(centers, groups, count);
  const props = PD.getDefaultValues(Spheres.Params);
  const values = Spheres.Utils.createValuesSimple(spheres, props, ColorNames.orange, 1);
  const state = Spheres.Utils.createRenderableState(props);
  return createRenderObject('spheres', values, state, -1);
}

function createMeshRenderObject(vertexCount: number, offset: number) {
  const vertices = new Float32Array(vertexCount * 3);
  const normals = new Float32Array(vertexCount * 3);
  const groups = new Float32Array(vertexCount);
  for (let i = 0; i < vertexCount; ++i) {
    vertices[i * 3] = i + offset;
    normals[i * 3 + 2] = 1;
  }
  const triangleCount = Math.floor(vertexCount / 3);
  const indices = new Uint32Array(triangleCount * 3);
  for (let i = 0, il = triangleCount * 3; i < il; ++i) indices[i] = i;
  const mesh = Mesh.create(vertices, indices, normals, groups, vertexCount, triangleCount);
  const props = PD.getDefaultValues(Mesh.Params);
  const values = Mesh.Utils.createValuesSimple(mesh, props, ColorNames.orange, 1);
  const state = Mesh.Utils.createRenderableState(props);
  return createRenderObject('mesh', values, state, -1);
}

function createMembers() {
  const a = createSpheresRenderObject(new Float32Array([0, 0, 0, 1, 0, 0, 0, 1, 0]), new Float32Array([0, 1, 2]), 3);
  const b = createSpheresRenderObject(new Float32Array([2, 0, 0, 0, 2, 0]), new Float32Array([0, 1]), 2);
  return [a, b];
}

describe('merged', () => {
  it('bounds the members drawing something: an empty member not at its origin', () => {
    const far = () => createSpheresRenderObject(new Float32Array([100, 0, 0, 101, 0, 0]), new Float32Array([0, 1]), 2);
    const empty = createSpheresRenderObject(new Float32Array(0), new Float32Array(0), 0);
    expect(empty.values.drawCount.ref.value).toBe(0);
    const merged = createMergedRenderObject([far(), empty, far()]);
    const sphere = merged.values.boundingSphere.ref.value;
    expect(sphere.center[0]).toBeGreaterThan(99);
    expect(sphere.radius).toBeLessThan(2);
  });

  it("merges members colored by a palette, the palette the shared values'", () => {
    const [a, b] = createMembers();
    for (const o of [a, b]) ValueCell.update((o.values as unknown as MergeableValues).dUsePalette!, true);
    expect(canMergeValues('spheres', [a.values as unknown as MergeableValues, b.values as any])).toBe(true);
    const merged = createMergedRenderObject([a, b]);
    const mv = merged.values as unknown as MergeableValues;
    expect(mv.dUsePalette!.ref.value).toBe(true);
    expect(mv.tPalette).toBe((a.values as unknown as MergeableValues).tPalette);
    expect(mv.uPaletteDomain).toBe((a.values as unknown as MergeableValues).uPaletteDomain);
  });

  it("checks the layers' types unless asked not to (layers the merged object sets again its own way)", () => {
    const [a] = createMembers();
    const v = a.values as unknown as MergeableValues;
    ValueCell.update(v.dOverpaint, true);
    ValueCell.update(v.overpaintType, 'volumeInstance');
    expect(canMergeValues('spheres', [v])).toBe(false);
    expect(canMergeValues('spheres', [v], false)).toBe(true);
  });

  it('values', () => {
    const [a, b] = createMembers();
    expect(canMergeRenderObjects([a, b])).toBe(true);

    const merged = createMergedRenderObject([a, b]);
    const mv = merged.values;

    expect(mv.drawCount.ref.value).toBe(5 * 6);
    // one data vertex per sphere ('vertex' color and size), six drawn
    expect(mv.uVertexCount.ref.value).toBe(5);
    expect(mv.instanceCount.ref.value).toBe(2);
    expect(mv.uInstanceCount.ref.value).toBe(2);

    const { segments } = merged.merged;
    expect(segments.count).toBe(2);
    expect(segments.instanceBases).toEqual([0, 1]);
    expect(segments.instanceCounts).toEqual([1, 1]);
    expect(segments.vertexBases).toEqual([0, 3]);

    expect(Array.from(mv.aInstance.ref.value.subarray(0, 2))).toEqual([0, 1]);

    const g0 = a.values.uGroupCount.ref.value;
    const g1 = b.values.uGroupCount.ref.value;
    // per instance its segment's group count and offsets and the segment's index; per segment its sphere, level scales,
    // strides and first item (`SegmentTexels` elements)
    const aSegment = (mv as any).aSegment.ref.value as Float32Array;
    expect(Array.from(aSegment.subarray(0, 8))).toEqual([g0, 0, 0, 0, g1, 1 * g0 - 1 * g1, g0, 1]);
    const segment = (mv as any).tSegment.ref.value.array as Float32Array;
    const sphere = b.values.invariantBoundingSphere.ref.value;
    const o = SegmentTexels * 4;
    expect(Array.from(segment.subarray(o, o + 4))).toEqual(
      Array.from(Float32Array.of(...sphere.center, sphere.radius)),
    );
    // no levels: every item (stride 1); the second member's first item 3 (its high and low 16 bits)
    expect(Array.from(segment.subarray(o + 8, o + 16))).toEqual([1, 1, 1, 1, 0, 3, 0, 0]);

    // second member's positions start at texel 3 in the merged position texture
    const pg = SegmentedArray.flatten(mv.tPositionGroup.ref.value.array as unknown as SegmentedArray<Float32Array>);
    expect(items(pg.subarray(3 * 4, 3 * 4 + 4), 1)).toEqual([2, 0, 0, 0]);
    expect(items(pg.subarray(4 * 4, 4 * 4 + 4), 1)).toEqual([0, 2, 0, 1]);

    expect((mv as any).dSegmented.ref.value).toBe(true);
  });

  it("keeps its members' object colors: a texel each, read by the segment; not merged with colors per instance", () => {
    const [a, b] = createMembers();
    createObjectColor(Color(0xff0000), a.values as any);
    createObjectColor(Color(0x00ff00), b.values as any);
    expect(canMergeRenderObjects([a, b])).toBe(true);
    const merged = createMergedRenderObject([a, b]);
    const mv = merged.values as any;
    expect(mv.colorType.ref.value).toBe('object');
    expect(mv.dColorType.ref.value).toBe('texture');
    // the segment's texel (`ColorTypeIndex`)
    expect(mv.uColorIndexType.ref.value).toBe(5);
    expect(mv.tColor.ref.value.count).toBe(2);
    expect(Array.from(SegmentedArray.flatten(mv.tColor.ref.value.array).subarray(0, 8))).toEqual([
      255, 0, 0, 255, 0, 255, 0, 255,
    ]);
    // a member's color changed: its texel
    createObjectColor(Color(0x0000ff), b.values as any);
    merged.merged.sync();
    expect(Array.from(SegmentedArray.flatten(mv.tColor.ref.value.array).subarray(4, 8))).toEqual([0, 0, 255, 255]);
    // colors per instance: the same define, another granularity
    const [c] = createMembers();
    createTextureColor(createDataArray(1, 4, Uint8Array), 'instance', c.values as any);
    expect(canMergeRenderObjects([a, c])).toBe(false);
  });

  it("packs the members' second color tables by site type, a member's first entry in its segment data", () => {
    const [a, b] = createMembers();
    const table = (count: number, first: number) => {
      const array = new Uint8Array(count * 4).map((_, i) => first + i);
      return { array, count };
    };
    for (const [m, t] of [
      [a, table(2, 10)],
      [b, table(3, 50)],
    ] as const) {
      const v = m.values as any;
      ValueCell.update(v.dSiteColor, 'siteType');
      ValueCell.update(v.tSiteColor, t);
    }
    expect(canMergeRenderObjects([a, b])).toBe(true);
    const mv = createMergedRenderObject([a, b]).values as any;
    const packed = mv.tSiteColor.ref.value;
    expect(packed.count).toBe(5);
    expect(Array.from(SegmentedArray.flatten(packed.array).subarray(0, 20))).toEqual([
      ...Array.from({ length: 8 }, (_, i) => 10 + i),
      ...Array.from({ length: 12 }, (_, i) => 50 + i),
    ]);
    // the size table's first entry, then the second color's (none, and the first member's 2)
    const segment = mv.tSegment.ref.value.array as Float32Array;
    expect(segment[(0 * SegmentTexels + 1) * 4 + 11]).toBe(0);
    expect(segment[(1 * SegmentTexels + 1) * 4 + 11]).toBe(2);
    // a member without the second color does not merge with them
    const [c] = createMembers();
    expect(canMergeRenderObjects([a, c])).toBe(false);
  });

  it('sync', () => {
    const [a, b] = createMembers();
    const merged = createMergedRenderObject([a, b]);

    // member texture update flows into the merged values: the member's layer is a part, its version advanced
    const layer = merged.values.tPositionGroup;
    const segmented = layer.ref.value.array as unknown as SegmentedArray<Float32Array>;
    const before = segmented.parts.map((p) => p.version);
    const version = layer.ref.version;
    const pgB = b.values.tPositionGroup.ref.value;
    pgB.array.set([3, 0, 0, 0], 0);
    ValueCell.update(b.values.tPositionGroup, pgB);

    merged.merged.sync();
    expect(layer.ref.value.array).toBe(segmented);
    expect(layer.ref.version).toBe(version + 1);
    expect(segmented.parts.map((p) => p.version)).toEqual([before[0], before[1] + 1]);
    expect(SegmentedArray.changedRanges(segmented, { array: segmented, versions: before })).toEqual([3 * 4, 5 * 4]);
    const pg = SegmentedArray.flatten(segmented);
    expect(Array.from(pg.subarray(3 * 4, 3 * 4 + 4))).toEqual([3, 0, 0, 0]);
  });

  it('builds a merged object of records alone as of its members, their changes followed', () => {
    const [a, b] = createMembers();
    const [c, d] = createMembers();
    const byMembers = createMergedRenderObject([a, b]);
    const records = [c, d].map((m) => pickMergeRecord('spheres', m.values as unknown as MergeableValues));
    // a record its cells of what the merged object reads alone, not a member's values
    expect(Object.keys(records[0]).length).toBeLessThan(Object.keys(c.values).length / 3);
    const byRecords = createMergedRecordsRenderObject('spheres', c.values, records, c.state, -1);
    expect(byRecords.members).toEqual([]);
    expect(byRecords.merged.records).toBe(records);

    const same = () => {
      const x = byMembers.values as any;
      const y = byRecords.values as any;
      expect(byRecords.merged.segments).toEqual(byMembers.merged.segments);
      for (const k of ['drawCount', 'uVertexCount', 'instanceCount', 'uGroupCount', 'boundingSphere']) {
        expect(y[k].ref.value).toEqual(x[k].ref.value);
      }
      for (const k of ['aInstance', 'aSegment']) expect(Array.from(y[k].ref.value)).toEqual(Array.from(x[k].ref.value));
      expect(Array.from(y.tSegment.ref.value.array)).toEqual(Array.from(x.tSegment.ref.value.array));
      for (const k of ['tPositionGroup', 'tColor']) {
        const flat = (v: any) => Array.from(SegmentedArray.flatten(v[k].ref.value.array as SegmentedArray));
        expect(flat(y)).toEqual(flat(x));
      }
    };
    same();

    // a record's cell changed as a member's: read by the next sync alike
    for (const m of [b, d]) {
      const pg = m.values.tPositionGroup.ref.value;
      pg.array.set([3, 0, 0, 0], 0);
      ValueCell.update(m.values.tPositionGroup, pg);
    }
    byMembers.merged.sync();
    byRecords.merged.sync();
    same();
  });

  it("keeps its records' writes to a full read: a sync of the observed cells between a write and its update", () => {
    const [c, d] = createMembers();
    const records = [c, d].map((m) => pickMergeRecord('spheres', m.values as unknown as MergeableValues));
    const recordsVersion = ValueCell.create(0);
    const { merged, values } = createMergedRecordsRenderObject('spheres', c.values, records, c.state, -1, {
      recordsVersion,
    });
    // as the merged object's update reads it (`MergedRenderable.update`)
    const update = (observedOnly = false) => merged.sync(observedOnly || !!merged.recordsUnchanged!());
    update();
    expect(merged.recordsUnchanged!()).toBe(true);

    // a record written (its owner's version advanced), then a mark's sync of the observed cells alone
    const pg = d.values.tPositionGroup.ref.value;
    pg.array.set([3, 0, 0, 0], 0);
    ValueCell.update(d.values.tPositionGroup, pg);
    ValueCell.update(recordsVersion, 1);
    update(true);
    expect(merged.recordsUnchanged!()).toBe(false);

    // the update after: the write read
    update();
    expect(merged.recordsUnchanged!()).toBe(true);
    const flat = SegmentedArray.flatten(
      values.tPositionGroup.ref.value.array as unknown as SegmentedArray<Float32Array>,
    );
    expect(Array.from(flat.subarray(3 * 4, 3 * 4 + 3))).toEqual([3, 0, 0]);
  });

  it("a merged object moves as one: its own transient transform drawn, its members' not read", () => {
    const [a, b] = createMembers();
    const { values, merged } = createMergedRenderObject([a, b]);
    const versions = () =>
      ['aTransform', 'tSegment', 'boundingSphere'].map((k) => (values as any)[k].ref.version as number);
    const before = versions();
    const cullVersion = merged.cullVersion.ref.version;
    const placed = Sphere3D.clone(values.boundingSphere.ref.value);

    // a member's own: nothing of the merged object follows it
    ValueCell.update(b.values.uObjectTransform, Mat4.fromTranslation(Mat4(), [100, 0, 0]));
    expect(merged.sync()).toBe(false);
    expect(versions()).toEqual(before);
    expect(merged.cullVersion.ref.version).toBe(cullVersion);

    // the object's: the drawn bounds move, the instances, segments and placed bounds stay
    ValueCell.update(values.uObjectTransform, Mat4.fromTranslation(Mat4(), [0, 50, 0]));
    merged.sync();
    expect(versions()).toEqual(before);
    const drawn = getDrawnBoundingSphere(values);
    expect(drawn.center[1]).toBeCloseTo(placed.center[1] + 50, 4);
    expect(drawn.radius).toBeCloseTo(placed.radius, 4);
    expect(values.boundingSphere.ref.value.center[1]).toBeCloseTo(placed.center[1], 4);
  });

  it('sync advances the cull version where what culling reads of the members changed', () => {
    const [a, b] = createMembers();
    const { merged } = createMergedRenderObject([a, b]);
    const version = () => merged.cullVersion.ref.version;

    let v = version();
    merged.sync();
    expect(version()).toBe(v);

    // a mark: drawn data, not culled
    ValueCell.update(b.values.markerStatus, 1);
    merged.sync();
    expect(version()).toBe(v);

    // a member's sphere: the segment culled by it
    ValueCell.update(b.values.boundingSphere, { ...b.values.boundingSphere.ref.value });
    merged.sync();
    expect(version()).toBeGreaterThan(v);

    // a member's transforms: its instance grid with them
    v = version();
    ValueCell.update(a.values.aTransform, a.values.aTransform.ref.value);
    merged.sync();
    expect(version()).toBeGreaterThan(v);
  });

  it('points', () => {
    const create = (centers: number[], groups: number[]) => {
      const points = Points.create(new Float32Array(centers), new Float32Array(groups), groups.length);
      const props = PD.getDefaultValues(Points.Params);
      const values = Points.Utils.createValuesSimple(points, props, ColorNames.orange, 1);
      return createRenderObject('points', values, Points.Utils.createRenderableState(props), -1);
    };
    const a = create([0, 0, 0, 1, 0, 0, 0, 1, 0], [0, 1, 2]);
    const b = create([2, 0, 0, 0, 2, 0], [0, 1]);
    expect(canMergeRenderObjects([a, b])).toBe(true);

    const merged = createMergedRenderObject([a, b]);
    const mv = merged.values;

    // quads of 6 drawn vertices, one data vertex per point
    expect(mv.drawCount.ref.value).toBe(5 * 6);
    expect(mv.uVertexCount.ref.value).toBe(5);
    expect(merged.merged.segments.vertexBases).toEqual([0, 3]);
    expect(merged.merged.segments.drawBases).toEqual([0, 3 * 6]);

    const pg = SegmentedArray.flatten(mv.tPositionGroup.ref.value.array as unknown as SegmentedArray<Float32Array>);
    expect(items(pg.subarray(3 * 4, 3 * 4 + 4), 1)).toEqual([2, 0, 0, 0]);
    expect(items(pg.subarray(4 * 4, 4 * 4 + 4), 1)).toEqual([0, 2, 0, 1]);
  });

  it('cylinders colored per vertex', () => {
    // colored by the x of their start and end, per data vertex: a start and an end per cylinder
    const create = (xs: [number, number][]) => {
      const data = new Float32Array(xs.length * CylinderStride);
      xs.forEach(([start, end], i) => {
        data[i * CylinderStride] = start;
        data[i * CylinderStride + 3] = i; // group
        data[i * CylinderStride + 4] = end;
        data[i * CylinderStride + 7] = 1; // radius scale
      });
      const cylinders = Cylinders.create(data, xs.length);
      const props = PD.getDefaultValues(Cylinders.Params);
      const s = BaseGeometry.createSimple();
      const color = { ...s.theme.color, granularity: 'vertex' as const, color: (l: any) => Color(l.position[0]) };
      const locationIt = LocationIterator(xs.length, 1, 1, () => NullLocation);
      const values = Cylinders.Utils.createValues(cylinders, s.transform, locationIt, { ...s.theme, color }, props);
      return createRenderObject('cylinders', values, Cylinders.Utils.createRenderableState(props), -1);
    };
    const a = create([
      [1, 2],
      [3, 4],
    ]);
    const b = create([[5, 6]]);
    expect(a.values.colorType.ref.value).toBe('vertex');
    expect(canMergeRenderObjects([a, b])).toBe(true);

    const merged = createMergedRenderObject([a, b]);
    const mv = merged.values;
    expect(mv.uVertexCount.ref.value).toBe(3 * 2);
    expect(merged.merged.segments.vertexBases).toEqual([0, 4]);
    // the colors per data vertex, the second member's at its vertex base (the blue byte is the x)
    const colors = SegmentedArray.flatten(mv.tColor.ref.value.array as unknown as SegmentedArray<Uint8Array>);
    expect(Array.from({ length: 6 }, (_, i) => colors[i * 4 + 2])).toEqual([1, 2, 3, 4, 5, 6]);
    // three data texels per cylinder
    const data = SegmentedArray.flatten(mv.tCylinder.ref.value.array as unknown as SegmentedArray<Float32Array>);
    expect(data.length).toBe(3 * 3 * 4);
    expect(data[2 * CylinderStride + 4]).toBe(6);
  });

  it('mesh', () => {
    const a = createMeshRenderObject(6, 0);
    const b = createMeshRenderObject(9, 100);
    expect(canMergeRenderObjects([a, b])).toBe(true);

    const merged = createMergedRenderObject([a, b]);
    const mv = merged.values;

    expect(mv.drawCount.ref.value).toBe(6 + 9);
    expect(mv.uVertexCount.ref.value).toBe(6 + 9);
    expect(mv.instanceCount.ref.value).toBe(2);

    const { segments } = merged.merged;
    expect(segments.vertexBases).toEqual([0, 6]);
    expect(segments.elementOffsets).toEqual([0, 6 * 4]);

    // the members' arrays as parts, not copied: the second member's indices biased to absolute vertex ids
    const segmentedElements = mv.elements.ref.value as unknown as SegmentedArray<Uint32Array>;
    expect(segmentedElements.parts[1].array).toBe(b.values.elements.ref.value);
    const elements = SegmentedArray.flatten(segmentedElements);
    expect(Array.from(elements.subarray(0, 6))).toEqual([0, 1, 2, 3, 4, 5]);
    expect(Array.from(elements.subarray(6, 15))).toEqual([6, 7, 8, 9, 10, 11, 12, 13, 14]);

    // vertex positions at the members' vertex bases
    const segmentedPositions = mv.aPosition.ref.value as unknown as SegmentedArray<Float32Array>;
    expect(segmentedPositions.parts[0].array).toBe(a.values.aPosition.ref.value);
    const positions = SegmentedArray.flatten(segmentedPositions);
    expect(positions[0]).toBe(0);
    expect(positions[6 * 3]).toBe(100);

    // a member's change: its part only, read from the member again
    const before = segmentedPositions.parts.map((p) => p.version);
    const moved = b.values.aPosition.ref.value.slice();
    moved[0] = 200;
    ValueCell.update(b.values.aPosition, moved);
    merged.merged.sync();
    expect(mv.aPosition.ref.value).toBe(segmentedPositions);
    expect(segmentedPositions.parts.map((p) => p.version)).toEqual([before[0], before[1] + 1]);
    expect(SegmentedArray.flatten(segmentedPositions)[6 * 3]).toBe(200);
  });

  it('syncs the members whose observed cells changed alone, as a full sync does', () => {
    const create = () => {
      const members = Array.from({ length: 6 }, (_, i) =>
        createSpheresRenderObject(new Float32Array([i, 0, 0, i, 1, 0]), new Float32Array([0, 1]), 2),
      );
      const records = members.map((m) => pickMergeRecord('spheres', m.values as unknown as MergeableValues));
      return createMergedRecordsRenderObject('spheres', members[0].values, records, members[0].state, -1);
    };
    const observed = create();
    const full = create();
    const markers = (o: ReturnType<typeof create>) => {
      const v = o.values as any;
      const array = v.tMarker.ref.value.array;
      return {
        layer: Array.from(SegmentedArray.is(array) ? SegmentedArray.flatten(array) : array),
        average: v.markerAverage.ref.value,
        status: v.markerStatus.ref.value,
        marked: o.merged.getMarkedMembers!(),
      };
    };
    // as `Visual.mark` writes a member's markers: its layer, its average and status; a uniform status's layer not
    // marked changed (its uniform draws it)
    const mark = (i: number, layer: number[], average: number, status: number) => {
      for (const o of [observed, full]) {
        const r = o.merged.records[i] as MergeRecord;
        const marker = r.tMarker!.ref.value as { array: Uint8Array };
        marker.array.set(layer);
        if (status === -1) ValueCell.update(r.tMarker!, marker);
        ValueCell.update(r.markerAverage, average);
        ValueCell.update(r.markerStatus, status);
      }
      observed.merged.sync(true);
      full.merged.sync();
      expect(markers(observed)).toEqual(markers(full));
    };

    mark(2, [1, 1], 1, 1);
    mark(4, [1, 0], 0.5, -1);
    expect(markers(observed).marked).toEqual([2, 4]);
    expect(markers(observed).status).toBe(-1);
    mark(2, [0, 0], 0, 0);
    expect(markers(observed).marked).toEqual([4]);
    // none marked: the average 0, whatever the updates' rounding
    mark(4, [0, 0], 0, 0);
    expect(markers(observed)).toMatchObject({ average: 0, status: 0, marked: [] });
  });

  const ctx = tryGetGLContext(32, 32, { fragDepth: true });

  (ctx ? it : it.skip)(
    'updates a merged object by its observed cells alone where only they changed',
    { tags: ['medium'] },
    () => {
      const ctx = getGLContext(32, 32);
      const scene = Scene.create(ctx);
      const [a, b] = createMembers();
      const merged = createMergedRenderObject([a, b]);
      scene.add(merged);
      scene.commit();
      scene.updateChanged(false);
      const renderable = scene.get(merged)!;
      const update = vi.spyOn(renderable, 'update');

      // a mark: the members' observed cells alone
      ValueCell.update(b.values.markerStatus, 1);
      scene.updateChanged(false);
      expect(update).toHaveBeenCalledExactlyOnceWith(true);
      expect(merged.values.markerStatus.ref.value).toBe(-1);

      // a representation's update: all of it read
      update.mockClear();
      scene.update([merged], false);
      expect(update).toHaveBeenCalledExactlyOnceWith();
      ctx.destroy();
    },
  );

  (ctx ? it : it.skip)('scene', { tags: ['medium'] }, async () => {
    const ctx = getGLContext(32, 32);
    const scene = Scene.create(ctx);
    const [a, b] = createMembers();
    const merged = createMergedRenderObject([a, b]);
    scene.add(merged);
    const mergedMesh = createMergedRenderObject([createMeshRenderObject(6, 0), createMeshRenderObject(9, 100)]);
    scene.add(mergedMesh);
    setDebugMode(true);
    expect(() => scene.commit()).not.toThrow();
    expect(() => scene.update(undefined, false)).not.toThrow();
    setDebugMode(false);
    ctx.destroy();
  });
});
