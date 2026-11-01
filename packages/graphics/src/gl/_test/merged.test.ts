/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { createRenderObject, createMergedRenderObject, canMergeRenderObjects } from '../render-object.js';
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
import { NullLocation } from '@molstar/model/model/location';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';
import { Cylinders, CylinderStride } from '@molstar/graphics/geo/geometry/cylinders/cylinders';
import { LocationIterator } from '@molstar/graphics/geo/util/location-iterator';

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
    // per instance its segment's group count and offsets and the segment's index; per segment its sphere, level scales
    const aSegment = (mv as any).aSegment.ref.value as Float32Array;
    expect(Array.from(aSegment.subarray(0, 8))).toEqual([g0, 0, 0, 0, g1, 1 * g0 - 1 * g1, g0, 1]);
    const segment = (mv as any).tSegment.ref.value.array as Float32Array;
    const sphere = b.values.invariantBoundingSphere.ref.value;
    expect(Array.from(segment.subarray(8, 12))).toEqual(Array.from(Float32Array.of(...sphere.center, sphere.radius)));

    // second member's positions start at texel 3 in the merged position texture
    const pg = SegmentedArray.flatten(mv.tPositionGroup.ref.value.array as unknown as SegmentedArray<Float32Array>);
    expect(Array.from(pg.subarray(3 * 4, 3 * 4 + 4))).toEqual([2, 0, 0, 0]);
    expect(Array.from(pg.subarray(4 * 4, 4 * 4 + 4))).toEqual([0, 2, 0, 1]);

    expect((mv as any).dSegmented.ref.value).toBe(true);
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
    expect(Array.from(pg.subarray(3 * 4, 3 * 4 + 4))).toEqual([2, 0, 0, 0]);
    expect(Array.from(pg.subarray(4 * 4, 4 * 4 + 4))).toEqual([0, 2, 0, 1]);
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

  const ctx = tryGetGLContext(32, 32, { fragDepth: true });

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
