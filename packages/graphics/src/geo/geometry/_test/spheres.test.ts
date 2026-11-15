/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { unpackRGBToInt } from '@molstar/core/util/number-packing';
import { BaseGeometry } from '../base.js';
import { Geometry } from '../geometry.js';
import { sizeDataFactor } from '../size-data.js';
import { createTransform } from '../transform-data.js';
import { describe, expect, it } from 'vitest';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Spheres } from '../spheres/spheres.js';
import { Points } from '../points/points.js';
import { getLodPixelLevels, getNestedStrides } from '../lod-levels.js';
import { getDataWords } from '../center-group.js';

/** `count` items of a position group's array: xyz and the group (an int32 as its bits) */
function items(array: Float32Array, count: number) {
  const words = getDataWords(array);
  return Array.from({ length: count * 4 }, (_, i) => (i % 4 === 3 ? words[i] : array[i]));
}

/** `n` spheres, sphere `i` at (i, 2i, 3i) of group `i` */
function line(n: number) {
  const centers = new Float32Array(n * 3);
  const groups = new Float32Array(n);
  for (let i = 0; i < n; i++) {
    centers.set([i, 2 * i, 3 * i], i * 3);
    groups[i] = i;
  }
  return Spheres.create(centers, groups, n);
}

const groupsOf = (s: Spheres) =>
  Array.from({ length: s.sphereCount }, (_, i) => s.positionGroup.ref.value.array[i * 4 + 3]);

const levels = (...strides: number[]) =>
  strides.map((stride) => ({ minDistance: 0, maxDistance: 0, stride, scaleBias: 1 }));

describe('Spheres', () => {
  it('keep center and group in one array, which createValues uses as it is (perf-audit-2 08/P11)', () => {
    const spheres = Spheres.create(Float32Array.of(0, 0, 0, 1, 2, 3, 4, 5, 6), Float32Array.of(0, 1, 2), 3);
    const { positionGroup } = spheres;
    expect(spheres.shaderData.positionGroup).toBe(positionGroup);
    const values = Spheres.Utils.createValuesSimple(spheres, PD.getDefaultValues(Spheres.Params), ColorNames.orange, 1);
    // no levels: not written again (no copy, the geometry's own array)
    expect(values.tPositionGroup).toBe(positionGroup);
    expect(positionGroup.ref.version).toBe(0);
    // xyz and group per sphere
    expect(items(values.tPositionGroup.ref.value.array, 3)).toEqual([0, 0, 0, 0, 1, 2, 3, 1, 4, 5, 6, 2]);
  });

  it('keep the items in their order with levels of detail: a level every stride-th', () => {
    const n = 23;
    const spheres = line(n);
    const version = spheres.positionGroup.ref.version;
    spheres.shaderData.update({ lodPixelSize: 3, sizeFactor: 1 });
    spheres.shaderData.update({ lodPixelSize: 3, sizeFactor: 0.5 });
    spheres.shaderData.update({ lodPixelSize: 0, sizeFactor: 1 });
    // not written: no reorder, the items as built
    expect(spheres.positionGroup.ref.version).toBe(version);
    expect(groupsOf(spheres)).toEqual(groupsOf(line(n)));
  });

  it('give the levels the quads they draw and their strides', () => {
    const spheres = line(25);
    const levelsOf = () => spheres.shaderData.lodLevels.ref.value as unknown as number[][];
    spheres.shaderData.update({ lodPixelSize: 3, sizeFactor: 1 });
    // finest first: every item, items 0, 10, 20 (stride 10), item 0 (strides 100, 1000); 6 vertices each
    expect(levelsOf().map((l) => l[2])).toEqual([150, 18, 6, 6]);
    expect(levelsOf().map((l) => l[4])).toEqual([1, 10, 100, 1000]);
    // at size factor 0.5: every 2nd item at the finest (13 of 25), item 0 at the others
    spheres.shaderData.update({ lodPixelSize: 3, sizeFactor: 0.5 });
    expect(levelsOf().map((l) => l[2])).toEqual([78, 6, 6, 6]);
    expect(levelsOf().map((l) => l[4])).toEqual([2, 80, 640, 5120]);
  });

  it('nest the strides: a coarser one rounded to a multiple of the finer', () => {
    expect(getNestedStrides(levels(2, 3), 1)).toEqual([2, 4]);
    expect(getNestedStrides(levels(1, 40, 150, 300), 1)).toEqual([1, 40, 160, 320]);
    // the size factor's adjustment nested too (1, 10, 100, 1000 at 1.2: 6, 61, 657 -> 6, 60, 660)
    const chain = getLodPixelLevels(3);
    expect(getNestedStrides(chain, 1.2)).toEqual([1, 6, 60, 660]);
    expect(getNestedStrides(chain, 1)).toEqual([1, 10, 100, 1000]);
    expect(getNestedStrides(chain, 0.5)).toEqual([2, 80, 640, 5120]);
  });

  it('map groups and bound the spheres, with levels as without', () => {
    const spheres = line(9);
    const sphere = { ...spheres.boundingSphere, center: [...spheres.boundingSphere.center] };
    spheres.shaderData.update({ lodPixelSize: 3, sizeFactor: 0.5 });
    const { indices, offsets } = spheres.groupMapping;
    const pg = spheres.positionGroup.ref.value.array;
    const words = getDataWords(pg);
    for (let g = 0; g < 9; g++) {
      expect(offsets[g + 1] - offsets[g]).toBe(1);
      // an offset into the position groups, at the group's center
      expect(words[indices[offsets[g]] + 3]).toBe(g);
      expect(pg[indices[offsets[g]]]).toBe(g);
    }
    expect(spheres.boundingSphere.radius).toBeCloseTo(sphere.radius, 5);

    // the same sphere, built with levels or without
    const reordered = line(9);
    reordered.shaderData.update({ lodPixelSize: 3, sizeFactor: 0.5 });
    expect(reordered.boundingSphere.radius).toBe(sphere.radius);
    expect(Array.from(reordered.boundingSphere.center)).toEqual(sphere.center);
  });
});

describe('Points', () => {
  it('keep center and group in one array, transformed in place', () => {
    const points = Points.create(Float32Array.of(0, 0, 0, 1, 2, 3), Float32Array.of(5, 6), 2);
    Points.transform(points, Float32Array.of(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 10, 20, 30, 1) as any);
    expect(items(points.positionGroup.ref.value.array, 2)).toEqual([10, 20, 30, 5, 11, 22, 33, 6]);
  });

  it('have levels of detail as spheres, by the strides of the chain (not adjusted by the size factor)', () => {
    const n = 25;
    const centers = new Float32Array(n * 3);
    const groups = Float32Array.from({ length: n }, (_, i) => i);
    const points = Points.create(centers, groups, n);
    const groupsOfPoints = () =>
      Array.from({ length: n }, (_, i) => getDataWords(points.positionGroup.ref.value.array)[i * 4 + 3]);

    const props = { ...PD.getDefaultValues(Points.Params), sizeFactor: 3, lodPixelSize: 3 };
    const values = Points.Utils.createValuesSimple(points, props, ColorNames.orange, 1);
    // the points in their order
    expect(groupsOfPoints()).toEqual(Array.from(groups));
    // vertex counts of the levels, finest first: every point, 3 of stride 10, point 0 of strides 100 and 1000
    expect((values.lodLevels.ref.value as unknown as number[][]).map((l) => l[2])).toEqual([150, 18, 6, 6]);
    expect((values.lodLevels.ref.value as unknown as number[][]).map((l) => l[4])).toEqual([1, 10, 100, 1000]);

    Points.Utils.updateValues(values, { ...props, lodPixelSize: 0 });
    expect(groupsOfPoints()).toEqual(Array.from(groups));
    expect(values.lodLevels.ref.value).toEqual([]);
  });
  // what is computed per vertex is read by item (a level draws every stride-th)
  it("colors and sizes per vertex in the items' order, with levels of detail", () => {
    const count = 20;
    const centers = new Float32Array(count * 3);
    for (let i = 0; i < count; ++i) centers[i * 3] = i;
    const groups = new Float32Array(Array.from({ length: count }, (_, i) => i));
    const spheres = Spheres.create(centers, groups, count);
    const props = { ...PD.getDefaultValues(Spheres.Params), lodPixelSize: 3 };
    const transform = createTransform(Mat4.identity() as unknown as Float32Array, 1, Sphere3D.create(Vec3(), 30), 0, 0);
    const simple = BaseGeometry.createSimple(0 as never, 1, transform);
    const byX = (l: { position: Vec3 }) => l.position[0];
    const theme = {
      ...simple.theme,
      color: { ...simple.theme.color, granularity: 'vertex', color: (l: unknown) => byX(l as { position: Vec3 }) },
      size: { ...simple.theme.size, granularity: 'vertex', size: (l: unknown) => 1 + byX(l as { position: Vec3 }) },
    } as never;
    const values = Spheres.Utils.createValues(spheres, simple.transform, simple.locationIterator, theme, props);
    const positions = values.tPositionGroup.ref.value.array;
    const colors = values.tColor.ref.value.array;
    const sizes = values.tSize.ref.value.array;
    expect((values.lodLevels.ref.value as unknown[]).length).toBe(4);
    for (let i = 0; i < count; ++i) {
      expect(positions[i * 4]).toBe(i);
      expect(colors[i * 4 + 2]).toBe(i);
      expect(unpackRGBToInt(sizes[i * 4], sizes[i * 4 + 1], sizes[i * 4 + 2]) / sizeDataFactor).toBeCloseTo(1 + i, 3);
    }
  });
});

describe('Spheres from quantized sites', () => {
  // centers multiples of 0.25 from the origin: exact quantized and as floats
  const origin = Vec3.create(1, -2, 0.5);
  const scale = 0.25;
  const n = 7;
  const data = new Int16Array(n * 4);
  for (let i = 0; i < n; ++i) data.set([i * 4, -i * 2, i * i, 3], i * 4);
  const centerOf = (i: number) => [origin[0] + i, origin[1] - i * 0.5, origin[2] + i * i * 0.25];

  it('read the centers from the array as it is, a sphere a group, as their float twin', () => {
    const sites = Spheres.fromSites({ data, origin, scale }, n);
    const floats = Spheres.create(
      Float32Array.from({ length: n * 3 }, (_, k) => centerOf(Math.floor(k / 3))[k % 3]),
      Float32Array.from({ length: n }, (_, i) => i),
      n,
    );
    expect(sites.sites!.data).toBe(data);
    expect(sites.boundingSphere.radius).toBeCloseTo(floats.boundingSphere.radius, 5);
    expect(Array.from(sites.boundingSphere.center)).toEqual(Array.from(floats.boundingSphere.center));
    expect(Geometry.getGroupCount(sites)).toBe(n);
    const { indices, offsets } = sites.groupMapping;
    for (let g = 0; g < n; ++g) expect(indices[offsets[g]]).toBe(g * 4);

    const values = Spheres.Utils.createValuesSimple(sites, {}, ColorNames.orange, 1);
    expect(values.dSites.ref.value).toBe(true);
    // the view itself, no copy
    expect(values.tSite.ref.value.array).toBe(data);
    expect(Array.from(values.uSiteOrigin.ref.value)).toEqual(Array.from(origin));
    expect(values.uSiteScale.ref.value).toBe(scale);
    const pg = Spheres.getPositionGroupArray(values);
    for (let i = 0; i < n; ++i) expect(items(pg.subarray(i * 4, i * 4 + 4), 1)).toEqual([...centerOf(i), i]);

    // the position iterator: the decoded centers
    const transform = createTransform(Mat4.identity() as unknown as Float32Array, 1, Sphere3D.create(Vec3(), 30), 0, 0);
    const it = Spheres.Utils.createPositionIterator(sites, transform);
    const p = it.move();
    expect(Array.from((p.location as { position: Vec3 }).position)).toEqual(centerOf(0));
  });

  it('take other sites in place, float spheres not', () => {
    const sites = Spheres.fromSites({ data, origin, scale }, n);
    const other = Spheres.fromSites({ data: data.subarray(0, 8), origin: Vec3.create(0, 0, 0), scale: 0.5 }, 2, sites);
    expect(other).toBe(sites);
    expect(other.sphereCount).toBe(2);
    expect(other.shaderData.uSiteScale.ref.value).toBe(0.5);
    // a float geometry from a quantized one: a new geometry, not the sites' placeholder
    const floats = Spheres.fromPositionGroup(Float32Array.of(1, 2, 3, 0), 1, sites);
    expect(floats).not.toBe(sites);
    expect(floats.sites).toBeUndefined();
    expect(floats.shaderData.dSites.ref.value).toBe(false);
  });
});
