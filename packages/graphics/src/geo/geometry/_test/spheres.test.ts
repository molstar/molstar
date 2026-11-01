/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Spheres } from '../spheres/spheres.js';
import { Points } from '../points/points.js';

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
  strides.map((stride) => ({ minDistance: 0, maxDistance: 0, overlap: 0, stride, scaleBias: 1 }));

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
    expect(Array.from(values.tPositionGroup.ref.value.array.subarray(0, 12))).toEqual([
      0, 0, 0, 0, 1, 2, 3, 1, 4, 5, 6, 2,
    ]);
  });

  it('order into levels of detail by the strides only, not by the levels set before', () => {
    const n = 23;
    const direct = line(n);
    direct.shaderData.update({ lodLevels: levels(1, 3), sizeFactor: 1 });

    const changed = line(n);
    changed.shaderData.update({ lodLevels: levels(1, 2, 5), sizeFactor: 1 });
    changed.shaderData.update({ lodLevels: levels(1, 4), sizeFactor: 1 });
    changed.shaderData.update({ lodLevels: levels(1, 3), sizeFactor: 1 });
    expect(groupsOf(changed)).toEqual(groupsOf(direct));

    // coarsest level first: every 3rd item, then the rest
    const coarse = Array.from({ length: Math.ceil(n / 3) }, (_, i) => i * 3);
    expect(groupsOf(direct)).toEqual([...coarse, ...groupsOf(line(n)).filter((g) => g % 3 !== 0)]);
    // the center stays with its group
    const pg = direct.positionGroup.ref.value.array;
    for (let i = 0; i < n; i++)
      expect(Array.from(pg.subarray(i * 4, i * 4 + 3))).toEqual([1, 2, 3].map((k) => k * pg[i * 4 + 3]));

    // without levels: back to the items' order
    direct.shaderData.update({ lodLevels: [], sizeFactor: 1 });
    expect(groupsOf(direct)).toEqual(groupsOf(line(n)));
  });

  it('give the levels the offsets of their vertices', () => {
    const spheres = line(10);
    spheres.shaderData.update({ lodLevels: levels(1, 4), sizeFactor: 1 });
    // level of stride 4: items 0, 4, 8; then the 7 others; 6 vertices each
    const counts = (spheres.shaderData.lodLevels.ref.value as unknown as number[][]).map((l) => l[3]);
    expect(counts).toEqual([60, 18]);
  });

  it('keep the items no level takes, after the levels', () => {
    const spheres = line(7);
    spheres.shaderData.update({ lodLevels: levels(2, 3), sizeFactor: 1 });
    // stride 3: 0, 3, 6; stride 2: 2, 4; no level: 1, 5
    expect(groupsOf(spheres)).toEqual([0, 3, 6, 2, 4, 1, 5]);
    spheres.shaderData.update({ lodLevels: [], sizeFactor: 1 });
    expect(groupsOf(spheres)).toEqual([0, 1, 2, 3, 4, 5, 6]);
  });

  it('map groups and bound the spheres in the current order', () => {
    const spheres = line(9);
    const sphere = { ...spheres.boundingSphere, center: [...spheres.boundingSphere.center] };
    spheres.shaderData.update({ lodLevels: levels(1, 3), sizeFactor: 1 });
    const { indices, offsets } = spheres.groupMapping;
    const pg = spheres.positionGroup.ref.value.array;
    for (let g = 0; g < 9; g++) {
      expect(offsets[g + 1] - offsets[g]).toBe(1);
      // an offset into the position groups, at the group's center
      expect(pg[indices[offsets[g]] + 3]).toBe(g);
      expect(pg[indices[offsets[g]]]).toBe(g);
    }
    expect(spheres.boundingSphere.radius).toBeCloseTo(sphere.radius, 5);

    // computed after a reorder: as in the items' order (the same sphere, built with levels or without)
    const reordered = line(9);
    reordered.shaderData.update({ lodLevels: levels(1, 3), sizeFactor: 1 });
    expect(reordered.boundingSphere.radius).toBe(sphere.radius);
    expect(Array.from(reordered.boundingSphere.center)).toEqual(sphere.center);
  });
});

describe('Points', () => {
  it('keep center and group in one array, transformed in place', () => {
    const points = Points.create(Float32Array.of(0, 0, 0, 1, 2, 3), Float32Array.of(5, 6), 2);
    Points.transform(points, Float32Array.of(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 10, 20, 30, 1) as any);
    expect(Array.from(points.positionGroup.ref.value.array.subarray(0, 8))).toEqual([10, 20, 30, 5, 11, 22, 33, 6]);
  });

  it('order into levels of detail as spheres, by the strides as given (not by the size factor)', () => {
    const n = 10;
    const centers = new Float32Array(n * 3);
    const groups = Float32Array.from({ length: n }, (_, i) => i);
    const points = Points.create(centers, groups, n);
    const groupsOfPoints = () => Array.from({ length: n }, (_, i) => points.positionGroup.ref.value.array[i * 4 + 3]);

    const props = { ...PD.getDefaultValues(Points.Params), sizeFactor: 3, lodLevels: levels(1, 4) };
    const values = Points.Utils.createValuesSimple(points, props, ColorNames.orange, 1);
    expect(groupsOfPoints()).toEqual([0, 4, 8, 1, 2, 3, 5, 6, 7, 9]);
    // vertex counts of the levels: 3 points of stride 4, then the other 7
    expect((values.lodLevels.ref.value as unknown as number[][]).map((l) => l[3])).toEqual([60, 18]);

    // levels changed by props: ordered again, back to the points' order without
    Points.Utils.updateValues(values, { ...props, lodLevels: [] });
    expect(groupsOfPoints()).toEqual(
      groups
        .map((g) => g)
        .values()
        .toArray(),
    );
    expect(values.lodLevels.ref.value).toEqual([]);
  });
});
