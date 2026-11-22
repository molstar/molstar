/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Mat3, Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { GaussianSplats } from '../gaussian-splats/gaussian-splats.js';
import {
  getMaxSplatScale,
  getSplatCovariance,
  getSplatOpacity,
  getSplatRotation,
  getSplatScale,
  setSplat,
  SplatLogScaleFactor,
  SplatStride,
} from '../gaussian-splats/splat-data.js';
import { int32ToFloatBits, PositionGroupStride } from '../center-group.js';
import { calculateInvariantBoundingSphere } from '@molstar/graphics/gl/renderable/util';
import { Sphere3D } from '@molstar/core/math/geometry';
import { ShapeGroup } from '@molstar/model/model/shape/shape';
import { OrderedSet } from '@molstar/core/data/int';
import { Color } from '@molstar/core/util/color';
import * as Shape from '@molstar/graphics/geo/shape/shape';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { Texture } from '@molstar/graphics/gl/texture';

function splatData(...splats: [Vec3, Quat, number][]) {
  const data = new Int16Array(splats.length * SplatStride);
  splats.forEach(([s, q, o], i) =>
    setSplat(data, i, Math.log(s[0]), Math.log(s[1]), Math.log(s[2]), q[0], q[1], q[2], q[3], o),
  );
  return data;
}

function expectCovariance(c: Mat3, expected: number[], digits = 3) {
  for (let i = 0; i < 9; ++i) expect(c[i]).toBeCloseTo(expected[i], digits);
}

describe('splat data', () => {
  it('keep scales in relative steps of 0.1%, the rotation and the opacity', () => {
    const q = Quat.normalize(Quat(), Quat.create(0.1, -0.7, 0.3, 0.5));
    const scales = [1e-6, 0.37, 1, 42, 1e5];
    const data = new Int16Array(scales.length * SplatStride);
    scales.forEach((s, i) =>
      setSplat(data, i, Math.log(s), Math.log(2 * s), Math.log(s / 3), q[0], q[1], q[2], q[3], i / 4),
    );
    const step = Math.exp(0.5 / SplatLogScaleFactor) - 1;
    for (let i = 0; i < scales.length; ++i) {
      const s = getSplatScale(Vec3(), data, i);
      expect(Math.abs(s[0] / scales[i] - 1)).toBeLessThanOrEqual(step);
      expect(Math.abs(s[1] / (2 * scales[i]) - 1)).toBeLessThanOrEqual(step);
      expect(Math.abs(s[2] / (scales[i] / 3) - 1)).toBeLessThanOrEqual(step);
      const r = getSplatRotation(Quat(), data, i);
      for (let k = 0; k < 4; ++k) expect(r[k]).toBeCloseTo(q[k], 4);
      expect(getSplatOpacity(data, i)).toBeCloseTo(i / 4, 4);
    }
  });

  it('normalize rotations, w non-negative (-q the same rotation), zero the identity', () => {
    const data = new Int16Array(3 * SplatStride);
    setSplat(data, 0, 0, 0, 0, 0, 0, 2, -2, 1);
    setSplat(data, 1, 0, 0, 0, 0, 0, 0, 0, 1);
    setSplat(data, 2, 0, 0, 0, 0, 0, 0, -1, 1);
    const r = Quat();
    expect(Array.from(getSplatRotation(r, data, 0)).map((v) => +v.toFixed(4))).toEqual([0, 0, -0.7071, 0.7071]);
    expect(Array.from(getSplatRotation(r, data, 1))).toEqual([0, 0, 0, 1]);
    // w = 1 the largest uint16, stored as its bits
    expect(Array.from(getSplatRotation(r, data, 2))).toEqual([0, 0, 0, 1]);
  });

  it('clamp scales out of range and opacities out of 0 to 1', () => {
    const data = new Int16Array(SplatStride);
    setSplat(data, 0, Math.log(0), 1000, Number.NaN, 0, 0, 0, 1, Number.NaN);
    const s = getSplatScale(Vec3(), data, 0);
    expect(s[0]).toBeCloseTo(Math.exp(-32767 / SplatLogScaleFactor), 20);
    expect(s[1]).toBeCloseTo(Math.exp(32767 / SplatLogScaleFactor), -10);
    expect(s[2]).toBe(s[0]);
    expect(getSplatOpacity(data, 0)).toBe(0);
    setSplat(data, 0, 0, 0, 0, 0, 0, 0, 1, 2);
    expect(getSplatOpacity(data, 0)).toBe(1);
  });

  it('build the covariance R S Sᵀ Rᵀ', () => {
    const c = Mat3();
    const identity = Quat.identity();
    expectCovariance(
      getSplatCovariance(c, splatData([Vec3.create(1, 1, 1), identity, 1]), 0),
      [1, 0, 0, 0, 1, 0, 0, 0, 1],
    );
    expectCovariance(
      getSplatCovariance(c, splatData([Vec3.create(2, 3, 4), identity, 1]), 0),
      [4, 0, 0, 0, 9, 0, 0, 0, 16],
      1,
    );
    // 90° about z: the x and y variances swapped, the trace kept, symmetric
    const z90 = Quat.setAxisAngle(Quat(), Vec3.create(0, 0, 1), Math.PI / 2);
    getSplatCovariance(c, splatData([Vec3.create(2, 3, 4), z90, 1]), 0);
    expectCovariance(c, [9, 0, 0, 0, 4, 0, 0, 0, 16], 1);
    const tilted = Quat.normalize(Quat(), Quat.create(0.3, 0.2, -0.4, 0.8));
    getSplatCovariance(c, splatData([Vec3.create(2, 3, 4), tilted, 1]), 0);
    expect(c[0] + c[4] + c[8]).toBeCloseTo(29, 1);
    expect(c[1]).toBeCloseTo(c[3], 6);
    expect(c[2]).toBeCloseTo(c[6], 6);
    expect(c[5]).toBeCloseTo(c[7], 6);
  });

  it('find the largest scale', () => {
    const identity = Quat.identity();
    const data = splatData([Vec3.create(1, 0.5, 2), identity, 1], [Vec3.create(0.1, 7, 3), identity, 1]);
    expect(getMaxSplatScale(data, 2)).toBeCloseTo(7, 2);
    expect(getMaxSplatScale(data, 1)).toBeCloseTo(2, 2);
    expect(getMaxSplatScale(data, 0)).toBe(0);
  });
});

describe('gaussian splats', () => {
  function create(count: number, scale = 1) {
    const { positionGroup, splatData } = GaussianSplats.createArrays(count);
    for (let i = 0; i < count; ++i) {
      positionGroup.set([i * 10, 0, 0, int32ToFloatBits(i)], i * 4);
      setSplat(splatData, i, Math.log(scale), Math.log(scale), Math.log(scale), 0, 0, 0, 1, 1);
    }
    return { positionGroup, splatData };
  }

  it('keep their arrays as the shaders read them, two rgba16i elements a splat', () => {
    const { positionGroup, splatData } = create(3);
    const splats = GaussianSplats.create(positionGroup, splatData, 3);
    expect(splats.splatCount).toBe(3);
    expect(splats.positionGroup.ref.value).toMatchObject({ array: positionGroup, count: 3 });
    expect(splats.splatData.ref.value).toMatchObject({ array: splatData, count: 6 });
    // splat i of group i, at offset 4i in the position groups
    expect(Array.from(splats.groupMapping.indices)).toEqual([0, 4, 8]);
  });

  it('pad the bounds by three standard deviations of the largest scale, again where the data changes', () => {
    const { positionGroup, splatData } = create(3, 0.5);
    const splats = GaussianSplats.create(positionGroup, splatData, 3);
    // the centers' sphere, as points and spheres bound theirs
    const centers = calculateInvariantBoundingSphere(positionGroup, 3, 1, PositionGroupStride).radius;
    expect(centers).toBeGreaterThanOrEqual(10);
    expect(splats.boundingSphere.radius).toBeCloseTo(centers + GaussianSplats.BoundsSigmas * 0.5, 2);
    expect(splats.hasBoundingSphere()).toBe(true);

    const larger = create(3, 2);
    GaussianSplats.create(larger.positionGroup, larger.splatData, 3, splats);
    expect(splats.hasBoundingSphere()).toBe(false);
    expect(splats.boundingSphere.radius).toBeCloseTo(centers + GaussianSplats.BoundsSigmas * 2, 2);
  });

  it('update in place, empty without splats', () => {
    const { positionGroup, splatData } = create(3);
    const splats = GaussianSplats.create(positionGroup, splatData, 3);
    const version = splats.splatData.ref.version;
    expect(GaussianSplats.createEmpty(splats)).toBe(splats);
    expect(splats.splatCount).toBe(0);
    expect(splats.splatData.ref.value.count).toBe(0);
    expect(splats.splatData.ref.version).toBeGreaterThan(version);
    expect(GaussianSplats.createEmpty().splatCount).toBe(0);
  });

  it('reject shape data shorter than the splats', () => {
    const { positionGroup } = create(2);
    expect(() => GaussianSplats.create(positionGroup, new Int16Array(SplatStride), 2)).toThrow(/splat data/);
  });

  it('bound a group of splats by its centers, padded by three standard deviations of its largest scale', () => {
    // two splats far apart: 0 at x 0 of scale 0.5, 1 at x 100 of scale 2
    const { positionGroup, splatData } = GaussianSplats.createArrays(2);
    positionGroup.set([0, 0, 0, int32ToFloatBits(0), 100, 0, 0, int32ToFloatBits(1)]);
    setSplat(splatData, 0, Math.log(0.5), Math.log(0.5), Math.log(0.5), 0, 0, 0, 1, 1);
    setSplat(splatData, 1, Math.log(2), Math.log(1), Math.log(1), 0, 0, 0, 1, 1);
    const shape = Shape.create(
      'splats',
      {},
      GaussianSplats.create(positionGroup, splatData, 2),
      () => Color(0),
      () => 1,
      () => 'splat',
    );
    const group0 = Shape.getGroupBoundingSphere(
      ShapeGroup.Loci(shape, [{ ids: OrderedSet.ofSingleton(0), instance: 0 }]),
    );
    expect(Array.from(group0.center)).toEqual([0, 0, 0]);
    expect(group0.radius).toBeCloseTo(GaussianSplats.BoundsSigmas * 0.5, 2);
    const group1 = Shape.getGroupBoundingSphere(
      ShapeGroup.Loci(shape, [{ ids: OrderedSet.ofSingleton(1), instance: 0 }]),
    );
    expect(Array.from(group1.center)).toEqual([100, 0, 0]);
    expect(group1.radius).toBeCloseTo(GaussianSplats.BoundsSigmas * 2, 2);
  });

  // SOG read in place: the textures the geometry's (destroyed when replaced), the theme the shader's colors
  it('of SOG images: the sog encoding, its textures owned, the colors direct, the bounds by the codebook', () => {
    const destroyed: string[] = [];
    const texture = (name: string, width = 2, height = 2) =>
      ({ getWidth: () => width, getHeight: () => height, destroy: () => destroyed.push(name) }) as unknown as Texture;
    const graphics = {} as GraphicsContext;
    const other = {} as GraphicsContext;
    const input = (tag: string): GaussianSplats.SogInput => ({
      graphics,
      count: 2,
      positionGroup: new Float32Array([0, 0, 0, int32ToFloatBits(0), 10, 0, 0, int32ToFloatBits(1)]),
      centersBounds: Sphere3D.create(Vec3.create(5, 0, 0), 8),
      textures: {
        meansLow: texture(`${tag}.meansLow`),
        meansHigh: texture(`${tag}.meansHigh`),
        scales: texture(`${tag}.scales`),
        quats: texture(`${tag}.quats`),
        sh0: texture(`${tag}.sh0`),
        shLabels: texture(`${tag}.labels`),
        shCentroids: texture(`${tag}.centroids`, 192, 1),
      },
      meansMin: Vec3.create(-1, -1, -1),
      meansMax: Vec3.create(1, 1, 1),
      codebook: new Float32Array(256 * 4),
      version: 3,
      sh: { degree: 1, coefficients: 3, count: 2 },
      maxScale: 2,
    });
    const splats = GaussianSplats.createSog(input('a'));
    expect(splats.encoding.ref.value).toBe('sog');
    expect(splats.splatCount).toBe(2);
    expect(splats.shDegree.ref.value).toBe(1);
    expect(Array.from(splats.sog.shDim.ref.value)).toEqual([192, 1]);
    // the GPU's centers a placeholder (the shader decodes the means), the CPU's whole
    expect(splats.gpuPositionGroup.ref.value.count).toBe(1);
    expect(splats.positionGroup.ref.value.count).toBe(2);
    // the centers' bounds padded by the codebook's largest scale
    const centers = calculateInvariantBoundingSphere(splats.positionGroup.ref.value.array, 2, 1, PositionGroupStride);
    expect(splats.boundingSphere.radius).toBeCloseTo(centers.radius + GaussianSplats.BoundsSigmas * 2, 5);
    expect(GaussianSplats.Utils.isCompatible!(splats, graphics)).toBe(true);
    expect(GaussianSplats.Utils.isCompatible!(splats, other)).toBe(false);
    const shape = Shape.create(
      'sog',
      {},
      splats,
      () => Color(0xffffff),
      () => 1,
      () => '',
    );
    expect(Shape.getTheme(shape).color.granularity).toBe('direct');

    // replaced: the previous textures destroyed, the new ones kept
    GaussianSplats.createSog(input('b'), splats);
    expect(destroyed.sort()).toEqual(
      ['a.centroids', 'a.labels', 'a.meansHigh', 'a.meansLow', 'a.quats', 'a.scales', 'a.sh0'].sort(),
    );
    destroyed.length = 0;
    // packed again: released
    const { positionGroup, splatData } = GaussianSplats.createArrays(1);
    GaussianSplats.create(positionGroup, splatData, 1, splats);
    expect(splats.encoding.ref.value).toBe('packed');
    expect(destroyed).toHaveLength(7);
    expect(GaussianSplats.Utils.isCompatible!(splats, other)).toBe(true);
    expect(Shape.getTheme(shape).color.granularity).not.toBe('direct');
    // disposed: nothing more to destroy
    GaussianSplats.dispose(splats);
    expect(destroyed).toHaveLength(7);
  });

  // the GPU decodes its own centers (WebGPU): the CPU's read after the first frame, the bounds the means' till then
  it("of SOG images, the centers read later: the groups mapped, the bounds the means' ranges until they are", () => {
    const texture = () => ({ getWidth: () => 2, getHeight: () => 2, destroy: () => {} }) as unknown as Texture;
    const splats = GaussianSplats.createSog({
      graphics: {} as GraphicsContext,
      count: 2,
      centersBounds: Sphere3D.create(Vec3.create(1, 2, 3), 20),
      textures: { meansLow: texture(), meansHigh: texture(), scales: texture(), quats: texture(), sh0: texture() },
      meansMin: Vec3.create(-1, -1, -1),
      meansMax: Vec3.create(1, 1, 1),
      codebook: new Float32Array(256 * 4),
      version: 3,
      maxScale: 1,
    });
    expect(splats.sog.centersPending).toBe(true);
    expect(Array.from(splats.boundingSphere.center)).toEqual([1, 2, 3]);
    expect(splats.boundingSphere.radius).toBeCloseTo(20 + GaussianSplats.BoundsSigmas, 5);
    // the groups each splat's index already
    const { indices, offsets } = splats.groupMapping;
    expect(Array.from(offsets.slice(0, 3))).toEqual([0, 1, 2]);
    expect(indices[1] / PositionGroupStride).toBe(1);
    // a splat's loci bounds the whole's
    const shape = Shape.create(
      'sog',
      {},
      splats,
      () => Color(0xffffff),
      () => 1,
      () => '',
    );
    const loci = ShapeGroup.Loci(shape, [{ ids: OrderedSet.ofSingleton(1), instance: 0 }]);
    expect(Shape.getGroupBoundingSphere(loci).radius).toBeCloseTo(splats.boundingSphere.radius, 5);

    GaussianSplats.setSogCenters(
      splats,
      new Float32Array([0, 0, 0, int32ToFloatBits(0), 4, 0, 0, int32ToFloatBits(1)]),
    );
    expect(splats.sog.centersPending).toBe(false);
    // the centers' own bounds now
    const centers = calculateInvariantBoundingSphere(splats.positionGroup.ref.value.array, 2, 1, PositionGroupStride);
    expect(splats.boundingSphere.center[0]).toBeCloseTo(centers.center[0], 5);
    expect(Array.from(Shape.getGroupBoundingSphere(loci).center)).toEqual([4, 0, 0]);
  });
});
