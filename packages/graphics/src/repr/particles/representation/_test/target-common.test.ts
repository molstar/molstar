/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { OrderedSet, SortedArray } from '@molstar/core/data/int';
import { CustomProperties } from '@molstar/model/model/custom-property';
import {
  getParticleTransforms,
  Particle,
  type ParticleList,
  ParticleValues,
} from '@molstar/model/model/particles/particle-list';
import { InstanceTransforms } from '@molstar/core/math/geometry/instance-transforms';
import { calculateTransformDataBoundingSphere } from '@molstar/graphics/gl/renderable/util';
import { getObjectTransforms } from '@molstar/graphics/gl/transform-parts';
import { Sphere3D } from '@molstar/core/math/geometry';
import { Mat4, Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import {
  createTargetParticleTransform,
  getTargetCellSize,
  isTargetCulledOnGpu,
  isTargetOccluderOff,
} from '../target/common.js';
import {
  moveTransformData,
  type TransformData,
  updateTransformData,
} from '@molstar/graphics/geo/geometry/transform-data';

function createTestParticles(): ParticleList {
  return {
    count: 3,
    keys: new Int32Array([0, 1, 2]),
    targets: new Int32Array(3),
    targetInfo: new Map([[0, {}]]),
    coordinates: new Float32Array([1, 2, 3, 4, 5, 6, 7, 8, 9]),
    radii: new Float32Array([2, 3, 4]),
    getParticleLabel: (index) => `${index}`,
    sourceData: { kind: 'test', name: 'test', data: {} },
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
  };
}

/** the transforms in the instances' order: `transform` kept as the grid's one copy, its rows' instances `aInstance` */
function instanceOrder(data: TransformData) {
  const rows = data.transform.ref.value;
  const n = data.instanceCount.ref.value;
  if (rows !== data.aTransform.ref.value) return rows.subarray(0, n * 16);
  const out = new Float32Array(n * 16);
  const instances = data.aInstance.ref.value;
  for (let i = 0; i < n; ++i) out.set(rows.subarray(i * 16, i * 16 + 16), instances[i] * 16);
  return out;
}

const sphereAtOrigin = Sphere3D.create(Vec3.create(0, 0, 0), 1);
const sphereOffCenter = Sphere3D.create(Vec3.create(10, 20, 30), 1);

// a target spread over the scene: its grid's cells grown to hold a few instances each, not one each
describe('isTargetOccluderOff', () => {
  it('turns the occluder grid off as asked and with clip objects (its rays would march what they cut away)', () => {
    const plane = { type: 'plane' as const, invert: false } as never;
    const clip = (n: number) => ({ variant: 'pixel' as const, objects: Array.from({ length: n }, () => plane) });
    expect(isTargetOccluderOff({ occluder: 'auto', clip: clip(0) })).toBe(false);
    expect(isTargetOccluderOff({ occluder: 'compute', clip: clip(0) })).toBe(false);
    expect(isTargetOccluderOff({ occluder: 'off', clip: clip(0) })).toBe(true);
    expect(isTargetOccluderOff({ occluder: 'auto', clip: clip(1) })).toBe(true);
  });
});

describe('getTargetCellSize', () => {
  const box = (side: number) => [Vec3.create(0, 0, 0), Vec3.create(side, side, side)] as const;
  it('keeps the cell size of dense and of small targets', () => {
    // 10,000 instances in 400 Å: 1,250 a 200 Å cell
    expect(getTargetCellSize(200, 10_000, ...box(400))).toBe(200);
    // a few instances, however spread
    expect(getTargetCellSize(200, 32, ...box(40_000))).toBe(200);
  });
  it('grows the cells of sparse targets to about 32 instances each', () => {
    // 5,000 instances over 40,000 Å: 0.6 a 200 Å cell
    const size = getTargetCellSize(200, 5_000, ...box(40_000));
    expect(size).toBeGreaterThan(200);
    expect((5_000 * size ** 3) / 40_000 ** 3).toBeCloseTo(32, 5);
  });
  it('counts a flat box as one cell deep', () => {
    const size = getTargetCellSize(200, 5_000, Vec3.create(0, 0, 0), Vec3.create(40_000, 40_000, 0));
    expect((5_000 * size ** 3) / (40_000 * 40_000 * 200)).toBeCloseTo(32, 5);
  });
});

// the merge's rule ahead: only where the GPU can cull and the representation merges
describe('isTargetCulledOnGpu', () => {
  const particles = { count: 10_000_000, targetInfo: new Map([[0, {}]]) } as unknown as ParticleList;
  it('is false without a GPU culler or without merging', () => {
    expect(isTargetCulledOnGpu(particles, { gpuCulling: true, renderMerged: 'auto' } as never, undefined)).toBe(false);
    expect(isTargetCulledOnGpu(particles, { gpuCulling: true, renderMerged: false } as never, undefined)).toBe(false);
  });
});

describe('createTargetParticleTransform', () => {
  it('places each particle (contiguous indices, origin-centered geometry)', () => {
    const particles = createTestParticles();
    const data = createTargetParticleTransform(particles, OrderedSet.ofBounds(0, 3), sphereAtOrigin, 100, 1000, false);
    const t = instanceOrder(data);

    expect(data.instanceCount.ref.value).toBe(3);
    for (let i = 0; i < 3; ++i) {
      // identity basis
      expect(t[i * 16 + 0]).toBe(1);
      expect(t[i * 16 + 5]).toBe(1);
      expect(t[i * 16 + 10]).toBe(1);
      expect(t[i * 16 + 15]).toBe(1);
      // translation = particle position
      expect(t[i * 16 + 12]).toBe(particles.coordinates[i * 3 + 0]);
      expect(t[i * 16 + 13]).toBe(particles.coordinates[i * 3 + 1]);
      expect(t[i * 16 + 14]).toBe(particles.coordinates[i * 3 + 2]);
    }
  });

  it('subtracts the geometry center from the translation', () => {
    const particles = createTestParticles();
    const data = createTargetParticleTransform(particles, OrderedSet.ofBounds(0, 3), sphereOffCenter, 100, 1000, false);
    const t = instanceOrder(data);

    for (let i = 0; i < 3; ++i) {
      expect(t[i * 16 + 12]).toBe(particles.coordinates[i * 3 + 0] - 10);
      expect(t[i * 16 + 13]).toBe(particles.coordinates[i * 3 + 1] - 20);
      expect(t[i * 16 + 14]).toBe(particles.coordinates[i * 3 + 2] - 30);
    }
  });

  it('subtracts the given center instead of the geometry center', () => {
    const particles = createTestParticles();
    const center = Vec3.create(1, 2, 3);
    for (const indices of [OrderedSet.ofBounds(0, 3), SortedArray.ofSortedArray([0, 1, 2])]) {
      const data = createTargetParticleTransform(
        particles,
        indices,
        sphereOffCenter,
        100,
        1000,
        false,
        undefined,
        center,
      );
      const t = instanceOrder(data);
      for (let i = 0; i < 3; ++i) {
        expect(t[i * 16 + 12]).toBe(particles.coordinates[i * 3 + 0] - 1);
        expect(t[i * 16 + 13]).toBe(particles.coordinates[i * 3 + 1] - 2);
        expect(t[i * 16 + 14]).toBe(particles.coordinates[i * 3 + 2] - 3);
      }
    }
  });

  it('gathers non-contiguous indices', () => {
    const particles = createTestParticles();
    const data = createTargetParticleTransform(
      particles,
      OrderedSet.ofSortedArray([0, 2]),
      sphereOffCenter,
      100,
      1000,
      false,
    );
    const t = instanceOrder(data);

    expect(data.instanceCount.ref.value).toBe(2);
    expect(t[12]).toBe(1 - 10);
    expect(t[16 + 12]).toBe(7 - 10);
  });

  it('scales by the particle radius', () => {
    const particles = createTestParticles();
    const data = createTargetParticleTransform(particles, OrderedSet.ofBounds(0, 3), sphereOffCenter, 100, 1000, true);
    const t = instanceOrder(data);

    for (let i = 0; i < 3; ++i) {
      const s = particles.radii![i];
      expect(t[i * 16 + 0]).toBe(s);
      expect(t[i * 16 + 5]).toBe(s);
      expect(t[i * 16 + 10]).toBe(s);
      // translation accounts for the scaled basis
      expect(t[i * 16 + 12]).toBe(particles.coordinates[i * 3 + 0] - s * 10);
      expect(t[i * 16 + 13]).toBe(particles.coordinates[i * 3 + 1] - s * 20);
      expect(t[i * 16 + 14]).toBe(particles.coordinates[i * 3 + 2] - s * 30);
    }
  });

  it('matches between the bulk and gather paths', () => {
    const particles = createTestParticles();
    const bulk = createTargetParticleTransform(particles, OrderedSet.ofBounds(0, 3), sphereOffCenter, 100, 1000, false);
    // a raw SortedArray is not normalized to an Interval, forcing the per-instance gather path
    const gather = createTargetParticleTransform(
      particles,
      SortedArray.ofSortedArray([0, 1, 2]),
      sphereOffCenter,
      100,
      1000,
      false,
    );

    expect(Array.from(instanceOrder(gather))).toEqual(Array.from(instanceOrder(bulk)));
  });

  it('keeps the transforms once, in the grid order, through updates, moves and a matrix', () => {
    // particles spread over several grid cells: the grid order not the particles'
    const n = 40;
    const particles: ParticleList = {
      ...createTestParticles(),
      count: n,
      keys: Int32Array.from({ length: n }, (_, i) => i),
      targets: new Int32Array(n),
      coordinates: Float32Array.from({ length: n * 3 }, (_, i) => ((i * 7919) % 97) * 3 - 150),
      radii: new Float32Array(n).fill(1),
    };
    const expected = (center: Vec3) => {
      const out = new Float32Array(n * 16);
      for (let i = 0; i < n; ++i) {
        Mat4.toArray(Mat4.identity(), out, i * 16);
        for (let k = 0; k < 3; ++k) out[i * 16 + 12 + k] = particles.coordinates[i * 3 + k] - center[k];
      }
      return Array.from(out);
    };
    const data = createTargetParticleTransform(particles, OrderedSet.ofBounds(0, n), sphereOffCenter, 20, 4, false);
    // one copy: the grid's, which is also what is drawn
    expect(data.transform.ref.value).toBe(data.aTransform.ref.value);
    expect(data.transform.ref.value).toBe(data.instanceGrid.ref.value.cellTransform);
    expect(Array.from(data.aInstance.ref.value)).not.toEqual(Array.from({ length: n }, (_, i) => i));
    expect(Array.from(instanceOrder(data))).toEqual(expected(sphereOffCenter.center));

    // updated: new transforms, not written into the grid's (still read by what was built from it)
    const grid = data.instanceGrid.ref.value;
    const before = Array.from(grid.cellTransform);
    createTargetParticleTransform(particles, OrderedSet.ofBounds(0, n), sphereAtOrigin, 20, 4, false, data);
    expect(Array.from(grid.cellTransform)).toEqual(before);
    expect(data.transform.ref.value).toBe(data.aTransform.ref.value);
    expect(Array.from(instanceOrder(data))).toEqual(expected(sphereAtOrigin.center));

    // moved: the grid's moved, the same one copy
    expect(moveTransformData(data, Mat4.fromTranslation(Mat4(), Vec3.create(5, 0, 0)))).toBe(true);
    expect(data.transform.ref.value).toBe(data.aTransform.ref.value);
    const moved = expected(Vec3.create(-5, 0, 0));
    expect(Array.from(instanceOrder(data))).toEqual(moved);

    // a matrix: the instances' order restored for the rebuild, the transforms as before under it; one copy again
    // once the matrix is gone
    Mat4.fromTranslation(data.matrix.ref.value, Vec3.create(0, 7, 0));
    updateTransformData(data, sphereAtOrigin, 20, 4);
    expect(data.transform.ref.value).not.toBe(data.aTransform.ref.value);
    expect(Array.from(data.transform.ref.value.subarray(0, n * 16))).toEqual(moved);
    const drawn = data.aTransform.ref.value;
    for (let i = 0; i < n; ++i) {
      const id = data.aInstance.ref.value[i];
      expect(drawn[i * 16 + 13]).toBe(moved[id * 16 + 13] + 7);
    }
    Mat4.setIdentity(data.matrix.ref.value);
    updateTransformData(data, sphereAtOrigin, 20, 4, true);
    expect(data.transform.ref.value).toBe(data.aTransform.ref.value);
    expect(Array.from(instanceOrder(data))).toEqual(moved);
  });

  it('rotates, scales and centers as the list transforms, without keeping a copy of them', () => {
    const rotated = (): ParticleList => {
      const q = Quat.setAxisAngle(Quat(), Vec3.normalize(Vec3(), Vec3.create(1, 2, 3)), 0.7);
      return {
        ...createTestParticles(),
        rotations: new Float32Array([0, 0, 0, 1, ...q, ...Quat.conjugate(Quat(), q)]),
      };
    };
    for (const scaleByRadius of [false, true]) {
      for (const indices of [OrderedSet.ofBounds(0, 3), SortedArray.ofSortedArray([0, 2])]) {
        const particles = rotated();
        const center = Vec3.create(10, 20, 30);
        const data = createTargetParticleTransform(particles, indices, sphereOffCenter, 100, 1000, scaleByRadius);
        // the transforms a target used before: the list's, scaled, then the center subtracted
        const t = instanceOrder(data);
        const reference = getParticleTransforms(rotated());
        for (let i = 0; i < OrderedSet.size(indices); ++i) {
          const p = OrderedSet.getAt(indices, i);
          const s = scaleByRadius ? particles.radii![p] : 1;
          const m = Mat4.fromArray(Mat4(), reference, p * 16);
          for (let j = 0; j < 12; ++j) m[j] *= s;
          for (let j = 0; j < 12; ++j) expect(t[i * 16 + j]).toBeCloseTo(m[j], 5);
          for (let k = 0; k < 3; ++k) {
            const moved = m[k] * center[0] + m[4 + k] * center[1] + m[8 + k] * center[2];
            expect(t[i * 16 + 12 + k]).toBeCloseTo(m[12 + k] - moved, 4);
          }
          expect(t[i * 16 + 15]).toBe(1);
        }
        expect(particles._propertyData['particle-transforms']).toBeUndefined();
      }
    }
  });

  it("sizes a pack's target's cells by its particles' bounds, no pass over them", () => {
    // two targets as runs (a pack's), a radius each; positions spread, rotated
    const count = 40;
    const coordinates = new Float32Array(count * 3);
    const rotations = new Float32Array(count * 4);
    for (let i = 0; i < count; ++i) {
      coordinates.set([Math.sin(i) * 300, Math.cos(i * 1.3) * 200, (i % 7) * 40], i * 3);
      const q = Quat.setAxisAngle(Quat(), Vec3.normalize(Vec3(), Vec3.create(1, i, 2)), i * 0.4);
      rotations.set(q, i * 4);
    }
    const offsets = Int32Array.from([0, 25, count]);
    const list = (runs: boolean): ParticleList => ({
      ...createTestParticles(),
      count,
      coordinates,
      rotations,
      targets: runs
        ? ParticleValues.ofRuns(offsets, Int32Array.from([0, 1]))
        : new Int32Array(count).fill(0, 0, 25).fill(1, 25),
      targetInfo: new Map([
        [0, {}],
        [1, {}],
      ]),
      radii: ParticleValues.ofRuns(offsets, Float32Array.from([2, 3])),
      customProperties: new CustomProperties(),
      _propertyData: Object.create(null),
    });
    const indices = OrderedSet.ofBounds(0, 25);
    const invariant = Sphere3D.create(Vec3.create(1, -2, 0.5), 4);
    const writes = vi.spyOn(InstanceTransforms, 'writePoints');
    for (const center of [Vec3.create(1, -2, 0.5), Vec3.create(4, 0, -3)]) {
      const particles = list(true);
      writes.mockClear();
      // poses (spheres' geometry), culled on the GPU (the cell size by the box), scaled by the radius: no pass, the grid
      // pending (built where a CPU cull needs it), a few instances' sphere by the boundary helper (a point at a time)
      const data = createTargetParticleTransform(
        particles,
        indices,
        invariant,
        100,
        1000,
        true,
        undefined,
        center,
        true,
        true,
      );
      const sphere = calculateTransformDataBoundingSphere(invariant, data);
      expect(writes).not.toHaveBeenCalled();
      expect(data.instanceGrid.ref.value.pending).toBeDefined();
      // every instance's center within the sphere less the invariant radius
      const t = getObjectTransforms(data);
      expect(t.kind).toBe('poses');
      const point = Vec3();
      for (let i = 0; i < 25; ++i) {
        InstanceTransforms.getPoint(t, i, invariant.center, point);
        expect(Vec3.distance(point, sphere.center)).toBeLessThanOrEqual(sphere.radius - invariant.radius + 1e-3);
      }
      // the placed point the invariant center: the cell size as from the instances' box
      if (Vec3.exactEquals(center, invariant.center)) {
        writes.mockClear();
        const computed = createTargetParticleTransform(
          list(false),
          indices,
          invariant,
          100,
          1000,
          true,
          undefined,
          center,
          true,
          true,
        );
        // without runs no bounds: the box's pass
        expect(writes).toHaveBeenCalledTimes(1);
        expect(data.instanceGrid.ref.value.cellSize).toBe(computed.instanceGrid.ref.value.cellSize);
      }
    }
    writes.mockRestore();
  });

  it("keeps a GPU-culled target's pending grid and cell size where its particles are posed otherwise", () => {
    const count = 40;
    const coordinates = new Float32Array(count * 3).map((_, k) => (k % 3 === 0 ? (k / 3) * 30 : (k % 3) * 5));
    const offsets = Int32Array.from([0, count]);
    const list: ParticleList = {
      ...createTestParticles(),
      count,
      coordinates,
      targets: ParticleValues.ofRuns(offsets, Int32Array.from([0])),
      radii: undefined,
      customProperties: new CustomProperties(),
      _propertyData: Object.create(null),
    };
    const indices = OrderedSet.ofBounds(0, count);
    const invariant = Sphere3D.create(Vec3(), 2);
    const create = (particles: ParticleList, data?: TransformData, posed = false) =>
      createTargetParticleTransform(particles, indices, invariant, 100, 1000, false, data, Vec3(), true, true, posed);
    const data = create(list);
    const grid = data.instanceGrid.ref.value;
    const cellSize = data.instanceCellSize.ref.value;
    expect(grid.pending).toBeDefined();
    // posed otherwise, spread wider (another cell size by the bounds): the grid and its size kept
    const spread = new Float32Array(coordinates).map((v, k) => (k % 3 === 0 ? v * 3 : v));
    create(Particle.withPoses(list, spread), data, true);
    expect(data.instanceGrid.ref.value).toBe(grid);
    expect(data.instanceCellSize.ref.value).toBe(cellSize);
    // not posed (other particles): sized by their bounds again
    create(
      { ...list, coordinates: spread, customProperties: new CustomProperties(), _propertyData: Object.create(null) },
      data,
    );
    expect(data.instanceCellSize.ref.value).not.toBe(cellSize);
    expect(data.instanceGrid.ref.value.pending).toBeDefined();
  });
});
