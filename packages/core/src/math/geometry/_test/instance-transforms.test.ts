/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Mat4 } from '../../linear-algebra/3d/mat4.js';
import { Quat } from '../../linear-algebra/3d/quat.js';
import { Vec3 } from '../../linear-algebra/3d/vec3.js';
import { GridPositions } from '../grid-positions.js';
import { InstanceTransforms, type InstancePoses } from '../instance-transforms.js';
import { calcInstanceGrid } from '../instance-grid.js';
import { Sphere3D } from '../primitives/sphere3d.js';
import { fillSerial } from '../../../util/array.js';

/** A deterministic generator, so that a failure repeats */
function random(seed: number) {
  let s = seed;
  return () => {
    s = (s * 1664525 + 1013904223) >>> 0;
    return s / 4294967296;
  };
}

/** Poses of `elements` particles, `count` instances of them in a shuffled order, with scales, a center and a matrix */
function createPoses(
  elements: number,
  count: number,
  options: { scales?: boolean; matrix?: boolean; range?: boolean },
) {
  const rnd = random(42);
  const positions = new Float32Array(elements * 3);
  const rotations = new Float32Array(elements * 4);
  const scales = options.scales ? new Float32Array(elements) : undefined;
  const q = Quat();
  for (let e = 0; e < elements; ++e) {
    for (let k = 0; k < 3; ++k) positions[e * 3 + k] = (rnd() - 0.5) * 1000;
    Quat.normalize(q, Quat.create(rnd() - 0.5, rnd() - 0.5, rnd() - 0.5, rnd() - 0.5));
    rotations.set(q, e * 4);
    if (scales) scales[e] = 0.5 + rnd() * 3;
  }
  const order = Array.from({ length: count }, (_, i) => i).sort(() => rnd() - 0.5);
  const indices = Array.from({ length: count }, (_, i) => (i * 7) % elements);
  const matrix = options.matrix
    ? Mat4.mul(
        Mat4(),
        Mat4.fromTranslation(Mat4(), Vec3.create(10, -20, 30)),
        Mat4.mul(Mat4(), Mat4.fromRotation(Mat4(), 0.7, Vec3.create(1, 2, 3)), Mat4.fromUniformScaling(Mat4(), 1.5)),
      )
    : undefined;
  const poses: InstancePoses = {
    kind: 'poses',
    count,
    positions,
    rotations,
    scales,
    order,
    ...(options.range ? { first: 3 } : { indices }),
    center: Vec3.create(4, -2, 7),
    matrix,
  };
  return poses;
}

/** The matrices of `poses` built by `Mat4`'s operations: T(position) R(q) S(scale) T(-center), then the matrix */
function referenceMatrices(t: InstancePoses) {
  const out = new Float32Array(t.count * 16);
  const m = Mat4();
  for (let i = 0; i < t.count; ++i) {
    const id = t.order ? t.order[i] : i;
    const e = t.indices ? t.indices[id] : (t.first ?? 0) + id;
    const q = Quat.create(t.rotations[e * 4], t.rotations[e * 4 + 1], t.rotations[e * 4 + 2], t.rotations[e * 4 + 3]);
    const s = t.scales ? t.scales[e] : 1;
    Mat4.fromTranslation(m, Vec3.fromArray(Vec3(), t.positions, e * 3));
    Mat4.mul(m, m, Mat4.fromQuat(Mat4(), q));
    Mat4.mul(m, m, Mat4.fromUniformScaling(Mat4(), s));
    Mat4.mul(m, m, Mat4.fromTranslation(Mat4(), Vec3.negate(Vec3(), t.center)));
    if (t.matrix) Mat4.mul(m, t.matrix, m);
    Mat4.toArray(m, out, i * 16);
  }
  return InstanceTransforms.ofMatrices(out, t.count);
}

// the readers' kernels give the same for poses as for their matrices
describe('InstanceTransforms', () => {
  const variants = [
    { name: 'indices', scales: false, matrix: false, range: false },
    { name: 'scales', scales: true, matrix: false, range: false },
    { name: 'a matrix', scales: true, matrix: true, range: false },
    { name: 'a range', scales: false, matrix: true, range: true },
  ];
  const local = Vec3.create(1.5, -3, 2.25);
  for (const v of variants) {
    it(`places poses as their matrices do, with ${v.name}`, () => {
      const poses = createPoses(50, 40, v);
      const matrices = referenceMatrices(poses);
      const a = Vec3();
      const b = Vec3();
      const m = Mat4();
      for (let i = 0; i < poses.count; ++i) {
        expect(
          Mat4.areEqual(
            InstanceTransforms.getMatrix(poses, i, m),
            Mat4.fromArray(Mat4(), matrices.matrices, i * 16),
            1e-3,
          ),
        ).toBe(true);
        InstanceTransforms.getPoint(poses, i, local, a);
        InstanceTransforms.getPoint(matrices, i, local, b);
        expect(Vec3.distance(a, b)).toBeLessThan(1e-3);
        InstanceTransforms.getDirection(poses, i, local, a);
        InstanceTransforms.getDirection(matrices, i, local, b);
        expect(Vec3.distance(a, b)).toBeLessThan(1e-5);
        expect(Vec3.magnitude(a)).toBeCloseTo(1, 5);
        expect(InstanceTransforms.getMaxScale(poses, i)).toBeCloseTo(InstanceTransforms.getMaxScale(matrices, i), 4);
      }
      expect(InstanceTransforms.getMaxScaleOfAll(poses)).toBeCloseTo(InstanceTransforms.getMaxScaleOfAll(matrices), 4);

      const [xa, ya, za, xb, yb, zb] = Array.from({ length: 6 }, () => new Float32Array(poses.count));
      const [minA, maxA, minB, maxB] = [Vec3(), Vec3(), Vec3(), Vec3()];
      InstanceTransforms.writePoints(poses, local, xa, ya, za, minA, maxA);
      InstanceTransforms.writePoints(matrices, local, xb, yb, zb, minB, maxB);
      for (let i = 0; i < poses.count; ++i) {
        expect(Math.abs(xa[i] - xb[i]) + Math.abs(ya[i] - yb[i]) + Math.abs(za[i] - zb[i])).toBeLessThan(3e-3);
      }
      expect(Vec3.distance(minA, minB)).toBeLessThan(3e-3);
      expect(Vec3.distance(maxA, maxB)).toBeLessThan(3e-3);
      const from = Vec3.create(5, 5, 5);
      const da = InstanceTransforms.getMaxDistanceSq(poses, local, from);
      const db = InstanceTransforms.getMaxDistanceSq(matrices, local, from);
      expect(Math.abs(Math.sqrt(da) - Math.sqrt(db))).toBeLessThan(3e-3);
    });
  }

  it('reads snorm16 rotations (`Int16Array`) as the unit quaternions they quantize, as they are', () => {
    const float = createPoses(50, 40, { scales: true, matrix: true, range: false });
    const r = float.rotations as Float32Array;
    const snorm = Int16Array.from(r, (v) => Math.round(v * 32767));
    // the quaternions the integers stand for, normalized
    const unit = new Float32Array(r.length);
    for (let o = 0; o < r.length; o += 4) {
      const n = Math.hypot(snorm[o], snorm[o + 1], snorm[o + 2], snorm[o + 3]);
      for (let k = 0; k < 4; ++k) unit[o + k] = snorm[o + k] / n;
    }
    const a = InstanceTransforms.ofPoses({ ...float, rotations: snorm });
    const b = InstanceTransforms.ofPoses({ ...float, rotations: unit });
    const local = Vec3.create(1.5, -3, 2.25);
    const [pa, pb, ma, mb] = [Vec3(), Vec3(), Mat4(), Mat4()];
    for (let i = 0; i < a.count; ++i) {
      InstanceTransforms.getPoint(a, i, local, pa);
      InstanceTransforms.getPoint(b, i, local, pb);
      expect(Vec3.distance(pa, pb)).toBeLessThan(1e-3);
      InstanceTransforms.getDirection(a, i, local, pa);
      InstanceTransforms.getDirection(b, i, local, pb);
      expect(Vec3.distance(pa, pb)).toBeLessThan(1e-5);
      expect(Mat4.areEqual(InstanceTransforms.getMatrix(a, i, ma), InstanceTransforms.getMatrix(b, i, mb), 1e-3)).toBe(
        true,
      );
    }
    const ra = InstanceTransforms.writeAffineRows(a, 0, a.count, new Float32Array(a.count * 12));
    const rb = InstanceTransforms.writeAffineRows(b, 0, b.count, new Float32Array(b.count * 12));
    for (let k = 0; k < ra.length; ++k) expect(Math.abs(ra[k] - rb[k])).toBeLessThan(1e-3);
    const [min, max, minB, maxB] = [Vec3(), Vec3(), Vec3(), Vec3()];
    InstanceTransforms.writePoints(a, local, undefined, undefined, undefined, min, max);
    InstanceTransforms.writePoints(b, local, undefined, undefined, undefined, minB, maxB);
    expect(Vec3.distance(min, minB)).toBeLessThan(1e-3);
    expect(Vec3.distance(max, maxB)).toBeLessThan(1e-3);
    const from = Vec3.create(5, 5, 5);
    const da = Math.sqrt(InstanceTransforms.getMaxDistanceSq(a, local, from));
    const db = Math.sqrt(InstanceTransforms.getMaxDistanceSq(b, local, from));
    expect(Math.abs(da - db)).toBeLessThan(1e-3);
  });

  it('reads positions on a grid (`GridPositions`) as the float32 positions they decode to, as they are', () => {
    const float = createPoses(50, 40, { scales: true, matrix: true, range: false });
    const grid = GridPositions.fromFloat(float.positions as Float32Array);
    const a = InstanceTransforms.ofPoses({ ...float, positions: grid });
    const b = InstanceTransforms.ofPoses({ ...float, positions: GridPositions.toFloat(grid) });
    const local = Vec3.create(1.5, -3, 2.25);
    const [pa, pb, ma, mb] = [Vec3(), Vec3(), Mat4(), Mat4()];
    for (let i = 0; i < a.count; ++i) {
      InstanceTransforms.getPoint(a, i, local, pa);
      InstanceTransforms.getPoint(b, i, local, pb);
      expect(Vec3.distance(pa, pb)).toBeLessThan(1e-3);
      expect(Mat4.areEqual(InstanceTransforms.getMatrix(a, i, ma), InstanceTransforms.getMatrix(b, i, mb), 1e-3)).toBe(
        true,
      );
    }
    const ra = InstanceTransforms.writeAffineRows(a, 0, a.count, new Float32Array(a.count * 12));
    const rb = InstanceTransforms.writeAffineRows(b, 0, b.count, new Float32Array(b.count * 12));
    for (let k = 0; k < ra.length; ++k) expect(Math.abs(ra[k] - rb[k])).toBeLessThan(1e-3);
    const [min, max, minB, maxB] = [Vec3(), Vec3(), Vec3(), Vec3()];
    InstanceTransforms.writePoints(a, local, undefined, undefined, undefined, min, max);
    InstanceTransforms.writePoints(b, local, undefined, undefined, undefined, minB, maxB);
    expect(Vec3.distance(min, minB)).toBeLessThan(1e-3);
    expect(Vec3.distance(max, maxB)).toBeLessThan(1e-3);
    const from = Vec3.create(5, 5, 5);
    const da = Math.sqrt(InstanceTransforms.getMaxDistanceSq(a, local, from));
    const db = Math.sqrt(InstanceTransforms.getMaxDistanceSq(b, local, from));
    expect(Math.abs(da - db)).toBeLessThan(1e-3);
  });

  it('writes the rows of poses as those of their matrices', () => {
    const poses = createPoses(50, 40, { scales: true, matrix: true, range: false });
    const matrices = referenceMatrices(poses);
    const a = InstanceTransforms.writeAffineRows(poses, 5, 30, new Float32Array(30 * 16), 4, 16);
    const b = InstanceTransforms.writeAffineRows(matrices, 5, 30, new Float32Array(30 * 16), 4, 16);
    for (let k = 0; k < a.length; ++k) expect(Math.abs(a[k] - b[k])).toBeLessThan(1e-3);
  });

  it('scales every element by `scale` where poses have no `scales`, as by `scales` of that value', () => {
    const base = createPoses(50, 40, { matrix: true, range: false });
    const uniform = InstanceTransforms.ofPoses({ ...base, scales: undefined, scale: 2.5 });
    const perElement = InstanceTransforms.ofPoses({ ...base, scales: new Float32Array(50).fill(2.5) });
    const rows = (t: InstancePoses) => InstanceTransforms.writeAffineRows(t, 0, 40, new Float32Array(40 * 12));
    expect(Array.from(rows(uniform))).toEqual(Array.from(rows(perElement)));
    expect(InstanceTransforms.getMaxScaleOfAll(uniform)).toBeCloseTo(
      InstanceTransforms.getMaxScaleOfAll(perElement),
      6,
    );
    const points = (t: InstancePoses) => {
      const x = new Float32Array(40);
      const y = new Float32Array(40);
      const z = new Float32Array(40);
      InstanceTransforms.writePoints(t, Vec3.create(1, 2, 3), x, y, z, Vec3(), Vec3());
      return [...x, ...y, ...z];
    };
    expect(points(uniform)).toEqual(points(perElement));
    const p = Vec3();
    const q = Vec3();
    for (let i = 0; i < 40; ++i) {
      InstanceTransforms.getPoint(uniform, i, Vec3.create(1, 0, 0), p);
      InstanceTransforms.getPoint(perElement, i, Vec3.create(1, 0, 0), q);
      expect(Vec3.distance(p, q)).toBeLessThan(1e-4);
    }
  });

  it('places poses without rotations and with scales not above 0 as translations and a scale of 1', () => {
    const t: InstancePoses = {
      kind: 'poses',
      count: 2,
      positions: Float32Array.of(1, 2, 3, 10, 20, 30),
      scales: Float32Array.of(2, 0),
      center: Vec3.create(1, 1, 1),
    };
    const out = Vec3();
    expect(Array.from(InstanceTransforms.getPoint(t, 0, Vec3.create(2, 1, 1), out))).toEqual([3, 2, 3]);
    expect(Array.from(InstanceTransforms.getPoint(t, 1, Vec3.create(2, 1, 1), out))).toEqual([11, 20, 30]);
    expect(InstanceTransforms.getMaxScaleOfAll(t)).toBe(2);
    const m = InstanceTransforms.getMatrix(t, 0, Mat4());
    expect(
      Mat4.areEqual(
        m,
        Mat4.mul(Mat4(), Mat4.fromTranslation(Mat4(), Vec3.create(-1, 0, 1)), Mat4.fromUniformScaling(Mat4(), 2)),
        1e-6,
      ),
    ).toBe(true);
  });

  it('builds the grid of poses as that of their matrices, without copying them', () => {
    const poses = { ...createPoses(500, 400, { scales: true, matrix: true, range: false }), order: undefined };
    const matrices = referenceMatrices(poses);
    const sphere = Sphere3D.create(Vec3.create(4, -2, 7), 3);
    const instance = fillSerial(new Float32Array(poses.count));
    const a = calcInstanceGrid(
      { instanceCount: poses.count, instance, transform: new Float32Array(0), poses, invariantBoundingSphere: sphere },
      50,
      8,
    );
    const b = calcInstanceGrid(
      { instanceCount: poses.count, instance, transform: matrices.matrices, invariantBoundingSphere: sphere },
      50,
      8,
    );
    expect(a.cellTransform.length).toBe(0);
    expect(Array.from(a.cellInstance)).toEqual(Array.from(b.cellInstance));
    expect(Array.from(a.cellOffsets)).toEqual(Array.from(b.cellOffsets));
    for (let k = 0; k < a.cellSpheres.length; ++k)
      expect(Math.abs(a.cellSpheres[k] - b.cellSpheres[k])).toBeLessThan(1e-2);
    // drawn in the grid's order: its matrices those of the poses ordered by it
    const ordered: InstancePoses = { ...poses, order: a.cellInstance };
    const m = Mat4();
    for (let i = 0; i < poses.count; ++i) {
      expect(
        Mat4.areEqual(
          InstanceTransforms.getMatrix(ordered, i, m),
          Mat4.fromArray(Mat4(), b.cellTransform, i * 16),
          1e-3,
        ),
      ).toBe(true);
    }
  });

  it('reads matrices from an offset', () => {
    const m = Mat4.fromTranslation(Mat4(), Vec3.create(1, 2, 3));
    const arr = new Float32Array(48);
    Mat4.toArray(m, arr, 16);
    const t = InstanceTransforms.ofMatrices(arr, 1, 16);
    expect(Vec3.equals(InstanceTransforms.getPoint(t, 0, Vec3(), Vec3()), Vec3.create(1, 2, 3))).toBe(true);
  });
});
