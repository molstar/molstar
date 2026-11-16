/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { InstanceTransforms } from '@molstar/core/math/geometry/instance-transforms';
import { Mat4, Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { CulledData, CulledSegmentVec4s } from '../culled-data.js';
import type { TransformPart } from '../transform-parts.js';
import { GridPositions } from '@molstar/core/math/geometry/grid-positions';

/** `data` written by `layout` */
function writeAll(layout: ReturnType<typeof CulledData.layout>) {
  const bytes = new ArrayBuffer(layout.vec4s * 16);
  const written = new Uint8Array(bytes);
  layout.write((d, o) => {
    expect(o + d.byteLength).toBeLessThanOrEqual(bytes.byteLength);
    written.set(new Uint8Array(d.buffer, d.byteOffset, d.byteLength), o);
  });
  return { f: new Float32Array(bytes), u: new Uint32Array(bytes) };
}

/** `CulledData.wgsl`'s readers on the CPU: instance `i`'s rows, `aInstance` and segment */
function read(data: { f: Float32Array; u: Uint32Array }, i: number) {
  const { f, u } = data;
  const vec4 = (k: number) => Array.from(f.subarray(k * 4, k * 4 + 4));
  const [n, , slotFirst] = u.subarray(0, 4);
  const q = Math.ceil(n / 4);
  const word = u[4 + i];
  const o = slotFirst + (word & 0x7fffffff);
  const segment = f[(1 + 2 * q) * 4 + i];
  const instance = u[(1 + q) * 4 + i];
  const groups = vec4(1 + 3 * q + CulledSegmentVec4s * segment).slice(0, 3);
  const pose = (word & 0x80000000) !== 0;
  if (!pose) return { rows: [vec4(o), vec4(o + 1), vec4(o + 2)].flat(), instance, segment, groups, pose };
  const [tx, ty, tz, sc] = vec4(o);
  const [x, y, z, w] = vec4(o + 1);
  const rows = [
    [sc * (1 - 2 * (y * y + z * z)), sc * 2 * (x * y - w * z), sc * 2 * (x * z + w * y), tx],
    [sc * 2 * (x * y + w * z), sc * (1 - 2 * (x * x + z * z)), sc * 2 * (y * z - w * x), ty],
    [sc * 2 * (x * z - w * y), sc * 2 * (y * z + w * x), sc * (1 - 2 * (x * x + y * y)), tz],
  ].flat();
  return { rows, instance, segment, groups, pose };
}

function createList(count: number) {
  const positions = new Float32Array(count * 3);
  const rotations = new Float32Array(count * 4);
  const scales = new Float32Array(count);
  const q = Quat();
  for (let e = 0; e < count; ++e) {
    positions.set([e * 1.5, -e, e * 0.25], e * 3);
    Quat.setAxisAngle(q, Vec3.normalize(Vec3(), Vec3.create(1, e, 2)), e * 0.9);
    Quat.toArray(q, rotations, e * 4);
    scales[e] = e === 3 ? 0 : 0.5 + e * 0.1;
  }
  return { positions, rotations, scales };
}

describe('CulledData', () => {
  it('reads the rows of matrices and of poses, by index, a pose an id of its part', () => {
    const list = createList(20);
    const snorm = Int16Array.from(list.rotations, (v) => Math.round(v * 32767));
    const center = Vec3.create(0.5, -1, 2);
    // a similarity, folded into the poses: rotated, scaled uniformly, moved
    const similar = Mat4.mul(
      Mat4(),
      Mat4.fromTranslation(Mat4(), Vec3.create(3, -2, 1)),
      Mat4.mul(Mat4(), Mat4.fromRotation(Mat4(), 1.1, Vec3.create(0, 1, 0)), Mat4.fromUniformScaling(Mat4(), 2)),
    );
    // not one: its poses' instances on rows
    const sheared = Mat4.fromScaling(Mat4(), Vec3.create(1, 2, 3));
    const matrices = new Float32Array(3 * 16);
    for (let i = 0; i < 3; ++i) {
      Mat4.toArray(Mat4.fromScaling(Mat4(), Vec3.create(1 + i, 2, 0.5)), matrices, i * 16);
      matrices[i * 16 + 12] = i;
    }
    // segments of matrices; of poses by a range (in the grid's order, with the similarity); by indices, twice (one
    // source); without rotations or scales; with the shear; an empty one between
    const transforms = [
      InstanceTransforms.ofMatrices(matrices, 3),
      InstanceTransforms.ofPoses({ ...list, count: 4, first: 2, order: [2, 0, 3, 1], center, matrix: similar }),
      InstanceTransforms.ofPoses({ ...list, count: 3, indices: [17, 4, 9], center }),
      // positions on a grid (a particle pack's), decoded where written
      InstanceTransforms.ofPoses({ positions: GridPositions.fromFloat(list.positions), count: 2, first: 10, center }),
      InstanceTransforms.ofPoses({ ...list, count: 2, indices: [6, 5], center }),
      // snorm16 rotations (a particle pack's), normalized where written
      InstanceTransforms.ofPoses({ ...list, rotations: snorm, count: 2, indices: [7, 8], center }),
      InstanceTransforms.ofPoses({ ...list, count: 2, first: 0, center, matrix: sheared }),
    ];
    const firsts = [0, 3, 7, 10, 12, 16, 14];
    const parts: TransformPart[] = transforms.map((t, k) => ({ transforms: t, first: firsts[k] }));
    const n = 18;
    const segments = {
      count: 8,
      instanceBases: [0, 3, 7, 7, 10, 12, 14, 16],
      instanceCounts: [3, 4, 0, 3, 2, 2, 2, 2],
    };
    const partSegments = [0, 1, 3, 4, 5, 7, 6];
    const groups = new Float32Array(8 * CulledSegmentVec4s * 4).map((_, k) => k + 1);
    const aInstance = new Uint32Array(n).map((_, k) => 100 + k);
    for (const block of [2, 5, 1 << 18]) {
      const layout = CulledData.layout(parts, aInstance, segments, groups, n, block);
      // header, words, aInstance, segment indices, segments' data (`CulledSegmentVec4s` each); rows of the matrices and the sheared poses; two vec4s an id
      // of the other poses' parts
      expect(layout.vec4s).toBe(1 + 3 * 5 + CulledSegmentVec4s * 8 + 3 * (3 + 2) + 2 * (4 + 3 + 2 + 2 + 2));
      const data = writeAll(layout);
      const expected = new Float32Array(12);
      for (let k = 0; k < parts.length; ++k) {
        const { transforms: t, first } = parts[k];
        const segment = partSegments[k];
        for (let i = 0; i < t.count; ++i) {
          InstanceTransforms.writeAffineRows(t, i, 1, expected);
          const r = read(data, first + i);
          expect(r.pose).toBe(k > 0 && k < 6);
          for (let j = 0; j < 12; ++j) expect(r.rows[j]).toBeCloseTo(expected[j], 4);
          expect(r.instance).toBe(100 + first + i);
          expect(r.segment).toBe(segment);
          const g = segment * CulledSegmentVec4s * 4;
          expect(r.groups).toEqual(Array.from(groups.subarray(g, g + 3)));
        }
      }
    }
  });

  it('writes the same words with its runs applied as without them', () => {
    const list = createList(20);
    const center = Vec3.create(0.5, -1, 2);
    const matrices = new Float32Array(3 * 16);
    for (let i = 0; i < 3; ++i) Mat4.toArray(Mat4.fromTranslation(Mat4(), Vec3.create(i, 2, -i)), matrices, i * 16);
    const transforms = [
      InstanceTransforms.ofMatrices(matrices, 3),
      // in an order of its own: its slot words written by the CPU
      InstanceTransforms.ofPoses({ ...list, count: 4, first: 2, order: [2, 0, 3, 1], center }),
      InstanceTransforms.ofPoses({ ...list, count: 3, indices: [17, 4, 9], center }),
      InstanceTransforms.ofPoses({ ...list, count: 5, first: 10, center }),
    ];
    const firsts = [0, 3, 7, 10];
    const parts: TransformPart[] = transforms.map((t, k) => ({ transforms: t, first: firsts[k] }));
    const n = 15;
    const segments = { count: 5, instanceBases: [0, 3, 7, 7, 10], instanceCounts: [3, 4, 0, 3, 5] };
    // serial in the first two segments (one run), not in the fourth, serial again in the last (from another start)
    const aInstance = new Uint32Array([0, 1, 2, 3, 4, 5, 6, 9, 7, 8, 50, 51, 52, 53, 54]);
    for (const block of [2, 1 << 18]) {
      const layout = CulledData.layout(parts, aInstance, segments, undefined, n, block);
      const expected = writeAll(layout).u;
      const bytes = new ArrayBuffer(layout.vec4s * 16);
      const written = new Uint8Array(bytes);
      layout.write((d, o) => written.set(new Uint8Array(d.buffer, d.byteOffset, d.byteLength), o), { runs: false });
      const u = new Uint32Array(bytes);
      for (const r of layout.runs) {
        for (let i = 0; i < r.count; ++i) u[r.word + i] = ((r.start + r.step * i) | r.bits) >>> 0;
      }
      expect(Array.from(u)).toEqual(Array.from(expected));
      // the ordered part's slot words and the fourth segment's aInstance not runs; aInstance's first runs merged
      const q = Math.ceil(n / 4);
      expect(layout.runs.some((r) => r.word <= 4 + 3 && 4 + 3 < r.word + r.count)).toBe(false);
      expect(layout.runs.some((r) => r.word <= (1 + q) * 4 + 7 && (1 + q) * 4 + 7 < r.word + r.count)).toBe(false);
      expect(layout.runs.find((r) => r.word === (1 + q) * 4)?.count).toBe(7);
    }
  });

  it('writes a header alone for no instances', () => {
    const layout = CulledData.layout(
      [],
      new Float32Array(0),
      { count: 0, instanceBases: [], instanceCounts: [] },
      undefined,
      0,
      16,
    );
    expect(layout.vec4s).toBe(1);
    expect(Array.from(writeAll(layout).u)).toEqual([0, 0, 1, 0]);
  });
});
