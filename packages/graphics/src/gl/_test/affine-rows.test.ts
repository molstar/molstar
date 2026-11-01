/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { ValueCell } from '@molstar/core/util';
import { affineRowsCell, packAffineRows } from '../affine-rows.js';

/** `m` applied to `p` by the rows of `packAffineRows`, as the shaders' `affineTransform` does */
function applyRows(rows: Float32Array, o: number, p: Vec3) {
  const r = (i: number) =>
    rows[o + i * 4] * p[0] + rows[o + i * 4 + 1] * p[1] + rows[o + i * 4 + 2] * p[2] + rows[o + i * 4 + 3];
  return Vec3.create(r(0), r(1), r(2));
}

describe('affine rows', () => {
  const a = Mat4.mul(
    Mat4(),
    Mat4.fromTranslation(Mat4(), Vec3.create(1, -2, 3)),
    Mat4.fromRotation(Mat4(), 0.7, Vec3.create(1, 1, 0)),
  );
  const b = Mat4.fromScaling(Mat4(), Vec3.create(2, 3, 4));
  const transforms = new Float32Array(32);
  Mat4.toArray(a, transforms, 0);
  Mat4.toArray(b, transforms, 16);

  it('packs the first three rows of each matrix', () => {
    const rows = packAffineRows(transforms, 2, new Float32Array(24));
    const p = Vec3.create(0.5, -1.5, 2);
    for (const [m, o] of [
      [a, 0],
      [b, 12],
    ] as const) {
      const expected = Vec3.transformMat4(Vec3(), p, m);
      const actual = applyRows(rows, o, p);
      for (let i = 0; i < 3; ++i) expect(actual[i]).toBeCloseTo(expected[i], 5);
    }
  });

  it('repacks when the cell changes, passes other values through', () => {
    const cell = ValueCell.create(transforms);
    const rows = affineRowsCell(cell);
    expect(rows.ref.value.length).toBe(24);
    expect(rows.ref.version).toBe(cell.ref.version);
    const first = rows.ref.value[3];
    expect(first).toBeCloseTo(1);

    const moved = transforms.slice();
    moved[12] = 5; // the first matrix's x translation
    ValueCell.update(cell, moved);
    expect(rows.ref.version).toBe(cell.ref.version);
    expect(rows.ref.value[3]).toBe(5);

    const other = { buffer: 'written on the GPU' };
    ValueCell.update(cell, other as unknown as Float32Array);
    expect(rows.ref.value).toBe(other);
  });
});
