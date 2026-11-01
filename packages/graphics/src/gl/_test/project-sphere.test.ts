/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Mat4, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { projectSphere } from '../project-sphere.js';

describe('projectSphere', () => {
  it('bounds a sphere above and right of the center near the top right, y from the top', () => {
    const projection = Mat4.perspective(Mat4(), -1, 1, 1, -1, 1, 100);
    const rect = projectSphere(Vec4(), Vec3.create(5, 5, -20), 1, projection);
    // the center at (0.625, 0.375): x right, y from the top
    expect(rect[0]).toBeLessThan(0.625);
    expect(rect[2]).toBeGreaterThan(0.625);
    expect(rect[1]).toBeLessThan(0.375);
    expect(rect[3]).toBeGreaterThan(0.375);
    expect(rect[3] - rect[1]).toBeLessThan(0.1);
  });

  it('bounds the same sphere alike for an orthographic projection', () => {
    const projection = Mat4.ortho(Mat4(), -20, 20, 20, -20, 1, 100);
    const rect = projectSphere(Vec4(), Vec3.create(5, 5, -20), 1, projection);
    expect(Array.from(rect).map((v) => +v.toFixed(3))).toEqual([0.6, 0.35, 0.65, 0.4]);
  });
});
