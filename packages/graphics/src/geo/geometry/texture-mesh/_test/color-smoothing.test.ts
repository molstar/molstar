/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Box3D } from '@molstar/core/math/geometry';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { getFittingResolution, getGridDim } from '../color-smoothing.js';

/** A square atlas of all cells, as a stand-in for the layouts */
const atlasSide = (gridDim: Vec3) => Math.ceil(Math.sqrt(gridDim[0] * gridDim[1] * gridDim[2]));

describe('color smoothing grid', () => {
  const box = Box3D.create(Vec3.create(0, 0, 0), Vec3.create(500, 500, 500));

  it('keeps the resolution when the atlas fits', () => {
    expect(getFittingResolution(box, 2, 8192, atlasSide)).toBe(2);
  });

  it('coarsens the resolution until the atlas fits', () => {
    const resolution = getFittingResolution(box, 1, 8192, atlasSide);
    expect(resolution).toBeGreaterThan(1);
    expect(atlasSide(getGridDim(box, resolution))).toBeLessThanOrEqual(8192);
    // the step before did not fit
    expect(atlasSide(getGridDim(box, resolution / 1.1))).toBeGreaterThan(8192);
  });
});
