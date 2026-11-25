/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { expect } from 'vitest';
import { Quat } from '@molstar/core/math/linear-algebra';
import { GaussianSplat, type GaussianSplatData } from '../data.js';

/** test helpers shared by the gaussian splat format tests, not a test file itself */

/** The splats of `data`, each decoded */
export function decode(data: GaussianSplatData) {
  return Array.from({ length: data.count }, (_, i) => {
    const s = GaussianSplat();
    data.getSplat(i, s);
    return s;
  });
}

/** The same rotation: `q` or `-q`, both normalized */
export function expectRotation(s: GaussianSplat, expected: Quat, digits = 3) {
  const a = Quat.normalize(Quat(), Quat.create(s.qx, s.qy, s.qz, s.qw));
  const e = Quat.normalize(Quat(), expected);
  const sign = Quat.dot(a, e) < 0 ? -1 : 1;
  for (let k = 0; k < 4; ++k) expect(a[k] * sign).toBeCloseTo(e[k], digits);
}
