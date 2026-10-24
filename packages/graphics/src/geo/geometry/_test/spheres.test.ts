/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Spheres } from '../spheres/spheres.js';

describe('Spheres', () => {
  it('fill their position-group texture once, in createValues (perf-audit-2 08/P11)', () => {
    const spheres = Spheres.create(Float32Array.of(0, 0, 0, 1, 2, 3, 4, 5, 6), Float32Array.of(0, 1, 2), 3);
    const { positionGroup } = spheres.shaderData;
    // not filled by `create` with the default LOD, then again by `createValues`
    expect(positionGroup.ref.version).toBe(0);
    const values = Spheres.Utils.createValuesSimple(spheres, PD.getDefaultValues(Spheres.Params), ColorNames.orange, 1);
    expect(positionGroup.ref.version).toBe(1);
    // xyz and group per sphere
    expect(Array.from(values.tPositionGroup.ref.value.array.subarray(0, 12))).toEqual([
      0, 0, 0, 0, 1, 2, 3, 1, 4, 5, 6, 2,
    ]);
  });
});
