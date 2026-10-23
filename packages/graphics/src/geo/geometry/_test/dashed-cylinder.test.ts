/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { addFixedCountDashedCylinder, type BasicCylinderProps } from '../mesh/builder/cylinder.js';
import { MeshBuilder } from '../mesh/mesh-builder.js';

describe('addFixedCountDashedCylinder', () => {
  it("leaves the caller's props as they were (perf-audit-2 08/B7)", () => {
    const props: BasicCylinderProps = { radiusTop: 0.2, radiusBottom: 0.2, topCap: true, bottomCap: true };
    const state = MeshBuilder.createState();
    // an odd count without a stub cap: the last dash has no top cap
    addFixedCountDashedCylinder(state, Vec3.create(0, 0, 0), Vec3.create(0, 0, 4), 1, 3, false, props);
    const first = state.vertices.elementCount;
    expect(props.topCap).toBe(true);

    // the same props again (MirroredAromatic's second row): the same dashes
    addFixedCountDashedCylinder(state, Vec3.create(1, 0, 0), Vec3.create(1, 0, 4), 1, 3, false, props);
    expect(state.vertices.elementCount).toBe(2 * first);
  });
});
