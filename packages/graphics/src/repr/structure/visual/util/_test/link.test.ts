/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CylindersBuilder } from '@molstar/graphics/geo/geometry/cylinders/cylinders-builder';
import type { VisualContext } from '../../../../visual.js';
import { createLinkCylinderImpostors, LinkCylinderParams, type LinkBuilderProps } from '../link.js';

describe('link cylinder impostors', () => {
  it('allocate buffers of the given cylinder count, kept without a copy (perf-audit-2 07/P11)', () => {
    const linkCount = 100;
    const builderProps: LinkBuilderProps = {
      linkCount,
      position: (a, b, i) => {
        Vec3.set(a, i, 0, 0);
        Vec3.set(b, i, 1, 0);
      },
      radius: () => 0.2,
      // one impostor per solid link
      cylinderCount: () => linkCount,
    };
    const create = vi.spyOn(CylindersBuilder, 'create');
    const f32 = vi.spyOn(globalThis, 'Float32Array');
    const { cylinders } = createLinkCylinderImpostors(
      {} as VisualContext,
      builderProps,
      PD.getDefaultValues(LinkCylinderParams),
    );
    expect(create.mock.calls[0][0]).toBe(linkCount);
    expect(cylinders.cylinderCount).toBe(linkCount);
    // six vertices a cylinder, one chunk per attribute
    expect(cylinders.groupBuffer.ref.value.length).toBe(linkCount * 6);
    const sized = f32.mock.calls.map((c) => c[0]).filter((n) => typeof n === 'number' && n > 0);
    expect(Math.max(...sized)).toBe(linkCount * 6 * 3);
    // starts, ends, groups, scales, caps, color modes once each, and the vertex mappings
    expect(sized.length).toBe(7);
  });
});
