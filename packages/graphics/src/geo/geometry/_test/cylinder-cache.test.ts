/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { MeshBuilder } from '../mesh/mesh-builder.js';
import { addCylinder, MaxCachedCylinders } from '../mesh/builder/cylinder.js';

// counts the cylinder primitives built
const built = vi.hoisted(() => ({ count: 0 }));
vi.mock('../../primitive/cylinder.js', async (importOriginal) => {
  const original = await importOriginal<typeof import('../../primitive/cylinder.js')>();
  return { ...original, Cylinder: (props: any) => (built.count++, original.Cylinder(props)) };
});

describe('cylinder primitive cache', () => {
  it('keeps a bounded number of primitives (perf-audit-2 08/P9)', () => {
    const state = MeshBuilder.createState(1024, 1024);
    const add = (radius: number) =>
      addCylinder(state, Vec3.create(0, 0, 0), Vec3.create(0, 1, 0), 1, {
        radiusTop: radius,
        radiusBottom: radius,
        radialSegments: 8,
      });

    built.count = 0;
    add(0.123);
    add(0.123);
    expect(built.count).toBe(1);

    // continuous radii, as a gizmo drag makes: the first one is evicted
    for (let i = 1; i <= MaxCachedCylinders; i++) add(1 + i * 0.001);
    expect(built.count).toBe(1 + MaxCachedCylinders);
    add(0.123);
    expect(built.count).toBe(2 + MaxCachedCylinders);
  });
});
