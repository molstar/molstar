/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CylindersBuilder } from '@molstar/graphics/geo/geometry/cylinders/cylinders-builder';
import { perVertexCylinders, perVertexLines } from '../../../../../geo/geometry/_test/item-vertices.js';
import type { VisualContext } from '../../../../visual.js';
import {
  createLinkCylinderImpostors,
  createLinkLines,
  LinkCylinderParams,
  LinkLineParams,
  LinkStyle,
  type LinkBuilderProps,
} from '../link.js';

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
    // three texels a cylinder, in one array allocated once
    expect(cylinders.cylinderData.ref.value.array.length).toBe(linkCount * 12);
    const sized = f32.mock.calls.map((c) => c[0]).filter((n) => typeof n === 'number' && n > 0);
    expect(sized).toEqual([linkCount * 12]);
  });
});

describe('disk links', () => {
  // two edges: a 10 Å link along x from the origin, and one along y from (0, 5, 0)
  const ends = [
    [0, 0, 0, 10, 0, 0],
    [0, 5, 0, 0, 15, 0],
  ];
  const builderProps: LinkBuilderProps = {
    linkCount: 2,
    position: (a, b, i) => {
      Vec3.set(a, ends[i][0], ends[i][1], ends[i][2]);
      Vec3.set(b, ends[i][3], ends[i][4], ends[i][5]);
    },
    radius: () => 0.2,
    style: () => LinkStyle.Disk,
  };
  // the mesh draws a disk from 0.475 to 0.5 of the way from a to b
  const expected = ends.map(([x, y, z, u, v, w]) => [
    [x + 0.475 * (u - x), y + 0.475 * (v - y), z + 0.475 * (w - z)],
    [x + 0.5 * (u - x), y + 0.5 * (v - y), z + 0.5 * (w - z)],
  ]);
  const endsOf = (start: Float32Array, end: Float32Array, stride: number) =>
    [0, 1].map((k) => [
      Array.from(start.subarray(k * stride, k * stride + 3)),
      Array.from(end.subarray(k * stride, k * stride + 3)),
    ]);

  it('impostors start and end on their own edge, as the mesh (perf-audit-2 07/B6)', () => {
    const { cylinders } = createLinkCylinderImpostors(
      {} as VisualContext,
      builderProps,
      PD.getDefaultValues(LinkCylinderParams),
    );
    const { start, end } = perVertexCylinders(cylinders.cylinderData.ref.value.array, cylinders.cylinderCount);
    const got = endsOf(start, end, 6 * 3);
    for (let k = 0; k < 2; k++)
      for (let e = 0; e < 2; e++) for (let c = 0; c < 3; c++) expect(got[k][e][c]).toBeCloseTo(expected[k][e][c], 5);
  });

  it('lines likewise', () => {
    const { lines } = createLinkLines({} as VisualContext, builderProps, PD.getDefaultValues(LinkLineParams));
    const { start, end } = perVertexLines(lines.lineData.ref.value.array, lines.lineCount);
    const got = endsOf(start, end, 4 * 3);
    for (let k = 0; k < 2; k++)
      for (let e = 0; e < 2; e++) for (let c = 0; c < 3; c++) expect(got[k][e][c]).toBeCloseTo(expected[k][e][c], 5);
  });
});
