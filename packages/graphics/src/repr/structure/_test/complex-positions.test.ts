/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { RuntimeContext } from '@molstar/core/task';
import { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { EmptyLoci } from '@molstar/model/model/loci';
import { NullLocation } from '@molstar/model/model/location';
import type { PositionLocation } from '@molstar/model/model/position-location';
import { Structure } from '@molstar/model/model/structure';
import type { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { MeshBuilder } from '@molstar/graphics/geo/geometry/mesh/mesh-builder';
import { LocationIterator } from '@molstar/graphics/geo/util/location-iterator';
import type { ColorTheme } from '@molstar/graphics/theme/color';
import { SizeTheme } from '@molstar/graphics/theme/size';
import type { Theme } from '@molstar/graphics/theme/theme';
import { ComplexMeshParams, ComplexMeshVisual } from '../complex-visual.js';

/** A color from the vertex's x */
const colorOfX = (x: number) => Color(Math.round((x + 100) * 50) & 0xffffff);

const Params = { ...ComplexMeshParams, count: PD.Numeric(1, { min: 1, max: 10, step: 1 }) };

/** One triangle per `count`, at x = 10 i, in a mesh made anew each time (no reuse of the old one) */
function TriangleVisual(built: Mesh[]) {
  return ComplexMeshVisual<typeof Params>(
    {
      defaultProps: PD.getDefaultValues(Params),
      createGeometry: (_ctx, _structure, _theme, props) => {
        const state = MeshBuilder.createState(16, 16);
        for (let i = 0; i < props.count; i++) {
          MeshBuilder.addTriangle(
            state,
            Vec3.create(10 * i, 0, 0),
            Vec3.create(10 * i + 1, 0, 0),
            Vec3.create(10 * i, 1, 0),
          );
        }
        const mesh = MeshBuilder.getMesh(state);
        built.push(mesh);
        return mesh;
      },
      createLocationIterator: () => LocationIterator(1, 1, 1, () => NullLocation),
      getLoci: () => EmptyLoci,
      eachLocation: () => false,
      setUpdateState: (state, newProps, currentProps) => {
        state.createGeometry = newProps.count !== currentProps.count;
      },
    },
    0,
  );
}

describe('complex visual positions', () => {
  it('a vertex color theme reads the new geometry after it changes (perf-audit-2 07/B7)', async () => {
    const theme: Theme = {
      color: {
        factory: () => theme.color,
        granularity: 'vertex',
        color: (l) => colorOfX((l as PositionLocation).position[0]),
        props: {},
        description: '',
      } as ColorTheme<any>,
      size: SizeTheme.Empty,
    };
    const built: Mesh[] = [];
    const visual = TriangleVisual(built);
    const ctx = { runtime: RuntimeContext.Synchronous };
    await visual.createOrUpdate(ctx, theme, { count: 1 }, Structure.Empty);
    await visual.createOrUpdate(ctx, theme, { count: 3 }, Structure.Empty);

    const { values } = visual.renderObject!;
    const vertexCount = values.uVertexCount.ref.value;
    expect(vertexCount).toBe(9);
    // the positions of the geometry the update built
    const positions = built[built.length - 1].vertexBuffer.ref.value;
    const colors = values.tColor.ref.value.array as Uint8Array;
    for (let v = 0; v < vertexCount; v++) {
      expect([colors[v * 4], colors[v * 4 + 1], colors[v * 4 + 2]], `vertex ${v}`).toEqual(
        Color.toRgb(colorOfX(positions[v * 3])),
      );
    }
  });
});
