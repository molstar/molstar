/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { ValueCell } from '@molstar/core/util';
import { Color } from '@molstar/core/util/color';
import { Shape } from '@molstar/model/model/shape';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { MeshBuilder } from '@molstar/graphics/geo/geometry/mesh/mesh-builder';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { ShapeRepresentation } from '../representation.js';

/** A triangle whose vertices belong to two groups, drawn with varying groups */
function getVaryingShape() {
  const state = MeshBuilder.createState(16, 16);
  MeshBuilder.addTriangle(state, Vec3.create(0, 0, 0), Vec3.create(1, 0, 0), Vec3.create(0, 1, 0));
  const mesh = MeshBuilder.getMesh(state);
  mesh.groupBuffer.ref.value[2] = 1;
  ValueCell.update(mesh.groupBuffer, mesh.groupBuffer.ref.value);
  ValueCell.updateIfChanged(mesh.varyingGroup, true);
  return Shape.create(
    'triangle',
    {},
    mesh,
    () => Color(0xff0000),
    () => 1,
    () => '',
    2,
  );
}

const graphics = (flatInterpolation: boolean) =>
  ({ capabilities: { flatInterpolation } }) as unknown as GraphicsContext;

describe('shape representation, geometry relying on capabilities', () => {
  it('makes varying groups uniform without flat interpolation', async () => {
    const shape = getVaryingShape();
    const repr = ShapeRepresentation((_, s: Shape<Mesh>) => s, Mesh.Utils, { graphics: graphics(false) });
    await repr.createOrUpdate({}, shape).run();
    const values = repr.renderObjects[0].values;
    expect(values.dVaryingGroup.ref.value).toBe(false);
    // the triangle split at the group border
    expect(values.uVertexCount.ref.value).toBeGreaterThan(3);
    expect(values.drawCount.ref.value).toBeGreaterThan(3);
  });

  it('keeps varying groups with flat interpolation, changes them after `setGraphics`', async () => {
    const shape = getVaryingShape();
    const repr = ShapeRepresentation((_, s: Shape<Mesh>) => s, Mesh.Utils, { graphics: graphics(true) });
    await repr.createOrUpdate({}, shape).run();
    const values = repr.renderObjects[0].values;
    expect(values.dVaryingGroup.ref.value).toBe(true);
    expect(values.uVertexCount.ref.value).toBe(3);

    repr.setGraphics!(graphics(false));
    await repr.createOrUpdate().run();
    expect(values.dVaryingGroup.ref.value).toBe(false);
    expect(values.uVertexCount.ref.value).toBeGreaterThan(3);
  });
});
