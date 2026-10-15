/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ColorNames } from '@molstar/core/util/color/names';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import type { PluginContext } from '@molstar/plugin/context';
import { Task } from '@molstar/core/task';
import { ShapeRepresentation } from '@molstar/graphics/repr/shape/representation';
import type { Box3D } from '@molstar/core/math/geometry';
import type { Color } from '@molstar/core/util/color';
import { getBoxMesh } from '@molstar/plugin/state/transforms/shape/box';
import * as ShapeUtils from '@molstar/graphics/geo/shape/shape';
import { StateTransformer } from '@molstar/core/state';

export { StructureBoundingBox3D };
type StructureBoundingBox3D = typeof StructureBoundingBox3D;
const StructureBoundingBox3D = PluginStateTransform.BuiltIn({
  name: 'structure-bounding-box-3d',
  display: 'Bounding Box',
  from: SO.Molecule.Structure,
  to: SO.Shape.Representation3D,
  params: {
    radius: PD.Numeric(0.05, { min: 0.01, max: 4, step: 0.01 }, { isEssential: true }),
    color: PD.Color(ColorNames.red, { isEssential: true }),
    ...Mesh.Params,
  },
})({
  canAutoUpdate() {
    return true;
  },
  apply({ a, params }, plugin: PluginContext) {
    return Task.create('Bounding Box', async (ctx) => {
      const repr = ShapeRepresentation((_, data: { box: Box3D; radius: number; color: Color }, __, shape) => {
        const mesh = getBoxMesh(data.box, data.radius, shape?.geometry);
        return ShapeUtils.create(
          'Bouding Box',
          data,
          mesh,
          () => data.color,
          () => 1,
          () => 'Bounding Box',
        );
      }, Mesh.Utils);
      await repr
        .createOrUpdate(params, { box: a.data.boundary.box, radius: params.radius, color: params.color })
        .runInContext(ctx);
      return new SO.Shape.Representation3D({ repr, sourceData: a.data }, { label: `Bounding Box` });
    });
  },
  update({ a, b, oldParams, newParams }, plugin: PluginContext) {
    return Task.create('Bounding Box', async (ctx) => {
      await b.data.repr
        .createOrUpdate(newParams, { box: a.data.boundary.box, radius: newParams.radius, color: newParams.color })
        .runInContext(ctx);
      b.data.sourceData = a.data;
      return StateTransformer.UpdateResult.Updated;
    });
  },
});
