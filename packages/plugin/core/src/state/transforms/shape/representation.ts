/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import type { PluginContext } from '@molstar/plugin/context';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';
import { Task } from '@molstar/core/task';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ShapeRepresentation } from '@molstar/graphics/repr/shape/representation';
import { StateTransformer } from '@molstar/core/state';

export { ShapeRepresentation3D };
type ShapeRepresentation3D = typeof ShapeRepresentation3D;
const ShapeRepresentation3D = PluginStateTransform.BuiltIn({
  name: 'shape-representation-3d',
  display: '3D Representation',
  from: SO.Shape.Provider,
  to: SO.Shape.Representation3D,
  params: (a, ctx: PluginContext) => {
    return a ? a.data.params : BaseGeometry.Params;
  },
})({
  canAutoUpdate() {
    return true;
  },
  apply({ a, params }) {
    return Task.create('Shape Representation', async (ctx) => {
      const props = { ...PD.getDefaultValues(a.data.params), ...params };
      const repr = ShapeRepresentation(a.data.getShape, a.data.geometryUtils);
      await repr.createOrUpdate(props, a.data.data).runInContext(ctx);
      return new SO.Shape.Representation3D({ repr, sourceData: a.data }, { label: a.data.label });
    });
  },
  update({ a, b, newParams }) {
    return Task.create('Shape Representation', async (ctx) => {
      const props = { ...b.data.repr.props, ...newParams };
      await b.data.repr.createOrUpdate(props, a.data.data).runInContext(ctx);
      b.data.sourceData = a.data;
      return StateTransformer.UpdateResult.Updated;
    });
  },
});
