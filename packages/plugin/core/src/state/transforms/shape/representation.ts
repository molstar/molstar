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
import { OrderedSet, SortedArray } from '@molstar/core/data/int';
import type { ShapeProvider } from '@molstar/graphics/geo/shape/provider';
import { Transparency } from '@molstar/graphics/theme/transparency';
import { Shape, ShapeGroup } from '@molstar/model/model/shape';
import { getShapeSignatures, registerTransformSignatures } from '@molstar/plugin/state/precompile';

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
  apply({ a, params }, plugin: PluginContext) {
    return Task.create('Shape Representation', async (ctx) => {
      const props = { ...PD.getDefaultValues(a.data.params), ...params };
      const repr = ShapeRepresentation(a.data.getShape, a.data.geometryUtils, { graphics: plugin.canvas3d?.graphics });
      await repr.createOrUpdate(props, a.data.data).runInContext(ctx);
      applyTransparency(repr, a.data);
      return new SO.Shape.Representation3D({ repr, sourceData: a.data }, { label: a.data.label });
    });
  },
  update({ a, b, newParams }) {
    return Task.create('Shape Representation', async (ctx) => {
      const props = { ...b.data.repr.props, ...newParams };
      await b.data.repr.createOrUpdate(props, a.data.data).runInContext(ctx);
      applyTransparency(b.data.repr, a.data);
      b.data.sourceData = a.data;
      return StateTransformer.UpdateResult.Updated;
    });
  },
});

/** The provider's per-group transparency (`ShapeProvider.transparency`) as a layer on the shape's groups */
/** The representations whose transparency is their provider's (what a provider without one takes away again) */
const ProviderTransparency = new WeakSet<object>();

function applyTransparency(repr: ShapeRepresentation<any, any, any>, provider: ShapeProvider<any, any, any>) {
  const of = provider.transparency;
  if (!of) {
    // a provider that had one and has not anymore: its layers go (another's transparency is left as it is)
    if (ProviderTransparency.delete(repr)) {
      repr.setState({ transparency: Transparency.Empty });
    }
    return;
  }
  ProviderTransparency.add(repr);
  const all = repr.getAllLoci()[0];
  if (!Shape.isLoci(all)) return;
  const { shape } = all;
  const byValue = new Map<number, number[]>();
  for (let g = 0; g < shape.groupCount; g++) {
    const value = of(g);
    if (!(value > 0)) continue;
    const ids = byValue.get(value);
    if (ids) ids.push(g);
    else byValue.set(value, [g]);
  }
  const instances = shape.transforms.length;
  const layers = [...byValue].map(([value, ids]) => {
    const set = OrderedSet.ofSortedArray(SortedArray.ofSortedArray(ids));
    const groups = Array.from({ length: instances }, (_, instance) => ({ ids: set, instance }));
    return { loci: ShapeGroup.Loci(shape, groups), value };
  });
  repr.setState({ transparency: { kind: 'group-loci', layers } as Transparency });
}

registerTransformSignatures(ShapeRepresentation3D, (params, ctx) => {
  const parent = ctx.tree.transforms.get(ctx.transform.parent);
  return getShapeSignatures(parent, params, ctx);
});
