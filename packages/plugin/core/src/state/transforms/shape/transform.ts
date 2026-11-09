/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import type { RuntimeContext } from '@molstar/core/task';
import { Task } from '@molstar/core/task';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ShapeProvider } from '@molstar/graphics/geo/shape/provider';
import { calculateTransformBoundingSphere } from '@molstar/graphics/gl/renderable/util';
import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import {
  getTransformFromParams,
  TransformParam,
  transformParamsNeedCentroid,
} from '@molstar/plugin/state/transforms/helpers';

/**
 * The center of a shape as its provider makes it with its default props: of the bounds of all its instances, as a
 * representation of it bounds them drawn (`calculateTransformBoundingSphere`)
 */
export async function getShapeProviderCenter(ctx: RuntimeContext, provider: ShapeProvider<any, any, any>) {
  const shape = await provider.getShape(ctx, provider.data, PD.getDefaultValues(provider.params));
  const sphere = shape.geometry.boundingSphere;
  const { transforms } = shape;
  if (transforms.length === 0) return Vec3.clone(sphere.center);
  const packed = new Float32Array(transforms.length * 16);
  for (let i = 0, il = transforms.length; i < il; ++i) Mat4.toArray(transforms[i], packed, i * 16);
  return Vec3.clone(calculateTransformBoundingSphere(sphere, packed, transforms.length, 0).center);
}

/** the transform a provider made by `transformShapeProvider` moves its shapes by, composed over nested ones */
const ProviderTransforms = new WeakMap<ShapeProvider<any, any, any>, Mat4>();

/** The transform `provider` moves its shapes by (`transformShapeProvider`), undefined if not moved */
export function getShapeProviderTransform(provider: ShapeProvider<any, any, any>): Mat4 | undefined {
  return ProviderTransforms.get(provider);
}

/** Whether `a` and `b` move their shapes alike (`getShapeProviderTransform`): neither moved, or by the same transform */
export function isSameShapeProviderTransform(a: ShapeProvider<any, any, any>, b: ShapeProvider<any, any, any>) {
  const ta = ProviderTransforms.get(a);
  const tb = ProviderTransforms.get(b);
  if (!ta || !tb) return ta === tb;
  return Mat4.areEqual(ta, tb, 0);
}

/** `provider` with every shape it makes moved by `transform` (its instances' transforms premultiplied) */
export function transformShapeProvider(
  provider: ShapeProvider<any, any, any>,
  transform: Mat4,
): ShapeProvider<any, any, any> {
  const { getShape } = provider;
  const moved: ShapeProvider<any, any, any> = {
    ...provider,
    getShape: async (ctx, data, props, shape) => {
      const made = await getShape(ctx, data, props, shape);
      return { ...made, transforms: made.transforms.map((t) => Mat4.mul(Mat4(), transform, t)) };
    },
  };
  const inner = ProviderTransforms.get(provider);
  ProviderTransforms.set(moved, inner ? Mat4.mul(Mat4(), transform, inner) : Mat4.copy(Mat4(), transform));
  return moved;
}

export { TransformShapeProvider };
type TransformShapeProvider = typeof TransformShapeProvider;
const TransformShapeProvider = PluginStateTransform.BuiltIn({
  name: 'transform-shape-provider',
  display: { name: 'Transform Shape', description: 'Move and turn a shape as a whole.' },
  isDecorator: true,
  from: SO.Shape.Provider,
  to: SO.Shape.Provider,
  params: {
    transform: TransformParam,
  },
})({
  canAutoUpdate({ newParams }) {
    return newParams.transform.name !== 'matrix';
  },
  apply({ a, params }) {
    return Task.create('Transform Shape', async (ctx) => {
      // as TransformStructureConformation: a turn about the shape's center
      const center = transformParamsNeedCentroid(params.transform)
        ? await getShapeProviderCenter(ctx, a.data)
        : Vec3.unit;
      const transform = getTransformFromParams(params.transform, center);
      return new SO.Shape.Provider(transformShapeProvider(a.data, transform), {
        label: a.label,
        description: `${a.description ?? ''} [Transformed]`.trim(),
      });
    });
  },
});
