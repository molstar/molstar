/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '@molstar/core/math/linear-algebra';
import { Particle } from '@molstar/model/model/particles/particle-list';
import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import {
  getTransformFromParams,
  TransformParam,
  transformParamsNeedCentroid,
} from '@molstar/plugin/state/transforms/helpers';

export { TransformParticleList };
type TransformParticleList = typeof TransformParticleList;
const TransformParticleList = PluginStateTransform.BuiltIn({
  name: 'transform-particle-list',
  display: { name: 'Transform Particles', description: 'Move and turn a particle list as a whole.' },
  isDecorator: true,
  from: SO.Particle.List,
  to: SO.Particle.List,
  params: {
    transform: TransformParam,
  },
})({
  canAutoUpdate({ newParams }) {
    return newParams.transform.name !== 'matrix';
  },
  apply({ a, params }) {
    // as TransformStructureConformation: a turn about the list's center
    const center = transformParamsNeedCentroid(params.transform)
      ? Particle.getBoundary(a.data).sphere.center
      : Vec3.unit;
    const transform = getTransformFromParams(params.transform, center);
    return new SO.Particle.List(Particle.transform(a.data, transform), {
      label: a.label,
      description: `${a.description ?? ''} [Transformed]`.trim(),
    });
  },
});
