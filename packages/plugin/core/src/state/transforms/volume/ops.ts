/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Yakov Pechersky <ffxen158@gmail.com>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import type { PluginContext } from '@molstar/plugin/context';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { StateSelection, type StateTransform, StateTransformer } from '@molstar/core/state';
import { Task, type RuntimeContext } from '@molstar/core/task';
import { Volume, Grid } from '@molstar/model/model/volume';
import {
  TransformParam,
  transformParamsNeedCentroid,
  getTransformFromParams,
} from '@molstar/plugin/state/transforms/helpers';
import { Vec3, Mat4 } from '@molstar/core/math/linear-algebra';

export { AssignColorVolume };
type AssignColorVolume = typeof AssignColorVolume;
const AssignColorVolume = PluginStateTransform.BuiltIn({
  name: 'assign-color-volume',
  display: { name: 'Assign Color Volume', description: 'Assigns another volume to be available for coloring.' },
  from: SO.Volume.Data,
  to: SO.Volume.Data,
  isDecorator: true,
  params(a, plugin: PluginContext) {
    if (!a) return { ref: PD.Text() };
    const cells = plugin.state.data.select(
      StateSelection.Generators.root
        .subtree()
        .ofType(SO.Volume.Data)
        .filter((cell) => !!cell.obj && !cell.obj?.data.colorVolume && cell.obj !== a),
    );
    if (cells.length === 0) return { ref: PD.Text('', { isHidden: true }) };
    return {
      ref: PD.Select(
        cells[0].transform.ref,
        cells.map((c) => [c.transform.ref, c.obj!.label]),
      ),
    };
  },
})({
  apply({ a, params, dependencies }) {
    return Task.create('Assign Color Volume', async (ctx) => {
      if (!dependencies || !dependencies[params.ref]) {
        throw new Error('Dependency not available.');
      }
      const colorVolume = dependencies[params.ref].data as Volume;
      const volume: Volume = {
        ...a.data,
        colorVolume,
        _localPropertyData: Object.create(null),
      };
      const props = { label: a.label, description: 'Volume + Colors' };
      return new SO.Volume.Data(volume, props);
    });
  },
  getDependencies: ({ ref }) => (ref ? [ref as StateTransform.Ref] : []),
});

export { VolumeTransform };
type VolumeTransform = typeof VolumeTransform;
const VolumeTransform = PluginStateTransform.BuiltIn({
  name: 'volume-transform',
  display: { name: 'Transform Volume' },
  isDecorator: true,
  from: SO.Volume.Data,
  to: SO.Volume.Data,
  params: {
    transform: TransformParam,
  },
})({
  canAutoUpdate() {
    return false;
  },
  apply({ a, params }) {
    // similar to StateTransforms.Model.TransformStructureConformation;
    const center = transformParamsNeedCentroid(params.transform)
      ? Grid.getBoundingSphere(a.data.grid).center
      : Vec3.unit;
    const transform = getTransformFromParams(params.transform, center);
    const gridTransform = {
      kind: 'matrix' as const,
      matrix: Mat4.mul(Mat4(), transform, Grid.getGridToCartesianTransform(a.data.grid)),
    };
    return new SO.Volume.Data(
      {
        ...a.data,
        grid: {
          ...a.data.grid,
          transform: gridTransform,
        },
        _localPropertyData: Object.create(null),
      },
      {
        label: a.label,
        description: `${a.description} [Transformed]`,
      },
    );
  },
});

export { VolumeInstances };
type VolumeInstances = typeof VolumeInstances;
const VolumeInstances = PluginStateTransform.BuiltIn({
  name: 'volume-instances',
  display: { name: 'Volume Instances' },
  isDecorator: true,
  from: SO.Volume.Data,
  to: SO.Volume.Data,
  params: {
    mode: PD.Select('transforms', PD.arrayToOptions(['transforms', 'periodicRange'])),
    transforms: PD.ObjectList({ transform: TransformParam }, () => 'Transform', {
      hideIf: (c) => c.mode !== 'transforms',
    }),
    periodicRange: PD.Group(
      {
        min: PD.Vec3(
          Vec3.create(0, 0, 0),
          { min: -10, max: 10, step: 1 },
          { label: 'Min', description: 'Inclusive lower bound of the translation range (x, y, z)' },
        ),
        max: PD.Vec3(
          Vec3.create(1, 1, 1),
          { min: -10, max: 10, step: 1 },
          { label: 'Max', description: 'Exclusive upper bound of the translation range (x, y, z)' },
        ),
      },
      { hideIf: (c) => c.mode !== 'periodicRange', isFlat: true },
    ),
  },
})({
  canAutoUpdate() {
    return true;
  },
  apply({ a, params }) {
    let instances: { transform: Mat4 }[] = [];
    if (params.mode === 'transforms') {
      const center = params.transforms.some((t) => transformParamsNeedCentroid(t.transform))
        ? Grid.getBoundingSphere(a.data.grid).center
        : Vec3.unit;
      instances = params.transforms.map((t) => ({ transform: getTransformFromParams(t.transform, center) }));
    } else if (params.mode === 'periodicRange' && Volume.isPeriodic(a.data)) {
      const dims = a.data.grid.cells.space.dimensions;
      const gridToCartn = Grid.getGridToCartesianTransform(a.data.grid);
      const { min, max } = params.periodicRange;

      const [minA, minB, minC] = min;
      const [maxA, maxB, maxC] = max;

      const t = Vec3();
      for (let ia = minA; ia < maxA; ia++) {
        for (let ib = minB; ib < maxB; ib++) {
          for (let ic = minC; ic < maxC; ic++) {
            Vec3.set(t, dims[0] * ia, dims[1] * ib, dims[2] * ic);
            Vec3.transformMat4(t, t, gridToCartn);
            instances.push({ transform: Mat4.fromTranslation(Mat4(), t) });
          }
        }
      }
    }

    if (!instances.length) {
      return a;
    }
    return new SO.Volume.Data(
      {
        ...a.data,
        instances,
        _localPropertyData: Object.create(null),
      },
      {
        label: a.label,
        description: `${a.description} [Instanced]`,
      },
    );
  },
});

export { CustomVolumeProperties };
type CustomVolumeProperties = typeof CustomVolumeProperties;
const CustomVolumeProperties = PluginStateTransform.BuiltIn({
  name: 'custom-volume-properties',
  display: { name: 'Custom Volume Properties' },
  isDecorator: true,
  from: SO.Volume.Data,
  to: SO.Volume.Data,
  params: (a, ctx: PluginContext) => {
    return ctx.customVolumeProperties.getParams(a?.data);
  },
})({
  apply({ a, params }, ctx: PluginContext) {
    return Task.create('Custom Volume Props', async (taskCtx) => {
      await attachVolumeProps(a.data, ctx, taskCtx, params);
      return new SO.Volume.Data(a.data, { label: a.label, description: a.description });
    });
  },
  update({ a, b, oldParams, newParams }, ctx: PluginContext) {
    return Task.create('Custom Volume Props', async (taskCtx) => {
      b.data = a.data;
      b.label = a.label;
      b.description = a.description;
      for (const name of oldParams.autoAttach) {
        const property = ctx.customVolumeProperties.get(name);
        if (!property) continue;
        a.data.customProperties.reference(property.descriptor, false);
      }
      await attachVolumeProps(a.data, ctx, taskCtx, newParams);
      return StateTransformer.UpdateResult.Updated;
    });
  },
  dispose({ b }) {
    b?.data.customProperties.dispose();
  },
});

async function attachVolumeProps(
  volume: Volume,
  ctx: PluginContext,
  taskCtx: RuntimeContext,
  params: ReturnType<CustomVolumeProperties['createDefaultParams']>,
) {
  const propertyCtx = { runtime: taskCtx, assetManager: ctx.managers.asset, errorContext: ctx.errorContext };
  const { autoAttach, properties } = params;
  for (const name of Object.keys(properties)) {
    const property = ctx.customVolumeProperties.get(name)!;
    const props = properties[name];
    if (autoAttach.includes(name) || property.isHidden) {
      try {
        await property.attach(propertyCtx, volume, props, true);
      } catch (e) {
        ctx.log.warn(`Error attaching volume prop '${name}': ${e}`);
      }
    } else {
      property.set(volume, props);
    }
  }
}
