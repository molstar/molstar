/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type StateTransform, StateTransformer } from '@molstar/core/state';
import { Task, type RuntimeContext } from '@molstar/core/task';
import { type Coordinates, Structure, Model, type Frame, type Trajectory } from '@molstar/model/model/structure';
import {
  PeriodicParams,
  periodicFrame,
  periodicParamsFor,
  periodicTrajectory,
  type PeriodicProps,
} from '@molstar/plugin/state/helpers/periodic';
import { parseXsc } from '@molstar/io/reader/xsc/parser';
import { Asset } from '@molstar/core/util/assets';
import { coordinatesWithXsc } from '@molstar/model/formats/structure/namd-xsc';
import { getTrajectory } from './trajectory-helpers.js';
import { RootStructureDefinition } from '@molstar/plugin/state/helpers/root-structure';
import { BondProviderRegistry } from '@molstar/plugin/state/bonds/registry';
import type { PluginContext } from '@molstar/plugin/context';
import { deepEqual } from '@molstar/core/util';
import {
  TransformParam,
  transformParamsNeedCentroid,
  getTransformFromParams,
} from '@molstar/plugin/state/transforms/helpers';
import { Vec3 } from '@molstar/core/math/linear-algebra';

export { CoordinatesWithBoxes };
type CoordinatesWithBoxes = typeof CoordinatesWithBoxes;
/**
 * Coordinates with their frames' boxes from a NAMD extended system file (`.xsc`, `.xst`: a text data cell), which a DCD
 * lacks the origin and orientation of: each frame takes the box of its step.
 */
const CoordinatesWithBoxes = PluginStateTransform.BuiltIn({
  name: 'coordinates-with-boxes',
  display: {
    name: 'Coordinates with Boxes',
    description: "Give coordinates their frames' periodic boxes from a NAMD extended system file (.xsc, .xst).",
  },
  from: SO.Molecule.Coordinates,
  to: SO.Molecule.Coordinates,
  params: {
    /** The parsed file's cell (`ParseXsc`: a dropped or downloaded file) */
    boxesRef: PD.Text('', { isHidden: true }),
    /** Or its address, fetched here (MolViewSpec's `custom.molstar_boxes`: no data cell of its own) */
    url: PD.Text('', { isHidden: true }),
  },
})({
  getDependencies: ({ boxesRef }: { boxesRef: string }) => (boxesRef ? [boxesRef as StateTransform.Ref] : []),
  apply({ a, params, dependencies, cache }, plugin: PluginContext) {
    return Task.create('Coordinates with boxes', async (ctx) => {
      let text: string;
      if (params.boxesRef) {
        // a parsed file (`ParseXsc`: a dropped or downloaded file)
        return new SO.Molecule.Coordinates(
          coordinatesWithXsc(a.data, (dependencies![params.boxesRef] as SO.Format.Xsc).data),
          { label: a.label, description: 'Coordinates with boxes' },
        );
      } else {
        const asset = await plugin.managers.asset
          .resolve(Asset.getUrlAsset(plugin.managers.asset, params.url), 'string')
          .runInContext(ctx);
        (cache as { asset?: Asset.Wrapper }).asset = asset;
        text = asset.data as string;
      }
      const parsed = await parseXsc(text).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      return new SO.Molecule.Coordinates(coordinatesWithXsc(a.data, parsed.result), {
        label: a.label,
        description: 'Coordinates with boxes',
      });
    });
  },
  dispose({ cache }) {
    (cache as { asset?: Asset.Wrapper } | undefined)?.asset?.dispose();
  },
});

export { TrajectoryFromModelAndCoordinates };
type TrajectoryFromModelAndCoordinates = typeof TrajectoryFromModelAndCoordinates;
const TrajectoryFromModelAndCoordinates = PluginStateTransform.BuiltIn({
  name: 'trajectory-from-model-and-coordinates',
  display: {
    name: 'Trajectory from Topology & Coordinates',
    description: 'Create a trajectory from existing model/topology and coordinates.',
  },
  from: SO.Root,
  to: SO.Molecule.Trajectory,
  params: {
    modelRef: PD.Text('', { isHidden: true }),
    coordinatesRef: PD.Text('', { isHidden: true }),
  },
})({
  getDependencies: ({ modelRef, coordinatesRef }: { modelRef: string; coordinatesRef: string }) => {
    const deps: StateTransform.Ref[] = [];
    if (modelRef) deps.push(modelRef as StateTransform.Ref);
    if (coordinatesRef) deps.push(coordinatesRef as StateTransform.Ref);
    return deps;
  },
  apply({ params, dependencies }) {
    return Task.create('Create trajectory from model/topology and coordinates', async (ctx) => {
      const coordinates = dependencies![params.coordinatesRef].data as Coordinates;
      const trajectory = await getTrajectory(ctx, dependencies![params.modelRef], coordinates);
      const props = {
        label: 'Trajectory',
        description: `${trajectory.frameCount} model${trajectory.frameCount === 1 ? '' : 's'}`,
      };
      return new SO.Molecule.Trajectory(trajectory, props);
    });
  },
});

const plus1 = (v: number) => v + 1,
  minus1 = (v: number) => v - 1;

async function modelAtIndex(trajectory: Trajectory, index: number, ctx: RuntimeContext, periodic: PeriodicProps) {
  let modelIndex = Math.round(index) % trajectory.frameCount;
  if (modelIndex < 0) modelIndex += trajectory.frameCount;
  // the same model for the same recent frame (cached): an unchanged frame stays `Unchanged`
  const model = await periodicFrame(trajectory, modelIndex, periodic, ctx);
  const label = `Model ${modelIndex + 1}`;
  const description = trajectory.frameCount === 1 ? undefined : `of ${trajectory.frameCount}`;
  return { model, label, description };
}

export { ModelFromTrajectory };
type ModelFromTrajectory = typeof ModelFromTrajectory;
const ModelFromTrajectory = PluginStateTransform.BuiltIn({
  name: 'model-from-trajectory',
  display: { name: 'Molecular Model', description: 'Create a molecular model from specified index in a trajectory.' },
  from: SO.Molecule.Trajectory,
  to: SO.Molecule.Model,
  params: (a) => {
    if (!a) {
      return {
        modelIndex: PD.Numeric(0, {}, { description: 'Zero-based index of the model', immediateUpdate: true }),
        ...PeriodicParams,
      };
    }
    return {
      modelIndex: PD.Converted(
        plus1,
        minus1,
        PD.Numeric(
          1,
          { min: 1, max: a.data.frameCount, step: 1 },
          { description: 'Model Index', immediateUpdate: true },
        ),
      ),
      ...periodicParamsFor(a.data),
    };
  },
})({
  isApplicable: (a) => a.data.frameCount > 0,
  apply({ a, params, cache }) {
    return Task.create('Model from Trajectory', async (ctx) => {
      const { model, label, description } = await modelAtIndex(a.data, params.modelIndex, ctx, params);
      (cache as { trajectory?: Trajectory }).trajectory = a.data;
      return new SO.Molecule.Model(model, { label, description });
    });
  },
  update({ a, b, newParams, cache }) {
    // Another trajectory (the data read again): recreated, which disposes the old model's custom properties.
    if ((cache as { trajectory?: Trajectory }).trajectory !== a.data) return StateTransformer.UpdateResult.Recreate;
    // A new frame updates the cell in place: recreating it disposed the old frame's custom properties, which the
    // frames of the same topology share (`{ ...previous }`), on every frame of a playing trajectory.
    return Task.create('Model from Trajectory', async (ctx) => {
      const { model, label, description } = await modelAtIndex(a.data, newParams.modelIndex, ctx, newParams);
      if (model === b.data) return StateTransformer.UpdateResult.Unchanged;
      // Frames of the same topology share them (`{ ...previous }`): only the old frame's own assets go (a dynamic
      // property's, kept by the model id); a frame of another topology has its own properties, all of which go.
      const props = b.data.customProperties;
      if (model.customProperties !== props) props.dispose();
      else for (const desc of props.all) props.assets(desc, undefined, b.data.id);
      b.data = model;
      b.label = label;
      b.description = description;
      return StateTransformer.UpdateResult.Updated;
    });
  },
  interpolate(a, b, t) {
    const modelIndex = t >= 1 ? b.modelIndex : a.modelIndex + Math.floor((b.modelIndex - a.modelIndex + 1) * t);
    return { modelIndex, periodic: b.periodic, periodicCenter: b.periodicCenter };
  },
  dispose({ b }) {
    b?.data.customProperties.dispose();
  },
});

export { StructureFromTrajectory };
type StructureFromTrajectory = typeof StructureFromTrajectory;
const StructureFromTrajectory = PluginStateTransform.BuiltIn({
  name: 'structure-from-trajectory',
  display: { name: 'Structure from Trajectory', description: 'Create a molecular structure from a trajectory.' },
  from: SO.Molecule.Trajectory,
  to: SO.Molecule.Structure,
  params(a, plugin: PluginContext | undefined) {
    // the same choice as `StructureFromModel`'s, for the frames' representative model (bond-providers plan §4)
    const registry = plugin?.bondProviders ?? new BondProviderRegistry();
    // the frames as `ModelFromTrajectory` shows them: whole by default
    return { bonds: PD.Optional(registry.getParam(a?.data.representative)), ...periodicParamsFor(a?.data) };
  },
})({
  apply({ a, params }, plugin: PluginContext) {
    return Task.create('Build Structure', async (ctx) => {
      const bondProvider = plugin.bondProviders.resolve(params?.bonds, a.data.representative);
      const s = await Structure.ofTrajectory(periodicTrajectory(a.data, params ?? {}), ctx, { bondProvider });
      const props = { label: 'Ensemble', description: Structure.elementDescription(s) };
      return new SO.Molecule.Structure(s, props);
    });
  },
  dispose({ b }) {
    b?.data.customPropertyDescriptors.dispose();
  },
});

export { StructureFromModel };
type StructureFromModel = typeof StructureFromModel;
const StructureFromModel = PluginStateTransform.BuiltIn({
  name: 'structure-from-model',
  display: {
    name: 'Structure',
    description: 'Create a molecular structure (model, assembly, or symmetry) from the specified model.',
  },
  from: SO.Molecule.Model,
  to: SO.Molecule.Structure,
  params(a, plugin: PluginContext | undefined) {
    return RootStructureDefinition.getParams(a && a.data, undefined, plugin);
  },
})({
  canAutoUpdate({ oldParams, newParams }) {
    return RootStructureDefinition.canAutoUpdate(oldParams.type, newParams.type);
  },
  apply({ a, params }, plugin: PluginContext) {
    return Task.create('Build Structure', async (ctx) => {
      return RootStructureDefinition.create(plugin, ctx, a.data, params && params.type);
    });
  },
  update: ({ a, b, oldParams, newParams }) => {
    if (!deepEqual(oldParams, newParams)) return StateTransformer.UpdateResult.Recreate;
    if (b.data.model === a.data) return StateTransformer.UpdateResult.Unchanged;
    if (!Model.areHierarchiesEqual(a.data, b.data.model)) return StateTransformer.UpdateResult.Recreate;

    b.data = b.data.remapModel(a.data);

    return StateTransformer.UpdateResult.Updated;
  },
  dispose({ b }) {
    b?.data.customPropertyDescriptors.dispose();
  },
});

export { TransformStructureConformation };
type TransformStructureConformation = typeof TransformStructureConformation;
const TransformStructureConformation = PluginStateTransform.BuiltIn({
  name: 'transform-structure-conformation',
  display: { name: 'Transform Conformation' },
  isDecorator: true,
  from: SO.Molecule.Structure,
  to: SO.Molecule.Structure,
  params: {
    transform: TransformParam,
  },
})({
  canAutoUpdate({ newParams }) {
    return newParams.transform.name !== 'matrix';
  },
  apply({ a, params }) {
    const center = transformParamsNeedCentroid(params.transform) ? a.data.boundary.sphere.center : Vec3.unit;
    const transform = getTransformFromParams(params.transform, center);
    const s = Structure.transform(a.data, transform);
    return new SO.Molecule.Structure(s, { label: a.label, description: `${a.description} [Transformed]` });
  },
  dispose({ b }) {
    b?.data.customPropertyDescriptors.dispose();
  },
  // interpolate(src, tar, t) {
  //     // TODO: optimize
  //     const u = Mat4.fromRotation(Mat4(), Math.PI / 180 * src.angle, Vec3.normalize(Vec3(), src.axis));
  //     Mat4.setTranslation(u, src.translation);
  //     const v = Mat4.fromRotation(Mat4(), Math.PI / 180 * tar.angle, Vec3.normalize(Vec3(), tar.axis));
  //     Mat4.setTranslation(v, tar.translation);
  //     const m = SymmetryOperator.slerp(Mat4(), u, v, t);
  //     const rot = Mat4.getRotation(Quat.zero(), m);
  //     const axis = Vec3();
  //     const angle = Quat.getAxisAngle(axis, rot);
  //     const translation = Mat4.getTranslation(Vec3(), m);
  //     return { axis, angle, translation };
  // }
});

export { StructureInstances };
type StructureInstances = typeof StructureInstances;
const StructureInstances = PluginStateTransform.BuiltIn({
  name: 'structure-instances',
  display: { name: 'Structure Instances' },
  isDecorator: true,
  from: SO.Molecule.Structure,
  to: SO.Molecule.Structure,
  params: {
    transforms: PD.ObjectList({ transform: TransformParam }, () => 'Transform'),
  },
})({
  canAutoUpdate() {
    return true;
  },
  apply({ a, params }) {
    const center = params.transforms.some((t) => transformParamsNeedCentroid(t.transform))
      ? a.data.boundary.sphere.center
      : Vec3.unit;
    const instances = params.transforms.map((t) => getTransformFromParams(t.transform, center));
    if (!instances.length) {
      return a;
    }

    const s = Structure.instances(a.data, instances);
    return new SO.Molecule.Structure(s, { label: a.label, description: `${a.description} [Instanced]` });
  },
  dispose({ b }) {
    b?.data.customPropertyDescriptors.dispose();
  },
});

export { ModelWithCoordinates };
type ModelWithCoordinates = typeof ModelWithCoordinates;
const ModelWithCoordinates = PluginStateTransform.BuiltIn({
  name: 'model-with-coordinates',
  display: { name: 'Model With Coordinates', description: 'Updates the current model with provided coordinate frame' },
  from: SO.Molecule.Model,
  to: SO.Molecule.Model,
  params: {
    frameIndex: PD.Optional(PD.Numeric(0, undefined, { isHidden: true })),
    frameCount: PD.Optional(PD.Numeric(1, undefined, { isHidden: true })),
    atomicCoordinateFrame: PD.Optional(PD.Value<Frame | undefined>(undefined, { isHidden: true })),
  },
  isDecorator: true,
})({
  apply({ a, params, cache }) {
    (cache as { source?: Model }).source = a.data;
    if (!params.atomicCoordinateFrame) {
      return a;
    }
    const model: Model = {
      ...a.data,
      atomicConformation: Model.getAtomicConformationFromFrame(a.data, params.atomicCoordinateFrame),
    };
    Model.TrajectoryInfo.set(model, { index: params.frameIndex ?? 0, size: params.frameCount ?? 1 });
    return new SO.Molecule.Model(model, { label: a.label, description: a.description });
  },
  update: ({ a, b, oldParams, newParams, cache }) => {
    // the parent model updated in place (a trajectory frame) is a new source, also with the frame as it was
    const sourceChanged = (cache as { source?: Model }).source !== a.data;
    if (!sourceChanged && oldParams.atomicCoordinateFrame === newParams.atomicCoordinateFrame) {
      return StateTransformer.UpdateResult.Unchanged;
    }
    (cache as { source?: Model }).source = a.data;
    if (!newParams.atomicCoordinateFrame) {
      b.data = a.data;
    } else {
      b.data = {
        ...a.data,
        atomicConformation: Model.getAtomicConformationFromFrame(a.data, newParams.atomicCoordinateFrame),
      };
    }
    Model.TrajectoryInfo.set(b.data, { index: newParams.frameIndex ?? 0, size: newParams.frameCount ?? 1 });
    return StateTransformer.UpdateResult.Updated;
  },
});

export { CustomModelProperties };
type CustomModelProperties = typeof CustomModelProperties;
const CustomModelProperties = PluginStateTransform.BuiltIn({
  name: 'custom-model-properties',
  display: { name: 'Custom Model Properties' },
  isDecorator: true,
  from: SO.Molecule.Model,
  to: SO.Molecule.Model,
  params: (a, ctx: PluginContext) => {
    return ctx.customModelProperties.getParams(a?.data);
  },
})({
  apply({ a, params }, ctx: PluginContext) {
    return Task.create('Custom Props', async (taskCtx) => {
      await attachModelProps(a.data, ctx, taskCtx, params);
      return new SO.Molecule.Model(a.data, { label: a.label, description: a.description });
    });
  },
  update({ a, b, oldParams, newParams }, ctx: PluginContext) {
    return Task.create('Custom Props', async (taskCtx) => {
      // the references are on the model they were added to (a new frame of another topology has its own properties)
      for (const name of oldParams.autoAttach) {
        const property = ctx.customModelProperties.get(name);
        if (!property) continue;
        b.data.customProperties.reference(property.descriptor, false);
      }
      b.data = a.data;
      b.label = a.label;
      b.description = a.description;
      await attachModelProps(a.data, ctx, taskCtx, newParams);
      return StateTransformer.UpdateResult.Updated;
    });
  },
  dispose({ b }) {
    b?.data.customProperties.dispose();
  },
});

async function attachModelProps(
  model: Model,
  ctx: PluginContext,
  taskCtx: RuntimeContext,
  params: ReturnType<CustomModelProperties['createDefaultParams']>,
) {
  const propertyCtx = { runtime: taskCtx, assetManager: ctx.managers.asset, errorContext: ctx.errorContext };
  const { autoAttach, properties } = params;
  for (const name of Object.keys(properties)) {
    const property = ctx.customModelProperties.get(name)!;
    const props = properties[name];
    if (autoAttach.includes(name) || property.isHidden) {
      try {
        await property.attach(propertyCtx, model, props, true);
      } catch (e) {
        ctx.log.warn(`Error attaching model prop '${name}': ${e}`);
      }
    } else {
      property.set(model, props);
    }
  }
}

export { CustomStructureProperties };
type CustomStructureProperties = typeof CustomStructureProperties;
const CustomStructureProperties = PluginStateTransform.BuiltIn({
  name: 'custom-structure-properties',
  display: { name: 'Custom Structure Properties' },
  isDecorator: true,
  from: SO.Molecule.Structure,
  to: SO.Molecule.Structure,
  params: (a, ctx: PluginContext) => {
    return ctx.customStructureProperties.getParams(a?.data.root);
  },
})({
  apply({ a, params }, ctx: PluginContext) {
    return Task.create('Custom Props', async (taskCtx) => {
      await attachStructureProps(a.data.root, ctx, taskCtx, params);
      return new SO.Molecule.Structure(a.data, { label: a.label, description: a.description });
    });
  },
  update({ a, b, oldParams, newParams }, ctx: PluginContext) {
    return Task.create('Custom Props', async (taskCtx) => {
      if (a.data !== b.data) {
        // A new structure (a trajectory frame): the old one's properties go with it, and the new one's are attached
        // below. A recreate did the same, and also recreated everything below.
        b.data.customPropertyDescriptors.dispose();
      } else {
        for (const name of oldParams.autoAttach) {
          const property = ctx.customStructureProperties.get(name);
          if (!property) continue;
          a.data.customPropertyDescriptors.reference(property.descriptor, false);
        }
      }
      b.data = a.data;
      b.label = a.label;
      b.description = a.description;
      await attachStructureProps(a.data.root, ctx, taskCtx, newParams);
      return StateTransformer.UpdateResult.Updated;
    });
  },
  dispose({ b }) {
    b?.data.customPropertyDescriptors.dispose();
  },
});

async function attachStructureProps(
  structure: Structure,
  ctx: PluginContext,
  taskCtx: RuntimeContext,
  params: ReturnType<CustomStructureProperties['createDefaultParams']>,
) {
  const propertyCtx = { runtime: taskCtx, assetManager: ctx.managers.asset, errorContext: ctx.errorContext };
  const { autoAttach, properties } = params;
  for (const name of Object.keys(properties)) {
    const property = ctx.customStructureProperties.get(name)!;
    const props = properties[name];
    if (autoAttach.includes(name) || property.isHidden) {
      try {
        await property.attach(propertyCtx, structure, props, true);
      } catch (e) {
        ctx.log.warn(`Error attaching structure prop '${name}': ${e}`);
      }
    } else {
      property.set(structure, props);
    }
  }
}
