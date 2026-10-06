/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Neli Fonseca <neli@ebi.ac.uk>
 * @author Adam Midlik <midlik@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 * @author Yakov Pechersky <ffxen158@gmail.com>
 * @author Aliaksei Chareshneu <chareshneu.tech@gmail.com>
 */

import { PluginStateTransform, PluginStateObject as SO, type PluginStateObject } from '@molstar/plugin/state/objects';
import { Task } from '@molstar/core/task';
import { parseCube } from '@molstar/io/reader/cube/parser';
import { trajectoryFromCube } from '@molstar/model/formats/structure/cube';
import { trajectoryProps } from '@molstar/plugin/state/formats/trajectory/helpers';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { volumeFromCube } from '@molstar/model/formats/volume/cube';
import { DataFormatProvider, applyTransformerRaw, rawDataObject } from '@molstar/plugin/state/formats/provider';
import { VolumeFormatCategory } from './category.js';
import { type VolumeFormatParams, tryObtainRecommendedIsoValue, type VolumeData } from './provider.js';
import type { PluginContext } from '@molstar/plugin/context';
import type { StateObjectSelector } from '@molstar/core/state';
import { Volume } from '@molstar/model/model/volume';
import { createVolumeRepresentationParams } from '@molstar/plugin/state/helpers/volume-representation-params';
import { ColorNames } from '@molstar/core/util/color/names';
import { objectForEach } from '@molstar/core/util/object';
import { CustomVolumeProperties } from '@molstar/plugin/state/transforms/volume/ops';
import {
  ModelFromTrajectory,
  CustomModelProperties,
  StructureFromModel,
  CustomStructureProperties,
} from '@molstar/plugin/state/transforms/structure/hierarchy';
import { VolumeRepresentation3D } from '@molstar/plugin/state/transforms/volume/representation';

export { ParseCube };
type ParseCube = typeof ParseCube;
const ParseCube = PluginStateTransform.BuiltIn({
  name: 'parse-cube',
  display: { name: 'Parse Cube', description: 'Parse Cube from String data' },
  from: SO.Data.String,
  to: SO.Format.Cube,
})({
  apply({ a }) {
    return Task.create('Parse Cube', async (ctx) => {
      const parsed = await parseCube(a.data, a.label).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      return new SO.Format.Cube(parsed.result);
    });
  },
});

export { TrajectoryFromCube };
type TrajectoryFromCube = typeof TrajectoryFromCube;
const TrajectoryFromCube = PluginStateTransform.BuiltIn({
  name: 'trajectory-from-cube',
  display: { name: 'Parse Cube', description: 'Parse Cube file to create a trajectory.' },
  from: SO.Format.Cube,
  to: SO.Molecule.Trajectory,
})({
  apply({ a }) {
    return Task.create('Parse MOL', async (ctx) => {
      const models = await trajectoryFromCube(a.data).runInContext(ctx);
      const props = trajectoryProps(models);
      return new SO.Molecule.Trajectory(models, props);
    });
  },
});

export { VolumeFromCube };
type VolumeFromCube = typeof VolumeFromCube;
const VolumeFromCube = PluginStateTransform.BuiltIn({
  name: 'volume-from-cube',
  display: { name: 'Volume from Cube', description: 'Create Volume from Cube data' },
  from: SO.Format.Cube,
  to: SO.Volume.Data,
  params(a) {
    const dataIndex = a
      ? PD.Select(
          0,
          a.data.header.dataSetIds.map((id, i) => [i, `${id}`] as const),
        )
      : PD.Numeric(0);
    return {
      dataIndex,
      entryId: PD.Text(''),
    };
  },
})({
  apply({ a, params }) {
    return Task.create('Create volume from Cube', async (ctx) => {
      const volume = await volumeFromCube(a.data, { ...params, label: a.data.name || a.label }).runInContext(ctx);
      const props = {
        label: volume.label || 'Volume',
        description: `Volume ${a.data.header.dim[0]}\u00D7${a.data.header.dim[1]}\u00D7${a.data.header.dim[2]}`,
      };
      return new SO.Volume.Data(volume, props);
    });
  },
  dispose({ b }) {
    b?.data.customProperties.dispose();
  },
});

export const CubeProvider = DataFormatProvider({
  name: 'cube',
  label: 'Cube',
  description: 'Cube',
  category: VolumeFormatCategory,
  stringExtensions: ['cub', 'cube'],
  parse: async (plugin, data, params?: VolumeFormatParams) => {
    const format = plugin
      .build()
      .to(data)
      .apply(ParseCube, {}, { state: { isGhost: true } });

    const volume = format.apply(VolumeFromCube, { entryId: params?.entryId }).apply(CustomVolumeProperties);
    const structure = format
      .apply(TrajectoryFromCube, void 0, { state: { isGhost: true } })
      .apply(ModelFromTrajectory)
      .apply(CustomModelProperties)
      .apply(StructureFromModel)
      .apply(CustomStructureProperties);

    await format.commit({ revertOnError: true });
    await tryObtainRecommendedIsoValue(plugin, volume.selector.data);

    return { format: format.selector, volume: volume.selector, structure: structure.selector };
  },
  parseRaw: async (plugin, ctx, data, params?: VolumeFormatParams) => {
    const format = await applyTransformerRaw(plugin, ctx, ParseCube, rawDataObject(data));
    const volume = await applyTransformerRaw(plugin, ctx, VolumeFromCube, format, {
      entryId: params?.entryId,
    });
    return { volume: volume.data };
  },
  visuals: async (
    plugin: PluginContext,
    data: { volume: VolumeData; structure: StateObjectSelector<PluginStateObject.Molecule.Structure> },
  ) => {
    const surfaces = plugin.build();

    const volumeReprs: StateObjectSelector<PluginStateObject.Volume.Representation3D>[] = [];
    const volumeData = data.volume.cell?.obj?.data;
    if (volumeData && Volume.isOrbitals(volumeData)) {
      const volumePos = surfaces.to(data.volume).apply(
        VolumeRepresentation3D,
        createVolumeRepresentationParams(plugin, volumeData, {
          type: 'isosurface',
          typeParams: { isoValue: Volume.IsoValue.relative(1), alpha: 0.4 },
          color: 'uniform',
          colorParams: { value: ColorNames.blue },
        }),
      );
      const volumeNeg = surfaces.to(data.volume).apply(
        VolumeRepresentation3D,
        createVolumeRepresentationParams(plugin, volumeData, {
          type: 'isosurface',
          typeParams: { isoValue: Volume.IsoValue.relative(-1), alpha: 0.4 },
          color: 'uniform',
          colorParams: { value: ColorNames.red },
        }),
      );
      volumeReprs.push(volumePos.selector, volumeNeg.selector);
    } else {
      const volume = surfaces.to(data.volume).apply(
        VolumeRepresentation3D,
        createVolumeRepresentationParams(plugin, volumeData, {
          type: 'isosurface',
          typeParams: { isoValue: Volume.IsoValue.relative(2), alpha: 0.4 },
          color: 'uniform',
          colorParams: { value: ColorNames.grey },
        }),
      );
      volumeReprs.push(volume.selector);
    }

    const structure = await plugin.builders.structure.representation.applyPreset(data.structure, 'auto');
    await surfaces.commit();

    const structureReprs: StateObjectSelector<PluginStateObject.Molecule.Structure.Representation3D>[] = [];
    objectForEach(structure?.representations as any, (r: any) => {
      if (r) structureReprs.push(r);
    });

    return [...volumeReprs, ...structureReprs];
  },
});
