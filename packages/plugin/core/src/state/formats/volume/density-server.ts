/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Yakov Pechersky <ffxen158@gmail.com>
 * @author Aliaksei Chareshneu <chareshneu.tech@gmail.com>
 */

import { PluginStateTransform, PluginStateObject as SO, type PluginStateObject } from '@molstar/plugin/state/objects';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Task } from '@molstar/core/task';
import { CIF } from '@molstar/io/reader/cif';
import { volumeFromDensityServerData } from '@molstar/model/formats/volume/density-server';
import {
  DataFormatProvider,
  guessCifVariant,
  applyTransformerRaw,
  rawDataObject,
} from '@molstar/plugin/state/formats/provider';
import { VolumeFormatCategory } from './category.js';
import { type VolumeData, tryObtainRecommendedIsoValue, tryGetRecomendedIsoValue } from './provider.js';
import { Volume } from '@molstar/model/model/volume';
import type { StateObjectSelector } from '@molstar/core/state';
import { VolumeRepresentation3DHelpers } from '@molstar/plugin/state/transforms/volume/representation-helpers';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParseCif } from '@molstar/plugin/state/formats/cif';
import { CustomVolumeProperties } from '@molstar/plugin/state/transforms/volume/ops';
import { VolumeRepresentation3D } from '@molstar/plugin/state/transforms/volume/representation';
import { IsosurfaceRepresentationProvider } from '@molstar/graphics/repr/volume/isosurface';
import { UniformColorThemeProvider } from '@molstar/graphics/theme/color/uniform';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { Isosurface } from '@molstar/plugin/registry/volume/isosurface';

export { VolumeFromDensityServerCif };
type VolumeFromDensityServerCif = typeof VolumeFromDensityServerCif;
const VolumeFromDensityServerCif = PluginStateTransform.BuiltIn({
  name: 'volume-from-density-server-cif',
  display: {
    name: 'Volume from density-server CIF',
    description: 'Identify and create all separate models in the specified CIF data block',
  },
  from: SO.Format.Cif,
  to: SO.Volume.Data,
  params(a) {
    if (!a) {
      return {
        blockHeader: PD.Optional(
          PD.Text(void 0, {
            description: 'Header of the block to parse. If none is specifed, the 1st data block in the file is used.',
          }),
        ),
        entryId: PD.Text(''),
      };
    }
    const blocks = a.data.blocks.slice(1); // zero block contains query meta-data
    return {
      blockHeader: PD.Optional(
        PD.Select(
          blocks[0] && blocks[0].header,
          blocks.map((b) => [b.header, b.header] as [string, string]),
          { description: 'Header of the block to parse' },
        ),
      ),
      entryId: PD.Text(''),
    };
  },
})({
  isApplicable: (a) => a.data.blocks.length > 0,
  apply({ a, params }) {
    return Task.create('Parse density-server CIF', async (ctx) => {
      const header = params.blockHeader || a.data.blocks[1].header; // zero block contains query meta-data
      const block = a.data.blocks.find((b) => b.header === header);
      if (!block) throw new Error(`Data block '${[header]}' not found.`);
      const densityServerCif = CIF.schema.densityServer(block);
      const volume = await volumeFromDensityServerData(densityServerCif, { entryId: params.entryId }).runInContext(ctx);
      const [x, y, z] = volume.grid.cells.space.dimensions;
      const props = {
        label: params.entryId ?? densityServerCif.volume_data_3d_info.name.value(0),
        description: `Volume ${x}\u00D7${y}\u00D7${z}`,
      };
      return new SO.Volume.Data(volume, props);
    });
  },
  dispose({ b }) {
    b?.data.customProperties.dispose();
  },
});

type DsCifParams = { entryId?: string | string[] };

export const DscifProvider = DataFormatProvider({
  name: 'dscif',
  label: 'DensityServer CIF',
  description: 'DensityServer CIF',
  category: VolumeFormatCategory,
  stringExtensions: ['cif'],
  binaryExtensions: ['bcif'],
  isApplicable: (info, data) => {
    return guessCifVariant(info, data) === 'dscif' ? true : false;
  },
  parse: async (plugin, data, params?: DsCifParams) => {
    const cifCell = await plugin.build().to(data).apply(ParseCif).commit();
    const b = plugin.build().to(cifCell);
    const blocks = cifCell.obj!.data.blocks;

    if (blocks.length === 0) throw new Error('no data blocks');

    const volumes: VolumeData[] = [];
    let i = 0;
    for (const block of blocks) {
      // Skip "server" data block.
      if (block.header.toUpperCase() === 'SERVER') continue;

      const entryId = Array.isArray(params?.entryId) ? params?.entryId[i] : params?.entryId;
      if (block.categories['volume_data_3d_info']?.rowCount > 0) {
        const volume = b
          .apply(VolumeFromDensityServerCif, { blockHeader: block.header, entryId })
          .apply(CustomVolumeProperties);
        volumes.push(volume.selector);
        i++;
      }
    }

    await b.commit();
    for (const v of volumes) await tryObtainRecommendedIsoValue(plugin, v.data);

    return { volumes };
  },
  parseRaw: async (plugin, ctx, data, params?: DsCifParams) => {
    const cif = await applyTransformerRaw(plugin, ctx, ParseCif, rawDataObject(data));
    const blocks = cif.data.blocks;
    if (blocks.length === 0) throw new Error('no data blocks');

    const volumes: Volume[] = [];
    for (const block of blocks) {
      if (block.header.toUpperCase() === 'SERVER') continue;
      if (!(block.categories['volume_data_3d_info']?.rowCount > 0)) continue;

      const entryId = Array.isArray(params?.entryId) ? params?.entryId[volumes.length] : params?.entryId;
      const volume = await applyTransformerRaw(plugin, ctx, VolumeFromDensityServerCif, cif, {
        blockHeader: block.header,
        entryId,
      });
      volumes.push(volume.data);
    }

    return { volumes };
  },
  visuals: async (plugin, data: { volumes: StateObjectSelector<PluginStateObject.Volume.Data>[] }) => {
    const { volumes } = data;
    const tree = plugin.build();
    const visuals: StateObjectSelector<PluginStateObject.Volume.Representation3D>[] = [];

    if (volumes.length > 0) {
      const isoValue = (volumes[0].data && tryGetRecomendedIsoValue(volumes[0].data)) || Volume.IsoValue.relative(1.5);

      visuals[0] = tree
        .to(volumes[0])
        .apply(
          VolumeRepresentation3D,
          VolumeRepresentation3DHelpers.getDefaultParamsStatic(
            plugin,
            IsosurfaceRepresentationProvider,
            { isoValue, alpha: 1 },
            UniformColorThemeProvider,
            { value: ColorNames.teal },
          ),
        ).selector;
    }

    if (volumes.length > 1) {
      const posParams = VolumeRepresentation3DHelpers.getDefaultParamsStatic(
        plugin,
        IsosurfaceRepresentationProvider,
        { isoValue: Volume.IsoValue.relative(3), alpha: 0.3 },
        UniformColorThemeProvider,
        { value: ColorNames.green },
      );
      const negParams = VolumeRepresentation3DHelpers.getDefaultParamsStatic(
        plugin,
        IsosurfaceRepresentationProvider,
        { isoValue: Volume.IsoValue.relative(-3), alpha: 0.3 },
        UniformColorThemeProvider,
        { value: ColorNames.red },
      );
      visuals[visuals.length] = tree.to(volumes[1]).apply(VolumeRepresentation3D, posParams).selector;
      visuals[visuals.length] = tree.to(volumes[1]).apply(VolumeRepresentation3D, negParams).selector;
    }

    await tree.commit();

    return visuals;
  },
});

/** The Dscif data format and the representations its visuals use. */
export const Dscif: PluginRegistryEntry = {
  formats: [DscifProvider],
  ...Isosurface,
};
