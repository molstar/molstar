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
import { volumeFromSegmentationData } from '@molstar/model/formats/volume/segmentation';
import {
  DataFormatProvider,
  guessCifVariant,
  applyTransformerRaw,
  rawDataObject,
} from '@molstar/plugin/state/formats/provider';
import { VolumeFormatCategory } from './category.js';
import type { VolumeData } from './provider.js';
import { Volume } from '@molstar/model/model/volume';
import type { StateObjectSelector } from '@molstar/core/state';
import { VolumeRepresentation3DHelpers } from '@molstar/plugin/state/transforms/volume/representation-helpers';
import { ParseCif } from '@molstar/plugin/state/formats/cif';
import { CustomVolumeProperties } from '@molstar/plugin/state/transforms/volume/ops';
import { VolumeRepresentation3D } from '@molstar/plugin/state/transforms/volume/representation';
import { SegmentRepresentationProvider } from '@molstar/graphics/repr/volume/segment';
import { VolumeSegmentColorThemeProvider } from '@molstar/graphics/theme/color/volume-segment';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { Segment } from '@molstar/plugin/registry/volume/segment';

export { VolumeFromSegmentationCif };
type VolumeFromSegmentationCif = typeof VolumeFromSegmentationCif;
const VolumeFromSegmentationCif = PluginStateTransform.BuiltIn({
  name: 'volume-from-segmentation-cif',
  display: { name: 'Volume from Segmentation CIF' },
  from: SO.Format.Cif,
  to: SO.Volume.Data,
  params(a) {
    const blocks = a?.data.blocks.slice(1);
    const blockHeaderParam = blocks
      ? PD.Optional(
          PD.Select(
            blocks[0] && blocks[0].header,
            blocks.map((b) => [b.header, b.header] as [string, string]),
            { description: 'Header of the block to parse' },
          ),
        )
      : PD.Optional(
          PD.Text(void 0, {
            description: 'Header of the block to parse. If none is specifed, the 1st data block in the file is used.',
          }),
        );
    return {
      blockHeader: blockHeaderParam,
      segmentLabels: PD.ObjectList({ id: PD.Numeric(-1), label: PD.Text('') }, (s) => `${s.id} = ${s.label}`, {
        description: 'Mapping of segment IDs to segment labels',
      }),
      ownerId: PD.Text('', { isHidden: true, description: 'Reference to the object which manages this volume' }),
    };
  },
})({
  isApplicable: (a) => a.data.blocks.length > 0,
  apply({ a, params }) {
    return Task.create('Parse segmentation CIF', async (ctx) => {
      const header = params.blockHeader || a.data.blocks[1].header; // zero block contains query meta-data
      const block = a.data.blocks.find((b) => b.header === header);
      if (!block) throw new Error(`Data block '${[header]}' not found.`);
      const segmentationCif = CIF.schema.segmentation(block);
      const segmentLabels: { [id: number]: string } = {};
      for (const segment of params.segmentLabels) segmentLabels[segment.id] = segment.label;
      const volume = await volumeFromSegmentationData(segmentationCif, {
        segmentLabels,
        ownerId: params.ownerId,
      }).runInContext(ctx);
      const [x, y, z] = volume.grid.cells.space.dimensions;
      const props = {
        label: segmentationCif.volume_data_3d_info.name.value(0),
        description: `Segmentation ${x}\u00D7${y}\u00D7${z}`,
      };
      return new SO.Volume.Data(volume, props);
    });
  },
  dispose({ b }) {
    b?.data.customProperties.dispose();
  },
});

export const SegcifProvider = DataFormatProvider({
  name: 'segcif',
  label: 'Segmentation CIF',
  description: 'Segmentation CIF',
  category: VolumeFormatCategory,
  stringExtensions: ['cif'],
  binaryExtensions: ['bcif'],
  isApplicable: (info, data) => {
    return guessCifVariant(info, data) === 'segcif' ? true : false;
  },
  parse: async (plugin, data) => {
    const cifCell = await plugin.build().to(data).apply(ParseCif).commit();
    const b = plugin.build().to(cifCell);
    const blocks = cifCell.obj!.data.blocks;

    if (blocks.length === 0) throw new Error('no data blocks');

    const volumes: VolumeData[] = [];
    for (const block of blocks) {
      // Skip "server" data block.
      if (block.header.toUpperCase() === 'SERVER') continue;

      if (block.categories['volume_data_3d_info']?.rowCount > 0) {
        const volume = b.apply(VolumeFromSegmentationCif, { blockHeader: block.header }).apply(CustomVolumeProperties);
        volumes.push(volume.selector);
      }
    }

    await b.commit();

    return { volumes };
  },
  parseRaw: async (plugin, ctx, data) => {
    const cif = await applyTransformerRaw(plugin, ctx, ParseCif, rawDataObject(data));
    const blocks = cif.data.blocks;
    if (blocks.length === 0) throw new Error('no data blocks');

    const volumes: Volume[] = [];
    for (const block of blocks) {
      if (block.header.toUpperCase() === 'SERVER') continue;
      if (!(block.categories['volume_data_3d_info']?.rowCount > 0)) continue;

      const volume = await applyTransformerRaw(plugin, ctx, VolumeFromSegmentationCif, cif, {
        blockHeader: block.header,
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
      const segmentation = Volume.Segmentation.get(volumes[0].data!);
      if (segmentation) {
        visuals[visuals.length] = tree
          .to(volumes[0])
          .apply(
            VolumeRepresentation3D,
            VolumeRepresentation3DHelpers.getDefaultParams(
              plugin,
              SegmentRepresentationProvider,
              volumes[0].data!,
              { alpha: 1, instanceGranularity: true },
              VolumeSegmentColorThemeProvider,
              {},
            ),
          ).selector;
      }
    }

    await tree.commit();

    return visuals;
  },
});

/** The Segcif data format and the representations its visuals use. */
export const Segcif: PluginRegistryEntry = {
  formats: [SegcifProvider],
  ...Segment,
};
