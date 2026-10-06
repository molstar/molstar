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
import {
  type StructureFactorMapType,
  volumeFromStructureFactors,
} from '@molstar/model/formats/volume/structure-factors';
import { Task } from '@molstar/core/task';
import { CIF } from '@molstar/io/reader/cif';
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
import { Color } from '@molstar/core/util/color/color';
import { ParseCif } from '@molstar/plugin/state/formats/cif';
import { CustomVolumeProperties } from '@molstar/plugin/state/transforms/volume/ops';
import { VolumeRepresentation3D } from '@molstar/plugin/state/transforms/volume/representation';
import { IsosurfaceRepresentationProvider } from '@molstar/graphics/repr/volume/isosurface';
import { UniformColorThemeProvider } from '@molstar/graphics/theme/color/uniform';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { Isosurface } from '@molstar/plugin/registry/volume/isosurface';

export { VolumeFromStructureFactorsCif };
type VolumeFromStructureFactorsCif = typeof VolumeFromStructureFactorsCif;
const VolumeFromStructureFactorsCif = PluginStateTransform.BuiltIn({
  name: 'volume-from-structure-factors-cif',
  display: {
    name: 'Volume from Structure Factors CIF',
    description:
      'Compute electron density map from structure factor reflection data (pdbx_FWT/pdbx_PHWT or pdbx_DELFWT/pdbx_DELPHWT)',
  },
  from: SO.Format.Cif,
  to: SO.Volume.Data,
  params(a) {
    const blockOptions: [string, string][] = a
      ? a.data.blocks.map((b) => [b.header, b.header] as [string, string])
      : [];
    return {
      blockHeader: PD.Optional(
        blockOptions.length > 0
          ? PD.Select(blockOptions[0][0], blockOptions, { description: 'Header of the data block to use' })
          : PD.Text(void 0, { description: 'Header of the block to parse. Defaults to first block.' }),
      ),
      entryId: PD.Text(''),
      mapType: PD.Select<StructureFactorMapType>('2fo-fc', [
        ['2fo-fc', '2Fo-Fc (pdbx_FWT/PHWT)'],
        ['fo-fc', 'Fo-Fc (pdbx_DELFWT/DELPHWT)'],
      ]),
    };
  },
})({
  isApplicable: (a) => a.data.blocks.length > 0,
  apply({ a, params }) {
    return Task.create('Parse structure factors CIF', async (ctx) => {
      const header = params.blockHeader || a.data.blocks[0].header;
      const block = a.data.blocks.find((b) => b.header === header);
      if (!block) throw new Error(`Data block '${header}' not found.`);
      const sfCif = CIF.schema.SF(block);
      const volume = await volumeFromStructureFactors(sfCif, {
        entryId: params.entryId,
        mapType: params.mapType,
      }).runInContext(ctx);
      const [x, y, z] = volume.grid.cells.space.dimensions;
      const label = params.entryId || header;
      const mapLabel = params.mapType === 'fo-fc' ? 'Fo-Fc' : '2Fo-Fc';
      const props = { label, description: `${mapLabel} Volume ${x}\u00D7${y}\u00D7${z}` };
      return new SO.Volume.Data(volume, props);
    });
  },
  dispose({ b }) {
    b?.data.customProperties.dispose();
  },
});

type SfcifParams = { entryId?: string };

export const SfcifProvider = DataFormatProvider({
  name: 'sfcif',
  label: 'Structure Factors CIF',
  description: 'mmCIF structure factor file with reflection data (pdbx_FWT/pdbx_PHWT map coefficients)',
  category: VolumeFormatCategory,
  stringExtensions: ['cif'],
  binaryExtensions: ['bcif'],
  isApplicable: (info, data) => {
    return guessCifVariant(info, data) === 'sfcif' ? true : false;
  },
  parse: async (plugin, data, params?: SfcifParams) => {
    const cifCell = await plugin.build().to(data).apply(ParseCif).commit();
    const b = plugin.build().to(cifCell);
    const blocks = cifCell.obj!.data.blocks;

    if (blocks.length === 0) throw new Error('no data blocks');

    const volumes: { '2fofc': VolumeData[]; fofc: VolumeData[] } = { '2fofc': [], fofc: [] };
    for (const block of blocks) {
      const reflnCat = block.categories['refln'];
      if (!(reflnCat?.rowCount > 0)) continue;

      const has2fofcCoeffs =
        (reflnCat.getField('pdbx_FWT')?.rowCount ?? 0) > 0 && (reflnCat.getField('pdbx_PHWT')?.rowCount ?? 0) > 0;
      const hasFofcCoeffs =
        (reflnCat.getField('pdbx_DELFWT')?.rowCount ?? 0) > 0 && (reflnCat.getField('pdbx_DELPHWT')?.rowCount ?? 0) > 0;
      if (!has2fofcCoeffs && !hasFofcCoeffs) continue;

      if (has2fofcCoeffs) {
        // 2Fo-Fc map
        const vol2FoFc = b
          .apply(VolumeFromStructureFactorsCif, {
            blockHeader: block.header,
            entryId: params?.entryId,
            mapType: '2fo-fc',
          })
          .apply(CustomVolumeProperties);
        volumes['2fofc'].push(vol2FoFc.selector);
      }

      if (hasFofcCoeffs) {
        // Fo-Fc difference map
        const volFoFc = b
          .apply(VolumeFromStructureFactorsCif, {
            blockHeader: block.header,
            entryId: params?.entryId,
            mapType: 'fo-fc',
          })
          .apply(CustomVolumeProperties);
        volumes['fofc'].push(volFoFc.selector);
      }
    }

    await b.commit();

    return { volumes };
  },
  parseRaw: async (plugin, ctx, data, params?: SfcifParams) => {
    const cif = await applyTransformerRaw(plugin, ctx, ParseCif, rawDataObject(data));
    const blocks = cif.data.blocks;
    if (blocks.length === 0) throw new Error('no data blocks');

    const twoFoFc: Volume[] = [];
    const foFc: Volume[] = [];
    for (const block of blocks) {
      const reflnCat = block.categories['refln'];
      if (!(reflnCat?.rowCount > 0)) continue;

      if ((reflnCat.getField('pdbx_FWT')?.rowCount ?? 0) > 0 && (reflnCat.getField('pdbx_PHWT')?.rowCount ?? 0) > 0) {
        const volume = await applyTransformerRaw(plugin, ctx, VolumeFromStructureFactorsCif, cif, {
          blockHeader: block.header,
          entryId: params?.entryId,
          mapType: '2fo-fc',
        });
        twoFoFc.push(volume.data);
      }

      if (
        (reflnCat.getField('pdbx_DELFWT')?.rowCount ?? 0) > 0 &&
        (reflnCat.getField('pdbx_DELPHWT')?.rowCount ?? 0) > 0
      ) {
        const volume = await applyTransformerRaw(plugin, ctx, VolumeFromStructureFactorsCif, cif, {
          blockHeader: block.header,
          entryId: params?.entryId,
          mapType: 'fo-fc',
        });
        foFc.push(volume.data);
      }
    }

    return { volumes: [...twoFoFc, ...foFc] };
  },
  visuals: async (plugin, data: { volumes: { '2fofc': VolumeData[]; fofc: VolumeData[] } }) => {
    const { volumes } = data;
    const tree = plugin.build();
    const visuals: StateObjectSelector<PluginStateObject.Volume.Representation3D>[] = [];

    // 2Fo-Fc map: teal solid isosurface at 2σ
    if (volumes['2fofc'].length > 0) {
      const isoValue = Volume.IsoValue.relative(2);
      visuals.push(
        tree
          .to(volumes['2fofc'][0])
          .apply(
            VolumeRepresentation3D,
            VolumeRepresentation3DHelpers.getDefaultParamsStatic(
              plugin,
              IsosurfaceRepresentationProvider,
              { isoValue, alpha: 1 },
              UniformColorThemeProvider,
              { value: Color(0x3362b2) },
            ),
          ).selector,
      );
    }

    // Fo-Fc map: green positive / red negative at ±3σ
    if (volumes['fofc'].length > 0) {
      const posParams = VolumeRepresentation3DHelpers.getDefaultParamsStatic(
        plugin,
        IsosurfaceRepresentationProvider,
        { isoValue: Volume.IsoValue.relative(3), alpha: 0.3 },
        UniformColorThemeProvider,
        { value: Color(0x33bb33) },
      );
      const negParams = VolumeRepresentation3DHelpers.getDefaultParamsStatic(
        plugin,
        IsosurfaceRepresentationProvider,
        { isoValue: Volume.IsoValue.relative(-3), alpha: 0.3 },
        UniformColorThemeProvider,
        { value: Color(0xbb3333) },
      );
      visuals.push(tree.to(volumes['fofc'][0]).apply(VolumeRepresentation3D, posParams).selector);
      visuals.push(tree.to(volumes['fofc'][0]).apply(VolumeRepresentation3D, negParams).selector);
    }

    await tree.commit();

    return visuals;
  },
});

/** The Sfcif data format and the representations its visuals use. */
export const Sfcif: PluginRegistryEntry = {
  formats: [SfcifProvider],
  ...Isosurface,
};
