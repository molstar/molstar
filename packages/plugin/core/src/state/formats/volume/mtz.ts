/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Neli Fonseca <neli@ebi.ac.uk>
 * @author Yakov Pechersky <ffxen158@gmail.com>
 * @author Aliaksei Chareshneu <chareshneu.tech@gmail.com>
 */

import { PluginStateTransform, PluginStateObject as SO, type PluginStateObject } from '@molstar/plugin/state/objects';
import { Task } from '@molstar/core/task';
import { parseMtz } from '@molstar/io/reader/mtz/parser';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { volumeFromMtz, detectMtzColumnPairs } from '@molstar/model/formats/volume/mtz';
import { DataFormatProvider, applyTransformerRaw, rawDataObject } from '@molstar/plugin/state/formats/provider';
import { VolumeFormatCategory } from './category.js';
import type { VolumeData } from './provider.js';
import { Volume } from '@molstar/model/model/volume';
import type { StateObjectSelector } from '@molstar/core/state';
import { VolumeRepresentation3DHelpers } from '@molstar/plugin/state/transforms/volume/representation-helpers';
import { Color } from '@molstar/core/util/color/color';
import { CustomVolumeProperties } from '@molstar/plugin/state/transforms/volume/ops';
import { VolumeRepresentation3D } from '@molstar/plugin/state/transforms/volume/representation';

export { ParseMtz };
type ParseMtz = typeof ParseMtz;
const ParseMtz = PluginStateTransform.BuiltIn({
  name: 'parse-mtz',
  display: { name: 'Parse MTZ', description: 'Parse CCP4 MTZ reflection file from Binary data' },
  from: [SO.Data.Binary],
  to: SO.Format.Mtz,
})({
  apply({ a }) {
    return Task.create('Parse MTZ', async (ctx) => {
      const parsed = await parseMtz(a.data, a.label).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      return new SO.Format.Mtz(parsed.result);
    });
  },
});

export { VolumeFromMtz };
type VolumeFromMtz = typeof VolumeFromMtz;
const VolumeFromMtz = PluginStateTransform.BuiltIn({
  name: 'volume-from-mtz',
  display: { name: 'Volume from MTZ', description: 'Compute electron density map from MTZ reflection data' },
  from: SO.Format.Mtz,
  to: SO.Volume.Data,
  params(a) {
    if (!a) {
      return {
        ampLabel: PD.Text('FWT', { description: 'Column label of the structure-factor amplitude (type F)' }),
        phiLabel: PD.Text('PHWT', { description: 'Column label of the phase in degrees (type P)' }),
        entryId: PD.Text(''),
      };
    }
    const fCols = a.data.header.columns
      .filter((c) => c.type === 'F')
      .map((c) => [c.label, c.label] as [string, string]);
    const pCols = a.data.header.columns
      .filter((c) => c.type === 'P')
      .map((c) => [c.label, c.label] as [string, string]);
    return {
      ampLabel:
        fCols.length > 0
          ? PD.Select(fCols[0][0], fCols, { description: 'Structure-factor amplitude column (type F)' })
          : PD.Text('FWT', { description: 'Structure-factor amplitude column label (type F)' }),
      phiLabel:
        pCols.length > 0
          ? PD.Select(pCols[0][0], pCols, { description: 'Phase column in degrees (type P)' })
          : PD.Text('PHWT', { description: 'Phase column label in degrees (type P)' }),
      entryId: PD.Text(''),
    };
  },
})({
  apply({ a, params }) {
    return Task.create('Compute volume from MTZ', async (ctx) => {
      const volume = await volumeFromMtz(a.data, {
        ampLabel: params.ampLabel,
        phiLabel: params.phiLabel,
        entryId: params.entryId || undefined,
        label: params.entryId || a.data.name || undefined,
      }).runInContext(ctx);
      const [x, y, z] = volume.grid.cells.space.dimensions;
      const label = params.entryId || a.data.name || 'Volume';
      const props = { label, description: `Volume ${x}\u00D7${y}\u00D7${z}` };
      return new SO.Volume.Data(volume, props);
    });
  },
  dispose({ b }) {
    b?.data.customProperties.dispose();
  },
});

type MtzParams = { entryId?: string };

export const MtzProvider = DataFormatProvider({
  label: 'MTZ',
  description: 'CCP4 MTZ reflection data file with amplitude and phase columns',
  category: VolumeFormatCategory,
  binaryExtensions: ['mtz'],
  isApplicable: (info, data) => {
    // Check magic bytes "MTZ"
    if (!(data instanceof Uint8Array) || data.length < 4) return false;
    return data[0] === 0x4d && data[1] === 0x54 && data[2] === 0x5a;
  },
  parse: async (plugin, data, params?: MtzParams) => {
    const mtzCell = await plugin.build().to(data).apply(ParseMtz).commit();
    const b = plugin.build().to(mtzCell);
    const mtz = mtzCell.obj!.data;
    const pairs = detectMtzColumnPairs(mtz.header);

    if (pairs.length === 0) {
      throw new Error(
        'MTZ file does not contain any recognised amplitude+phase column pairs (FWT/PHWT, DELFWT/DELPHWT, 2FOFCWT/PH2FOFCWT, FOFCWT/PHFOFCWT).',
      );
    }

    const volumes: { '2fofc': VolumeData[]; fofc: VolumeData[] } = { '2fofc': [], fofc: [] };
    for (const pair of pairs) {
      const vol = b
        .apply(VolumeFromMtz, {
          ampLabel: pair.ampLabel,
          phiLabel: pair.phiLabel,
          entryId: params?.entryId,
        })
        .apply(CustomVolumeProperties);
      if (pair.label === '2fo-fc') volumes['2fofc'].push(vol.selector);
      else volumes['fofc'].push(vol.selector);
    }

    await b.commit();
    return { volumes };
  },
  parseRaw: async (plugin, ctx, data, params?: MtzParams) => {
    const mtz = await applyTransformerRaw(plugin, ctx, ParseMtz, rawDataObject(data));
    const pairs = detectMtzColumnPairs(mtz.data.header);

    if (pairs.length === 0) {
      throw new Error(
        'MTZ file does not contain any recognised amplitude+phase column pairs (FWT/PHWT, DELFWT/DELPHWT, 2FOFCWT/PH2FOFCWT, FOFCWT/PHFOFCWT).',
      );
    }

    const twoFoFc: Volume[] = [];
    const foFc: Volume[] = [];
    for (const pair of pairs) {
      const volume = await applyTransformerRaw(plugin, ctx, VolumeFromMtz, mtz, {
        ampLabel: pair.ampLabel,
        phiLabel: pair.phiLabel,
        entryId: params?.entryId,
      });
      if (pair.label === '2fo-fc') twoFoFc.push(volume.data);
      else foFc.push(volume.data);
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
              'isosurface',
              { isoValue, alpha: 1 },
              'uniform',
              { value: Color(0x3362b2) },
            ),
          ).selector,
      );
    }

    // Fo-Fc map: green positive / red negative at ±3σ
    if (volumes['fofc'].length > 0) {
      const posParams = VolumeRepresentation3DHelpers.getDefaultParamsStatic(
        plugin,
        'isosurface',
        { isoValue: Volume.IsoValue.relative(3), alpha: 0.3 },
        'uniform',
        { value: Color(0x33bb33) },
      );
      const negParams = VolumeRepresentation3DHelpers.getDefaultParamsStatic(
        plugin,
        'isosurface',
        { isoValue: Volume.IsoValue.relative(-3), alpha: 0.3 },
        'uniform',
        { value: Color(0xbb3333) },
      );
      visuals.push(tree.to(volumes['fofc'][0]).apply(VolumeRepresentation3D, posParams).selector);
      visuals.push(tree.to(volumes['fofc'][0]).apply(VolumeRepresentation3D, negParams).selector);
    }

    await tree.commit();
    return visuals;
  },
});
