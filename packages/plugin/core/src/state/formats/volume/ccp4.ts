/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Neli Fonseca <neli@ebi.ac.uk>
 * @author Yakov Pechersky <ffxen158@gmail.com>
 * @author Aliaksei Chareshneu <chareshneu.tech@gmail.com>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { Task } from '@molstar/core/task';
import * as CCP4 from '@molstar/io/reader/ccp4/parser';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { volumeFromCcp4 } from '@molstar/model/formats/volume/ccp4';
import { DataFormatProvider, applyTransformerRaw, rawDataObject } from '@molstar/plugin/state/formats/provider';
import { VolumeFormatCategory } from './category.js';
import { type VolumeFormatParams, tryObtainRecommendedIsoValue, defaultVisuals } from './provider.js';
import { CustomVolumeProperties } from '@molstar/plugin/state/transforms/volume/ops';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { Isosurface } from '@molstar/plugin/registry/volume/isosurface';

export { ParseCcp4 };
type ParseCcp4 = typeof ParseCcp4;
const ParseCcp4 = PluginStateTransform.BuiltIn({
  name: 'parse-ccp4',
  display: { name: 'Parse CCP4/MRC/MAP', description: 'Parse CCP4/MRC/MAP from Binary data' },
  from: [SO.Data.Binary],
  to: SO.Format.Ccp4,
})({
  apply({ a }) {
    return Task.create('Parse CCP4/MRC/MAP', async (ctx) => {
      const parsed = await CCP4.parse(a.data, a.label).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      return new SO.Format.Ccp4(parsed.result);
    });
  },
});

export { VolumeFromCcp4 };
type VolumeFromCcp4 = typeof VolumeFromCcp4;
const VolumeFromCcp4 = PluginStateTransform.BuiltIn({
  name: 'volume-from-ccp4',
  display: { name: 'Volume from CCP4/MRC/MAP', description: 'Create Volume from CCP4/MRC/MAP data' },
  from: SO.Format.Ccp4,
  to: SO.Volume.Data,
  params(a) {
    return {
      voxelSize: PD.Vec3(Vec3.create(1, 1, 1)),
      offset: PD.Vec3(Vec3.create(0, 0, 0)),
      entryId: PD.Text(''),
    };
  },
})({
  apply({ a, params }) {
    return Task.create('Create volume from CCP4/MRC/MAP', async (ctx) => {
      const volume = await volumeFromCcp4(a.data, { ...params, label: a.data.name || a.label }).runInContext(ctx);
      const props = {
        label: volume.label || 'Volume',
        description: `Volume ${a.data.header.NX}\u00D7${a.data.header.NX}\u00D7${a.data.header.NX}`,
      };
      return new SO.Volume.Data(volume, props);
    });
  },
  dispose({ b }) {
    b?.data.customProperties.dispose();
  },
});

export const Ccp4Provider = DataFormatProvider({
  name: 'ccp4',
  label: 'CCP4/MRC/MAP',
  description: 'CCP4/MRC/MAP',
  category: VolumeFormatCategory,
  binaryExtensions: ['ccp4', 'mrc', 'map'],
  parse: async (plugin, data, params?: VolumeFormatParams) => {
    const format = plugin
      .build()
      .to(data)
      .apply(ParseCcp4, {}, { state: { isGhost: true } });

    const volume = format.apply(VolumeFromCcp4, { entryId: params?.entryId }).apply(CustomVolumeProperties);

    await format.commit({ revertOnError: true });
    await tryObtainRecommendedIsoValue(plugin, volume.selector.data);

    return { format: format.selector, volume: volume.selector };
  },
  parseRaw: async (plugin, ctx, data, params?: VolumeFormatParams) => {
    const format = await applyTransformerRaw(plugin, ctx, ParseCcp4, rawDataObject(data));
    const volume = await applyTransformerRaw(plugin, ctx, VolumeFromCcp4, format, {
      entryId: params?.entryId,
    });
    return { volume: volume.data };
  },
  visuals: defaultVisuals,
});

/** The Ccp4 data format with its actions and the representations its visuals use. */
export const Ccp4: PluginRegistryEntry = {
  formats: [Ccp4Provider],
  actions: [ParseCcp4, VolumeFromCcp4],
  ...Isosurface,
};
