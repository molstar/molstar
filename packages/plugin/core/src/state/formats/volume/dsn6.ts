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
import * as DSN6 from '@molstar/io/reader/dsn6/parser';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { volumeFromDsn6 } from '@molstar/model/formats/volume/dsn6';
import { DataFormatProvider, applyTransformerRaw, rawDataObject } from '@molstar/plugin/state/formats/provider';
import { VolumeFormatCategory } from './category.js';
import { type VolumeFormatParams, tryObtainRecommendedIsoValue, defaultVisuals } from './provider.js';
import { CustomVolumeProperties } from '@molstar/plugin/state/transforms/volume/ops';

export { ParseDsn6 };
type ParseDsn6 = typeof ParseDsn6;
const ParseDsn6 = PluginStateTransform.BuiltIn({
  name: 'parse-dsn6',
  display: { name: 'Parse DSN6/BRIX', description: 'Parse CCP4/BRIX from Binary data' },
  from: [SO.Data.Binary],
  to: SO.Format.Dsn6,
})({
  apply({ a }) {
    return Task.create('Parse DSN6/BRIX', async (ctx) => {
      const parsed = await DSN6.parse(a.data, a.label).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      return new SO.Format.Dsn6(parsed.result);
    });
  },
});

export { VolumeFromDsn6 };
type VolumeFromDsn6 = typeof VolumeFromDsn6;
const VolumeFromDsn6 = PluginStateTransform.BuiltIn({
  name: 'volume-from-dsn6',
  display: { name: 'Volume from DSN6/BRIX', description: 'Create Volume from DSN6/BRIX data' },
  from: SO.Format.Dsn6,
  to: SO.Volume.Data,
  params(a) {
    return {
      voxelSize: PD.Vec3(Vec3.create(1, 1, 1)),
      entryId: PD.Text(''),
    };
  },
})({
  apply({ a, params }) {
    return Task.create('Create volume from DSN6/BRIX', async (ctx) => {
      const volume = await volumeFromDsn6(a.data, { ...params, label: a.data.name || a.label }).runInContext(ctx);
      const props = {
        label: volume.label || 'Volume',
        description: `Volume ${a.data.header.xExtent}\u00D7${a.data.header.yExtent}\u00D7${a.data.header.zExtent}`,
      };
      return new SO.Volume.Data(volume, props);
    });
  },
  dispose({ b }) {
    b?.data.customProperties.dispose();
  },
});

export const Dsn6Provider = DataFormatProvider({
  name: 'dsn6',
  label: 'DSN6/BRIX',
  description: 'DSN6/BRIX',
  category: VolumeFormatCategory,
  binaryExtensions: ['dsn6', 'brix'],
  parse: async (plugin, data, params?: VolumeFormatParams) => {
    const format = plugin
      .build()
      .to(data)
      .apply(ParseDsn6, {}, { state: { isGhost: true } });

    const volume = format.apply(VolumeFromDsn6, { entryId: params?.entryId }).apply(CustomVolumeProperties);

    await format.commit({ revertOnError: true });
    await tryObtainRecommendedIsoValue(plugin, volume.selector.data);

    return { format: format.selector, volume: volume.selector };
  },
  parseRaw: async (plugin, ctx, data, params?: VolumeFormatParams) => {
    const format = await applyTransformerRaw(plugin, ctx, ParseDsn6, rawDataObject(data));
    const volume = await applyTransformerRaw(plugin, ctx, VolumeFromDsn6, format, {
      entryId: params?.entryId,
    });
    return { volume: volume.data };
  },
  visuals: defaultVisuals,
});
