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
import { parseDx } from '@molstar/io/reader/dx/parser';
import { volumeFromDx } from '@molstar/model/formats/volume/dx';
import { DataFormatProvider, applyTransformerRaw, rawDataObject } from '@molstar/plugin/state/formats/provider';
import { VolumeFormatCategory } from './category.js';
import { type VolumeFormatParams, tryObtainRecommendedIsoValue, defaultVisuals } from './provider.js';
import { CustomVolumeProperties } from '@molstar/plugin/state/transforms/volume/ops';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { Isosurface } from '@molstar/plugin/registry/volume/isosurface';

export { ParseDx };
type ParseDx = typeof ParseDx;
const ParseDx = PluginStateTransform.BuiltIn({
  name: 'parse-dx',
  display: { name: 'Parse DX', description: 'Parse DX from Binary/String data' },
  from: [SO.Data.Binary, SO.Data.String],
  to: SO.Format.Dx,
})({
  apply({ a }) {
    return Task.create('Parse DX', async (ctx) => {
      const parsed = await parseDx(a.data, a.label).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      return new SO.Format.Dx(parsed.result);
    });
  },
});

export { VolumeFromDx };
type VolumeFromDx = typeof VolumeFromDx;
const VolumeFromDx = PluginStateTransform.BuiltIn({
  name: 'volume-from-dx',
  display: { name: 'Parse DX', description: 'Create volume from DX data.' },
  from: SO.Format.Dx,
  to: SO.Volume.Data,
})({
  apply({ a }) {
    return Task.create('Parse DX', async (ctx) => {
      const volume = await volumeFromDx(a.data, { label: a.data.name || a.label }).runInContext(ctx);
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

export const DxProvider = DataFormatProvider({
  name: 'dx',
  label: 'DX',
  description: 'DX',
  category: VolumeFormatCategory,
  stringExtensions: ['dx'],
  binaryExtensions: ['dxbin'],
  parse: async (plugin, data, params?: VolumeFormatParams) => {
    const format = plugin
      .build()
      .to(data)
      .apply(ParseDx, {}, { state: { isGhost: true } });

    const volume = format.apply(VolumeFromDx, { entryId: params?.entryId }).apply(CustomVolumeProperties);

    await volume.commit({ revertOnError: true });
    await tryObtainRecommendedIsoValue(plugin, volume.selector.data);

    return { volume: volume.selector };
  },
  parseRaw: async (plugin, ctx, data, params?: VolumeFormatParams) => {
    const format = await applyTransformerRaw(plugin, ctx, ParseDx, rawDataObject(data));
    const volume = await applyTransformerRaw(plugin, ctx, VolumeFromDx, format, {
      entryId: params?.entryId,
    });
    return { volume: volume.data };
  },
  visuals: defaultVisuals,
});

/** The Dx data format with its actions and the representations its visuals use. */
export const Dx: PluginRegistryEntry = {
  formats: [DxProvider],
  actions: [VolumeFromDx],
  ...Isosurface,
};
