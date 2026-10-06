/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Aliaksei Chareshneu <chareshneu.tech@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';
import { Volume } from '@molstar/model/model/volume';
import { Task } from '@molstar/core/task';
import { getContourLevelEmdb } from '@molstar/plugin/behavior/dynamic/volume-streaming/util';
import { RecommendedIsoValue } from '@molstar/model/formats/volume/property';
import type { StateObjectSelector } from '@molstar/core/state';
import type { PluginStateObject } from '@molstar/plugin/state/objects';
import { createVolumeRepresentationParams } from '@molstar/plugin/state/helpers/volume-representation-params';
import { VolumeRepresentation3D } from '@molstar/plugin/state/transforms/volume/representation';

export type VolumeFormatParams = { entryId?: string };

export async function tryObtainRecommendedIsoValue(plugin: PluginContext, volume?: Volume) {
  if (!volume) return;

  const { entryId } = volume;
  if (!entryId || !entryId.toLowerCase().startsWith('emd')) return;

  return plugin.runTask(
    Task.create('Try Set Recommended IsoValue', async (ctx) => {
      try {
        const absIsoLevel = await getContourLevelEmdb(plugin, ctx, entryId);
        RecommendedIsoValue.Provider.set(volume, Volume.IsoValue.absolute(absIsoLevel));
      } catch (e) {
        console.warn(e);
      }
    }),
  );
}

export function tryGetRecomendedIsoValue(volume: Volume) {
  const recommendedIsoValue = RecommendedIsoValue.Provider.get(volume);
  if (!recommendedIsoValue) return;
  if (recommendedIsoValue.kind === 'relative') return recommendedIsoValue;

  return Volume.adjustedIsoValue(volume, recommendedIsoValue.absoluteValue, 'absolute');
}

export type VolumeData = StateObjectSelector<PluginStateObject.Volume.Data>;

export async function defaultVisuals(plugin: PluginContext, data: { volume: VolumeData }) {
  const typeParams: { isoValue?: Volume.IsoValue } = {};
  const isoValue = data.volume.data && tryGetRecomendedIsoValue(data.volume.data);
  if (isoValue) typeParams.isoValue = isoValue;

  const visual = plugin
    .build()
    .to(data.volume)
    .apply(
      VolumeRepresentation3D,
      createVolumeRepresentationParams(plugin, data.volume.data, {
        type: 'isosurface',
        typeParams,
      }),
    );
  return [await visual.commit()];
}
