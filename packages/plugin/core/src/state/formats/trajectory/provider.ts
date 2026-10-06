/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import type { StateObjectRef, StateObjectSelector, StateTransformer } from '@molstar/core/state';
import type { PluginStateObject } from '@molstar/plugin/state/objects';
import type { ParseBlob } from '@molstar/plugin/state/formats/cif';
import { type DataFormatProvider, applyTransformerRaw, rawDataObject } from '@molstar/plugin/state/formats/provider';
import type { PluginContext } from '@molstar/plugin/context';
import { PluginConfig } from '@molstar/plugin/config';

export interface TrajectoryFormatProvider<
  P extends { trajectoryTags?: string | string[] } = { trajectoryTags?: string | string[] },
  R extends { trajectory: StateObjectRef<PluginStateObject.Molecule.Trajectory> } = {
    trajectory: StateObjectRef<PluginStateObject.Molecule.Trajectory>;
  },
  Id extends string = string,
> extends DataFormatProvider<P, R, any, any, Id> {
  /**
   * Parses a blob of several data entries into a single trajectory. Used by `builders.structure.parseTrajectory(blob)`
   * and implemented by the mmCIF format only.
   */
  parseBlob?(
    plugin: PluginContext,
    blob: StateObjectRef<PluginStateObject.Data.Blob>,
    params: StateTransformer.Params<typeof ParseBlob>,
  ): Promise<StateObjectSelector<PluginStateObject.Molecule.Trajectory>>;
}

/** Identity helper that type checks a trajectory provider and keeps its `name` a literal type. */
export function TrajectoryFormatProvider<const T extends TrajectoryFormatProvider>(provider: T): T {
  return provider;
}

export function defaultVisuals(
  plugin: PluginContext,
  data: { trajectory: StateObjectRef<PluginStateObject.Molecule.Trajectory> },
) {
  return plugin.builders.structure.hierarchy.applyPreset(
    data.trajectory,
    plugin.config.get(PluginConfig.Structure.DefaultHierarchyPreset) ?? '',
  );
}

export function directTrajectory<P extends {}>(
  transformer: StateTransformer<
    PluginStateObject.Data.String | PluginStateObject.Data.Binary,
    PluginStateObject.Molecule.Trajectory,
    P
  >,
  transformerParams?: P,
): Pick<TrajectoryFormatProvider, 'parse' | 'parseRaw'> {
  return {
    parse: async (plugin, data, params) => {
      const state = plugin.state.data;
      const trajectory = await state
        .build()
        .to(data)
        .apply(transformer, transformerParams, { tags: params?.trajectoryTags })
        .commit({ revertOnError: true });
      return { trajectory };
    },
    parseRaw: async (plugin, ctx, data) => {
      const trajectory = await applyTransformerRaw(plugin, ctx, transformer, rawDataObject(data), transformerParams);
      return { trajectory: trajectory.data };
    },
  };
}
