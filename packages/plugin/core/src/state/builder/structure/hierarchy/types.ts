/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { PresetProvider } from '../../preset-provider.js';
import type { PluginStateObject } from '../../../objects.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { StateTransformer } from '@molstar/core/state';
import { CustomModelProperties, CustomStructureProperties } from '@molstar/plugin/state/transforms/structure/hierarchy';
import type { PluginContext } from '@molstar/plugin/context';
import type { BuiltInStructureRepresentationPresetAlias } from '../presets/catalog.js';

export interface TrajectoryHierarchyPresetProvider<
  P = any,
  S = {},
  Id extends string = string,
  Alias extends string = string,
> extends PresetProvider<PluginStateObject.Molecule.Trajectory, P, S, Id, Alias> {}
export function TrajectoryHierarchyPresetProvider<P, S, const Id extends string, const Alias extends string = string>(
  preset: TrajectoryHierarchyPresetProvider<P, S, Id, Alias>,
) {
  return preset;
}
export namespace TrajectoryHierarchyPresetProvider {
  export type Params<P extends TrajectoryHierarchyPresetProvider> =
    P extends TrajectoryHierarchyPresetProvider<infer T> ? T : never;
  export type State<P extends TrajectoryHierarchyPresetProvider> =
    P extends TrajectoryHierarchyPresetProvider<infer _, infer S> ? S : never;

  export const CommonParams = (a: PluginStateObject.Molecule.Trajectory | undefined, plugin: PluginContext) => ({
    modelProperties: PD.Optional(PD.Group(StateTransformer.getParamDefinition(CustomModelProperties, void 0, plugin))),
    structureProperties: PD.Optional(
      PD.Group(StateTransformer.getParamDefinition(CustomStructureProperties, void 0, plugin)),
    ),
    representationPreset: PD.Optional(PD.Text<BuiltInStructureRepresentationPresetAlias>('auto' as const)),
  });
}
