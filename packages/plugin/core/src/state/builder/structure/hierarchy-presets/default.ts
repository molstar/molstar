/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { StateTransformer } from '@molstar/core/state';
import { ModelFromTrajectory } from '@molstar/plugin/state/transforms/structure/hierarchy';
import type { PluginStateObject } from '../../../objects.js';
import { RootStructureDefinition } from '../../../helpers/root-structure.js';
import type { PluginContext } from '@molstar/plugin/context';
import { StructureRepresentationPresetProvider } from '../representation-presets/types.js';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { TrajectoryHierarchyPresetProvider } from './types.js';

const CommonParams = TrajectoryHierarchyPresetProvider.CommonParams;

const DefaultParams = (a: PluginStateObject.Molecule.Trajectory | undefined, plugin: PluginContext) => ({
  model: PD.Optional(PD.Group(StateTransformer.getParamDefinition(ModelFromTrajectory, a, plugin))),
  showUnitcell: PD.Optional(PD.Boolean(false)),
  structure: PD.Optional(RootStructureDefinition.getParams(void 0, 'assembly').type),
  representationPresetParams: PD.Optional(PD.Group(StructureRepresentationPresetProvider.CommonParams)),
  ...CommonParams(a, plugin),
});

export const DefaultHierarchyPreset = TrajectoryHierarchyPresetProvider({
  id: 'preset-trajectory-default',
  alias: 'default',
  display: {
    name: 'Default (Assembly)',
    group: 'Preset',
    description: 'Shows the first assembly or, if that is unavailable, the first model.',
  },
  isApplicable: (o) => {
    return true;
  },
  params: DefaultParams,
  async apply(trajectory, params, plugin) {
    const builder = plugin.builders.structure;

    const model = await builder.createModel(trajectory, params.model);
    const modelProperties = await builder.insertModelProperties(model, params.modelProperties);

    const structure = await builder.createStructure(modelProperties || model, params.structure);
    const structureProperties = await builder.insertStructureProperties(structure, params.structureProperties);

    const unitcell =
      params.showUnitcell === void 0 || !!params.showUnitcell
        ? await builder.tryCreateUnitcell(modelProperties, undefined, { isHidden: true })
        : void 0;
    const representationPreset = TrajectoryHierarchyPresetProvider.getRepresentationPreset(plugin, params);
    const representation = await plugin.builders.structure.representation.applyPreset(
      structureProperties,
      representationPreset,
      params.representationPresetParams,
    );

    return {
      model,
      modelProperties,
      unitcell,
      structure,
      structureProperties,
      representation,
    };
  },
});

/** The default hierarchy preset. The representation preset it applies comes from config, so the entry has none. */
export const DefaultHierarchyPresetEntry: PluginRegistryEntry = {
  structure: { presets: { hierarchy: [DefaultHierarchyPreset] } },
};
