/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { StateObjectRef } from '@molstar/core/state';
import type { PluginStateObject } from '../../../objects.js';
import type { PluginContext } from '@molstar/plugin/context';
import { getStructureQuality } from '@molstar/graphics/repr/util';
import { PluginConfig } from '@molstar/plugin/config';
import { StructureRepresentationPresetProvider } from '../representation-presets/types.js';
import { AutoPreset } from '../representation-presets/auto.js';
import { DefaultHierarchyPreset } from './default.js';
import { TrajectoryHierarchyPresetProvider } from './types.js';

const CommonParams = TrajectoryHierarchyPresetProvider.CommonParams;

const AllModelsParams = (a: PluginStateObject.Molecule.Trajectory | undefined, plugin: PluginContext) => ({
  useDefaultIfSingleModel: PD.Optional(PD.Boolean(false)),
  representationPresetParams: PD.Optional(PD.Group(StructureRepresentationPresetProvider.CommonParams)),
  ...CommonParams(a, plugin),
});

export const AllModelsHierarchyPreset = TrajectoryHierarchyPresetProvider({
  id: 'preset-trajectory-all-models',
  alias: 'all-models',
  display: {
    name: 'All Models',
    group: 'Preset',
    description: 'Shows all models; colored by trajectory-index.',
  },
  isApplicable: (o) => {
    return o.data.frameCount > 1;
  },
  params: AllModelsParams,
  async apply(trajectory, params, plugin) {
    const tr = StateObjectRef.resolveAndCheck(plugin.state.data, trajectory)?.obj?.data;
    if (!tr) return {};

    if (tr.frameCount === 1 && params.useDefaultIfSingleModel) {
      return DefaultHierarchyPreset.apply(trajectory, params as any, plugin);
    }

    const builder = plugin.builders.structure;

    const models = [],
      structures = [];

    for (let i = 0; i < tr.frameCount; i++) {
      const model = await builder.createModel(trajectory, { modelIndex: i });
      const modelProperties = await builder.insertModelProperties(model, params.modelProperties, { isCollapsed: true });
      const structure = await builder.createStructure(modelProperties || model, { name: 'model', params: {} });
      const structureProperties = await builder.insertStructureProperties(structure, params.structureProperties);

      models.push(model);
      structures.push(structure);

      const quality = structure.obj
        ? getStructureQuality(structure.obj.data, { elementCountFactor: tr.frameCount })
        : 'medium';
      const representationPreset =
        params.representationPreset ||
        plugin.config.get(PluginConfig.Structure.DefaultRepresentationPreset) ||
        AutoPreset.id;
      await builder.representation.applyPreset(structureProperties, representationPreset, {
        theme: { globalName: 'trajectory-index' },
        quality,
      });
    }

    return { models, structures };
  },
});
