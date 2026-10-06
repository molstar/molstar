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
import type { PluginContext } from '@molstar/plugin/context';
import { Model } from '@molstar/model/model/structure';
import { OperatorNameColorThemeProvider } from '@molstar/graphics/theme/color/operator-name';
import { ElementSymbolColorThemeProvider } from '@molstar/graphics/theme/color/element-symbol';
import { mergeRegistryEntries } from '@molstar/plugin/registry/merge';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { TrajectoryHierarchyPresetProvider } from './types.js';

const CommonParams = TrajectoryHierarchyPresetProvider.CommonParams;

const CrystalContactsParams = (a: PluginStateObject.Molecule.Trajectory | undefined, plugin: PluginContext) => ({
  model: PD.Optional(PD.Group(StateTransformer.getParamDefinition(ModelFromTrajectory, a, plugin))),
  ...CommonParams(a, plugin),
});

export const CrystalContactsHierarchyPreset = TrajectoryHierarchyPresetProvider({
  id: 'preset-trajectory-crystal-contacts',
  alias: 'crystalContacts',
  display: {
    name: 'Crystal Contacts',
    group: 'Preset',
    description: 'Showsasymetric unit and chains from neighbours within 5 \u212B, i.e., symmetry mates.',
  },
  isApplicable: (o) => {
    return Model.hasCrystalSymmetry(o.data.representative);
  },
  params: CrystalContactsParams,
  async apply(trajectory, params, plugin) {
    const builder = plugin.builders.structure;

    const model = await builder.createModel(trajectory, params.model);
    const modelProperties = await builder.insertModelProperties(model, params.modelProperties);

    const structure = await builder.createStructure(modelProperties || model, {
      name: 'symmetry-mates',
      params: { radius: 5 },
    });
    const structureProperties = await builder.insertStructureProperties(structure, params.structureProperties);

    const unitcell = await builder.tryCreateUnitcell(modelProperties, undefined, { isHidden: true });
    const representationPreset = TrajectoryHierarchyPresetProvider.getRepresentationPreset(plugin, params);
    const representation = await plugin.builders.structure.representation.applyPreset(
      structureProperties,
      representationPreset,
      {
        theme: {
          globalName: 'operator-name',
          carbonColor: 'operator-name',
          focus: {
            name: 'element-symbol',
            params: { carbonColor: { name: 'operator-name', params: OperatorNameColorThemeProvider.defaultValues } },
          },
        },
      },
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

/** The crystal contacts hierarchy preset with the `operator-name` and `element-symbol` color themes it applies. */
export const CrystalContactsHierarchyPresetEntry: PluginRegistryEntry = mergeRegistryEntries(
  { structure: { presets: { hierarchy: [CrystalContactsHierarchyPreset] } } },
  { structure: { themes: { color: [OperatorNameColorThemeProvider, ElementSymbolColorThemeProvider] } } },
);
