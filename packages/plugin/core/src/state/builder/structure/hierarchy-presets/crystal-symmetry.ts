/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type StateObjectRef, StateTransformer } from '@molstar/core/state';
import { ModelFromTrajectory } from '@molstar/plugin/state/transforms/structure/hierarchy';
import type { PluginStateObject } from '../../../objects.js';
import type { PluginContext } from '@molstar/plugin/context';
import type { Vec3 } from '@molstar/core/math/linear-algebra';
import { TrajectoryHierarchyPresetProvider } from './types.js';

const CommonParams = TrajectoryHierarchyPresetProvider.CommonParams;

export const CrystalSymmetryParams = (a: PluginStateObject.Molecule.Trajectory | undefined, plugin: PluginContext) => ({
  model: PD.Optional(PD.Group(StateTransformer.getParamDefinition(ModelFromTrajectory, a, plugin))),
  ...CommonParams(a, plugin),
});

export async function applyCrystalSymmetry(
  props: { ijkMin: Vec3; ijkMax: Vec3; theme?: string },
  trajectory: StateObjectRef<PluginStateObject.Molecule.Trajectory>,
  params: PD.ValuesFor<ReturnType<typeof CrystalSymmetryParams>>,
  plugin: PluginContext,
) {
  const builder = plugin.builders.structure;

  const model = await builder.createModel(trajectory, params.model);
  const modelProperties = await builder.insertModelProperties(model, params.modelProperties);

  const structure = await builder.createStructure(modelProperties || model, {
    name: 'symmetry',
    params: props,
  });
  const structureProperties = await builder.insertStructureProperties(structure, params.structureProperties);

  const unitcell = await builder.tryCreateUnitcell(modelProperties, undefined, { isHidden: false });
  const representationPreset = TrajectoryHierarchyPresetProvider.getRepresentationPreset(plugin, params);
  const representation = await plugin.builders.structure.representation.applyPreset(
    structureProperties,
    representationPreset,
    { theme: { globalName: props.theme } },
  );

  return {
    model,
    modelProperties,
    unitcell,
    structure,
    structureProperties,
    representation,
  };
}
