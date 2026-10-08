/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { PluginContext } from '../../mol-plugin/context';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { StateTransformer } from '../../mol-state';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { RootStructureDefinition } from '../../mol-plugin-state/helpers/root-structure';
import { PresetStructureRepresentations, StructureRepresentationPresetProvider } from '../../mol-plugin-state/builder/structure/representation-preset';
import { TrajectoryHierarchyPresetProvider } from '../../mol-plugin-state/builder/structure/hierarchy-preset';
import { PluginConfig } from '../../mol-plugin/config';
import { ModelBondProvider } from '../../mol-model/structure/structure/unit/bonds/bond-provider';
import { BondOrderProviderName } from './provider';
import { BondOrdersMode } from './perceiver';

const BondOrdersParams = (a: PluginStateObject.Molecule.Trajectory | undefined, plugin: PluginContext) => ({
    model: PD.Optional(PD.Group(StateTransformer.getParamDefinition(StateTransforms.Model.ModelFromTrajectory, a, plugin))),
    showUnitcell: PD.Optional(PD.Boolean(false)),
    structure: PD.Optional(RootStructureDefinition.getParams(void 0, 'assembly').type),
    representationPresetParams: PD.Optional(PD.Group(StructureRepresentationPresetProvider.CommonParams)),
    bondOrdersMode: PD.Select<BondOrdersMode>('auto', [
        ['model', 'Model'],
        ['auto', 'Auto'],
        ['forceCompute', 'Force Compute'],
    ]),
    ...TrajectoryHierarchyPresetProvider.CommonParams(a, plugin)
});

export const BondOrdersTrajectoryPreset = TrajectoryHierarchyPresetProvider({
    id: 'preset-trajectory-bond-orders',
    display: {
        name: 'Bond Orders', group: 'Preset',
        description: 'Default assembly hierarchy with perceived missing bond orders exposed through unit.bonds.'
    },
    isApplicable: () => true,
    params: BondOrdersParams,
    async apply(trajectory, params, plugin) {
        const builder = plugin.builders.structure;

        const model = await builder.createModel(trajectory, params.model);
        const defaultModelProperties = PD.getDefaultValues(plugin.customModelProperties.getParams(model.obj?.data));
        const requestedModelProperties = params.modelProperties ?? defaultModelProperties;
        const modelProperties = await builder.insertModelProperties(model, {
            ...requestedModelProperties,
            autoAttach: Array.from(new Set([
                ...requestedModelProperties.autoAttach,
                ModelBondProvider.Descriptor.name,
            ])),
            properties: {
                ...requestedModelProperties.properties,
                [ModelBondProvider.Descriptor.name]: {
                    provider: {
                        name: BondOrderProviderName,
                        params: { mode: params.bondOrdersMode },
                    },
                },
            },
        });
        const structure = await builder.createStructure(modelProperties || model, params.structure);
        const structureProperties = await builder.insertStructureProperties(structure, params.structureProperties);

        const unitcell = params.showUnitcell === void 0 || !!params.showUnitcell ? await builder.tryCreateUnitcell(modelProperties, undefined, { isHidden: true }) : void 0;
        const representationPreset = params.representationPreset || plugin.config.get(PluginConfig.Structure.DefaultRepresentationPreset) || PresetStructureRepresentations.auto.id;
        const representation = await plugin.builders.structure.representation.applyPreset(structureProperties, representationPreset, params.representationPresetParams);

        return {
            model,
            modelProperties,
            unitcell,
            structure,
            structureProperties,
            representation
        };
    }
});
