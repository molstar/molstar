/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PresetProvider } from '../preset-provider';
import { PluginStateObject } from '../../objects';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { StateObjectRef, StateTransformer } from '../../../mol-state';
import { StateTransforms } from '../../transforms';
import { RootStructureDefinition } from '../../helpers/root-structure';
import { PresetStructureReprentations } from './representation-preset';
import { PluginContext } from '../../../mol-plugin/context';
import { isProductionMode } from '../../../mol-util/debug';
import { Task } from '../../../mol-task';

export interface TrajectoryHierarchyPresetProvider<P = any, S = {}> extends PresetProvider<PluginStateObject.Molecule.Trajectory, P, S> { }
export namespace TrajectoryHierarchyPresetProvider {
    export type Params<P extends TrajectoryHierarchyPresetProvider> = P extends TrajectoryHierarchyPresetProvider<infer T> ? T : never;
    export type State<P extends TrajectoryHierarchyPresetProvider> = P extends TrajectoryHierarchyPresetProvider<infer _, infer S> ? S : never;
}
export function TrajectoryHierarchyPresetProvider<P, S>(preset: TrajectoryHierarchyPresetProvider<P, S>) { return preset; }

const CommonParams = (a: PluginStateObject.Molecule.Trajectory | undefined, plugin: PluginContext) => ({
    modelProperties: PD.Optional(PD.Group(StateTransformer.getParamDefinition(StateTransforms.Model.CustomModelProperties, void 0, plugin))),
    structureProperties: PD.Optional(PD.Group(StateTransformer.getParamDefinition(StateTransforms.Model.CustomStructureProperties, void 0, plugin))),
    representationPreset: PD.Optional(PD.Text<keyof PresetStructureReprentations>('auto' as const)),
})

const FirstModelParams = (a: PluginStateObject.Molecule.Trajectory | undefined, plugin: PluginContext) =>  ({
    model: PD.Optional(PD.Group(StateTransformer.getParamDefinition(StateTransforms.Model.ModelFromTrajectory, a, plugin))),
    showUnitcell: PD.Optional(PD.Boolean(true)),
    structure: PD.Optional(RootStructureDefinition.getParams(void 0, 'assembly').type),
    ...CommonParams(a, plugin)
});

const firstModel = TrajectoryHierarchyPresetProvider({
    id: 'preset-trajectory-first-model',
    display: { name: 'First Model', group: 'Preset' },
    params: FirstModelParams,
    async apply(trajectory, params, plugin) {
        const builder = plugin.builders.structure;

        const model = await builder.createModel(trajectory, params.model);
        const modelProperties = await builder.insertModelProperties(model, params.modelProperties);

        const structure = await builder.createStructure(modelProperties || model, params.structure);
        const structureProperties = await builder.insertStructureProperties(structure, params.structureProperties);

        const unitcell = params.showUnitcell === void 0 || !!params.showUnitcell ? await builder.tryCreateUnitcell(modelProperties, undefined, { isHidden: true }) : void 0;
        const representation =  await plugin.builders.structure.representation.applyPreset(structureProperties, params.representationPreset || 'auto');

        return {
            model: modelProperties,
            modelRoot: model,
            modelProperties,
            unitcell,
            structure: structureProperties,
            structureRoot: structure,
            structureProperties,
            representation
        };
    }
});

const allModels = TrajectoryHierarchyPresetProvider({
    id: 'preset-trajectory-all-models',
    display: { name: 'All Models', group: 'Preset' },
    params: CommonParams,
    async apply(trajectory, params, plugin) {
        const tr = StateObjectRef.resolveAndCheck(plugin.state.data, trajectory)?.obj?.data;
        if (!tr) return { };

        const builder = plugin.builders.structure;

        const models = [], structures = [];

        for (let i = 0; i < tr.length; i++) {
            const model = await builder.createModel(trajectory, { modelIndex: i }, { isCollapsed: true });
            const modelProperties = await builder.insertModelProperties(model, params.modelProperties);
            const structure = await builder.createStructure(modelProperties || model, { name: 'deposited', params: {} });
            const structureProperties = await builder.insertStructureProperties(structure, params.structureProperties);

            models.push(model);
            structures.push(structure);
            await builder.representation.applyPreset(structureProperties, params.representationPreset || 'auto', { globalThemeName: 'model-index' });
        }

        return { models, structures };
    }
});

export const PresetStructureTrajectoryHierarchy = {
    'first-model': firstModel,
    'all-models': allModels
};
export type PresetStructureTrajectoryHierarchy = typeof PresetStructureTrajectoryHierarchy;

// TODO: should there be a registry like for representations?

export function applyTrajectoryHierarchyPreset<K extends keyof PresetStructureTrajectoryHierarchy>(plugin: PluginContext, parent: StateObjectRef<PluginStateObject.Molecule.Trajectory>, preset: K, params?: Partial<TrajectoryHierarchyPresetProvider.Params<PresetStructureTrajectoryHierarchy[K]>>): Promise<TrajectoryHierarchyPresetProvider.State<PresetStructureTrajectoryHierarchy[K]>> | undefined
export function applyTrajectoryHierarchyPreset<P = any, S = {}>(plugin: PluginContext, parent: StateObjectRef<PluginStateObject.Molecule.Trajectory>, provider: TrajectoryHierarchyPresetProvider<P, S>, params?: P): Promise<S> | undefined
export function applyTrajectoryHierarchyPreset(plugin: PluginContext, parent: StateObjectRef, providerRef: string | TrajectoryHierarchyPresetProvider, params?: any): Promise<any> | undefined {
    const provider = typeof providerRef === 'string' ? (PresetStructureTrajectoryHierarchy as any)[providerRef] : providerRef;
    if (!provider) return;

    const state = plugin.state.data;
    const cell = StateObjectRef.resolveAndCheck(state, parent);
    if (!cell) {
        if (!isProductionMode) console.warn(`Applying hierarchy preset provider to bad cell.`);
        return;
    }
    const prms = { ...PD.getDefaultValues(provider.params(cell.obj!, plugin) as PD.Params), ...params };
    const task = Task.create(`${provider.display.name}`, () => provider.apply(cell, prms, plugin) as Promise<any>);
    return plugin.runTask(task);
}