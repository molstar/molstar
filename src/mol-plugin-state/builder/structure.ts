/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from '../../mol-plugin/context';
import { StateObjectRef, StateObjectSelector, StateTransformer } from '../../mol-state';
import { PluginStateObject as SO } from '../objects';
import { StateTransforms } from '../transforms';
import { RootStructureDefinition } from '../helpers/root-structure';
import { StructureComponentParams, StaticStructureComponentType } from '../helpers/structure-component';
import { BuildInTrajectoryFormat, TrajectoryFormatProvider } from '../formats/trajectory';
import { StructureRepresentationBuilder } from './structure/representation';
import { StructureSelectionQuery } from '../helpers/structure-selection-query';
import { Task } from '../../mol-task';

export type TrajectoryFormat = 'pdb' | 'cif' | 'gro' | '3dg'

export enum StructureBuilderTags {
    Trajectory = 'trajectory',
    Model = 'model',
    ModelProperties = 'model-properties',
    Structure = 'structure',
    StructureProperties = 'structure-properties',
    Component = 'structure-component'
}

export class StructureBuilder {
    private get dataState() {
        return this.plugin.state.dataState;
    }

    private async parseTrajectoryData(data: StateObjectRef<SO.Data.Binary | SO.Data.String>, format: BuildInTrajectoryFormat | TrajectoryFormatProvider) {
        const provider = typeof format === 'string' ? this.plugin.dataFormat.trajectory.get(format) : format;
        if (!provider) throw new Error(`'${format}' is not a supported data format.`);
        const { trajectory } = await provider.parse(this.plugin, data, { trajectoryTags: StructureBuilderTags.Trajectory });
        return trajectory;
    }

    private async parseTrajectoryBlob(data: StateObjectRef<SO.Data.Blob>, params: StateTransformer.Params<StateTransforms['Data']['ParseBlob']>) {
        const state = this.dataState;
        const trajectory = state.build().to(data)
            .apply(StateTransforms.Data.ParseBlob, params, { state: { isGhost: true } })
            .apply(StateTransforms.Model.TrajectoryFromBlob, void 0, { tags: StructureBuilderTags.Trajectory });        
        await this.plugin.runTask(this.dataState.updateTree(trajectory, { revertOnError: true }));
        return trajectory.selector;
    }

    readonly representation = new StructureRepresentationBuilder(this.plugin);

    async parseStructure(params: {
        data?: StateObjectRef<SO.Data.Binary | SO.Data.String>,
        dataFormat?: BuildInTrajectoryFormat | TrajectoryFormatProvider,
        blob?: StateObjectRef<SO.Data.Blob>
        blobParams?: StateTransformer.Params<StateTransforms['Data']['ParseBlob']>,
        model?: StateTransformer.Params<StateTransforms['Model']['ModelFromTrajectory']>,
        modelProperties?: boolean | StateTransformer.Params<StateTransforms['Model']['CustomModelProperties']>,
        structure?: RootStructureDefinition.Params,
        structureProperties?: boolean | StateTransformer.Params<StateTransforms['Model']['CustomStructureProperties']>
    }) {
        const trajectory = params.data 
            ? await this.parseTrajectory(params.data, params.dataFormat! || 'cif')
            : await this.parseTrajectoryBlob(params.blob!, params.blobParams!);
        
        const model = await this.createModel(trajectory, params.model);
        const modelProperties = !!params.modelProperties 
            ? await this.insertModelProperties(model, typeof params?.modelProperties !== 'boolean' ? params?.modelProperties : void 0) : void 0;
        
        const structure = await this.createStructure(modelProperties || model, params.structure);
        const structureProperties = !!params.structureProperties 
            ? await this.insertStructureProperties(structure, typeof params?.structureProperties !== 'boolean' ? params?.structureProperties : void 0) : void 0;
    
        return {
            trajectory,
            model: modelProperties || model,
            modelRoot: model,
            modelProperties,
            structure: structureProperties || structure,
            structureRoot: structure,
            structureProperties
        };
    }

    async parseTrajectory(data: StateObjectRef<SO.Data.Binary | SO.Data.String>, format: BuildInTrajectoryFormat | TrajectoryFormatProvider): Promise<StateObjectSelector<SO.Molecule.Trajectory>>
    async parseTrajectory(blob: StateObjectRef<SO.Data.Blob>, params: StateTransformer.Params<StateTransforms['Data']['ParseBlob']>): Promise<StateObjectSelector<SO.Molecule.Trajectory>>
    async parseTrajectory(data: StateObjectRef, params: any) {
        const cell = StateObjectRef.resolveAndCheck(this.dataState, data as StateObjectRef);
        if (!cell) throw new Error('Invalid data cell.');

        if (SO.Data.Blob.is(cell.obj)) {
            return this.parseTrajectoryBlob(data, params);
        } else {
            return this.parseTrajectoryData(data, params);
        }
    }

    async createModel(trajectory: StateObjectRef<SO.Molecule.Trajectory>, params?: StateTransformer.Params<StateTransforms['Model']['ModelFromTrajectory']>) {
        const state = this.dataState;
        const model = state.build().to(trajectory)
            .apply(StateTransforms.Model.ModelFromTrajectory, params || { modelIndex: 0 }, { tags: StructureBuilderTags.Model });

        await this.plugin.runTask(this.dataState.updateTree(model, { revertOnError: true }));
        return model.selector;
    }

    async insertModelProperties(model: StateObjectRef<SO.Molecule.Model>, params?: StateTransformer.Params<StateTransforms['Model']['CustomModelProperties']>) {
        const state = this.dataState;
        const props = state.build().to(model)
            .insert(StateTransforms.Model.CustomModelProperties, params, { tags: StructureBuilderTags.ModelProperties, isDecorator: true });
        await this.plugin.runTask(this.dataState.updateTree(props, { revertOnError: true }));
        return props.selector;
    }

    async createStructure(model: StateObjectRef<SO.Molecule.Model>, params?: RootStructureDefinition.Params) {
        const state = this.dataState;
        const structure = state.build().to(model)
            .apply(StateTransforms.Model.StructureFromModel, { type: params || { name: 'assembly', params: { } } }, { tags: StructureBuilderTags.Structure });        

        await this.plugin.runTask(this.dataState.updateTree(structure, { revertOnError: true }));
        return structure.selector;
    }

    async insertStructureProperties(structure: StateObjectRef<SO.Molecule.Structure>, params?: StateTransformer.Params<StateTransforms['Model']['CustomStructureProperties']>) {
        const state = this.dataState;
        const props = state.build().to(structure)
            .insert(StateTransforms.Model.CustomStructureProperties, params, { tags: StructureBuilderTags.StructureProperties, isDecorator: true });
        await this.plugin.runTask(this.dataState.updateTree(props, { revertOnError: true }));
        return props.selector;
    }

    /** returns undefined if the component is empty/null */
    async tryCreateComponent(structure: StateObjectRef<SO.Molecule.Structure>, params: StructureComponentParams, key: string, tags?: string[]): Promise<StateObjectRef<SO.Molecule.Structure> | undefined> {
        const state = this.dataState;

        const root = state.build().to(structure);

        const keyTag = `structure-component-${key}`;
        const component = root.applyOrUpdateTagged(keyTag, StateTransforms.Model.StructureComponent, params, { 
            tags: tags ? [...tags, StructureBuilderTags.Component, keyTag] : [StructureBuilderTags.Component, keyTag]
        });

        await this.plugin.runTask(this.dataState.updateTree(component));

        const selector = component.selector;

        if (!selector.isOk || selector.cell?.obj?.data.elementCount === 0) {
            const del = state.build().delete(selector.ref);
            await this.plugin.runTask(this.dataState.updateTree(del));
            return;
        }

        return selector;
    }

    tryCreateStaticComponent(params: { structure: StateObjectRef<SO.Molecule.Structure>, type: StaticStructureComponentType, key: string, label?: string, tags?: string[] }) { 
        return this.tryCreateComponent(params.structure, {
            type: { name: 'static', params: params.type },
            nullIfEmpty: true,
            label: ''
        }, params.key, params.tags);
    }

    tryCreateQueryComponent(params: { structure: StateObjectRef<SO.Molecule.Structure>, query: StructureSelectionQuery, key: string, label?: string, tags?: string[] }): Promise<StateObjectRef<SO.Molecule.Structure> | undefined> {
        return this.plugin.runTask(Task.create('Query Component', async taskCtx => {
            let { structure, query, key, label, tags } = params;        
            label = (label || '').trim();

            const structureData = StateObjectRef.resolveAndCheck(this.dataState, structure)?.obj?.data;

            if (!structureData) return;
    
            const transformParams: StructureComponentParams = query.referencesCurrent
                ? {
                    type: { name: 'bundle', params: await StructureSelectionQuery.getBundle(this.plugin, taskCtx, query, structureData) },
                    nullIfEmpty: true,
                    label: label || query.label
                } : {
                    type: { name: 'expression', params: query.expression },
                    nullIfEmpty: true,
                    label: label || query.label
                };

            if (query.ensureCustomProperties) {
                await query.ensureCustomProperties({ fetch: this.plugin.fetch, runtime: taskCtx }, structureData);
            }
            
            const state = this.dataState;
            const root = state.build().to(structure);
            const keyTag = `structure-component-${key}`;
            const component = root.applyOrUpdateTagged(keyTag, StateTransforms.Model.StructureComponent, transformParams, { 
                tags: tags ? [...tags, StructureBuilderTags.Component, keyTag] : [StructureBuilderTags.Component, keyTag]
            });
    
            await this.dataState.updateTree(component).runInContext(taskCtx);
    
            const selector = component.selector;
    
            if (!selector.isOk || selector.cell?.obj?.data.elementCount === 0) {
                const del = state.build().delete(selector.ref);
                await this.plugin.runTask(this.dataState.updateTree(del));
                return;
            }
    
            return selector; 
        }))
    }

    constructor(public plugin: PluginContext) {
    }
}