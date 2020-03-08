/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from '../../mol-plugin/context';
import { StateBuilder, StateObjectRef, StateObjectSelector, StateTransformer } from '../../mol-state';
import { PluginStateObject as SO } from '../objects';
import { StateTransforms } from '../transforms';
import { RootStructureDefinition } from '../helpers/root-structure';
import { StructureComponentParams } from '../helpers/structure-component';

type TrajectoryFormat = 'pdb' | 'cif' | 'gro' | '3dg'

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

    private async parseTrajectoryData(data: StateObjectRef<SO.Data.Binary | SO.Data.String>, format: TrajectoryFormat) {
        const state = this.dataState;
        const root = state.build().to(data);
        let parsed: StateBuilder.To<SO.Molecule.Trajectory>;
        const tag = { tags: StructureBuilderTags.Trajectory };
        switch (format) {
            case 'cif':
                parsed = root.apply(StateTransforms.Data.ParseCif, void 0, { state: { isGhost: true } })
                    .apply(StateTransforms.Model.TrajectoryFromMmCif, void 0, tag)
                break
            case 'pdb':
                parsed = root.apply(StateTransforms.Model.TrajectoryFromPDB, void 0, tag);
                break
            case 'gro':
                parsed = root.apply(StateTransforms.Model.TrajectoryFromGRO, void 0, tag);
                break
            case '3dg':
                parsed = root.apply(StateTransforms.Model.TrajectoryFrom3DG, void 0, tag);
                break
            default:
                throw new Error('unsupported format')
        }

        await this.plugin.runTask(this.dataState.updateTree(parsed, { revertOnError: true }));

        return parsed.selector;
    }

    private async parseTrajectoryBlob(data: StateObjectRef<SO.Data.Blob>, params: StateTransformer.Params<StateTransforms['Data']['ParseBlob']>) {
        const state = this.dataState;
        const trajectory = state.build().to(data)
            .apply(StateTransforms.Data.ParseBlob, params, { state: { isGhost: true } })
            .apply(StateTransforms.Model.TrajectoryFromBlob, void 0, { tags: StructureBuilderTags.Trajectory });        
        await this.plugin.runTask(this.dataState.updateTree(trajectory, { revertOnError: true }));
        return trajectory.selector;
    }

    async parseStructure(params: {
        data?: StateObjectRef<SO.Data.Binary | SO.Data.String>,
        dataFormat?: TrajectoryFormat,
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
            modelBase: model,
            modelProperties,
            structure: structureProperties || structure,
            structureBase: structure,
            structureProperties
        };
    }

    async parseTrajectory(data: StateObjectRef<SO.Data.Binary | SO.Data.String>, format: TrajectoryFormat): Promise<StateObjectSelector<SO.Molecule.Trajectory>>
    async parseTrajectory(blob: StateObjectRef<SO.Data.Blob>, params: StateTransformer.Params<StateTransforms['Data']['ParseBlob']>): Promise<StateObjectSelector<SO.Molecule.Trajectory>>
    async parseTrajectory(data: StateObjectRef, params: any) {
        // TODO: proper format support
        // needs to integrated to transforms directly because of blobs
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

        // const props = !!params?.properties 
        //     ? model.apply(StateTransforms.Model.CustomModelProperties, typeof params?.properties !== 'boolean' ? params?.properties : void 0, { tags: StructureBuilderTags.ModelProperties, isDecorator: true })
        //     : void 0;

        await this.plugin.runTask(this.dataState.updateTree(model, { revertOnError: true }));

        return model.selector;

        // const modelSelector = model.selector, propertiesSelector = props?.selector;

        // return { model: propertiesSelector || modelSelector, index: modelSelector, properties: propertiesSelector };
    }

    async insertModelProperties(model: StateObjectRef<SO.Molecule.Model>, params?: StateTransformer.Params<StateTransforms['Model']['CustomModelProperties']>) {
        const state = this.dataState;
        const props = state.build().to(model)
            .apply(StateTransforms.Model.CustomModelProperties, params, { tags: StructureBuilderTags.ModelProperties, isDecorator: true });
        await this.plugin.runTask(this.dataState.updateTree(props, { revertOnError: true }));
        return props.selector;
    }

    async createStructure(model: StateObjectRef<SO.Molecule.Model>, params?: RootStructureDefinition.Params) {
        const state = this.dataState;
        const structure = state.build().to(model)
            .apply(StateTransforms.Model.StructureFromModel, { type: params || { name: 'assembly', params: { } } }, { tags: StructureBuilderTags.Structure });        

        // const props = !!params?.properties 
        //     ? structure.apply(StateTransforms.Model.CustomStructureProperties, typeof params?.properties !== 'boolean' ? params?.properties : void 0, { tags: StructureBuilderTags.StructureProperties, isDecorator: true })
        //     : void 0;

        await this.plugin.runTask(this.dataState.updateTree(structure, { revertOnError: true }));

        return structure.selector;

        // const structureSelector = structure.selector, propertiesSelector = props?.selector;

        // return { structure: propertiesSelector || structureSelector, definition: structureSelector, properties: propertiesSelector };
    }

    async insertStructureProperties(structure: StateObjectRef<SO.Molecule.Structure>, params?: StateTransformer.Params<StateTransforms['Model']['CustomStructureProperties']>) {
        const state = this.dataState;
        const props = state.build().to(structure)
            .apply(StateTransforms.Model.CustomStructureProperties, params, { tags: StructureBuilderTags.StructureProperties, isDecorator: true });
        await this.plugin.runTask(this.dataState.updateTree(props, { revertOnError: true }));
        return props.selector;
    }

    /** returns undefined if the component is empty/null */
    async tryCreateComponent(structure: StateObjectRef<SO.Molecule.Structure>, params: StructureComponentParams, tag?: string): Promise<StateObjectRef<SO.Molecule.Structure> | undefined> {
        const state = this.dataState;

        const root = state.build().to(structure);
        let component: StateBuilder.To<SO.Molecule.Structure>;

        if (tag) {
            const typeTag = `structure-component-${tag}`;
            component = root.applyOrUpdateTagged(typeTag, StateTransforms.Model.StructureComponent, params, { tags: [StructureBuilderTags.Component, typeTag] });
        } else {
            component = root.apply(StateTransforms.Model.StructureComponent, params, { tags: StructureBuilderTags.Component });
        }

        await this.plugin.runTask(this.dataState.updateTree(component));

        const selector = component.selector;

        if (!selector.isOk || selector.cell?.obj?.data.elementCount === 0) {
            const del = state.build().delete(selector.ref);
            await this.plugin.runTask(this.dataState.updateTree(del));
            return;
        }

        return selector;
    }

    constructor(public plugin: PluginContext) {
    }
}