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

    async createModel(trajectory: StateObjectRef<SO.Molecule.Trajectory>, params?: {
        model?: StateTransformer.Params<StateTransforms['Model']['ModelFromTrajectory']>,
        properties?: boolean | StateTransformer.Params<StateTransforms['Model']['CustomModelProperties']>
    }) {
        const state = this.dataState;

        const model = state.build().to(trajectory)
            .apply(StateTransforms.Model.ModelFromTrajectory, params?.model || void 0, { tags: StructureBuilderTags.Model });

        const props = !!params?.properties 
            ? model.apply(StateTransforms.Model.CustomModelProperties, typeof params?.properties !== 'boolean' ? params?.properties : void 0, { tags: StructureBuilderTags.ModelProperties, isDecorator: true })
            : void 0;

        await this.plugin.runTask(this.dataState.updateTree(model, { revertOnError: true }));

        const modelSelector = model.selector, propertiesSelector = props?.selector;

        return { model: propertiesSelector || modelSelector, index: modelSelector, properties: propertiesSelector };
    }

    async createStructure(model: StateObjectRef<SO.Molecule.Model>, params?: {
        structure?: RootStructureDefinition.Params,
        properties?: boolean | StateTransformer.Params<StateTransforms['Model']['CustomStructureProperties']>
    }) {
        const state = this.dataState;
        const structure = state.build().to(model)
            .apply(StateTransforms.Model.StructureFromModel, { type: params?.structure || { name: 'assembly', params: { } } }, { tags: StructureBuilderTags.Structure });        

        const props = !!params?.properties 
            ? structure.apply(StateTransforms.Model.CustomStructureProperties, typeof params?.properties !== 'boolean' ? params?.properties : void 0, { tags: StructureBuilderTags.StructureProperties, isDecorator: true })
            : void 0;

        await this.plugin.runTask(this.dataState.updateTree(structure, { revertOnError: true }));

        const structureSelector = structure.selector, propertiesSelector = props?.selector;

        return { structure: propertiesSelector || structureSelector, definition: structureSelector, properties: propertiesSelector };
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