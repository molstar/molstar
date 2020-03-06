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

type TrajectoryFormat = 'pdb' | 'cif' | 'gro' | '3dg'

export enum StructureBuilderTags {
    Trajectory = 'trajectory',
    Model = 'model',
    ModelProperties = 'model-properties',
    Structure = 'structure'
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

    async createModel(trajectory: StateObjectRef<SO.Molecule.Trajectory>, params?: StateTransformer.Params<StateTransforms['Model']['ModelFromTrajectory']>, supportProps?: boolean) {
        const state = this.dataState;
        if (supportProps) {
            const model = state.build().to(trajectory)
                .apply(StateTransforms.Model.ModelFromTrajectory, params || { modelIndex: 0 })
                .apply(StateTransforms.Model.CustomModelProperties, void 0, { tags: [StructureBuilderTags.Model, StructureBuilderTags.ModelProperties] });
            await this.plugin.runTask(this.dataState.updateTree(model, { revertOnError: true }));
            return model.selector;
        } else {
            const model = state.build().to(trajectory)
                .apply(StateTransforms.Model.ModelFromTrajectory, params || { modelIndex: 0 }, { tags: StructureBuilderTags.Model });        
            await this.plugin.runTask(this.dataState.updateTree(model, { revertOnError: true }));
            return model.selector;
        }
    }

    async createStructure(model: StateObjectRef<SO.Molecule.Model>, params?: RootStructureDefinition.Params) {
        const state = this.dataState;
        const structure = state.build().to(model)
            .apply(StateTransforms.Model.StructureFromModel, { type: params || { name: 'assembly', params: { } } }, { tags: StructureBuilderTags.Structure });        
        await this.plugin.runTask(this.dataState.updateTree(structure, { revertOnError: true }));
        return structure.selector;
    }

    // TODO
    async makeComponent() {

    }

    constructor(public plugin: PluginContext) {
    }
}