/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateTransform } from '../base';
import { PluginStateObjects as SO } from '../objects';
import { Task } from 'mol-task';
import { Model, Format, Structure } from 'mol-model/structure';

export { CreateModelsFromMmCif }
namespace CreateModelsFromMmCif { export interface Params { blockHeader?: string } }
const CreateModelsFromMmCif = PluginStateTransform.Create<SO.Data.Cif, SO.Models, CreateModelsFromMmCif.Params>({
    name: 'create-models-from-mmcif',
    display: {
        name: 'Models from mmCIF',
        description: 'Identify and create all separate models in the specified CIF data block'
    },
    from: [SO.Data.Cif],
    to: [SO.Models],
    params: { default: a => ({ blockHeader: a.data.blocks[0].header }) },
    apply({ a, params }) {
        return Task.create('Parse mmCIF', async ctx => {
            const header = params.blockHeader || a.data.blocks[0].header;
            const block = a.data.blocks.find(b => b.header === header);
            if (!block) throw new Error(`Data block '${[header]}' not found.`);
            const models = await Model.create(Format.mmCIF(block)).runInContext(ctx);
            if (models.length === 0) throw new Error('No models found.');
            const label = models.length === 1 ? `${models[0].label}` : `${models[0].label} (${models.length} models)`;
            return new SO.Models({ label }, models);
        });
    }
});

export { CreateStructureFromModel }
namespace CreateStructureFromModel { export interface Params { modelIndex: number } }
const CreateStructureFromModel = PluginStateTransform.Create<SO.Models, SO.Structure, CreateStructureFromModel.Params>({
    name: 'structure-from-model',
    display: {
        name: 'Structure from Model',
        description: 'Create a molecular structure from the specified model.'
    },
    from: [SO.Models],
    to: [SO.Structure],
    params: { default: () => ({ modelIndex: 0 }) },
    apply({ a, params }) {
        if (params.modelIndex < 0 || params.modelIndex >= a.data.length) throw new Error(`Invalid modelIndex ${params.modelIndex}`);
        // TODO: make Structure.ofModel async?
        const s = Structure.ofModel(a.data[params.modelIndex]);
        return new SO.Structure({ label: `${a.data[params.modelIndex].label} (model ${s.models[0].modelNum})`, desctiption: s.elementCount === 1 ? '1 element' : `${s.elementCount} elements` }, s);
    }
});