/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateTransform } from '../base';
import { PluginStateObjects as SO } from '../objects';
import { Task } from 'mol-task';
import { Model, Format, Structure } from 'mol-model/structure';

export const CreateModelsFromMmCif = PluginStateTransform.Create<SO.Data.Cif, SO.Models, { blockHeader?: string }>({
    name: 'create-models-from-mmcif',
    from: [SO.Data.Cif],
    to: [SO.Models],
    defaultParams: a => ({ blockHeader: a.data.blocks[0].header }),
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

export const CreateStructureFromModel = PluginStateTransform.Create<SO.Models, SO.Structure, { modelIndex: number }>({
    name: 'structure-from-model',
    from: [SO.Models],
    to: [SO.Structure],
    defaultParams: () => ({ modelIndex: 0 }),
    apply({ a, params }) {
        if (params.modelIndex < 0 || params.modelIndex >= a.data.length) throw new Error(`Invalid modelIndex ${params.modelIndex}`);
        // TODO: make Structure.ofModel async?
        const s = Structure.ofModel(a.data[params.modelIndex]);
        return new SO.Structure({ label: `${a.data[params.modelIndex].label} (model ${s.models[0].modelNum})`, desctiption: s.elementCount === 1 ? '1 element' : `${s.elementCount} elements` }, s);
    }
});