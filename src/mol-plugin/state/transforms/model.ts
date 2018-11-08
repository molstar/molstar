/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateTransform } from '../base';
import { PluginStateObjects as SO } from '../objects';
import { Task } from 'mol-task';
import { Model, Format, Structure, ModelSymmetry, StructureSymmetry, QueryContext, StructureSelection } from 'mol-model/structure';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import Expression from 'mol-script/language/expression';
import { compile } from 'mol-script/runtime/query/compiler';
import { Mat4 } from 'mol-math/linear-algebra';

export { ParseModelsFromMmCif }
namespace ParseModelsFromMmCif { export interface Params { blockHeader?: string } }
const ParseModelsFromMmCif = PluginStateTransform.Create<SO.Data.Cif, SO.Molecule.Models, ParseModelsFromMmCif.Params>({
    name: 'parse-models-from-mmcif',
    display: {
        name: 'Models from mmCIF',
        description: 'Identify and create all separate models in the specified CIF data block'
    },
    from: [SO.Data.Cif],
    to: [SO.Molecule.Models],
    params: {
        default: a => ({ blockHeader: a.data.blocks[0].header }),
        controls(a) {
            const { blocks } = a.data;
            if (blocks.length === 0) return {};
            return {
                blockHeader: PD.Select('Header', 'Header of the block to parse', blocks[0].header, blocks.map(b => [b.header, b.header] as [string, string]))
            };
        }
    },
    apply({ a, params }) {
        return Task.create('Parse mmCIF', async ctx => {
            const header = params.blockHeader || a.data.blocks[0].header;
            const block = a.data.blocks.find(b => b.header === header);
            if (!block) throw new Error(`Data block '${[header]}' not found.`);
            const models = await Model.create(Format.mmCIF(block)).runInContext(ctx);
            if (models.length === 0) throw new Error('No models found.');
            const label = models.length === 1 ? `${models[0].label}` : `${models[0].label} (${models.length} models)`;
            return new SO.Molecule.Models({ label }, models);
        });
    }
});

export { CreateStructureFromModel }
namespace CreateStructureFromModel { export interface Params { modelIndex: number, transform3d?: Mat4 } }
const CreateStructureFromModel = PluginStateTransform.Create<SO.Molecule.Models, SO.Molecule.Structure, CreateStructureFromModel.Params>({
    name: 'create-structure-from-model',
    display: {
        name: 'Structure from Model',
        description: 'Create a molecular structure from the specified model.'
    },
    from: [SO.Molecule.Models],
    to: [SO.Molecule.Structure],
    params: {
        default: () => ({ modelIndex: 0 }),
        controls: a => ({ modelIndex: PD.Range('Model Index', 'Model Index', 0, 0, Math.max(0, a.data.length - 1), 1) })
    },
    isApplicable: a => a.data.length > 0,
    apply({ a, params }) {
        if (params.modelIndex < 0 || params.modelIndex >= a.data.length) throw new Error(`Invalid modelIndex ${params.modelIndex}`);
        let s = Structure.ofModel(a.data[params.modelIndex]);
        if (params.transform3d) s = Structure.transform(s, params.transform3d);
        return new SO.Molecule.Structure({ label: `Model ${s.models[0].modelNum}`, description: s.elementCount === 1 ? '1 element' : `${s.elementCount} elements` }, s);
    }
});


export { CreateStructureAssembly }
namespace CreateStructureAssembly { export interface Params { /** if not specified, use the 1st */ id?: string } }
const CreateStructureAssembly = PluginStateTransform.Create<SO.Molecule.Structure, SO.Molecule.Structure, CreateStructureAssembly.Params>({
    name: 'create-structure-assembly',
    display: {
        name: 'Structure Assembly',
        description: 'Create a molecular structure assembly.'
    },
    from: [SO.Molecule.Structure],
    to: [SO.Molecule.Structure],
    params: {
        default: () => ({ id: void 0 }),
        controls(a) {
            const { model } = a.data;
            const ids = model.symmetry.assemblies.map(a => [a.id, a.id] as [string, string]);
            return { id: PD.Select('Asm Id', 'Assembly Id', ids.length ? ids[0][0] : '', ids) };
        }
    },
    isApplicable: a => a.data.models.length === 1 && a.data.model.symmetry.assemblies.length > 0,
    apply({ a, params }) {
        return Task.create('Build Assembly', async ctx => {
            let id = params.id;
            const model = a.data.model;
            if (!id && model.symmetry.assemblies.length) id = model.symmetry.assemblies[0].id;
            const asm = ModelSymmetry.findAssembly(a.data.model, id || '');
            if (!asm) throw new Error(`Assembly '${id}' not found`);

            const s = await StructureSymmetry.buildAssembly(a.data, id!).runInContext(ctx);
            return new SO.Molecule.Structure({ label: `Assembly ${id}`, description: s.elementCount === 1 ? '1 element' : `${s.elementCount} elements` }, s);
        })
    }
});

export { CreateStructureSelection }
namespace CreateStructureSelection { export interface Params { query: Expression, label?: string } }
const CreateStructureSelection = PluginStateTransform.Create<SO.Molecule.Structure, SO.Molecule.Structure, CreateStructureSelection.Params>({
    name: 'create-structure-selection',
    display: {
        name: 'Structure Selection',
        description: 'Create a molecular structure from the specified model.'
    },
    from: [SO.Molecule.Structure],
    to: [SO.Molecule.Structure],
    apply({ a, params }) {
        // TODO: use cache, add "update"
        const compiled = compile<StructureSelection>(params.query);
        const result = compiled(new QueryContext(a.data));
        const s = StructureSelection.unionStructure(result);
        return new SO.Molecule.Structure({ label: `${params.label || 'Selection'}`, description: s.elementCount === 1 ? '1 element' : `${s.elementCount} elements` }, s);
    }
});
