/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateTransform } from '../objects';
import { PluginStateObject as SO } from '../objects';
import { Task } from 'mol-task';
import { Model, Format, Structure, ModelSymmetry, StructureSymmetry, QueryContext, StructureSelection as Sel, StructureQuery, Queries } from 'mol-model/structure';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import Expression from 'mol-script/language/expression';
import { compile } from 'mol-script/runtime/query/compiler';
import { MolScriptBuilder } from 'mol-script/language/builder';
import { StateObject } from 'mol-state';
import { PluginContext } from 'mol-plugin/context';

export { TrajectoryFromMmCif }
namespace TrajectoryFromMmCif { export interface Params { blockHeader?: string } }
const TrajectoryFromMmCif = PluginStateTransform.Create<SO.Data.Cif, SO.Molecule.Trajectory, TrajectoryFromMmCif.Params>({
    name: 'trajectory-from-mmcif',
    display: {
        name: 'Models from mmCIF',
        description: 'Identify and create all separate models in the specified CIF data block'
    },
    from: [SO.Data.Cif],
    to: [SO.Molecule.Trajectory],
    params(a) {
        const { blocks } = a.data;
        if (blocks.length === 0) return { };
        return {
            blockHeader: PD.Select(blocks[0].header, blocks.map(b => [b.header, b.header] as [string, string]), { description: 'Header of the block to parse' })
        };
    },
    isApplicable: a => a.data.blocks.length > 0,
    apply({ a, params }) {
        return Task.create('Parse mmCIF', async ctx => {
            const header = params.blockHeader || a.data.blocks[0].header;
            const block = a.data.blocks.find(b => b.header === header);
            if (!block) throw new Error(`Data block '${[header]}' not found.`);
            const models = await Model.create(Format.mmCIF(block)).runInContext(ctx);
            if (models.length === 0) throw new Error('No models found.');
            const label = { label: models[0].label, description: `${models.length} model${models.length === 1 ? '' : 's'}` };
            return new SO.Molecule.Trajectory(models, label);
        });
    }
});

export { ModelFromTrajectory }
const plus1 = (v: number) => v + 1, minus1 = (v: number) => v - 1;
namespace ModelFromTrajectory { export interface Params { modelIndex: number } }
const ModelFromTrajectory = PluginStateTransform.Create<SO.Molecule.Trajectory, SO.Molecule.Model, ModelFromTrajectory.Params>({
    name: 'model-from-trajectory',
    display: {
        name: 'Model from Trajectory',
        description: 'Create a molecular structure from the specified model.'
    },
    from: [SO.Molecule.Trajectory],
    to: [SO.Molecule.Model],
    params: a => ({ modelIndex: PD.Converted(plus1, minus1, PD.Numeric(1, { min: 1, max: a.data.length, step: 1 }, { description: 'Model Index' })) }),
    isApplicable: a => a.data.length > 0,
    apply({ a, params }) {
        if (params.modelIndex < 0 || params.modelIndex >= a.data.length) throw new Error(`Invalid modelIndex ${params.modelIndex}`);
        const model = a.data[params.modelIndex];
        const label = { label: `Model ${model.modelNum}` };
        return new SO.Molecule.Model(model, label);
    }
});

export { StructureFromModel }
namespace StructureFromModel { export interface Params { } }
const StructureFromModel = PluginStateTransform.Create<SO.Molecule.Model, SO.Molecule.Structure, StructureFromModel.Params>({
    name: 'structure-from-model',
    display: {
        name: 'Structure from Model',
        description: 'Create a molecular structure from the specified model.'
    },
    from: [SO.Molecule.Model],
    to: [SO.Molecule.Structure],
    apply({ a }) {
        let s = Structure.ofModel(a.data);
        const label = { label: a.data.label, description: s.elementCount === 1 ? '1 element' : `${s.elementCount} elements` };
        return new SO.Molecule.Structure(s, label);
    }
});

function structureDesc(s: Structure) {
    return s.elementCount === 1 ? '1 element' : `${s.elementCount} elements`;
}

export { StructureAssemblyFromModel }
namespace StructureAssemblyFromModel { export interface Params { /** if not specified, use the 1st */ id?: string } }
const StructureAssemblyFromModel = PluginStateTransform.Create<SO.Molecule.Model, SO.Molecule.Structure, StructureAssemblyFromModel.Params>({
    name: 'structure-assembly-from-model',
    display: {
        name: 'Structure Assembly',
        description: 'Create a molecular structure assembly.'
    },
    from: [SO.Molecule.Model],
    to: [SO.Molecule.Structure],
    params(a) {
        const model = a.data;
        const ids = model.symmetry.assemblies.map(a => [a.id, a.id] as [string, string]);
        return { id: PD.Select(ids.length ? ids[0][0] : '', ids, { label: 'Asm Id', description: 'Assembly Id' }) };
    },
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Build Assembly', async ctx => {
            let id = (params.id || '').trim();
            const model = a.data;
            if (!id && model.symmetry.assemblies.length) id = model.symmetry.assemblies[0].id;
            const asm = ModelSymmetry.findAssembly(model, id);
            if (id && !asm) throw new Error(`Assembly '${id}' not found`);

            const base = Structure.ofModel(model);
            if (!asm) {
                plugin.log.warn(`Model '${a.label}' has no assembly, returning default structure.`);
                const label = { label: a.data.label, description: structureDesc(base) };
                return new SO.Molecule.Structure(base, label);
            }

            const s = await StructureSymmetry.buildAssembly(base, id!).runInContext(ctx);
            const label = { label: `Assembly ${id}`, description: structureDesc(s) };
            return new SO.Molecule.Structure(s, label);
        })
    }
});

export { StructureSelection }
namespace StructureSelection { export interface Params { query: Expression, label?: string } }
const StructureSelection = PluginStateTransform.Create<SO.Molecule.Structure, SO.Molecule.Structure, StructureSelection.Params>({
    name: 'structure-selection',
    display: {
        name: 'Structure Selection',
        description: 'Create a molecular structure from the specified model.'
    },
    from: [SO.Molecule.Structure],
    to: [SO.Molecule.Structure],
    params: () => ({
        query: PD.Value<Expression>(MolScriptBuilder.struct.generator.all, { isHidden: true }),
        label: PD.Text('', { isOptional: true, isHidden: true })
    }),
    apply({ a, params }) {
        // TODO: use cache, add "update"
        const compiled = compile<Sel>(params.query);
        const result = compiled(new QueryContext(a.data));
        const s = Sel.unionStructure(result);
        const label = { label: `${params.label || 'Selection'}`, description: structureDesc(s) };
        return new SO.Molecule.Structure(s, label);
    }
});

export { StructureComplexElement }
namespace StructureComplexElement { export interface Params { type: 'atomic-sequence' | 'water' | 'atomic-het' | 'spheres' } }
const StructureComplexElement = PluginStateTransform.Create<SO.Molecule.Structure, SO.Molecule.Structure, StructureComplexElement.Params>({
    name: 'structure-complex-element',
    display: {
        name: 'Complex Element',
        description: 'Create a molecular structure from the specified model.'
    },
    from: [SO.Molecule.Structure],
    to: [SO.Molecule.Structure],
    params: () => ({ type: PD.Text('sequence', { isHidden: true }) }),
    apply({ a, params }) {
        // TODO: update function.

        let query: StructureQuery, label: string;
        switch (params.type) {
            case 'atomic-sequence': query = Queries.internal.atomicSequence(); label = 'Sequence'; break;
            case 'water': query = Queries.internal.water(); label = 'Water'; break;
            case 'atomic-het': query = Queries.internal.atomicHet(); label = 'HET Groups/Ligands'; break;
            case 'spheres': query = Queries.internal.spheres(); label = 'Coarse Spheres'; break;
            default: throw new Error(`${params.type} is a valid complex element.`);
        }

        const result = query(new QueryContext(a.data));
        const s = Sel.unionStructure(result);

        if (s.elementCount === 0) return StateObject.Null;
        return new SO.Molecule.Structure(s, { label, description: structureDesc(s) });
    }
});

