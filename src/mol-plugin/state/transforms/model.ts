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
type TrajectoryFromMmCif = typeof TrajectoryFromMmCif
const TrajectoryFromMmCif = PluginStateTransform.BuiltIn({
    name: 'trajectory-from-mmcif',
    from: SO.Format.Cif,
    to: SO.Molecule.Trajectory,
    params(a) {
        const { blocks } = a.data;
        return {
            blockHeader: PD.makeOptional(PD.Select(blocks[0] && blocks[0].header, blocks.map(b => [b.header, b.header] as [string, string]), { description: 'Header of the block to parse' }))
        };
    },
})({
    display: {
        name: 'Models from mmCIF',
        description: 'Identify and create all separate models in the specified CIF data block'
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
type ModelFromTrajectory = typeof ModelFromTrajectory
const ModelFromTrajectory = PluginStateTransform.BuiltIn({
    name: 'model-from-trajectory',
    from: SO.Molecule.Trajectory,
    to: SO.Molecule.Model,
    params: a => ({ modelIndex: PD.Converted(plus1, minus1, PD.Numeric(1, { min: 1, max: a.data.length, step: 1 }, { description: 'Model Index' })) })
})({
    display: {
        name: 'Model from Trajectory',
        description: 'Create a molecular structure from the specified model.'
    },
    isApplicable: a => a.data.length > 0,
    apply({ a, params }) {
        if (params.modelIndex < 0 || params.modelIndex >= a.data.length) throw new Error(`Invalid modelIndex ${params.modelIndex}`);
        const model = a.data[params.modelIndex];
        const label = { label: `Model ${model.modelNum}` };
        return new SO.Molecule.Model(model, label);
    }
});

export { StructureFromModel }
type StructureFromModel = typeof StructureFromModel
const StructureFromModel = PluginStateTransform.BuiltIn({
    name: 'structure-from-model',
    from: SO.Molecule.Model,
    to: SO.Molecule.Structure,
})({
    display: {
        name: 'Structure from Model',
        description: 'Create a molecular structure from the specified model.'
    },
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
type StructureAssemblyFromModel = typeof StructureAssemblyFromModel
const StructureAssemblyFromModel = PluginStateTransform.BuiltIn({
    name: 'structure-assembly-from-model',
    from: SO.Molecule.Model,
    to: SO.Molecule.Structure,
    params(a) {
        const model = a.data;
        const ids = model.symmetry.assemblies.map(a => [a.id, a.id] as [string, string]);
        return { id: PD.makeOptional(PD.Select(ids.length ? ids[0][0] : '', ids, { label: 'Asm Id', description: 'Assembly Id' })) };
    }
})({
    display: {
        name: 'Structure Assembly',
        description: 'Create a molecular structure assembly.'
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
type StructureSelection = typeof StructureSelection
const StructureSelection = PluginStateTransform.BuiltIn({
    name: 'structure-selection',
    from: SO.Molecule.Structure,
    to: SO.Molecule.Structure,
    params: () => ({
        query: PD.Value<Expression>(MolScriptBuilder.struct.generator.all, { isHidden: true }),
        label: PD.makeOptional(PD.Text('', { isHidden: true }))
    })
})({
    display: {
        name: 'Structure Selection',
        description: 'Create a molecular structure from the specified model.'
    },
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
namespace StructureComplexElement { export type Types = 'atomic-sequence' | 'water' | 'atomic-het' | 'spheres' }
type StructureComplexElement = typeof StructureComplexElement
const StructureComplexElement = PluginStateTransform.BuiltIn({
    name: 'structure-complex-element',
    from: SO.Molecule.Structure,
    to: SO.Molecule.Structure,
    params: () => ({ type: PD.Text<StructureComplexElement.Types>('atomic-sequence', { isHidden: true }) }),
})({
    display: {
        name: 'Complex Element',
        description: 'Create a molecular structure from the specified model.'
    },
    apply({ a, params }) {
        // TODO: update function.

        let query: StructureQuery, label: string;
        switch (params.type) {
            case 'atomic-sequence': query = Queries.internal.atomicSequence(); label = 'Sequence'; break;
            case 'water': query = Queries.internal.water(); label = 'Water'; break;
            case 'atomic-het': query = Queries.internal.atomicHet(); label = 'HET Groups/Ligands'; break;
            case 'spheres': query = Queries.internal.spheres(); label = 'Coarse Spheres'; break;
            default: throw new Error(`${params.type} is a not valid complex element.`);
        }

        const result = query(new QueryContext(a.data));
        const s = Sel.unionStructure(result);

        if (s.elementCount === 0) return StateObject.Null;
        return new SO.Molecule.Structure(s, { label, description: structureDesc(s) });
    }
});

