/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateTransform } from '../objects';
import { PluginStateObject as SO } from '../objects';
import { Task, RuntimeContext } from 'mol-task';
import { Model, Structure, ModelSymmetry, StructureSymmetry, QueryContext, StructureSelection as Sel, StructureQuery, Queries } from 'mol-model/structure';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import Expression from 'mol-script/language/expression';
import { compile } from 'mol-script/runtime/query/compiler';
import { MolScriptBuilder } from 'mol-script/language/builder';
import { StateObject } from 'mol-state';
import { PluginContext } from 'mol-plugin/context';
import { stringToWords } from 'mol-util/string';
import { volumeFromCcp4 } from 'mol-model-formats/volume/ccp4';
import { Vec3 } from 'mol-math/linear-algebra';
import CIF from 'mol-io/reader/cif';
import { volumeFromDsn6 } from 'mol-model-formats/volume/dsn6';
import { volumeFromDensityServerData } from 'mol-model-formats/volume/density-server';
import { trajectoryFromMmCIF } from 'mol-model-formats/structure/mmcif';
import { parsePDB } from 'mol-io/reader/pdb/parser';
import { trajectoryFromPDB } from 'mol-model-formats/structure/pdb';
import { Assembly } from 'mol-model/structure/model/properties/symmetry';

export { TrajectoryFromMmCif }
type TrajectoryFromMmCif = typeof TrajectoryFromMmCif
const TrajectoryFromMmCif = PluginStateTransform.BuiltIn({
    name: 'trajectory-from-mmcif',
    display: { name: 'Trajectory from mmCIF', description: 'Identify and create all separate models in the specified CIF data block' },
    from: SO.Format.Cif,
    to: SO.Molecule.Trajectory,
    params(a) {
        if (!a) {
            return {
                blockHeader: PD.makeOptional(PD.Text(void 0, { description: 'Header of the block to parse. If none is specifed, the 1st data block in the file is used.' }))
            };
        }
        const { blocks } = a.data;
        return {
            blockHeader: PD.makeOptional(PD.Select(blocks[0] && blocks[0].header, blocks.map(b => [b.header, b.header] as [string, string]), { description: 'Header of the block to parse' }))
        };
    }
})({
    isApplicable: a => a.data.blocks.length > 0,
    apply({ a, params }) {
        return Task.create('Parse mmCIF', async ctx => {
            const header = params.blockHeader || a.data.blocks[0].header;
            const block = a.data.blocks.find(b => b.header === header);
            if (!block) throw new Error(`Data block '${[header]}' not found.`);
            const models = await trajectoryFromMmCIF(block).runInContext(ctx);
            if (models.length === 0) throw new Error('No models found.');
            const props = { label: models[0].label, description: `${models.length} model${models.length === 1 ? '' : 's'}` };
            return new SO.Molecule.Trajectory(models, props);
        });
    }
});


export { TrajectoryFromPDB }
type TrajectoryFromPDB = typeof TrajectoryFromPDB
const TrajectoryFromPDB = PluginStateTransform.BuiltIn({
    name: 'trajectory-from-pdb',
    display: { name: 'Parse PDB', description: 'Parse PDB string and create trajectory.' },
    from: [SO.Data.String],
    to: SO.Molecule.Trajectory
})({
    apply({ a }) {
        return Task.create('Parse PDB', async ctx => {
            const parsed = await parsePDB(a.data).runInContext(ctx);
            if (parsed.isError) throw new Error(parsed.message);
            const models = await trajectoryFromPDB(parsed.result).runInContext(ctx);
            const props = { label: models[0].label, description: `${models.length} model${models.length === 1 ? '' : 's'}` };
            return new SO.Molecule.Trajectory(models, props);
        });
    }
});


export { ModelFromTrajectory }
const plus1 = (v: number) => v + 1, minus1 = (v: number) => v - 1;
type ModelFromTrajectory = typeof ModelFromTrajectory
const ModelFromTrajectory = PluginStateTransform.BuiltIn({
    name: 'model-from-trajectory',
    display: { name: 'Molecular Model', description: 'Create a molecular model from specified index in a trajectory.' },
    from: SO.Molecule.Trajectory,
    to: SO.Molecule.Model,
    params: a => {
        if (!a) {
            return { modelIndex: PD.Numeric(0, {}, { description: 'Zero-based index of the model' }) };
        }
        return { modelIndex: PD.Converted(plus1, minus1, PD.Numeric(1, { min: 1, max: a.data.length, step: 1 }, { description: 'Model Index' })) }
    }
})({
    isApplicable: a => a.data.length > 0,
    apply({ a, params }) {
        if (params.modelIndex < 0 || params.modelIndex >= a.data.length) throw new Error(`Invalid modelIndex ${params.modelIndex}`);
        const model = a.data[params.modelIndex];
        const props = a.data.length === 1
            ? { label: `${model.label}` }
            : { label: `${model.label}:${model.modelNum}`, description: `Model ${model.modelNum} of ${a.data.length}` };
        return new SO.Molecule.Model(model, props);
    }
});

export { StructureFromModel }
type StructureFromModel = typeof StructureFromModel
const StructureFromModel = PluginStateTransform.BuiltIn({
    name: 'structure-from-model',
    display: { name: 'Structure from Model', description: 'Create a molecular structure from the specified model.' },
    from: SO.Molecule.Model,
    to: SO.Molecule.Structure
})({
    apply({ a }) {
        let s = Structure.ofModel(a.data);
        const props = { label: a.data.label, description: s.elementCount === 1 ? '1 element' : `${s.elementCount} elements` };
        return new SO.Molecule.Structure(s, props);
    }
});

function structureDesc(s: Structure) {
    return s.elementCount === 1 ? '1 element' : `${s.elementCount} elements`;
}

export { StructureAssemblyFromModel }
type StructureAssemblyFromModel = typeof StructureAssemblyFromModel
const StructureAssemblyFromModel = PluginStateTransform.BuiltIn({
    name: 'structure-assembly-from-model',
    display: { name: 'Structure Assembly', description: 'Create a molecular structure assembly.' },
    from: SO.Molecule.Model,
    to: SO.Molecule.Structure,
    params(a) {
        if (!a) {
            return { id: PD.makeOptional(PD.Text('', { label: 'Assembly Id', description: 'Assembly Id. Value \'deposited\' can be used to specify deposited asymmetric unit.' })) };
        }
        const model = a.data;
        const ids = model.symmetry.assemblies.map(a => [a.id, `${a.id}: ${stringToWords(a.details)}`] as [string, string]);
        ids.push(['deposited', 'Deposited']);
        return { id: PD.makeOptional(PD.Select(ids[0][0], ids, { label: 'Asm Id', description: 'Assembly Id' })) };
    }
})({
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Build Assembly', async ctx => {
            const model = a.data;
            let id = params.id;
            let asm: Assembly | undefined = void 0;

            // if no id is specified, use the 1st assembly.
            if (!id && model.symmetry.assemblies.length !== 0) {
                id = model.symmetry.assemblies[0].id;
            }

            if (model.symmetry.assemblies.length === 0) {
                if (id !== 'deposited') {
                    plugin.log.warn(`Model '${a.label}' has no assembly, returning deposited structure.`);
                }
            } else {
                asm = ModelSymmetry.findAssembly(model, id || '');
                if (!asm) {
                    plugin.log.warn(`Model '${a.label}' has no assembly called '${id}', returning deposited structure.`);
                }
            }

            const base = Structure.ofModel(model);
            if (!asm) {
                const label = { label: a.data.label, description: structureDesc(base) };
                return new SO.Molecule.Structure(base, label);
            }

            id = asm.id;
            const s = await StructureSymmetry.buildAssembly(base, id!).runInContext(ctx);
            const props = { label: `Assembly ${id}`, description: structureDesc(s) };
            return new SO.Molecule.Structure(s, props);
        })
    }
});

export { StructureSymmetryFromModel }
type StructureSymmetryFromModel = typeof StructureSymmetryFromModel
const StructureSymmetryFromModel = PluginStateTransform.BuiltIn({
    name: 'structure-symmetry-from-model',
    display: { name: 'Structure Symmetry', description: 'Create a molecular structure symmetry.' },
    from: SO.Molecule.Model,
    to: SO.Molecule.Structure,
    params(a) {
        return {
            ijkMin: PD.Vec3(Vec3.create(-1, -1, -1), { label: 'Min IJK', fieldLabels: { x: 'I', y: 'J', z: 'K' } }),
            ijkMax: PD.Vec3(Vec3.create(1, 1, 1), { label: 'Max IJK', fieldLabels: { x: 'I', y: 'J', z: 'K' } })
        }
    }
})({
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Build Symmetry', async ctx => {
            const { ijkMin, ijkMax } = params
            const model = a.data;
            const base = Structure.ofModel(model);
            const s = await StructureSymmetry.buildSymmetryRange(base, ijkMin, ijkMax).runInContext(ctx);
            const props = { label: `Symmetry [${ijkMin}] to [${ijkMax}]`, description: structureDesc(s) };
            return new SO.Molecule.Structure(s, props);
        })
    }
});

export { StructureSelection }
type StructureSelection = typeof StructureSelection
const StructureSelection = PluginStateTransform.BuiltIn({
    name: 'structure-selection',
    display: { name: 'Structure Selection', description: 'Create a molecular structure from the specified query expression.' },
    from: SO.Molecule.Structure,
    to: SO.Molecule.Structure,
    params: {
        query: PD.Value<Expression>(MolScriptBuilder.struct.generator.all, { isHidden: true }),
        label: PD.makeOptional(PD.Text('', { isHidden: true }))
    }
})({
    apply({ a, params }) {
        // TODO: use cache, add "update"
        const compiled = compile<Sel>(params.query);
        const result = compiled(new QueryContext(a.data));
        const s = Sel.unionStructure(result);
        const props = { label: `${params.label || 'Selection'}`, description: structureDesc(s) };
        return new SO.Molecule.Structure(s, props);
    }
});

export { StructureComplexElement }
namespace StructureComplexElement {
    export type Types = 'atomic-sequence' | 'water' | 'atomic-het' | 'spheres'
}
const StructureComplexElementTypes: [StructureComplexElement.Types, StructureComplexElement.Types][] = ['atomic-sequence', 'water', 'atomic-het', 'spheres'].map(t => [t, t] as any);
type StructureComplexElement = typeof StructureComplexElement
const StructureComplexElement = PluginStateTransform.BuiltIn({
    name: 'structure-complex-element',
    display: { name: 'Complex Element', description: 'Create a molecular structure from the specified model.' },
    from: SO.Molecule.Structure,
    to: SO.Molecule.Structure,
    params: { type: PD.Select<StructureComplexElement.Types>('atomic-sequence', StructureComplexElementTypes, { isHidden: true }) }
})({
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

export { CustomModelProperties }
type CustomModelProperties = typeof CustomModelProperties
const CustomModelProperties = PluginStateTransform.BuiltIn({
    name: 'custom-model-properties',
    display: { name: 'Custom Model Properties' },
    from: SO.Molecule.Model,
    to: SO.Molecule.Model,
    params: (a, ctx: PluginContext) => {
        if (!a) return { properties: PD.MultiSelect([], [], { description: 'A list of property descriptor ids.' }) };
        return { properties: ctx.customModelProperties.getSelect(a.data) };
    }
})({
    apply({ a, params }, ctx: PluginContext) {
        return Task.create('Custom Props', async taskCtx => {
            await attachProps(a.data, ctx, taskCtx, params.properties);
            return new SO.Molecule.Model(a.data, { label: 'Props', description: `${params.properties.length} Selected` });
        });
    }
});
async function attachProps(model: Model, ctx: PluginContext, taskCtx: RuntimeContext, names: string[]) {
    for (const name of names) {
        const p = ctx.customModelProperties.get(name);
        await p.attach(model).runInContext(taskCtx);
    }
}

//

export { VolumeFromCcp4 }
type VolumeFromCcp4 = typeof VolumeFromCcp4
const VolumeFromCcp4 = PluginStateTransform.BuiltIn({
    name: 'volume-from-ccp4',
    display: { name: 'Volume from CCP4/MRC/MAP', description: 'Create Volume from CCP4/MRC/MAP data' },
    from: SO.Format.Ccp4,
    to: SO.Volume.Data,
    params(a) {
        return {
            voxelSize: PD.Vec3(Vec3.create(1, 1, 1))
        };
    }
})({
    apply({ a, params }) {
        return Task.create('Create volume from CCP4/MRC/MAP', async ctx => {
            const volume = await volumeFromCcp4(a.data, params).runInContext(ctx)
            const props = { label: 'Volume' };
            return new SO.Volume.Data(volume, props);
        });
    }
});

export { VolumeFromDsn6 }
type VolumeFromDsn6 = typeof VolumeFromDsn6
const VolumeFromDsn6 = PluginStateTransform.BuiltIn({
    name: 'volume-from-dsn6',
    display: { name: 'Volume from DSN6/BRIX', description: 'Create Volume from DSN6/BRIX data' },
    from: SO.Format.Dsn6,
    to: SO.Volume.Data,
    params(a) {
        return {
            voxelSize: PD.Vec3(Vec3.create(1, 1, 1))
        };
    }
})({
    apply({ a, params }) {
        return Task.create('Create volume from DSN6/BRIX', async ctx => {
            const volume = await volumeFromDsn6(a.data, params).runInContext(ctx)
            const props = { label: 'Volume' };
            return new SO.Volume.Data(volume, props);
        });
    }
});

export { VolumeFromDensityServerCif }
type VolumeFromDensityServerCif = typeof VolumeFromDensityServerCif
const VolumeFromDensityServerCif = PluginStateTransform.BuiltIn({
    name: 'volume-from-density-server-cif',
    display: { name: 'Volume from density-server CIF', description: 'Identify and create all separate models in the specified CIF data block' },
    from: SO.Format.Cif,
    to: SO.Volume.Data,
    params(a) {
        if (!a) {
            return {
                blockHeader: PD.makeOptional(PD.Text(void 0, { description: 'Header of the block to parse. If none is specifed, the 1st data block in the file is used.' }))
            };
        }
        const blocks = a.data.blocks.slice(1); // zero block contains query meta-data
        return {
            blockHeader: PD.makeOptional(PD.Select(blocks[0] && blocks[0].header, blocks.map(b => [b.header, b.header] as [string, string]), { description: 'Header of the block to parse' }))
        };
    }
})({
    isApplicable: a => a.data.blocks.length > 0,
    apply({ a, params }) {
        return Task.create('Parse density-server CIF', async ctx => {
            const header = params.blockHeader || a.data.blocks[1].header; // zero block contains query meta-data
            const block = a.data.blocks.find(b => b.header === header);
            if (!block) throw new Error(`Data block '${[header]}' not found.`);
            const densityServerCif = CIF.schema.densityServer(block)
            const volume = await volumeFromDensityServerData(densityServerCif).runInContext(ctx)
            const props = { label: densityServerCif.volume_data_3d_info.name.value(0), description: `${densityServerCif.volume_data_3d_info.name.value(0)}` };
            return new SO.Volume.Data(volume, props);
        });
    }
});