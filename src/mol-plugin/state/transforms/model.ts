/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { parsePDB } from 'mol-io/reader/pdb/parser';
import { Vec3, Mat4, Quat } from 'mol-math/linear-algebra';
import { trajectoryFromMmCIF } from 'mol-model-formats/structure/mmcif';
import { trajectoryFromPDB } from 'mol-model-formats/structure/pdb';
import { Model, ModelSymmetry, Queries, QueryContext, Structure, StructureQuery, StructureSelection as Sel, StructureSymmetry, QueryFn } from 'mol-model/structure';
import { Assembly } from 'mol-model/structure/model/properties/symmetry';
import { PluginContext } from 'mol-plugin/context';
import { MolScriptBuilder } from 'mol-script/language/builder';
import Expression from 'mol-script/language/expression';
import { compile } from 'mol-script/runtime/query/compiler';
import { StateObject, StateTransformer } from 'mol-state';
import { RuntimeContext, Task } from 'mol-task';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { stringToWords } from 'mol-util/string';
import { PluginStateObject as SO, PluginStateTransform } from '../objects';
import { trajectoryFromGRO } from 'mol-model-formats/structure/gro';
import { parseGRO } from 'mol-io/reader/gro/parser';
import { parseMolScript } from 'mol-script/language/parser';
import { transpileMolScript } from 'mol-script/script/mol-script/symbols';
import { shapeFromPly } from 'mol-model-formats/shape/ply';
import { SymmetryOperator } from 'mol-math/geometry';
import { ensureSecondaryStructure } from './helpers';

export { TrajectoryFromBlob };
export { TrajectoryFromMmCif };
export { TrajectoryFromPDB };
export { TrajectoryFromGRO };
export { ModelFromTrajectory };
export { StructureFromModel };
export { StructureAssemblyFromModel };
export { StructureSymmetryFromModel };
export { TransformStructureConformation }
export { StructureSelection };
export { UserStructureSelection };
export { StructureComplexElement };
export { CustomModelProperties };
export { CustomStructureProperties };

type TrajectoryFromBlob = typeof TrajectoryFromBlob
const TrajectoryFromBlob = PluginStateTransform.BuiltIn({
    name: 'trajectory-from-blob',
    display: { name: 'Parse Blob', description: 'Parse format blob into a single trajectory.' },
    from: SO.Format.Blob,
    to: SO.Molecule.Trajectory
})({
    apply({ a }) {
        return Task.create('Parse Format Blob', async ctx => {
            const models: Model[] = [];
            for (const e of a.data) {
                if (e.kind !== 'cif') continue;
                const block = e.data.blocks[0];
                const xs = await trajectoryFromMmCIF(block).runInContext(ctx);
                if (xs.length === 0) throw new Error('No models found.');
                for (const x of xs) models.push(x);
            }

            const props = { label: `Trajectory`, description: `${models.length} model${models.length === 1 ? '' : 's'}` };
            return new SO.Molecule.Trajectory(models, props);
        });
    }
});

type TrajectoryFromMmCif = typeof TrajectoryFromMmCif
const TrajectoryFromMmCif = PluginStateTransform.BuiltIn({
    name: 'trajectory-from-mmcif',
    display: { name: 'Trajectory from mmCIF', description: 'Identify and create all separate models in the specified CIF data block' },
    from: SO.Format.Cif,
    to: SO.Molecule.Trajectory,
    params(a) {
        if (!a) {
            return {
                blockHeader: PD.Optional(PD.Text(void 0, { description: 'Header of the block to parse. If none is specifed, the 1st data block in the file is used.' }))
            };
        }
        const { blocks } = a.data;
        return {
            blockHeader: PD.Optional(PD.Select(blocks[0] && blocks[0].header, blocks.map(b => [b.header, b.header] as [string, string]), { description: 'Header of the block to parse' }))
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

type TrajectoryFromPDB = typeof TrajectoryFromPDB
const TrajectoryFromPDB = PluginStateTransform.BuiltIn({
    name: 'trajectory-from-pdb',
    display: { name: 'Parse PDB', description: 'Parse PDB string and create trajectory.' },
    from: [SO.Data.String],
    to: SO.Molecule.Trajectory
})({
    apply({ a }) {
        return Task.create('Parse PDB', async ctx => {
            const parsed = await parsePDB(a.data, a.label).runInContext(ctx);
            if (parsed.isError) throw new Error(parsed.message);
            const models = await trajectoryFromPDB(parsed.result).runInContext(ctx);
            const props = { label: models[0].label, description: `${models.length} model${models.length === 1 ? '' : 's'}` };
            return new SO.Molecule.Trajectory(models, props);
        });
    }
});

type TrajectoryFromGRO = typeof TrajectoryFromGRO
const TrajectoryFromGRO = PluginStateTransform.BuiltIn({
    name: 'trajectory-from-gro',
    display: { name: 'Parse GRO', description: 'Parse GRO string and create trajectory.' },
    from: [SO.Data.String],
    to: SO.Molecule.Trajectory
})({
    apply({ a }) {
        return Task.create('Parse GRO', async ctx => {
            const parsed = await parseGRO(a.data).runInContext(ctx);
            if (parsed.isError) throw new Error(parsed.message);
            const models = await trajectoryFromGRO(parsed.result).runInContext(ctx);
            const props = { label: models[0].label, description: `${models.length} model${models.length === 1 ? '' : 's'}` };
            return new SO.Molecule.Trajectory(models, props);
        });
    }
});

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
        const label = a.data.length === 1 ? model.entry : `${model.entry}: ${model.modelNum}`
        const description = a.data.length === 1 ? undefined : `Model ${params.modelIndex + 1} of ${a.data.length}`
        return new SO.Molecule.Model(model, { label, description });
    }
});

type StructureFromModel = typeof StructureFromModel
const StructureFromModel = PluginStateTransform.BuiltIn({
    name: 'structure-from-model',
    display: { name: 'Structure from Model', description: 'Create a molecular structure from the specified model.' },
    from: SO.Molecule.Model,
    to: SO.Molecule.Structure
})({
    apply({ a }) {
        return Task.create('Build Structure', async ctx => {
            const s = Structure.ofModel(a.data);
            await ensureSecondaryStructure(s)
            const props = { label: a.data.label, description: s.elementCount === 1 ? '1 element' : `${s.elementCount} elements` };
            return new SO.Molecule.Structure(s, props);
        })
    }
});

function structureDesc(s: Structure) {
    return s.elementCount === 1 ? '1 element' : `${s.elementCount} elements`;
}

type StructureAssemblyFromModel = typeof StructureAssemblyFromModel
const StructureAssemblyFromModel = PluginStateTransform.BuiltIn({
    name: 'structure-assembly-from-model',
    display: { name: 'Structure Assembly', description: 'Create a molecular structure assembly.' },
    from: SO.Molecule.Model,
    to: SO.Molecule.Structure,
    params(a) {
        if (!a) {
            return { id: PD.Optional(PD.Text('', { label: 'Assembly Id', description: 'Assembly Id. Value \'deposited\' can be used to specify deposited asymmetric unit.' })) };
        }
        const model = a.data;
        const ids = model.symmetry.assemblies.map(a => [a.id, `${a.id}: ${stringToWords(a.details)}`] as [string, string]);
        ids.push(['deposited', 'Deposited']);
        return {
            id: PD.Optional(PD.Select(ids[0][0], ids, { label: 'Asm Id', description: 'Assembly Id' })) };
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
                    plugin.log.warn(`Model '${a.data.label}' has no assembly, returning deposited structure.`);
                }
            } else {
                asm = ModelSymmetry.findAssembly(model, id || '');
                if (!asm) {
                    plugin.log.warn(`Model '${a.data.label}' has no assembly called '${id}', returning deposited structure.`);
                }
            }

            const base = Structure.ofModel(model);
            if (!asm) {
                await ensureSecondaryStructure(base)
                const label = { label: 'Deposited', description: structureDesc(base) };
                return new SO.Molecule.Structure(base, label);
            }

            id = asm.id;
            const s = await StructureSymmetry.buildAssembly(base, id!).runInContext(ctx);
            await ensureSecondaryStructure(s)
            const props = { label: `Assembly ${id}`, description: structureDesc(s) };
            return new SO.Molecule.Structure(s, props);
        })
    }
});

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
            await ensureSecondaryStructure(s)
            const props = { label: `Symmetry [${ijkMin}] to [${ijkMax}]`, description: structureDesc(s) };
            return new SO.Molecule.Structure(s, props);
        })
    }
});

const _translation = Vec3.zero(), _m = Mat4.zero(), _n = Mat4.zero();
type TransformStructureConformation = typeof TransformStructureConformation
const TransformStructureConformation = PluginStateTransform.BuiltIn({
    name: 'transform-structure-conformation',
    display: { name: 'Transform Conformation' },
    from: SO.Molecule.Structure,
    to: SO.Molecule.Structure,
    params: {
        axis: PD.Vec3(Vec3.create(1, 0, 0)),
        angle: PD.Numeric(0, { min: -180, max: 180, step: 0.1 }),
        translation: PD.Vec3(Vec3.create(0, 0, 0)),
    }
})({
    canAutoUpdate() {
        return true;
    },
    apply({ a, params }) {
        // TODO: optimze

        const center = a.data.boundary.sphere.center;
        Mat4.fromTranslation(_m, Vec3.negate(_translation, center));
        Mat4.fromTranslation(_n, Vec3.add(_translation, center, params.translation));
        const rot = Mat4.fromRotation(Mat4.zero(), Math.PI / 180 * params.angle, Vec3.normalize(Vec3.zero(), params.axis));

        const m = Mat4.zero();
        Mat4.mul3(m, _n, rot, _m);

        const s = Structure.transform(a.data, m);
        const props = { label: `${a.label}`, description: `Transformed` };
        return new SO.Molecule.Structure(s, props);
    },
    interpolate(src, tar, t) {
        // TODO: optimize
        const u = Mat4.fromRotation(Mat4.zero(), Math.PI / 180 * src.angle, Vec3.normalize(Vec3.zero(), src.axis));
        Mat4.setTranslation(u, src.translation);
        const v = Mat4.fromRotation(Mat4.zero(), Math.PI / 180 * tar.angle, Vec3.normalize(Vec3.zero(), tar.axis));
        Mat4.setTranslation(v, tar.translation);
        const m = SymmetryOperator.slerp(Mat4.zero(), u, v, t);
        const rot = Mat4.getRotation(Quat.zero(), m);
        const axis = Vec3.zero();
        const angle = Quat.getAxisAngle(axis, rot);
        const translation = Mat4.getTranslation(Vec3.zero(), m);
        return { axis, angle, translation };
    }
});

type StructureSelection = typeof StructureSelection
const StructureSelection = PluginStateTransform.BuiltIn({
    name: 'structure-selection',
    display: { name: 'Structure Selection', description: 'Create a molecular structure from the specified query expression.' },
    from: SO.Molecule.Structure,
    to: SO.Molecule.Structure,
    params: {
        query: PD.Value<Expression>(MolScriptBuilder.struct.generator.all, { isHidden: true }),
        label: PD.Optional(PD.Text('', { isHidden: true }))
    }
})({
    apply({ a, params, cache }) {
        const compiled = compile<Sel>(params.query);
        (cache as { compiled: QueryFn<Sel> }).compiled = compiled;
        (cache as { source: Structure }).source = a.data;

        const result = compiled(new QueryContext(a.data));
        const s = Sel.unionStructure(result);
        if (s.elementCount === 0) return StateObject.Null;
        const props = { label: `${params.label || 'Selection'}`, description: structureDesc(s) };
        return new SO.Molecule.Structure(s, props);
    },
    update: ({ a, b, oldParams, newParams, cache }) => {
        if (oldParams.query !== newParams.query) return StateTransformer.UpdateResult.Recreate;

        if ((cache as { source: Structure }).source === a.data) {
            return StateTransformer.UpdateResult.Unchanged;
        }
        (cache as { source: Structure }).source = a.data;

        if (updateStructureFromQuery((cache as { compiled: QueryFn<Sel> }).compiled, a.data, b, newParams.label)) {
            return StateTransformer.UpdateResult.Updated;
        }
        return StateTransformer.UpdateResult.Null;
    }
});

type UserStructureSelection = typeof UserStructureSelection
const UserStructureSelection = PluginStateTransform.BuiltIn({
    name: 'user-structure-selection',
    display: { name: 'Structure Selection', description: 'Create a molecular structure from the specified query expression.' },
    from: SO.Molecule.Structure,
    to: SO.Molecule.Structure,
    params: {
        query: PD.ScriptExpression({ language: 'mol-script', expression: '(sel.atom.atom-groups :residue-test (= atom.resname ALA))' }),
        label: PD.Optional(PD.Text(''))
    }
})({
    apply({ a, params, cache }) {
        const parsed = parseMolScript(params.query.expression);
        if (parsed.length === 0) throw new Error('No query');
        const query = transpileMolScript(parsed[0]);
        const compiled = compile<Sel>(query);
        (cache as { compiled: QueryFn<Sel> }).compiled = compiled;
        (cache as { source: Structure }).source = a.data;
        const result = compiled(new QueryContext(a.data));
        const s = Sel.unionStructure(result);
        const props = { label: `${params.label || 'Selection'}`, description: structureDesc(s) };
        return new SO.Molecule.Structure(s, props);
    },
    update: ({ a, b, oldParams, newParams, cache }) => {
        if (oldParams.query.language !== newParams.query.language || oldParams.query.expression !== newParams.query.expression) {
            return StateTransformer.UpdateResult.Recreate;
        }

        if ((cache as { source: Structure }).source === a.data) {
            return StateTransformer.UpdateResult.Unchanged;
        }
        (cache as { source: Structure }).source = a.data;

        updateStructureFromQuery((cache as { compiled: QueryFn<Sel> }).compiled, a.data, b, newParams.label);
        return StateTransformer.UpdateResult.Updated;
    }
});

function updateStructureFromQuery(query: QueryFn<Sel>, src: Structure, obj: SO.Molecule.Structure, label?: string) {
    const result = query(new QueryContext(src));
    const s = Sel.unionStructure(result);
    if (s.elementCount === 0) {
        return false;
    }

    obj.label = `${label || 'Selection'}`;
    obj.description = structureDesc(s);
    obj.data = s;
    return true;
}

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
            await attachModelProps(a.data, ctx, taskCtx, params.properties);
            return new SO.Molecule.Model(a.data, { label: 'Model Props', description: `${params.properties.length} Selected` });
        });
    }
});
async function attachModelProps(model: Model, ctx: PluginContext, taskCtx: RuntimeContext, names: string[]) {
    for (const name of names) {
        const p = ctx.customModelProperties.get(name);
        await p.attach(model).runInContext(taskCtx);
    }
}

type CustomStructureProperties = typeof CustomStructureProperties
const CustomStructureProperties = PluginStateTransform.BuiltIn({
    name: 'custom-structure-properties',
    display: { name: 'Custom Structure Properties' },
    from: SO.Molecule.Structure,
    to: SO.Molecule.Structure,
    params: (a, ctx: PluginContext) => {
        if (!a) return { properties: PD.MultiSelect([], [], { description: 'A list of property descriptor ids.' }) };
        return { properties: ctx.customStructureProperties.getSelect(a.data) };
    }
})({
    apply({ a, params }, ctx: PluginContext) {
        return Task.create('Custom Props', async taskCtx => {
            await attachStructureProps(a.data, ctx, taskCtx, params.properties);
            return new SO.Molecule.Structure(a.data, { label: 'Structure Props', description: `${params.properties.length} Selected` });
        });
    }
});
async function attachStructureProps(structure: Structure, ctx: PluginContext, taskCtx: RuntimeContext, names: string[]) {
    for (const name of names) {
        const p = ctx.customStructureProperties.get(name);
        await p.attach(structure).runInContext(taskCtx);
    }
}

export { ShapeFromPly }
type ShapeFromPly = typeof ShapeFromPly
const ShapeFromPly = PluginStateTransform.BuiltIn({
    name: 'shape-from-ply',
    display: { name: 'Shape from PLY', description: 'Create Shape from PLY data' },
    from: SO.Format.Ply,
    to: SO.Shape.Provider,
    params(a) {
        return { };
    }
})({
    apply({ a, params }) {
        return Task.create('Create shape from PLY', async ctx => {
            const shape = await shapeFromPly(a.data, params).runInContext(ctx)
            const props = { label: 'Shape' };
            return new SO.Shape.Provider(shape, props);
        });
    }
});