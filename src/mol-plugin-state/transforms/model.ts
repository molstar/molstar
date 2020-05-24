/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { parse3DG } from '../../mol-io/reader/3dg/parser';
import { parseDcd } from '../../mol-io/reader/dcd/parser';
import { parseGRO } from '../../mol-io/reader/gro/parser';
import { parsePDB } from '../../mol-io/reader/pdb/parser';
import { Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { shapeFromPly } from '../../mol-model-formats/shape/ply';
import { trajectoryFrom3DG } from '../../mol-model-formats/structure/3dg';
import { coordinatesFromDcd } from '../../mol-model-formats/structure/dcd';
import { trajectoryFromGRO } from '../../mol-model-formats/structure/gro';
import { trajectoryFromMmCIF } from '../../mol-model-formats/structure/mmcif';
import { trajectoryFromPDB } from '../../mol-model-formats/structure/pdb';
import { topologyFromPsf } from '../../mol-model-formats/structure/psf';
import { Coordinates, Model, Queries, QueryContext, Structure, StructureElement, StructureQuery, StructureSelection as Sel, Topology } from '../../mol-model/structure';
import { PluginContext } from '../../mol-plugin/context';
import { MolScriptBuilder } from '../../mol-script/language/builder';
import Expression from '../../mol-script/language/expression';
import { Script } from '../../mol-script/script';
import { StateObject, StateTransformer } from '../../mol-state';
import { RuntimeContext, Task } from '../../mol-task';
import { deepEqual } from '../../mol-util';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { RootStructureDefinition } from '../helpers/root-structure';
import { createStructureComponent, StructureComponentParams, updateStructureComponent } from '../helpers/structure-component';
import { StructureQueryHelper } from '../helpers/structure-query';
import { StructureSelectionQueries } from '../helpers/structure-selection-query';
import { PluginStateObject as SO, PluginStateTransform } from '../objects';
import { parseMol } from '../../mol-io/reader/mol/parser';
import { trajectoryFromMol } from '../../mol-model-formats/structure/mol';
import { trajectoryFromCifCore } from '../../mol-model-formats/structure/cif-core';
import { trajectoryFromCube } from '../../mol-model-formats/structure/cube';
import { parseMol2 } from '../../mol-io/reader/mol2/parser';
import { trajectoryFromMol2 } from '../../mol-model-formats/structure/mol2';
import { parseXtc } from '../../mol-io/reader/xtc/parser';
import { coordinatesFromXtc } from '../../mol-model-formats/structure/xtc';

export { CoordinatesFromDcd };
export { CoordinatesFromXtc };
export { TopologyFromPsf };
export { TrajectoryFromModelAndCoordinates };
export { TrajectoryFromBlob };
export { TrajectoryFromMmCif };
export { TrajectoryFromPDB };
export { TrajectoryFromGRO };
export { TrajectoryFromMOL };
export { TrajectoryFromMOL2 };
export { TrajectoryFromCube };
export { TrajectoryFromCifCore };
export { TrajectoryFrom3DG };
export { ModelFromTrajectory };
export { StructureFromTrajectory };
export { StructureFromModel };
export { TransformStructureConformation };
export { StructureSelectionFromExpression };
export { MultiStructureSelectionFromExpression };
export { StructureSelectionFromScript };
export { StructureSelectionFromBundle };
export { StructureComplexElement };
export { StructureComponent };
export { CustomModelProperties };
export { CustomStructureProperties };
export { ShapeFromPly };

type CoordinatesFromDcd = typeof CoordinatesFromDcd
const CoordinatesFromDcd = PluginStateTransform.BuiltIn({
    name: 'coordinates-from-dcd',
    display: { name: 'Parse DCD', description: 'Parse DCD binary data.' },
    from: [SO.Data.Binary],
    to: SO.Molecule.Coordinates
})({
    apply({ a }) {
        return Task.create('Parse DCD', async ctx => {
            const parsed = await parseDcd(a.data).runInContext(ctx);
            if (parsed.isError) throw new Error(parsed.message);
            const coordinates = await coordinatesFromDcd(parsed.result).runInContext(ctx);
            return new SO.Molecule.Coordinates(coordinates, { label: a.label, description: 'Coordinates' });
        });
    }
});

type CoordinatesFromXtc = typeof CoordinatesFromDcd
const CoordinatesFromXtc = PluginStateTransform.BuiltIn({
    name: 'coordinates-from-xtc',
    display: { name: 'Parse XTC', description: 'Parse XTC binary data.' },
    from: [SO.Data.Binary],
    to: SO.Molecule.Coordinates
})({
    apply({ a }) {
        return Task.create('Parse XTC', async ctx => {
            const parsed = await parseXtc(a.data).runInContext(ctx);
            if (parsed.isError) throw new Error(parsed.message);
            const coordinates = await coordinatesFromXtc(parsed.result).runInContext(ctx);
            return new SO.Molecule.Coordinates(coordinates, { label: a.label, description: 'Coordinates' });
        });
    }
});

type TopologyFromPsf = typeof TopologyFromPsf
const TopologyFromPsf = PluginStateTransform.BuiltIn({
    name: 'topology-from-psf',
    display: { name: 'PSF Topology', description: 'Parse PSF string data.' },
    from: [SO.Format.Psf],
    to: SO.Molecule.Topology
})({
    apply({ a }) {
        return Task.create('Create Topology', async ctx => {
            const topology = await topologyFromPsf(a.data).runInContext(ctx);
            return new SO.Molecule.Topology(topology, { label: topology.label || a.label, description: 'Topology' });
        });
    }
});

async function getTrajectory(ctx: RuntimeContext, obj: StateObject, coordinates: Coordinates) {
    if (obj.type === SO.Molecule.Topology.type) {
        const topology = obj.data as Topology;
        return await Model.trajectoryFromTopologyAndCoordinates(topology, coordinates).runInContext(ctx);
    } else if (obj.type === SO.Molecule.Model.type) {
        const model = obj.data as Model;
        return Model.trajectoryFromModelAndCoordinates(model, coordinates);
    }
    throw new Error('no model/topology found');
}

type TrajectoryFromModelAndCoordinates = typeof TrajectoryFromModelAndCoordinates
const TrajectoryFromModelAndCoordinates = PluginStateTransform.BuiltIn({
    name: 'trajectory-from-model-and-coordinates',
    display: { name: 'Trajectory from Topology & Coordinates', description: 'Create a trajectory from existing model/topology and coordinates.' },
    from: SO.Root,
    to: SO.Molecule.Trajectory,
    params: {
        modelRef: PD.Text('', { isHidden: true }),
        coordinatesRef: PD.Text('', { isHidden: true }),
    }
})({
    apply({ params, dependencies }) {
        return Task.create('Create trajectory from model/topology and coordinates', async ctx => {
            const coordinates = dependencies![params.coordinatesRef].data as Coordinates;
            const trajectory = await getTrajectory(ctx, dependencies![params.modelRef], coordinates);
            const props = { label: 'Trajectory', description: `${trajectory.length} model${trajectory.length === 1 ? '' : 's'}` };
            return new SO.Molecule.Trajectory(trajectory, props);
        });
    }
});

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

            const props = { label: 'Trajectory', description: `${models.length} model${models.length === 1 ? '' : 's'}` };
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
            const props = { label: `${models[0].entry}`, description: `${models.length} model${models.length === 1 ? '' : 's'}` };
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
            const props = { label: `${models[0].entry}`, description: `${models.length} model${models.length === 1 ? '' : 's'}` };
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
            const props = { label: `${models[0].entry}`, description: `${models.length} model${models.length === 1 ? '' : 's'}` };
            return new SO.Molecule.Trajectory(models, props);
        });
    }
});

type TrajectoryFromMOL = typeof TrajectoryFromMOL
const TrajectoryFromMOL = PluginStateTransform.BuiltIn({
    name: 'trajectory-from-mol',
    display: { name: 'Parse MOL', description: 'Parse MOL string and create trajectory.' },
    from: [SO.Data.String],
    to: SO.Molecule.Trajectory
})({
    apply({ a }) {
        return Task.create('Parse MOL', async ctx => {
            const parsed = await parseMol(a.data).runInContext(ctx);
            if (parsed.isError) throw new Error(parsed.message);
            const models = await trajectoryFromMol(parsed.result).runInContext(ctx);
            const props = { label: `${models[0].entry}`, description: `${models.length} model${models.length === 1 ? '' : 's'}` };
            return new SO.Molecule.Trajectory(models, props);
        });
    }
});

type TrajectoryFromMOL2 = typeof TrajectoryFromMOL
const TrajectoryFromMOL2 = PluginStateTransform.BuiltIn({
    name: 'trajectory-from-mol2',
    display: { name: 'Parse MOL2', description: 'Parse MOL2 string and create trajectory.' },
    from: [SO.Data.String],
    to: SO.Molecule.Trajectory
})({
    apply({ a }) {
        return Task.create('Parse MOL2', async ctx => {
            const parsed = await parseMol2(a.data, a.label).runInContext(ctx);
            if (parsed.isError) throw new Error(parsed.message);
            const models = await trajectoryFromMol2(parsed.result).runInContext(ctx);
            const props = { label: `${models[0].entry}`, description: `${models.length} model${models.length === 1 ? '' : 's'}` };
            return new SO.Molecule.Trajectory(models, props);
        });
    }
});

type TrajectoryFromCube = typeof TrajectoryFromCube
const TrajectoryFromCube = PluginStateTransform.BuiltIn({
    name: 'trajectory-from-cube',
    display: { name: 'Parse Cube', description: 'Parse Cube file to create a trajectory.' },
    from: SO.Format.Cube,
    to: SO.Molecule.Trajectory
})({
    apply({ a }) {
        return Task.create('Parse MOL', async ctx => {
            const models = await trajectoryFromCube(a.data).runInContext(ctx);
            const props = { label: `${models[0].entry}`, description: `${models.length} model${models.length === 1 ? '' : 's'}` };
            return new SO.Molecule.Trajectory(models, props);
        });
    }
});

type TrajectoryFromCifCore = typeof TrajectoryFromCifCore
const TrajectoryFromCifCore = PluginStateTransform.BuiltIn({
    name: 'trajectory-from-cif-core',
    display: { name: 'Parse CIF Core', description: 'Identify and create all separate models in the specified CIF data block' },
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
    apply({ a, params }) {
        return Task.create('Parse CIF Core', async ctx => {
            const header = params.blockHeader || a.data.blocks[0].header;
            const block = a.data.blocks.find(b => b.header === header);
            if (!block) throw new Error(`Data block '${[header]}' not found.`);
            const models = await trajectoryFromCifCore(block).runInContext(ctx);
            if (models.length === 0) throw new Error('No models found.');
            const props = { label: `${models[0].entry}`, description: `${models.length} model${models.length === 1 ? '' : 's'}` };
            return new SO.Molecule.Trajectory(models, props);
        });
    }
});

type TrajectoryFrom3DG = typeof TrajectoryFrom3DG
const TrajectoryFrom3DG = PluginStateTransform.BuiltIn({
    name: 'trajectory-from-3dg',
    display: { name: 'Parse 3DG', description: 'Parse 3DG string and create trajectory.' },
    from: [SO.Data.String],
    to: SO.Molecule.Trajectory
})({
    apply({ a }) {
        return Task.create('Parse 3DG', async ctx => {
            const parsed = await parse3DG(a.data).runInContext(ctx);
            if (parsed.isError) throw new Error(parsed.message);
            const models = await trajectoryFrom3DG(parsed.result).runInContext(ctx);
            const props = { label: `${models[0].entry}`, description: `${models.length} model${models.length === 1 ? '' : 's'}` };
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
        return { modelIndex: PD.Converted(plus1, minus1, PD.Numeric(1, { min: 1, max: a.data.length, step: 1 }, { description: 'Model Index' })) };
    }
})({
    isApplicable: a => a.data.length > 0,
    apply({ a, params }) {
        let modelIndex = params.modelIndex % a.data.length;
        if (modelIndex < 0) modelIndex += a.data.length;
        const model = a.data[params.modelIndex];
        const label = `Model ${model.modelNum}`;
        const description = a.data.length === 1 ? undefined : `of ${a.data.length}`;
        return new SO.Molecule.Model(model, { label, description });
    },
    dispose({ b }) {
        b?.data.customProperties.dispose();
    }
});

type StructureFromTrajectory = typeof StructureFromTrajectory
const StructureFromTrajectory = PluginStateTransform.BuiltIn({
    name: 'structure-from-trajectory',
    display: { name: 'Structure from Trajectory', description: 'Create a molecular structure from a trajectory.' },
    from: SO.Molecule.Trajectory,
    to: SO.Molecule.Structure
})({
    apply({ a }) {
        return Task.create('Build Structure', async ctx => {
            const s = Structure.ofTrajectory(a.data);
            const props = { label: 'Ensemble', description: Structure.elementDescription(s) };
            return new SO.Molecule.Structure(s, props);
        });
    },
    dispose({ b }) {
        b?.data.customPropertyDescriptors.dispose();
    }
});

type StructureFromModel = typeof StructureFromModel
const StructureFromModel = PluginStateTransform.BuiltIn({
    name: 'structure-from-model',
    display: { name: 'Structure', description: 'Create a molecular structure (model, assembly, or symmetry) from the specified model.' },
    from: SO.Molecule.Model,
    to: SO.Molecule.Structure,
    params(a) { return RootStructureDefinition.getParams(a && a.data); }
})({
    canAutoUpdate({ oldParams, newParams }) {
        return RootStructureDefinition.canAutoUpdate(oldParams.type, newParams.type);
    },
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Build Structure', async ctx => {
            return RootStructureDefinition.create(plugin, ctx, a.data, params && params.type);
        });
    },
    update: ({ a, b, oldParams, newParams }) => {
        if (!b.data.models.includes(a.data)) return StateTransformer.UpdateResult.Recreate;
        if (!deepEqual(oldParams, newParams)) return StateTransformer.UpdateResult.Recreate;
        return StateTransformer.UpdateResult.Unchanged;
    },
    dispose({ b }) {
        b?.data.customPropertyDescriptors.dispose();
    }
});

const _translation = Vec3(), _m = Mat4(), _n = Mat4();

// type StructureCoordinateSystem = typeof StructureCoordinateSystem
// const StructureCoordinateSystem = PluginStateTransform.BuiltIn({
//     name: 'structure-coordinate-system',
//     display: { name: 'Coordinate System' },
//     isDecorator: true,
//     from: SO.Molecule.Structure,
//     to: SO.Molecule.Structure,
//     params: {
//         transform: PD.MappedStatic('components', {
//             components: PD.Group({
//                 axis: PD.Vec3(Vec3.create(1, 0, 0)),
//                 angle: PD.Numeric(0, { min: -180, max: 180, step: 0.1 }),
//                 translation: PD.Vec3(Vec3.create(0, 0, 0)),
//             }, { isFlat: true }),
//             matrix: PD.Group({
//                 data: PD.Mat4(Mat4.identity()),
//                 transpose: PD.Boolean(false)
//             }, { isFlat: true })
//         }, { label: 'Kind' })
//     }
// })({
//     canAutoUpdate({ newParams }) {
//         return newParams.transform.name === 'components';
//     },
//     apply({ a, params }) {
//         // TODO: optimze

//         const transform = Mat4();

//         if (params.transform.name === 'components') {
//             const { axis, angle, translation } = params.transform.params;
//             const center = a.data.boundary.sphere.center;
//             Mat4.fromTranslation(_m, Vec3.negate(_translation, center));
//             Mat4.fromTranslation(_n, Vec3.add(_translation, center, translation));
//             const rot = Mat4.fromRotation(Mat4(), Math.PI / 180 * angle, Vec3.normalize(Vec3(), axis));
//             Mat4.mul3(transform, _n, rot, _m);
//         } else if (params.transform.name === 'matrix') {
//             Mat4.copy(transform, params.transform.params.data);
//             if (params.transform.params.transpose) Mat4.transpose(transform, transform);
//         }

//         // TODO: compose with parent's coordinate system
//         a.data.coordinateSystem = SymmetryOperator.create('CS', transform);
//         return new SO.Molecule.Structure(a.data, { label: a.label, description: `${a.description} [Transformed]` });
//     }
// });

type TransformStructureConformation = typeof TransformStructureConformation
const TransformStructureConformation = PluginStateTransform.BuiltIn({
    name: 'transform-structure-conformation',
    display: { name: 'Transform Conformation' },
    isDecorator: true,
    from: SO.Molecule.Structure,
    to: SO.Molecule.Structure,
    params: {
        transform: PD.MappedStatic('components', {
            components: PD.Group({
                axis: PD.Vec3(Vec3.create(1, 0, 0)),
                angle: PD.Numeric(0, { min: -180, max: 180, step: 0.1 }),
                translation: PD.Vec3(Vec3.create(0, 0, 0)),
            }, { isFlat: true }),
            matrix: PD.Group({
                data: PD.Mat4(Mat4.identity()),
                transpose: PD.Boolean(false)
            }, { isFlat: true })
        }, { label: 'Kind' })
    }
})({
    canAutoUpdate({ newParams }) {
        return newParams.transform.name !== 'matrix';
    },
    apply({ a, params }) {
        // TODO: optimze
        // TODO: think of ways how to fast-track changes to this for animations

        const transform = Mat4();

        if (params.transform.name === 'components') {
            const { axis, angle, translation } = params.transform.params;
            const center = a.data.boundary.sphere.center;
            Mat4.fromTranslation(_m, Vec3.negate(_translation, center));
            Mat4.fromTranslation(_n, Vec3.add(_translation, center, translation));
            const rot = Mat4.fromRotation(Mat4(), Math.PI / 180 * angle, Vec3.normalize(Vec3(), axis));
            Mat4.mul3(transform, _n, rot, _m);
        } else if (params.transform.name === 'matrix') {
            Mat4.copy(transform, params.transform.params.data);
            if (params.transform.params.transpose) Mat4.transpose(transform, transform);
        }

        const s = Structure.transform(a.data, transform);
        return new SO.Molecule.Structure(s, { label: a.label, description: `${a.description} [Transformed]` });
    },
    dispose({ b }) {
        b?.data.customPropertyDescriptors.dispose();
    }
    // interpolate(src, tar, t) {
    //     // TODO: optimize
    //     const u = Mat4.fromRotation(Mat4(), Math.PI / 180 * src.angle, Vec3.normalize(Vec3(), src.axis));
    //     Mat4.setTranslation(u, src.translation);
    //     const v = Mat4.fromRotation(Mat4(), Math.PI / 180 * tar.angle, Vec3.normalize(Vec3(), tar.axis));
    //     Mat4.setTranslation(v, tar.translation);
    //     const m = SymmetryOperator.slerp(Mat4(), u, v, t);
    //     const rot = Mat4.getRotation(Quat.zero(), m);
    //     const axis = Vec3();
    //     const angle = Quat.getAxisAngle(axis, rot);
    //     const translation = Mat4.getTranslation(Vec3(), m);
    //     return { axis, angle, translation };
    // }
});

type StructureSelectionFromExpression = typeof StructureSelectionFromExpression
const StructureSelectionFromExpression = PluginStateTransform.BuiltIn({
    name: 'structure-selection-from-expression',
    display: { name: 'Selection', description: 'Create a molecular structure from the specified expression.' },
    from: SO.Molecule.Structure,
    to: SO.Molecule.Structure,
    params: {
        expression: PD.Value<Expression>(MolScriptBuilder.struct.generator.all, { isHidden: true }),
        label: PD.Optional(PD.Text('', { isHidden: true }))
    }
})({
    apply({ a, params, cache }) {
        const { selection, entry } = StructureQueryHelper.createAndRun(a.data, params.expression);
        (cache as any).entry = entry;

        if (Sel.isEmpty(selection)) return StateObject.Null;
        const s = Sel.unionStructure(selection);
        const props = { label: `${params.label || 'Selection'}`, description: Structure.elementDescription(s) };
        return new SO.Molecule.Structure(s, props);
    },
    update: ({ a, b, oldParams, newParams, cache }) => {
        if (oldParams.expression !== newParams.expression) return StateTransformer.UpdateResult.Recreate;

        const entry = (cache as { entry: StructureQueryHelper.CacheEntry }).entry;

        if (entry.currentStructure === a.data) {
            return StateTransformer.UpdateResult.Unchanged;
        }

        const selection = StructureQueryHelper.updateStructure(entry, a.data);
        if (Sel.isEmpty(selection)) return StateTransformer.UpdateResult.Null;

        StructureQueryHelper.updateStructureObject(b, selection, newParams.label);
        return StateTransformer.UpdateResult.Updated;
    },
    dispose({ b }) {
        b?.data.customPropertyDescriptors.dispose();
    }
});

type MultiStructureSelectionFromExpression = typeof MultiStructureSelectionFromExpression
const MultiStructureSelectionFromExpression = PluginStateTransform.BuiltIn({
    name: 'structure-multi-selection-from-expression',
    display: { name: 'Multi-structure Measurement Selection', description: 'Create selection object from multiple structures.' },
    from: SO.Root,
    to: SO.Molecule.Structure.Selections,
    params: {
        selections: PD.ObjectList({
            key: PD.Text(void 0, { description: 'A unique key.' }),
            ref: PD.Text(),
            groupId: PD.Optional(PD.Text()),
            expression: PD.Value<Expression>(MolScriptBuilder.struct.generator.empty)
        }, e => e.ref, { isHidden: true }),
        isTransitive: PD.Optional(PD.Boolean(false, { isHidden: true, description: 'Remap the selections from the original structure if structurally equivalent.' })),
        label: PD.Optional(PD.Text('', { isHidden: true }))
    }
})({
    apply({ params, cache, dependencies }) {
        const entries = new Map<string, StructureQueryHelper.CacheEntry>();

        const selections: SO.Molecule.Structure.SelectionEntry[] = [];
        let totalSize = 0;

        for (const sel of params.selections) {
            const { selection, entry } = StructureQueryHelper.createAndRun(dependencies![sel.ref].data as Structure, sel.expression);
            entries.set(sel.key, entry);
            const loci = Sel.toLociWithSourceUnits(selection);
            selections.push({ key: sel.key, loci, groupId: sel.groupId });
            totalSize += StructureElement.Loci.size(loci);
        }

        (cache as object as any).entries = entries;

        const props = { label: `${params.label || 'Multi-selection'}`, description: `${params.selections.length} source(s), ${totalSize} element(s) total` };
        return new SO.Molecule.Structure.Selections(selections, props);
    },
    update: ({ b, oldParams, newParams, cache, dependencies }) => {
        if (!!oldParams.isTransitive !== !!newParams.isTransitive) return StateTransformer.UpdateResult.Recreate;

        const cacheEntries = (cache as any).entries as Map<string, StructureQueryHelper.CacheEntry>;
        const entries = new Map<string, StructureQueryHelper.CacheEntry>();

        const current = new Map<string, SO.Molecule.Structure.SelectionEntry>();
        for (const e of b.data) current.set(e.key, e);

        let changed = false;
        let totalSize = 0;

        const selections: SO.Molecule.Structure.SelectionEntry[] = [];
        for (const sel of newParams.selections) {
            const structure = dependencies![sel.ref].data as Structure;

            let recreate = false;

            if (cacheEntries.has(sel.key)) {
                const entry = cacheEntries.get(sel.key)!;
                if (StructureQueryHelper.isUnchanged(entry, sel.expression, structure) && current.has(sel.key)) {
                    const loci = current.get(sel.key)!;
                    if (loci.groupId !== sel.groupId) {
                        loci.groupId = sel.groupId;
                        changed = true;
                    }
                    entries.set(sel.key, entry);
                    selections.push(loci);
                    totalSize += StructureElement.Loci.size(loci.loci);

                    continue;
                } if (entry.expression !== sel.expression) {
                    recreate = true;
                } else {
                    // TODO: properly support "transitive" queries. For that Structure.areUnitAndIndicesEqual needs to be fixed;
                    let update = false;

                    if (!!newParams.isTransitive) {
                        if (Structure.areUnitAndIndicesEqual(entry.originalStructure, structure)) {
                            const selection = StructureQueryHelper.run(entry, entry.originalStructure);
                            entry.currentStructure = structure;
                            entries.set(sel.key, entry);
                            const loci = StructureElement.Loci.remap(Sel.toLociWithSourceUnits(selection), structure);
                            selections.push({ key: sel.key, loci, groupId: sel.groupId });
                            totalSize += StructureElement.Loci.size(loci);
                            changed = true;
                        } else {
                            update = true;
                        }
                    } else {
                        update = true;
                    }

                    if (update) {
                        changed = true;
                        const selection = StructureQueryHelper.updateStructure(entry, structure);
                        entries.set(sel.key, entry);
                        const loci = Sel.toLociWithSourceUnits(selection);
                        selections.push({ key: sel.key, loci, groupId: sel.groupId });
                        totalSize += StructureElement.Loci.size(loci);
                    }
                }
            } else {
                recreate = true;
            }

            if (recreate) {
                changed = true;

                // create new selection
                const { selection, entry } = StructureQueryHelper.createAndRun(structure, sel.expression);
                entries.set(sel.key, entry);
                const loci = Sel.toLociWithSourceUnits(selection);
                selections.push({ key: sel.key, loci });
                totalSize += StructureElement.Loci.size(loci);
            }
        }

        if (!changed) return StateTransformer.UpdateResult.Unchanged;

        (cache as object as any).entries = entries;
        b.data = selections;
        b.label = `${newParams.label || 'Multi-selection'}`;
        b.description = `${selections.length} source(s), ${totalSize} element(s) total`;

        return StateTransformer.UpdateResult.Updated;
    }
});

type StructureSelectionFromScript = typeof StructureSelectionFromScript
const StructureSelectionFromScript = PluginStateTransform.BuiltIn({
    name: 'structure-selection-from-script',
    display: { name: 'Selection', description: 'Create a molecular structure from the specified script.' },
    from: SO.Molecule.Structure,
    to: SO.Molecule.Structure,
    params: {
        script: PD.Script({ language: 'mol-script', expression: '(sel.atom.atom-groups :residue-test (= atom.resname ALA))' }),
        label: PD.Optional(PD.Text(''))
    }
})({
    apply({ a, params, cache }) {
        const { selection, entry } = StructureQueryHelper.createAndRun(a.data, params.script);
        (cache as any).entry = entry;

        const s = Sel.unionStructure(selection);
        const props = { label: `${params.label || 'Selection'}`, description: Structure.elementDescription(s) };
        return new SO.Molecule.Structure(s, props);
    },
    update: ({ a, b, oldParams, newParams, cache }) => {
        if (!Script.areEqual(oldParams.script, newParams.script)) {
            return StateTransformer.UpdateResult.Recreate;
        }

        const entry = (cache as { entry: StructureQueryHelper.CacheEntry }).entry;

        if (entry.currentStructure === a.data) {
            return StateTransformer.UpdateResult.Unchanged;
        }

        const selection = StructureQueryHelper.updateStructure(entry, a.data);
        StructureQueryHelper.updateStructureObject(b, selection, newParams.label);
        return StateTransformer.UpdateResult.Updated;
    },
    dispose({ b }) {
        b?.data.customPropertyDescriptors.dispose();
    }
});

type StructureSelectionFromBundle = typeof StructureSelectionFromBundle
const StructureSelectionFromBundle = PluginStateTransform.BuiltIn({
    name: 'structure-selection-from-bundle',
    display: { name: 'Selection', description: 'Create a molecular structure from the specified structure-element bundle.' },
    from: SO.Molecule.Structure,
    to: SO.Molecule.Structure,
    params: {
        bundle: PD.Value<StructureElement.Bundle>(StructureElement.Bundle.Empty, { isHidden: true }),
        label: PD.Optional(PD.Text('', { isHidden: true }))
    }
})({
    apply({ a, params, cache }) {
        if (params.bundle.hash !== a.data.hashCode) {
            return StateObject.Null;
        }

        (cache as { source: Structure }).source = a.data;

        const s = StructureElement.Bundle.toStructure(params.bundle, a.data);
        if (s.elementCount === 0) return StateObject.Null;

        const props = { label: `${params.label || 'Selection'}`, description: Structure.elementDescription(s) };
        return new SO.Molecule.Structure(s, props);
    },
    update: ({ a, b, oldParams, newParams, cache }) => {
        if (!StructureElement.Bundle.areEqual(oldParams.bundle, newParams.bundle)) {
            return StateTransformer.UpdateResult.Recreate;
        }

        if (newParams.bundle.hash !== a.data.hashCode) {
            return StateTransformer.UpdateResult.Null;
        }

        if ((cache as { source: Structure }).source === a.data) {
            return StateTransformer.UpdateResult.Unchanged;
        }
        (cache as { source: Structure }).source = a.data;

        const s = StructureElement.Bundle.toStructure(newParams.bundle, a.data);
        if (s.elementCount === 0) return StateTransformer.UpdateResult.Null;

        b.label = `${newParams.label || 'Selection'}`;
        b.description = Structure.elementDescription(s);
        b.data = s;
        return StateTransformer.UpdateResult.Updated;
    },
    dispose({ b }) {
        b?.data.customPropertyDescriptors.dispose();
    }
});

export const StructureComplexElementTypes = {
    'polymer': 'polymer',

    'protein': 'protein',
    'nucleic': 'nucleic',
    'water': 'water',

    'branched': 'branched', // = carbs
    'ligand': 'ligand',
    'non-standard': 'non-standard',

    'coarse': 'coarse',

    // Legacy
    'atomic-sequence': 'atomic-sequence',
    'atomic-het': 'atomic-het',
    'spheres': 'spheres'
} as const;
export type StructureComplexElementTypes = keyof typeof StructureComplexElementTypes

const StructureComplexElementTypeTuples = PD.objectToOptions(StructureComplexElementTypes);

type StructureComplexElement = typeof StructureComplexElement
const StructureComplexElement = PluginStateTransform.BuiltIn({
    name: 'structure-complex-element',
    display: { name: 'Complex Element', description: 'Create a molecular structure from the specified model.' },
    from: SO.Molecule.Structure,
    to: SO.Molecule.Structure,
    params: { type: PD.Select<StructureComplexElementTypes>('atomic-sequence', StructureComplexElementTypeTuples, { isHidden: true }) }
})({
    apply({ a, params }) {
        // TODO: update function.

        let query: StructureQuery, label: string;
        switch (params.type) {
            case 'polymer': query = StructureSelectionQueries.polymer.query; label = 'Polymer'; break;

            case 'protein': query = StructureSelectionQueries.protein.query; label = 'Protein'; break;
            case 'nucleic': query = StructureSelectionQueries.nucleic.query; label = 'Nucleic'; break;
            case 'water': query = Queries.internal.water(); label = 'Water'; break;

            case 'branched': query = StructureSelectionQueries.branchedPlusConnected.query; label = 'Branched'; break;
            case 'ligand': query = StructureSelectionQueries.ligandPlusConnected.query; label = 'Ligand'; break;

            case 'non-standard': query = StructureSelectionQueries.nonStandardPolymer.query; label = 'Non-standard'; break;

            case 'coarse': query = StructureSelectionQueries.coarse.query; label = 'Coarse'; break;

            case 'atomic-sequence': query = Queries.internal.atomicSequence(); label = 'Sequence'; break;
            case 'atomic-het': query = Queries.internal.atomicHet(); label = 'HET Groups/Ligands'; break;
            case 'spheres': query = Queries.internal.spheres(); label = 'Coarse Spheres'; break;

            default: throw new Error(`${params.type} is a not valid complex element.`);
        }

        const result = query(new QueryContext(a.data));
        const s = Sel.unionStructure(result);

        if (s.elementCount === 0) return StateObject.Null;
        return new SO.Molecule.Structure(s, { label, description: Structure.elementDescription(s) });
    },
    dispose({ b }) {
        b?.data.customPropertyDescriptors.dispose();
    }
});

type StructureComponent = typeof StructureComponent
const StructureComponent = PluginStateTransform.BuiltIn({
    name: 'structure-component',
    display: { name: 'Component', description: 'A molecular structure component.' },
    from: SO.Molecule.Structure,
    to: SO.Molecule.Structure,
    params: StructureComponentParams
})({
    apply({ a, params, cache }) {
        return createStructureComponent(a.data, params, cache as any);
    },
    update: ({ a, b, oldParams, newParams, cache }) => {
        return updateStructureComponent(a.data, b, oldParams, newParams, cache as any);
    },
    dispose({ b }) {
        b?.data.customPropertyDescriptors.dispose();
    }
});

type CustomModelProperties = typeof CustomModelProperties
const CustomModelProperties = PluginStateTransform.BuiltIn({
    name: 'custom-model-properties',
    display: { name: 'Custom Model Properties' },
    isDecorator: true,
    from: SO.Molecule.Model,
    to: SO.Molecule.Model,
    params: (a, ctx: PluginContext) => {
        return ctx.customModelProperties.getParams(a?.data);
    }
})({
    apply({ a, params }, ctx: PluginContext) {
        return Task.create('Custom Props', async taskCtx => {
            await attachModelProps(a.data, ctx, taskCtx, params);
            return new SO.Molecule.Model(a.data, { label: a.label, description: a.description });
        });
    },
    update({ a, b, oldParams, newParams }, ctx: PluginContext) {
        return Task.create('Custom Props', async taskCtx => {
            b.data = a.data;
            b.label = a.label;
            b.description = a.description;
            for (const name of oldParams.autoAttach) {
                const property = ctx.customModelProperties.get(name);
                if (!property) continue;
                a.data.customProperties.reference(property.descriptor, false);
            }
            await attachModelProps(a.data, ctx, taskCtx, newParams);
            return StateTransformer.UpdateResult.Updated;
        });
    },
    dispose({ b }) {
        b?.data.customProperties.dispose();
    }
});
async function attachModelProps(model: Model, ctx: PluginContext, taskCtx: RuntimeContext, params: ReturnType<CustomModelProperties['createDefaultParams']>) {
    const propertyCtx = { runtime: taskCtx, assetManager: ctx.managers.asset };
    const { autoAttach, properties } = params;
    for (const name of Object.keys(properties)) {
        const property = ctx.customModelProperties.get(name);
        const props = properties[name];
        if (autoAttach.includes(name) || property.isHidden) {
            try {
                await property.attach(propertyCtx, model, props, true);
            } catch (e) {
                ctx.log.warn(`Error attaching model prop '${name}': ${e}`);
            }
        } else {
            property.set(model, props);
        }
    }
}

type CustomStructureProperties = typeof CustomStructureProperties
const CustomStructureProperties = PluginStateTransform.BuiltIn({
    name: 'custom-structure-properties',
    display: { name: 'Custom Structure Properties' },
    isDecorator: true,
    from: SO.Molecule.Structure,
    to: SO.Molecule.Structure,
    params: (a, ctx: PluginContext) => {
        return ctx.customStructureProperties.getParams(a?.data.root);
    }
})({
    apply({ a, params }, ctx: PluginContext) {
        return Task.create('Custom Props', async taskCtx => {
            await attachStructureProps(a.data.root, ctx, taskCtx, params);
            return new SO.Molecule.Structure(a.data, { label: a.label, description: a.description });
        });
    },
    update({ a, b, oldParams, newParams }, ctx: PluginContext) {
        return Task.create('Custom Props', async taskCtx => {
            b.data = a.data;
            b.label = a.label;
            b.description = a.description;
            for (const name of oldParams.autoAttach) {
                const property = ctx.customStructureProperties.get(name);
                if (!property) continue;
                a.data.customPropertyDescriptors.reference(property.descriptor, false);
            }
            await attachStructureProps(a.data.root, ctx, taskCtx, newParams);
            return StateTransformer.UpdateResult.Updated;
        });
    },
    dispose({ b }) {
        b?.data.customPropertyDescriptors.dispose();
    }
});
async function attachStructureProps(structure: Structure, ctx: PluginContext, taskCtx: RuntimeContext, params: ReturnType<CustomStructureProperties['createDefaultParams']>) {
    const propertyCtx = { runtime: taskCtx, assetManager: ctx.managers.asset };
    const { autoAttach, properties } = params;
    for (const name of Object.keys(properties)) {
        const property = ctx.customStructureProperties.get(name);
        const props = properties[name];
        if (autoAttach.includes(name) || property.isHidden) {
            try {
                await property.attach(propertyCtx, structure, props, true);
            } catch (e) {
                ctx.log.warn(`Error attaching structure prop '${name}': ${e}`);
            }
        } else {
            property.set(structure, props);
        }
    }
}

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
            const shape = await shapeFromPly(a.data, params).runInContext(ctx);
            const props = { label: 'Shape' };
            return new SO.Shape.Provider(shape, props);
        });
    }
});