/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { parseRelionStar, getRelionStarTomogramNames, getRelionStarMicrographNames } from '../../mol-io/reader/relion/star';
import { createParticleListFromCryoEtDataPortalNdjson } from '../../mol-model-formats/particles/ndjson';
import { createParticleListFromRelionStar } from '../../mol-model-formats/particles/star';
import { createParticleListFromDynamoTbl, getDynamoTblTomogramIds } from '../../mol-model-formats/particles/tbl';
import { createParticleListFromArtiatomiEm, getArtiatomiMotivelistTomogramIds } from '../../mol-model-formats/particles/em';
import { createSimulariumParticleTrajectory, getSimulariumAgentTypeNames, getSimulariumFrameCount } from '../../mol-model-formats/particles/simularium';
import { createParticleListFromMmcifAssembly, getAssemblyIdsFromMmcif, getAsymIdsFromMmcif } from '../../mol-model-formats/particles/mmcif';
import { PluginContext } from '../../mol-plugin/context';
import { StateTransform, StateTransformer } from '../../mol-state';
import { Task } from '../../mol-task';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Theme } from '../../mol-theme/theme';
import { PluginStateObject as SO, PluginStateTransform } from '../objects';
import { Structure, Trajectory } from '../../mol-model/structure';
import { ParticleList, Particle, ParticleTarget } from '../../mol-model/particles/particle-list';
import { buildTargetStructuresFromMapping } from '../../mol-model/particles/structure-mapping';
import { ShapeProvider } from '../../mol-model/shape/provider';

export { ParticleListFromRelionStar };
export { ParticleListFromDynamoTbl };
export { ParticleListFromCryoEtDataPortalNdjson };
export { ParticleListFromArtiatomiEm };
export { ParticleTrajectoryFromSimularium };
export { ParticleListFromTrajectory };
export { ParticleListFromMmcifAssembly };
export { ParticleListWithTargets };
export { ParticlesRepresentation3D };

type ParticleListFromRelionStar = typeof ParticleListFromRelionStar
const ParticleListFromRelionStar = PluginStateTransform.BuiltIn({
    name: 'particle-list-from-relion-star',
    display: { name: 'Particle List from RELION STAR', description: 'Create ParticleList from RELION STAR data.' },
    from: SO.Format.Cif,
    to: SO.Particle.List,
    params: a => {
        if (!a) {
            return {
                tomograms: PD.MultiSelect<string>([], [], { description: 'Empty selection includes all tomograms.' }),
                micrographs: PD.MultiSelect<string>([], [], { description: 'Empty selection includes all micrographs. Combined with the tomogram filter using AND.' }),
                pixelSize: PD.Optional(PD.Numeric(0, { min: 0, step: 0.001 }, { description: 'Override pixel size in Å/pixel for converting pixel-space coordinates to angstrom. Leave 0 to auto-detect from STAR optics/particle metadata.' })),
                particleRadius: PD.Numeric(0, { min: 0, step: 0.1 }, { description: 'Uniform particle radius in angstrom. Leave 0 to omit.' }),
            };
        }
        let tomoNames: string[] = [];
        let micrographNames: string[] = [];
        try {
            tomoNames = getRelionStarTomogramNames(a.data);
            micrographNames = getRelionStarMicrographNames(a.data);
        } catch {
            // ignore; apply will surface parse errors
        }
        const tomoOptions = tomoNames.map(n => [n, n] as [string, string]);
        const micrographOptions = micrographNames.map(n => [n, n] as [string, string]);
        const tomoDefault = tomoNames.length > 0 ? [tomoNames[0]] : [];
        const micrographDefault = micrographNames.length > 0 ? [micrographNames[0]] : [];
        return {
            tomograms: PD.MultiSelect<string>(tomoDefault, tomoOptions, { description: 'Empty selection includes all tomograms.' }),
            micrographs: PD.MultiSelect<string>(micrographDefault, micrographOptions, { description: 'Empty selection includes all micrographs. Combined with the tomogram filter using AND.' }),
            pixelSize: PD.Optional(PD.Numeric(0, { min: 0, step: 0.001 }, { description: 'Override pixel size in Å/pixel for converting pixel-space coordinates to angstrom. Leave 0 to auto-detect from STAR optics/particle metadata.' })),
            particleRadius: PD.Numeric(0, { min: 0, step: 0.1 }, { description: 'Uniform particle radius in angstrom. Leave 0 to omit.' }),
        };
    }
})({
    apply({ a, params }) {
        return Task.create('Create Particle List from RELION STAR', async ctx => {
            const relion = parseRelionStar(a.data);
            if (relion.isError) throw new Error(relion.message);

            const list = createParticleListFromRelionStar(relion.result, {
                tomograms: params.tomograms,
                micrographs: params.micrographs,
                pixelSize: params.pixelSize && params.pixelSize > 0 ? params.pixelSize : void 0,
                particleRadius: params.particleRadius > 0 ? params.particleRadius : void 0,
            });

            return new SO.Particle.List(list, { label: list.label || 'Particles', description: 'RELION Particle List' });
        });
    }
});

type ParticleListFromDynamoTbl = typeof ParticleListFromDynamoTbl
const ParticleListFromDynamoTbl = PluginStateTransform.BuiltIn({
    name: 'particle-list-from-dynamo-tbl',
    display: { name: 'Particle List from Dynamo TBL', description: 'Create ParticleList from Dynamo TBL data.' },
    from: SO.Format.DynamoTbl,
    to: SO.Particle.List,
    params: a => {
        if (!a) {
            return {
                tomos: PD.MultiSelect<string>([], [], { description: 'Empty selection includes all tomograms.' }),
                pixelSize: PD.Optional(PD.Numeric(0, { min: 0, step: 0.001 }, { description: 'Override pixel size in Å/pixel for converting pixel-space coordinates to angstrom. Leave 0 to auto-detect from the table’s `apix` field.' })),
                particleRadius: PD.Numeric(0, { min: 0, step: 0.1 }, { description: 'Uniform particle radius in angstrom. Leave 0 to omit.' }),
            };
        }
        const ids = getDynamoTblTomogramIds(a.data);
        const options = ids.map(id => [String(id), String(id)] as [string, string]);
        const defaultValue = ids.length > 0 ? [String(ids[0])] : [];
        return {
            tomos: PD.MultiSelect<string>(defaultValue, options, { description: 'Empty selection includes all tomograms.' }),
            pixelSize: PD.Optional(PD.Numeric(0, { min: 0, step: 0.001 }, { description: 'Override pixel size in Å/pixel for converting pixel-space coordinates to angstrom. Leave 0 to auto-detect from the table’s `apix` field.' })),
            particleRadius: PD.Numeric(0, { min: 0, step: 0.1 }, { description: 'Uniform particle radius in angstrom. Leave 0 to omit.' }),
        };
    }
})({
    apply({ a, params }) {
        return Task.create('Create Particle List from Dynamo TBL', async ctx => {
            const list = createParticleListFromDynamoTbl(a.data, {
                tomos: params.tomos.map(v => Number(v)),
                pixelSize: params.pixelSize && params.pixelSize > 0 ? params.pixelSize : void 0,
                particleRadius: params.particleRadius > 0 ? params.particleRadius : void 0,
            });
            return new SO.Particle.List(list, { label: list.label || 'Particles', description: 'Dynamo Particle List' });
        });
    }
});

type ParticleListFromCryoEtDataPortalNdjson = typeof ParticleListFromCryoEtDataPortalNdjson
const ParticleListFromCryoEtDataPortalNdjson = PluginStateTransform.BuiltIn({
    name: 'particle-list-from-cryoet-data-portal-ndjson',
    display: { name: 'Particle List from CryoET NDJSON', description: 'Create ParticleList from CryoET Data Portal NDJSON data.' },
    from: SO.Format.CryoEtDataPortalNdjson,
    to: SO.Particle.List,
    params: {
        pixelSize: PD.Numeric(1, { min: 0, step: 0.001 }, { description: 'Pixel size in Å/pixel used to convert pixel-space NDJSON coordinates to angstrom. Required because CryoET Data Portal NDJSON does not encode distance units.' }),
        type: PD.Optional(PD.Text('')),
        particleRadius: PD.Numeric(0, { min: 0, step: 0.1 }, { description: 'Uniform particle radius in angstrom. Leave 0 to omit.' }),
    }
})({
    apply({ a, params }) {
        return Task.create('Create Particle List from CryoET NDJSON', async ctx => {
            const list = createParticleListFromCryoEtDataPortalNdjson(a.data, {
                pixelSize: params.pixelSize,
                type: params.type || void 0,
                particleRadius: params.particleRadius > 0 ? params.particleRadius : void 0,
            });
            return new SO.Particle.List(list, { label: list.label || 'Particles', description: 'CryoET NDJSON Particle List' });
        });
    }
});

type ParticleListFromArtiatomiEm = typeof ParticleListFromArtiatomiEm
const ParticleListFromArtiatomiEm = PluginStateTransform.BuiltIn({
    name: 'particle-list-from-artiatomi-em',
    display: { name: 'Particle List from Artiatomi EM', description: 'Create ParticleList from Artiatomi EM motivelist data.' },
    from: SO.Format.ArtiatomiEm,
    to: SO.Particle.List,
    params: a => {
        if (!a) {
            return {
                tomos: PD.MultiSelect<string>([], [], { description: 'Empty selection includes all tomograms.' }),
                pixelSize: PD.Numeric(1, { min: 0, step: 0.001 }, { description: 'Pixel size in Å/pixel used to convert voxel-space coordinates to angstrom. Required because Artiatomi EM files do not encode distance units.' }),
                particleRadius: PD.Numeric(0, { min: 0, step: 0.1 }, { description: 'Uniform particle radius in angstrom. Leave 0 to omit.' }),
            };
        }
        const ids = getArtiatomiMotivelistTomogramIds(a.data);
        const options = ids.map(id => [String(id), String(id)] as [string, string]);
        const defaultValue = ids.length > 0 ? [String(ids[0])] : [];
        return {
            tomos: PD.MultiSelect<string>(defaultValue, options, { description: 'Empty selection includes all tomograms.' }),
            pixelSize: PD.Numeric(1, { min: 0, step: 0.001 }, { description: 'Pixel size in Å/pixel used to convert voxel-space coordinates to angstrom. Required because Artiatomi EM files do not encode distance units.' }),
            particleRadius: PD.Numeric(0, { min: 0, step: 0.1 }, { description: 'Uniform particle radius in angstrom. Leave 0 to omit.' }),
        };
    }
})({
    apply({ a, params }) {
        return Task.create('Create Particle List from Artiatomi EM', async () => {
            const list = createParticleListFromArtiatomiEm(a.data, {
                tomos: params.tomos.map(v => Number(v)),
                pixelSize: params.pixelSize,
                label: a.label,
                particleRadius: params.particleRadius > 0 ? params.particleRadius : void 0,
            });
            return new SO.Particle.List(list, { label: list.label || 'Particles', description: 'Artiatomi EM Particle List' });
        });
    }
});

type ParticleTrajectoryFromSimularium = typeof ParticleTrajectoryFromSimularium
const ParticleTrajectoryFromSimularium = PluginStateTransform.BuiltIn({
    name: 'particle-trajectory-from-simularium',
    display: { name: 'Particle Trajectory from Simularium', description: 'Create a ParticleTrajectory wrapping all frames of a Simularium file.' },
    from: SO.Format.Simularium,
    to: SO.Particle.Trajectory,
    params: a => {
        if (!a) {
            return {
                types: PD.MultiSelect<string>([], [], { description: 'Agent types to include. Empty selection includes all types.' }),
                scale: PD.Numeric(0, { min: 0, step: 0.001 }, { description: 'Spatial scale to angstrom. Leave 0 to auto-detect from the file spatial units.' }),
            };
        }
        const typeOptions = getSimulariumAgentTypeNames(a.data).map(t => [String(t.id), t.name] as [string, string]);
        return {
            types: PD.MultiSelect<string>([], typeOptions, { description: 'Agent types to include. Empty selection includes all types.' }),
            scale: PD.Numeric(0, { min: 0, step: 0.001 }, { description: 'Spatial scale to angstrom. Leave 0 to auto-detect from the file spatial units.' }),
        };
    }
})({
    apply({ a, params }) {
        const traj = createSimulariumParticleTrajectory(a.data, {
            scale: params.scale && params.scale > 0 ? params.scale : void 0,
            typeFilter: params.types.length > 0 ? params.types.map(v => Number(v)) : void 0,
        });
        const frameCount = getSimulariumFrameCount(a.data);
        return new SO.Particle.Trajectory(traj, { label: a.label, description: `${frameCount} frame${frameCount !== 1 ? 's' : ''}` });
    }
});

type ParticleListFromTrajectory = typeof ParticleListFromTrajectory
const ParticleListFromTrajectory = PluginStateTransform.BuiltIn({
    name: 'particle-list-from-trajectory',
    display: { name: 'Particle List from Trajectory', description: 'Extract a single frame from a ParticleTrajectory.' },
    from: SO.Particle.Trajectory,
    to: SO.Particle.List,
    params: a => ({
        frameIndex: PD.Numeric(0, { min: 0, max: a ? Math.max(0, a.data.frameCount - 1) : 0, step: 1 }, { description: 'Index of the trajectory frame to display.' }),
    })
})({
    apply({ a, params }) {
        const list = a.data.getFrameAtIndex(Math.max(0, Math.min(params.frameIndex, a.data.frameCount - 1)));
        return new SO.Particle.List(list, { label: list.label || 'Particles', description: `Frame ${params.frameIndex + 1} of ${a.data.frameCount}` });
    },
    // update({ a, b, oldParams, newParams }) {
    //     if (oldParams.frameIndex === newParams.frameIndex) return StateTransformer.UpdateResult.Unchanged;
    //     const list = a.data.getFrameAtIndex(Math.max(0, Math.min(newParams.frameIndex, a.data.frameCount - 1)));
    //     b.data = list;
    //     b.label = list.label || 'Particles';
    //     b.description = `Frame ${newParams.frameIndex + 1} of ${a.data.frameCount}`;
    //     return StateTransformer.UpdateResult.Updated;
    // }
});

type ParticleListFromMmcifAssembly = typeof ParticleListFromMmcifAssembly
const ParticleListFromMmcifAssembly = PluginStateTransform.BuiltIn({
    name: 'particle-list-from-mmcif-assembly',
    display: { name: 'Particle List from mmCIF Assembly', description: 'Create a ParticleList from _pdbx_struct_assembly_gen/_pdbx_struct_oper_list in a CIF file. Each expanded operator combination becomes one particle.' },
    from: SO.Format.Cif,
    to: SO.Particle.List,
    params: a => {
        if (!a) {
            return {
                assemblyId: PD.Text('1', { description: 'Assembly identifier from _pdbx_struct_assembly.id.' }),
                asymIds: PD.MultiSelect<string>([], [], { description: 'Individual chain IDs to include. Empty selection includes all chains.' }),
            };
        }
        let assemblyIds: string[] = [];
        let asymIds: string[] = [];
        try {
            assemblyIds = getAssemblyIdsFromMmcif(a.data);
            if (assemblyIds.length > 0) {
                asymIds = getAsymIdsFromMmcif(a.data, assemblyIds[0]);
            }
        } catch {
            // ignore; apply will surface parse errors
        }
        const assemblyOptions = assemblyIds.map(id => [id, id] as [string, string]);
        const assemblyDefault = assemblyIds[0] ?? '1';
        const asymOptions = asymIds.map(id => [id, id] as [string, string]);
        return {
            assemblyId: PD.Select<string>(assemblyDefault, assemblyOptions, { description: 'Assembly identifier from _pdbx_struct_assembly.id.' }),
            asymIds: PD.MultiSelect<string>([], asymOptions, { description: 'Individual chain IDs to include. Empty selection includes all chains.' }),
        };
    }
})({
    apply({ a, params }) {
        return Task.create('Create Particle List from mmCIF Assembly', async () => {
            const list = createParticleListFromMmcifAssembly(a.data, {
                assemblyId: params.assemblyId,
                asymIds: params.asymIds.length > 0 ? params.asymIds : undefined,
            });
            return new SO.Particle.List(list, { label: list.label || 'Particles', description: 'mmCIF Assembly Particle List' });
        });
    }
});

function getStructureRefOptions(ctx: PluginContext): [string, string][] {
    const out: [string, string][] = [];
    (ctx.state.data.cells as Map<string, any>).forEach((cell: any) => {
        if (cell.obj instanceof SO.Molecule.Structure) {
            out.push([cell.transform.ref, cell.obj.label ?? '<structure>']);
        }
    });
    return out;
}

function getTrajectoryRefOptions(ctx: PluginContext): [string, string][] {
    const out: [string, string][] = [];
    (ctx.state.data.cells as Map<string, any>).forEach((cell: any) => {
        if (cell.obj instanceof SO.Molecule.Trajectory) {
            out.push([cell.transform.ref, cell.obj.label ?? '<trajectory>']);
        }
    });
    return out;
}

function getShapeRefOptions(ctx: PluginContext): [string, string][] {
    const out: [string, string][] = [];
    (ctx.state.data.cells as Map<string, any>).forEach((cell: any) => {
        if (cell.obj instanceof SO.Shape.Provider) {
            out.push([cell.transform.ref, cell.obj.label ?? '<shape>']);
        }
    });
    return out;
}

type ParticleListWithTargets = typeof ParticleListWithTargets
const ParticleListWithTargets = PluginStateTransform.BuiltIn({
    name: 'particle-list-with-targets',
    display: { name: 'Particle List with Targets', description: 'Associate reference structures and shapes with a particle list for target-based representation.' },
    isDecorator: true,
    from: SO.Particle.List,
    to: SO.Particle.List,
    params(a) {
        const hasMappedTargets = !!(a?.data.targetMapping);
        const hasModelTargets = !!(a?.data.targetModels);
        return {
            // Used when ParticleList.targetModels is set: a trajectory whose models provide
            // one reference structure per target (petworld variant).
            trajectory: PD.ValueRef<Trajectory>(
                getTrajectoryRefOptions,
                (ref, getData) => getData(ref) as Trajectory,
                {
                    description: 'Trajectory whose models provide per-target reference structures using the particle list\'s targetModels mapping.',
                    isHidden: !hasModelTargets,
                }
            ),
            // Used when ParticleList.targetMapping is set: one parent structure that is
            // automatically split into per-target sub-structures by chain ID.
            structure: PD.ValueRef<Structure>(
                getStructureRefOptions,
                (ref, getData) => getData(ref) as Structure,
                {
                    description: 'Parent structure to automatically split into per-target sub-structures using the chain IDs in the particle list\'s targetMapping.',
                    isHidden: !hasMappedTargets,
                }
            ),
            // Used when ParticleList.targetMapping is NOT set: explicit per-target mapping.
            structures: PD.ObjectList({
                targetId: PD.Numeric(0, { min: 0, step: 1 }, { description: 'Target ID in the particle list this structure maps to (matches ParticleList.targets).' }),
                structure: PD.ValueRef<Structure>(
                    getStructureRefOptions,
                    (ref, getData) => getData(ref) as Structure,
                ),
            }, e => `Target ${e.targetId}`, {
                description: 'Reference structures mapped per particle target ID.',
                isHidden: hasMappedTargets || hasModelTargets,
            }),
            // Explicit per-target shape mapping (e.g. meshes from OBJ).
            shapes: PD.ObjectList({
                targetId: PD.Numeric(0, { min: 0, step: 1 }, { description: 'Target ID in the particle list this shape maps to (matches ParticleList.targets).' }),
                shape: PD.ValueRef<ShapeProvider<any, any, any>>(
                    getShapeRefOptions,
                    (ref, getData) => getData(ref) as ShapeProvider<any, any, any>,
                ),
            }, e => `Target ${e.targetId}`, {
                description: 'Reference shapes mapped per particle target ID.',
                isHidden: hasMappedTargets || hasModelTargets,
            }),
        };
    },
})({
    getDependencies: (params) => {
        const deps: StateTransform.Ref[] = [];
        if (params.trajectory?.ref) deps.push(params.trajectory.ref as StateTransform.Ref);
        if (params.structure?.ref) deps.push(params.structure.ref as StateTransform.Ref);
        for (const e of params.structures) {
            if (e.structure.ref) deps.push(e.structure.ref as StateTransform.Ref);
        }
        for (const e of params.shapes) {
            if (e.shape.ref) deps.push(e.shape.ref as StateTransform.Ref);
        }
        return deps;
    },
    apply({ a, params }) {
        return Task.create('Associate targets with particle list', async _ctx => {
            const targetMap = new Map<number, ParticleTarget>();

            // Structure targets: from trajectory models, a chain-split parent, or explicit refs.
            if (a.data.targetModels) {
                try {
                    const trajectory = params.trajectory.getValue();
                    for (const [targetId, modelIndex] of a.data.targetModels) {
                        const model = await Task.resolveInContext(trajectory.getFrameAtIndex(modelIndex), _ctx);
                        targetMap.set(targetId, { kind: 'structure', structure: Structure.ofModel(model) });
                    }
                } catch {
                    // leave structure targets empty on failure
                }
            } else if (a.data.targetMapping) {
                try {
                    const parentStructure = params.structure.getValue();
                    const map = buildTargetStructuresFromMapping(parentStructure, a.data.targetMapping);
                    for (const [targetId, structure] of map) targetMap.set(targetId, { kind: 'structure', structure });
                } catch {
                    // leave structure targets empty on failure
                }
            } else {
                for (const entry of params.structures) {
                    try {
                        const s = entry.structure.getValue();
                        if (s) targetMap.set(entry.targetId, { kind: 'structure', structure: s });
                    } catch {
                        // ref not resolved yet; skip
                    }
                }
            }

            // Shape targets: explicit per-target refs resolved to concrete shapes.
            for (const entry of params.shapes) {
                try {
                    const provider = entry.shape.getValue();
                    if (!provider) continue;
                    const shape = await provider.getShape(_ctx, provider.data, PD.getDefaultValues(provider.params));
                    targetMap.set(entry.targetId, { kind: 'shape', shape });
                } catch {
                    // ref not resolved yet; skip
                }
            }

            const particles: ParticleList = { ...a.data };
            Particle.setParticleTargets(particles, targetMap);
            return new SO.Particle.List(particles, { label: a.label, description: a.description });
        });
    },
});

type ParticlesRepresentation3D = typeof ParticlesRepresentation3D
const ParticlesRepresentation3D = PluginStateTransform.BuiltIn({
    name: 'particles-representation-3d',
    display: '3D Representation',
    from: SO.Particle.List,
    to: SO.Particle.Representation3D,
    params: (a, ctx: PluginContext) => {
        const { registry, themes: themeCtx } = ctx.representation.particles;
        const type = registry.get(registry.default.name);

        if (!a) {
            return {
                type: PD.Mapped<any>(
                    registry.default.name,
                    registry.types,
                    name => PD.Group<any>(registry.get(name).getParams(themeCtx, undefined as any))),
                colorTheme: PD.Mapped<any>(
                    type.defaultColorTheme.name,
                    themeCtx.colorThemeRegistry.types,
                    name => PD.Group<any>(themeCtx.colorThemeRegistry.get(name).getParams({}))
                ),
                sizeTheme: PD.Mapped<any>(
                    type.defaultSizeTheme.name,
                    themeCtx.sizeThemeRegistry.types,
                    name => PD.Group<any>(themeCtx.sizeThemeRegistry.get(name).getParams({}))
                )
            };
        }

        const dataCtx = { particles: a.data };
        return ({
            type: PD.Mapped<any>(
                registry.default.name,
                registry.getApplicableTypes(a.data),
                name => PD.Group<any>(registry.get(name).getParams(themeCtx, a.data))),
            colorTheme: PD.Mapped<any>(
                type.defaultColorTheme.name,
                themeCtx.colorThemeRegistry.getApplicableTypes(dataCtx),
                name => PD.Group<any>(themeCtx.colorThemeRegistry.get(name).getParams(dataCtx))
            ),
            sizeTheme: PD.Mapped<any>(
                type.defaultSizeTheme.name,
                themeCtx.sizeThemeRegistry.getApplicableTypes(dataCtx),
                name => PD.Group<any>(themeCtx.sizeThemeRegistry.get(name).getParams(dataCtx))
            )
        });
    }
})({
    canAutoUpdate({ oldParams, newParams }) {
        return oldParams.type.name === newParams.type.name;
    },
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Particles Representation', async ctx => {
            const themes = plugin.representation.particles.themes;
            const provider = plugin.representation.particles.registry.get(params.type.name);
            const repr = provider.factory({ webgl: plugin.canvas3d?.webgl, ...themes }, provider.getParams);
            repr.setTheme(Theme.create(themes, { particles: a.data }, params));
            const props = params.type.params || {};
            await repr.createOrUpdate(props, a.data).runInContext(ctx);
            return new SO.Particle.Representation3D({ repr, sourceData: a.data }, { label: provider.label });
        });
    },
    update({ a, b, oldParams, newParams }, plugin: PluginContext) {
        return Task.create('Particles Representation', async ctx => {
            if (newParams.type.name !== oldParams.type.name) return StateTransformer.UpdateResult.Recreate;

            const provider = plugin.representation.particles.registry.get(newParams.type.name);
            if (provider.mustRecreate?.(oldParams.type.params, newParams.type.params)) return StateTransformer.UpdateResult.Recreate;

            const themes = plugin.representation.particles.themes;
            b.data.repr.setTheme(Theme.create(themes, { particles: a.data }, newParams));
            const props = { ...b.data.repr.props, ...newParams.type.params };
            await b.data.repr.createOrUpdate(props, a.data).runInContext(ctx);
            b.data.sourceData = a.data;
            return StateTransformer.UpdateResult.Updated;
        });
    }
});


