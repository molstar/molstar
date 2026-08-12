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
import { createParticleListFromMmcifAssembly, getAssemblyIdsFromMmcif, getAsymIdsFromMmcif, MmcifVariant } from '../../mol-model-formats/particles/mmcif';
import { createSimulariumParticleTrajectory, getSimulariumAgentTypeNames, getSimulariumFrameCount } from '../../mol-model-formats/particles/simularium';
import { PluginContext } from '../../mol-plugin/context';
import { StateTransformer } from '../../mol-state';
import { Task } from '../../mol-task';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Theme } from '../../mol-theme/theme';
import { PluginStateObject as SO, PluginStateTransform } from '../objects';
import { Particle } from '../../mol-model/particles/particle-list';
import { StateObject } from '../../mol-state/object';
import { getUnitcellDataFromSymmetry, UnitcellParams, UnitcellRepresentation } from '../../mol-repr/shape/model/unitcell';
import { Cell } from '../../mol-math/geometry/spacegroup/cell';

export { ParticleListFromRelionStar };
export { ParticleListFromDynamoTbl };
export { ParticleListFromCryoEtDataPortalNdjson };
export { ParticleListFromArtiatomiEm };
export { ParticleListFromMmcifAssembly };
export { ParticleTrajectoryFromSimularium };
export { ParticleListFromTrajectory };
export { ParticlesRepresentation3D };
export { ParticleListUnitcell3D };

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

type ParticleListFromMmcifAssembly = typeof ParticleListFromMmcifAssembly
const ParticleListFromMmcifAssembly = PluginStateTransform.BuiltIn({
    name: 'particle-list-from-mmcif-assembly',
    display: { name: 'Particle List from mmCIF Assembly', description: 'Create ParticleList from a mmCIF (CellPack/PetWorld) assembly.' },
    from: SO.Format.Cif,
    to: SO.Particle.List,
    params: a => {
        const variant = PD.Select<MmcifVariant>('auto', [
            ['auto', 'Auto'],
            ['cellpack', 'CellPack'],
            ['standard', 'Standard'],
            ['petworld', 'PetWorld'],
        ], { description: 'mmCIF variant used to interpret entity names/compartments. Auto detects from the file header and categories.' });
        const label = PD.Optional(PD.Text(''));

        if (!a) {
            return {
                assemblyId: PD.Text('', { description: 'Assembly identifier from _pdbx_struct_assembly.id.' }),
                asymIds: PD.MultiSelect<string>([], [], { description: 'Empty selection includes all chains for the assembly.' }),
                variant,
                label,
            };
        }

        let assemblyIds: string[] = [];
        try {
            assemblyIds = getAssemblyIdsFromMmcif(a.data);
        } catch {
            // ignore; apply will surface parse errors
        }
        const assemblyOptions = assemblyIds.map(id => [id, id] as [string, string]);
        const defaultAssemblyId = assemblyIds.length > 0 ? assemblyIds[0] : '';

        let asymIds: string[] = [];
        try {
            asymIds = defaultAssemblyId ? getAsymIdsFromMmcif(a.data, defaultAssemblyId) : [];
        } catch {
            // ignore; apply will surface parse errors
        }
        const asymOptions = asymIds.map(id => [id, id] as [string, string]);

        return {
            assemblyId: PD.Select(defaultAssemblyId, assemblyOptions, { description: 'Assembly identifier from _pdbx_struct_assembly.id.' }),
            asymIds: PD.MultiSelect<string>([], asymOptions, { description: 'Empty selection includes all chains for the assembly.' }),
            variant,
            label,
        };
    }
})({
    apply({ a, params }) {
        return Task.create('Create Particle List from mmCIF Assembly', async ctx => {
            // Params default to the first assembly ID only when interactively edited; fall back
            // here too so auto-loading a file (which calls parse() with no params) still works.
            const assemblyId = params.assemblyId || getAssemblyIdsFromMmcif(a.data)[0];
            if (!assemblyId) {
                throw new Error('mmCIF file contains no _pdbx_struct_assembly_gen assemblies; cannot create a particle list.');
            }
            const list = await createParticleListFromMmcifAssembly(a.data, {
                assemblyId,
                asymIds: params.asymIds.length > 0 ? params.asymIds : void 0,
                variant: params.variant,
                label: params.label || a.label,
            }, ctx);
            return new SO.Particle.List(list, { label: list.label || 'Particles', description: 'mmCIF Particle List' });
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

const plus1 = (v: number) => v + 1, minus1 = (v: number) => v - 1;
type ParticleListFromTrajectory = typeof ParticleListFromTrajectory
const ParticleListFromTrajectory = PluginStateTransform.BuiltIn({
    name: 'particle-list-from-trajectory',
    display: { name: 'Particle List from Trajectory', description: 'Extract a single frame from a ParticleTrajectory.' },
    from: SO.Particle.Trajectory,
    to: SO.Particle.List,
    params: a => {
        if (!a) {
            return { frameIndex: PD.Numeric(0, {}, { description: 'Zero-based index of the frame', immediateUpdate: true }) };
        }
        return { frameIndex: PD.Converted(plus1, minus1, PD.Numeric(1, { min: 1, max: a.data.frameCount, step: 1 }, { description: 'Frame Index', immediateUpdate: true })) };
    }
})({
    isApplicable: a => a.data.frameCount > 0,
    apply({ a, params }) {
        return Task.create('Particle List from Trajectory', async ctx => {
            const idx = Math.max(0, Math.min(params.frameIndex, a.data.frameCount - 1));
            const list = await Task.resolveInContext(a.data.getFrameAtIndex(idx), ctx);
            return new SO.Particle.List(list, { label: list.label || 'Particles', description: `Frame ${params.frameIndex + 1} of ${a.data.frameCount}` });
        });
    }
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
            const colorThemeInfo = {
                help: (value: { name: string, params: {} }) => {
                    const { name, params } = value;
                    const p = themeCtx.colorThemeRegistry.get(name);
                    const ct = p.factory({}, params);
                    return { description: ct.description, legend: ct.legend };
                }
            };

            return {
                type: PD.Mapped<any>(
                    registry.default.name,
                    registry.types,
                    name => PD.Group<any>(registry.get(name).getParams(themeCtx, undefined as any))),
                colorTheme: PD.Mapped<any>(
                    type.defaultColorTheme.name,
                    themeCtx.colorThemeRegistry.types,
                    name => PD.Group<any>(themeCtx.colorThemeRegistry.get(name).getParams({})),
                    colorThemeInfo
                ),
                sizeTheme: PD.Mapped<any>(
                    type.defaultSizeTheme.name,
                    themeCtx.sizeThemeRegistry.types,
                    name => PD.Group<any>(themeCtx.sizeThemeRegistry.get(name).getParams({}))
                )
            };
        }

        const dataCtx = { particles: a.data };
        const colorThemeInfo = {
            help: (value: { name: string, params: {} }) => {
                const { name, params } = value;
                const p = themeCtx.colorThemeRegistry.get(name);
                const ct = p.factory(dataCtx, params);
                return { description: ct.description, legend: ct.legend };
            }
        };

        return ({
            type: PD.Mapped<any>(
                registry.default.name,
                registry.getApplicableTypes(a.data),
                name => PD.Group<any>(registry.get(name).getParams(themeCtx, a.data))),
            colorTheme: PD.Mapped<any>(
                type.defaultColorTheme.name,
                themeCtx.colorThemeRegistry.getApplicableTypes(dataCtx),
                name => PD.Group<any>(themeCtx.colorThemeRegistry.get(name).getParams(dataCtx)),
                colorThemeInfo
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

type ParticleListUnitcell3D = typeof ParticleListUnitcell3D
const ParticleListUnitcell3D = PluginStateTransform.BuiltIn({
    name: 'particle-list-unitcell-3d',
    display: 'Particle List Unit Cell',
    from: SO.Particle.List,
    to: SO.Shape.Representation3D,
    params: () => ({
        ...UnitcellParams,
    })
})({
    isApplicable: a => !!a.data.cell,
    canAutoUpdate({ oldParams, newParams }) {
        return true;
    },
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Particle List Unit Cell', async ctx => {
            const { cell } = a.data;
            if (!cell) return StateObject.Null;

            const center = Particle.getBoundary(a.data).sphere.center;
            const data = getUnitcellDataFromSymmetry(cell, center, params);
            const repr = UnitcellRepresentation({ webgl: plugin.canvas3d?.webgl, ...plugin.representation.structure.themes }, () => UnitcellParams);
            await repr.createOrUpdate(params, data).runInContext(ctx);
            return new SO.Shape.Representation3D({ repr, sourceData: data }, { label: 'Unit Cell', description: Cell.getLabel(cell) });
        });
    },
    update({ a, b, newParams }) {
        return Task.create('Particle List Unit Cell', async ctx => {
            const { cell } = a.data;
            if (!cell) return StateTransformer.UpdateResult.Null;

            const props = { ...b.data.repr.props, ...newParams };
            const center = Particle.getBoundary(a.data).sphere.center;
            const data = getUnitcellDataFromSymmetry(cell, center, props);
            await b.data.repr.createOrUpdate(props, data).runInContext(ctx);
            b.data.sourceData = data;
            b.description = Cell.getLabel(cell);
            return StateTransformer.UpdateResult.Updated;
        });
    }
});
