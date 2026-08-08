/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginContext } from '../../mol-plugin/context';
import { State, StateAction, StateObjectSelector } from '../../mol-state';
import { Task } from '../../mol-task';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { PluginStateObject } from '../objects';
import { StateTransforms } from '../transforms';
import { ModelFromTrajectory, StructureFromModel, ParticlesStructure } from '../transforms/model';
import { ParticlesVolume } from '../transforms/volume';
import { Asset } from '../../mol-util/assets';
import { getFileNameInfo } from '../../mol-util/file-info';
import { BuiltInTrajectoryFormat, TrajectoryFormatCategory } from '../formats/trajectory';
import { VolumeFormatCategory } from '../formats/volume';
import { ParticlesFormatCategory } from '../formats/particles';
import { createVolumeRepresentationParams } from '../helpers/volume-representation-params';
import { Volume } from '../../mol-model/volume';
import { Structure, Trajectory } from '../../mol-model/structure';
import { PluginConfig } from '../../mol-plugin/config';
import { getSimulariumAgentGeometries, getSimulariumFrameCount, getSimulariumSphereAgentTypes, SimulariumAgentGeometry } from '../../mol-model-formats/particles/simularium';

async function applyParticlesVolumeVisuals(plugin: PluginContext, volume: StateObjectSelector<PluginStateObject.Volume.Data>) {
    const typeParams: { isoValue?: Volume.IsoValue, instanceGranularity: boolean } = { instanceGranularity: true };
    if (volume.data) {
        typeParams.isoValue = Volume.IsoValue.relative(2);
    }
    const update = plugin.build().to(volume).apply(StateTransforms.Representation.VolumeRepresentation3D, createVolumeRepresentationParams(plugin, volume.data, {
        type: 'isosurface',
        typeParams,
    }));
    await update.commit();
}

const TargetDescription = 'The structure or volume to which the particles will be added. The target will be replicated at each particle position.';

export const LoadParticles = StateAction.build({
    display: { name: 'Load Particles', description: 'Load a particle list and a structure or volume from URL or file, replicating them at each particle position.' },
    from: PluginStateObject.Root,
    params(a, ctx: PluginContext) {
        const { options } = ctx.dataFormats;
        const structureOptions = options.filter(o => o[2] === TrajectoryFormatCategory);
        const volumeOptions = options.filter(o => o[2] === VolumeFormatCategory);
        const particlesOptions = options.filter(o => o[2] === ParticlesFormatCategory);

        const structureExts: string[] = [];
        const volumeExts: string[] = [];
        const particlesExts: string[] = [];
        for (const { provider } of ctx.dataFormats.list) {
            if (provider.category === TrajectoryFormatCategory) {
                if (provider.binaryExtensions) structureExts.push(...provider.binaryExtensions);
                if (provider.stringExtensions) structureExts.push(...provider.stringExtensions);
            } else if (provider.category === VolumeFormatCategory) {
                if (provider.binaryExtensions) volumeExts.push(...provider.binaryExtensions);
                if (provider.stringExtensions) volumeExts.push(...provider.stringExtensions);
            } else if (provider.category === ParticlesFormatCategory) {
                if (provider.binaryExtensions) particlesExts.push(...provider.binaryExtensions);
                if (provider.stringExtensions) particlesExts.push(...provider.stringExtensions);
            }
        }

        const targetFormatOptions = [...structureOptions, ...volumeOptions];
        const targetExts = [...structureExts, ...volumeExts];

        return {
            source: PD.MappedStatic('file', {
                url: PD.Group({
                    particles: PD.Group({
                        url: PD.Url('', { label: 'URL' }),
                        format: PD.Select(particlesOptions[0]?.[0] ?? '', particlesOptions),
                        isBinary: PD.Boolean(false),
                    }, { isExpanded: true }),
                    target: PD.Group({
                        url: PD.Url('', { label: 'URL', description: TargetDescription }),
                        format: PD.Select(targetFormatOptions[0]?.[0] ?? '', targetFormatOptions),
                        isBinary: PD.Boolean(false),
                    }, { isExpanded: true }),
                }, { isFlat: true }),
                file: PD.Group({
                    particles: PD.File({ accept: particlesExts.map(e => `.${e}`).join(','), label: 'Particles' }),
                    target: PD.File({ accept: targetExts.map(e => `.${e}`).join(','), label: 'Target', description: TargetDescription }),
                }, { isFlat: true }),
            }, { options: [['url', 'URL'], ['file', 'File']] as ['url' | 'file', string][] }),
        };
    }
})(({ params, state }, ctx: PluginContext) => Task.create('Load Particles', taskCtx => {
    return state.transaction(async () => {
        const s = params.source;

        if (s.name === 'file' && s.params.particles === null) {
            ctx.log.error('No particles file selected');
            return;
        }

        if (s.name === 'file' && s.params.target === null) {
            ctx.log.error('No target file selected');
            return;
        }

        const processUrl = async (url: string | Asset.Url, format: string, isBinary: boolean) => {
            const data = await ctx.builders.data.download({ url, isBinary });
            const provider = ctx.dataFormats.get(format);

            if (!provider) {
                ctx.log.warn(`LoadParticles: could not find data provider for '${format}'`);
                return;
            }

            return provider.parse(ctx, data);
        };

        const processFile = async (file: Asset.File | null) => {
            if (!file) throw new Error('No file selected');

            const info = getFileNameInfo(file.file?.name ?? '');
            const isBinary = ctx.dataFormats.binaryExtensions.has(info.ext);
            const { data } = await ctx.builders.data.readFile({ file, isBinary });
            const provider = ctx.dataFormats.auto(info, data.cell?.obj!);

            if (!provider) {
                ctx.log.warn(`LoadParticles: could not find data provider for '${info.ext}'`);
                await ctx.state.data.build().delete(data).commit();
                return;
            }

            return provider.parse(ctx, data);
        };

        try {
            const particlesParsed = s.name === 'url'
                ? await processUrl(s.params.particles.url, s.params.particles.format, s.params.particles.isBinary)
                : await processFile(s.params.particles);

            if (!particlesParsed || !('list' in particlesParsed)) {
                ctx.log.error('Expected a particles format for the particles input');
                return;
            }

            const targetParsed = s.name === 'url'
                ? await processUrl(s.params.target.url, s.params.target.format, s.params.target.isBinary)
                : await processFile(s.params.target);

            if (!targetParsed) {
                ctx.log.error('Could not parse target');
                return;
            }

            if ('trajectory' in targetParsed) {
                const model = await state.build().to(targetParsed.trajectory)
                    .apply(ModelFromTrajectory, { modelIndex: 0 })
                    .commit();
                const structure = await state.build().to(model)
                    .apply(StructureFromModel)
                    .commit();
                const structureTree = state.build().to(structure.ref)
                    .apply(ParticlesStructure, { particles: PD.Ref(particlesParsed.list.ref) });
                await state.updateTree(structureTree).runInContext(taskCtx);
                const structureProperties = await ctx.builders.structure.insertStructureProperties(structureTree.ref);
                await ctx.builders.structure.representation.applyPreset(structureProperties, 'mesoscale', {
                    theme: { globalName: 'chain-id' }
                });
            } else if ('volume' in targetParsed) {
                const volumeTree = state.build().to(targetParsed.volume.ref)
                    .apply(ParticlesVolume, { particles: PD.Ref(particlesParsed.list.ref) });
                await state.updateTree(volumeTree).runInContext(taskCtx);
                await applyParticlesVolumeVisuals(ctx, volumeTree.selector);
            } else {
                ctx.log.error('Expected a structure or volume format for the target input');
            }
        } catch (e) {
            console.error(e);
            ctx.log.error(`Error loading particles`);
        }
    }).runInContext(taskCtx);
}));

export const LoadMmcifAsParticles = StateAction.build({
    display: { name: 'Load mmCIF as Particles', description: 'Load an mmCIF file, extract assembly operators as a particle list, and display a structure-based particle representation.' },
    from: PluginStateObject.Root,
    params: {
        source: PD.MappedStatic('file', {
            file: PD.Group({
                file: PD.File({ accept: '.cif,.bcif,.mmcif', label: 'mmCIF File' }),
            }, { isFlat: true }),
            url: PD.Group({
                url: PD.Url(''),
                isBinary: PD.Boolean(false),
            }, { isFlat: true }),
        }, { options: [['url', 'URL'], ['file', 'File']] as ['url' | 'file', string][] }),
        assemblyId: PD.Text('1', { description: 'Assembly ID from _pdbx_struct_assembly.id.' }),
    }
})(({ params, state }, ctx: PluginContext) => Task.create('Load mmCIF as Particles', async taskCtx => {
    // 1. Download or read the mmCIF data.
    let data: StateObjectSelector;
    if (params.source.name === 'url') {
        data = await ctx.builders.data.download({ url: params.source.params.url, isBinary: params.source.params.isBinary });
    } else {
        if (!params.source.params.file) {
            ctx.log.error('No mmCIF file selected');
            return;
        }
        const isBinary = (params.source.params.file.file?.name ?? '').toLowerCase().endsWith('.bcif');
        const result = await ctx.builders.data.readFile({ file: params.source.params.file, isBinary });
        data = result.data;
    }

    // 2. Parse CIF (ghost node – not shown in the state tree UI).
    const cif = await state.build()
        .to(data)
        .apply(StateTransforms.Data.ParseCif, undefined, { state: { isGhost: true } })
        .commit({ revertOnError: true });

    // 3. Create particle list from mmCIF assembly.
    const particles = await state.build()
        .to(cif)
        .apply(StateTransforms.Particles.ParticleListFromMmcifAssembly, {
            assemblyId: params.assemblyId,
            asymIds: [],
        })
        .commit({ revertOnError: true });

    // 4. Build a structure from the same CIF data (trajectory → model → structure).
    const trajectory = await state.build()
        .to(cif)
        .apply(StateTransforms.Model.TrajectoryFromMmCif)
        .commit({ revertOnError: true });
    const model = await state.build()
        .to(trajectory)
        .apply(StateTransforms.Model.ModelFromTrajectory, { modelIndex: 0 })
        .commit({ revertOnError: true });
    const structure = await state.build()
        .to(model)
        .apply(StateTransforms.Model.StructureFromModel)
        .commit({ revertOnError: true });

    // 5. Associate reference structures with the particle list.
    //    - targetMapping variants use the model-0 `structure` (chain-split).
    //    - targetModels variants (petworld) use the `trajectory` (one structure per model).
    const decorated = await state.build()
        .to(particles)
        .apply(StateTransforms.Particles.ParticleListWithTargets, {
            trajectory: PD.Ref<Trajectory>(trajectory.ref),
            structure: PD.Ref<Structure>(structure.ref),
            structures: [],
            shapes: [],
        })
        .commit({ revertOnError: true });

    // 6. Show the ParticleTargetRepresentation.
    await state.build()
        .to(decorated)
        .apply(StateTransforms.Particles.ParticlesRepresentation3D, {
            type: { name: 'target', params: {} },
            colorTheme: { name: 'particle-entity', params: {} },
            sizeTheme: { name: 'uniform', params: { value: 1.6 } },
        })
        .commit({ revertOnError: true });
}));

export const AddParticles = StateAction.build({
    display: { name: 'Add Particles', description: 'Add particles to an existing structure or volume, replicating them at each particle position.' },
    from: PluginStateObject.Root,
    params(a, ctx: PluginContext) {
        const state = ctx.state.data;

        const particleLists = state.selectQ(q => q.rootsOfType(PluginStateObject.Particle.List));
        const particlesOptions = particleLists.map(p => [p.transform.ref, p.obj!.label]) as [string, string][];

        const structures = state.selectQ(q => q.rootsOfType(PluginStateObject.Molecule.Structure));
        const volumes = state.selectQ(q => q.rootsOfType(PluginStateObject.Volume.Data));

        const targetOptions = [
            ...structures.map(s => [s.transform.ref, s.obj!.label] as [string, string]),
            ...volumes.map(v => [v.transform.ref, v.obj!.label] as [string, string]),
        ];

        return {
            particles: PD.Select(particlesOptions.length ? particlesOptions[0][0] : '', particlesOptions),
            target: PD.Select(targetOptions.length ? targetOptions[0][0] : '', targetOptions, { description: TargetDescription }),
        };
    }
})(({ params, state }, ctx: PluginContext) => Task.create('Add Particles', async taskCtx => {
    const structures = state.selectQ(q => q.rootsOfType(PluginStateObject.Molecule.Structure));
    const isStructure = structures.some(s => s.transform.ref === params.target);

    if (isStructure) {
        const tree = state.build().to(params.target)
            .apply(ParticlesStructure, { particles: PD.Ref(params.particles) });

        await state.updateTree(tree).runInContext(taskCtx);
        const structureProperties = await ctx.builders.structure.insertStructureProperties(tree.ref);
        await ctx.builders.structure.representation.applyPreset(structureProperties, 'mesoscale', {
            theme: { globalName: 'chain-id' }
        });
    } else {
        const tree = state.build().to(params.target)
            .apply(ParticlesVolume, { particles: PD.Ref(params.particles) });

        await state.updateTree(tree).runInContext(taskCtx);
        await applyParticlesVolumeVisuals(ctx, tree.selector);
    }
}));

type SimulariumPdbSource = { url: string, isBinary: boolean, format: BuiltInTrajectoryFormat };

/**
 * Resolve a Simularium PDB geometry `url` into download parameters, mirroring the
 * `DownloadStructure` action: a value matching a 4-character PDB id is fetched from
 * the default PDB provider, otherwise the value is treated as a direct URL and the
 * format is inferred from its extension.
 */
function resolveSimulariumPdbSource(ctx: PluginContext, url: string): SimulariumPdbSource {
    const id = url.trim();
    if (/^[1-9][a-z0-9]{3}$/i.test(id)) {
        const provider = ctx.config.get(PluginConfig.Download.DefaultPdbProvider) || 'pdbe';
        const downloadUrl = provider === 'rcsb'
            ? `https://models.rcsb.org/${id.toUpperCase()}.bcif`
            : `https://www.ebi.ac.uk/pdbe/entry-files/download/${id.toLowerCase()}.bcif`;
        return { url: downloadUrl, isBinary: true, format: 'mmcif' };
    }
    switch (getFileNameInfo(id).ext) {
        case 'bcif': return { url: id, isBinary: true, format: 'mmcif' };
        case 'cif': case 'mmcif': return { url: id, isBinary: false, format: 'mmcif' };
        case 'pdb': case 'ent': return { url: id, isBinary: false, format: 'pdb' };
        default: return { url: id, isBinary: false, format: 'mmcif' };
    }
}

/** Load the per-agent-type PDB structures and OBJ meshes and instance them at each particle. */
async function loadSimulariumGeometries(ctx: PluginContext, state: State, simulariumRef: string, frameIndex: number, scale: number, pdbGeometries: SimulariumAgentGeometry[], objGeometries: SimulariumAgentGeometry[], sphereAgentTypes: { id: number, name: string }[], sphereFiles: Asset.File[]) {
    const geometries = [...pdbGeometries, ...objGeometries];

    // Resolve which sphere agent types actually have a matching file supplied.
    const matchedSphereTypes = sphereAgentTypes.filter(t =>
        sphereFiles.some(f => getFileNameInfo(f.file?.name ?? '').base === t.name)
    );

    // A single particle list holding the particles of every geometry agent type. Each
    // particle's target id is its Simularium type id, so the per-target structures and
    // shapes below are instanced on the matching particles.
    const trajNode = await state.build().to(simulariumRef)
        .apply(StateTransforms.Particles.ParticleTrajectoryFromSimularium, {
            scale,
            types: [...geometries, ...matchedSphereTypes].map(g => String(g.id)),
        })
        .commit({ revertOnError: true });

    const list = await state.build().to(trajNode.ref)
        .apply(StateTransforms.Particles.ParticleListFromTrajectory, {
            frameIndex,
        })
        .commit({ revertOnError: true });

    // Group node to collect all reference geometries, collapsed by default.
    const group = await state.build().toRoot()
        .apply(StateTransforms.Misc.CreateGroup, { label: 'Simularium Geometries' }, { state: { isCollapsed: true } })
        .commit({ revertOnError: true });

    // PDB structures (from a PDB id or direct URL, like DownloadStructure), mapped to
    // their agent type id (= particle target id).
    const structures: { targetId: number, structure: PD.Ref<Structure> }[] = [];
    for (const geometry of pdbGeometries) {
        const src = resolveSimulariumPdbSource(ctx, geometry.url);
        const data = await state.build().to(group.ref)
            .apply(StateTransforms.Data.Download, { url: Asset.Url(src.url), isBinary: src.isBinary, label: geometry.name }, { state: { isGhost: true } })
            .commit({ revertOnError: true });
        const trajectory = await ctx.builders.structure.parseTrajectory(data, src.format);
        const model = await ctx.builders.structure.createModel(trajectory);
        const structure = await ctx.builders.structure.createStructure(model, { name: 'model', params: {} });
        structures.push({ targetId: geometry.id, structure: PD.Ref<Structure>(structure.ref) });
    }

    // OBJ meshes loaded as shapes, mapped to their agent type id. The shape target scales
    // each mesh instance by the per-particle radius (see ParticleTargetRepresentation).
    const shapes: { targetId: number, shape: PD.Ref<any> }[] = [];
    for (const geometry of objGeometries) {
        const obj = await state.build().to(group.ref)
            .apply(StateTransforms.Data.Download, { url: Asset.Url(geometry.url), isBinary: false, label: geometry.name }, { state: { isGhost: true } })
            .apply(StateTransforms.Data.ParseObj, undefined, { state: { isGhost: true } })
            .commit({ revertOnError: true });
        const shape = await state.build().to(obj.ref)
            .apply(StateTransforms.Shape.ShapeFromObj, { label: geometry.name })
            .commit({ revertOnError: true });
        shapes.push({ targetId: geometry.id, shape: PD.Ref<any>(shape.ref) });
    }

    // SPHERE structures from user-supplied files, matched by file base name to agent type name.
    for (const file of sphereFiles) {
        const info = getFileNameInfo(file.file?.name ?? '');
        const match = sphereAgentTypes.find(t => t.name === info.base);
        if (!match) {
            ctx.log.warn(`LoadSimulariumGeometries: no SPHERE agent type matches file '${info.name}' (expected base name to match a typeMapping entry).`);
            continue;
        }
        const isBinary = ctx.dataFormats.binaryExtensions.has(info.ext);
        const data = await state.build().toRoot()
            .apply(StateTransforms.Data.ReadFile, { file, isBinary, label: info.name }, { state: { isGhost: true } })
            .commit({ revertOnError: true });
        const provider = ctx.dataFormats.auto(info, data.cell?.obj!);
        if (!provider) {
            ctx.log.warn(`LoadSimulariumGeometries: could not find a data provider for '${info.ext}'.`);
            continue;
        }
        const parsed = await provider.parse(ctx, data);
        if (!parsed || !('trajectory' in parsed)) {
            ctx.log.warn(`LoadSimulariumGeometries: expected a structure format for '${info.name}'.`);
            continue;
        }
        const model = await ctx.builders.structure.createModel(parsed.trajectory);
        const structure = await ctx.builders.structure.createStructure(model, { name: 'model', params: {} });
        structures.push({ targetId: match.id, structure: PD.Ref<Structure>(structure.ref) });
    }

    // Associate the reference structures and shapes with the particle list via the target
    // mapping, then render them instanced at each particle position.
    const decorated = await state.build().to(list.ref)
        .apply(StateTransforms.Particles.ParticleListWithTargets, {
            trajectory: PD.Ref<Trajectory>(''),
            structure: PD.Ref<Structure>(''),
            structures,
            shapes,
        })
        .commit({ revertOnError: true });
    await state.build().to(decorated.ref)
        .apply(StateTransforms.Particles.ParticlesRepresentation3D, {
            type: { name: 'target', params: {} },
            colorTheme: { name: 'particle-entity', params: {} },
            sizeTheme: { name: 'uniform', params: { value: 1.6 } },
        })
        .commit({ revertOnError: true });
}

export const LoadSimulariumGeometries = StateAction.build({
    display: { name: 'Load Simularium Geometries', description: 'Load the PDB structures and OBJ meshes referenced by agent types in a Simularium file and replicate them at each particle position. Optionally supply structure files for SPHERE-type agents, matched by file base name to the agent type name.' },
    from: PluginStateObject.Format.Simularium,
    params(a, ctx: PluginContext) {
        const frameCount = a ? getSimulariumFrameCount(a.data) : 1;
        const sphereTypes = a ? getSimulariumSphereAgentTypes(a.data) : [];

        const structureExts: string[] = [];
        if (ctx) {
            for (const { provider } of ctx.dataFormats.list) {
                if (provider.category === TrajectoryFormatCategory) {
                    if (provider.binaryExtensions) structureExts.push(...provider.binaryExtensions);
                    if (provider.stringExtensions) structureExts.push(...provider.stringExtensions);
                }
            }
        }
        const accept = structureExts.length > 0
            ? structureExts.map(e => `.${e}`).join(',')
            : '.cif,.bcif,.pdb,.ent,.mmcif';

        const sphereDescription = sphereTypes.length > 0
            ? `One file per SPHERE agent type; the file base name must match the agent type name. Known SPHERE types: ${sphereTypes.map(t => t.name).join(', ')}.`
            : 'One file per SPHERE agent type; the file base name must match the agent type name.';

        return {
            frameIndex: PD.Numeric(0, { min: 0, max: Math.max(0, frameCount - 1), step: 1 }, { description: 'Index of the trajectory frame to display.' }),
            scale: PD.Numeric(0, { min: 0, step: 0.001 }, { description: 'Spatial scale to angstrom. Leave 0 to auto-detect from the file spatial units.' }),
            sphereFiles: PD.FileList({ accept, label: 'SPHERE Structure Files', description: sphereDescription }),
        };
    }
})(({ a, ref, params, state }, ctx: PluginContext) => Task.create('Load Simularium Geometries', async () => {
    const file = a.data;
    const geometries = getSimulariumAgentGeometries(file);
    const sphereAgentTypes = getSimulariumSphereAgentTypes(file);
    const sphereFiles = params.sphereFiles ?? [];

    if (geometries.length === 0 && sphereFiles.length === 0) {
        ctx.log.warn('Simularium file has no PDB or OBJ agent geometries and no SPHERE structure files were provided.');
        return;
    }

    const pdbGeometries = geometries.filter(g => g.displayType === 'PDB');
    const objGeometries = geometries.filter(g => g.displayType === 'OBJ');

    // All PDB structures, OBJ shapes, and SPHERE structures share a single particle list
    // and one representation.
    try {
        await loadSimulariumGeometries(ctx, state, ref, params.frameIndex, params.scale, pdbGeometries, objGeometries, sphereAgentTypes, sphereFiles);
    } catch (e) {
        console.error(e);
        ctx.log.warn('Could not load the Simularium geometries.');
    }
}));
