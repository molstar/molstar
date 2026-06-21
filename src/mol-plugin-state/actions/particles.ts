/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginContext } from '../../mol-plugin/context';
import { StateAction, StateObjectSelector } from '../../mol-state';
import { Task } from '../../mol-task';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { PluginStateObject } from '../objects';
import { StateTransforms } from '../transforms';
import { ModelFromTrajectory, StructureFromModel, ParticlesStructure } from '../transforms/model';
import { ParticlesVolume } from '../transforms/volume';
import { Asset } from '../../mol-util/assets';
import { getFileNameInfo } from '../../mol-util/file-info';
import { TrajectoryFormatCategory } from '../formats/trajectory';
import { VolumeFormatCategory } from '../formats/volume';
import { ParticlesFormatCategory } from '../formats/particles';
import { createVolumeRepresentationParams } from '../helpers/volume-representation-params';
import { Volume } from '../../mol-model/volume';
import { Structure, Trajectory } from '../../mol-model/structure';

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
        .apply(StateTransforms.Particles.ParticleListWithStructures, {
            trajectory: PD.Ref<Trajectory>(trajectory.ref),
            structure: PD.Ref<Structure>(structure.ref),
            structures: [],
        })
        .commit({ revertOnError: true });

    // 6. Show the ParticlesStructureRepresentation.
    await state.build()
        .to(decorated)
        .apply(StateTransforms.Particles.ParticlesRepresentation3D, {
            type: { name: 'particles-structure', params: {} },
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
