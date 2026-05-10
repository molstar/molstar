/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginContext } from '../../mol-plugin/context';
import { StateAction, StateTransformer, StateSelection } from '../../mol-state';
import { Task } from '../../mol-task';
import { getFileNameInfo } from '../../mol-util/file-info';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { PluginStateObject } from '../objects';
import { Download } from '../transforms/data';
import { DataFormatProvider } from '../formats/provider';
import { Asset } from '../../mol-util/assets';
import { StateTransforms } from '../transforms';
import { assertUnreachable } from '../../mol-util/type-helpers';
import { VolumeFormatCategory } from '../formats/volume';
import { ParticlesFormatCategory } from '../formats/particles';
import { VolumeFromVolumeAndParticles } from '../transforms/volume';
import { createVolumeRepresentationParams } from '../helpers/volume-representation-params';
import { Volume } from '../../mol-model/volume';
import { StateObjectSelector } from '../../mol-state';

export type EmdbDownloadProvider = 'pdbe' | 'rcsb'

export { DownloadDensity };
type DownloadDensity = typeof DownloadDensity
const DownloadDensity = StateAction.build({
    from: PluginStateObject.Root,
    display: { name: 'Download Density', description: 'Load a density from the provided source and create its default visual.' },
    params: (a, ctx: PluginContext) => {
        const { options } = ctx.dataFormats;
        return {
            source: PD.MappedStatic('pdb-xray', {
                'pdb-xray': PD.Group({
                    provider: PD.Group({
                        id: PD.Text('1tqn', { label: 'Id' }),
                        server: PD.Select('pdbe', [['pdbe', 'PDBe']]),
                    }, { pivot: 'id' }),
                    type: PD.Select('2fofc', [['2fofc', '2Fo-Fc'], ['fofc', 'Fo-Fc']]),
                }, { isFlat: true }),
                'pdb-xray-ds': PD.Group({
                    provider: PD.Group({
                        id: PD.Text('1tqn', { label: 'Id' }),
                        server: PD.Select('pdbe', [['pdbe', 'PDBe'], ['rcsb', 'RCSB PDB']]),
                    }, { pivot: 'id' }),
                    detail: PD.Numeric(3, { min: 0, max: 6, step: 1 }, { label: 'Detail' }),
                }, { isFlat: true }),
                'pdb-emd-ds': PD.Group({
                    provider: PD.Group({
                        id: PD.Text('emd-8004', { label: 'Id' }),
                        server: PD.Select<EmdbDownloadProvider>('pdbe', [['pdbe', 'PDBe'], ['rcsb', 'RCSB PDB']]),
                    }, { pivot: 'id' }),
                    detail: PD.Numeric(3, { min: 0, max: 6, step: 1 }, { label: 'Detail' }),
                }, { isFlat: true }),
                'url': PD.Group({
                    url: PD.Url(''),
                    isBinary: PD.Boolean(false),
                    format: PD.Select('auto', options),
                }, { isFlat: true })
            }, {
                options: [
                    ['pdb-xray', 'PDB X-ray maps'],
                    ['pdb-emd-ds', 'PDB EMD Density Server'],
                    ['pdb-xray-ds', 'PDB X-ray Density Server'],
                    ['url', 'URL']
                ]
            })
        };
    }
})(({ params }, plugin: PluginContext) => Task.create('Download Density', async taskCtx => {
    const src = params.source;
    let downloadParams: StateTransformer.Params<Download>;
    let provider: DataFormatProvider | undefined;

    switch (src.name) {
        case 'url':
            downloadParams = src.params;
            break;
        case 'pdb-xray':
            downloadParams = {
                url: Asset.Url(src.params.type === '2fofc'
                    ? `https://www.ebi.ac.uk/pdbe/coordinates/files/${src.params.provider.id.toLowerCase()}.ccp4`
                    : `https://www.ebi.ac.uk/pdbe/coordinates/files/${src.params.provider.id.toLowerCase()}_diff.ccp4`),
                isBinary: true,
                label: `PDBe X-ray map: ${src.params.provider.id}`
            };
            break;
        case 'pdb-emd-ds':
            downloadParams = src.params.provider.server === 'pdbe' ? {
                url: Asset.Url(`https://www.ebi.ac.uk/pdbe/densities/emd/${src.params.provider.id.toLowerCase()}/cell?detail=${src.params.detail}`),
                isBinary: true,
                label: `PDBe EMD Density Server: ${src.params.provider.id}`
            } : {
                url: Asset.Url(`https://maps.rcsb.org/em/${src.params.provider.id.toLowerCase()}/cell?detail=${src.params.detail}`),
                isBinary: true,
                label: `RCSB PDB EMD Density Server: ${src.params.provider.id}`
            };
            break;
        case 'pdb-xray-ds':
            downloadParams = src.params.provider.server === 'pdbe' ? {
                url: Asset.Url(`https://www.ebi.ac.uk/pdbe/densities/x-ray/${src.params.provider.id.toLowerCase()}/cell?detail=${src.params.detail}`),
                isBinary: true,
                label: `PDBe X-ray Density Server: ${src.params.provider.id}`
            } : {
                url: Asset.Url(`https://maps.rcsb.org/x-ray/${src.params.provider.id.toLowerCase()}/cell?detail=${src.params.detail}`),
                isBinary: true,
                label: `RCSB PDB X-ray Density Server: ${src.params.provider.id}`
            };
            break;
        default: assertUnreachable(src);
    }

    const data = await plugin.builders.data.download(downloadParams);
    let entryId: string | undefined = undefined;

    switch (src.name) {
        case 'url':
            downloadParams = src.params;
            provider = src.params.format === 'auto' ? plugin.dataFormats.auto(getFileNameInfo(Asset.getUrl(downloadParams.url)), data.cell?.obj!) : plugin.dataFormats.get(src.params.format);
            break;
        case 'pdb-xray':
            entryId = src.params.provider.id;
            provider = plugin.dataFormats.get('ccp4');
            break;
        case 'pdb-emd-ds':
        case 'pdb-xray-ds':
            entryId = src.params.provider.id;
            provider = plugin.dataFormats.get('dscif');
            break;
        default: assertUnreachable(src);
    }

    if (!provider) {
        plugin.log.warn('DownloadDensity: Format provider not found.');
        return;
    }

    const volumes = await provider.parse(plugin, data, { entryId });
    await provider.visuals?.(plugin, volumes);
}));

export const AssignColorVolume = StateAction.build({
    display: { name: 'Assign Volume Colors', description: 'Assigns another volume to be available for coloring.' },
    from: PluginStateObject.Volume.Data,
    isApplicable(a) { return !a.data.colorVolume; },
    params(a, plugin: PluginContext) {
        const cells = plugin.state.data.select(StateSelection.Generators.root.subtree().ofType(PluginStateObject.Volume.Data).filter(cell => !!cell.obj && !cell.obj?.data.colorVolume && cell.obj !== a));
        if (cells.length === 0) return { ref: PD.Text('', { isHidden: true, label: 'Volume' }) };
        return { ref: PD.Select(cells[0].transform.ref, cells.map(c => [c.transform.ref, c.obj!.label]), { label: 'Volume' }) };
    }
})(({ ref, params, state }, plugin: PluginContext) => {
    return plugin.build().to(ref).apply(StateTransforms.Volume.AssignColorVolume, { ref: params.ref }, { dependsOn: [params.ref] }).commit();
});

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

export const AddParticlesVolume = StateAction.build({
    display: { name: 'Add Particles Volume', description: 'Replicate an existing volume at the positions/orientations of an existing particle list.' },
    from: PluginStateObject.Root,
    params(a, ctx: PluginContext) {
        const state = ctx.state.data;
        const volumes = state.selectQ(q => q.rootsOfType(PluginStateObject.Volume.Data));
        const volumeOptions = volumes.map(v => [v.transform.ref, v.obj!.label]) as [string, string][];
        const particles = state.selectQ(q => q.rootsOfType(PluginStateObject.Particle.List));
        const particleOptions = particles.map(p => [p.transform.ref, p.obj!.label]) as [string, string][];
        return {
            volume: PD.Select(volumeOptions.length ? volumeOptions[0][0] : '', volumeOptions),
            particles: PD.Select(particleOptions.length ? particleOptions[0][0] : '', particleOptions),
        };
    }
})(({ params, state }, ctx: PluginContext) => Task.create('Add Particles Volume', taskCtx => {
    return state.transaction(async () => {
        const dependsOn = [params.volume, params.particles];
        const tree = state.build().toRoot()
            .apply(VolumeFromVolumeAndParticles, {
                volumeRef: params.volume,
                particlesRef: params.particles,
            }, { dependsOn });

        await state.updateTree(tree).runInContext(taskCtx);
        await applyParticlesVolumeVisuals(ctx, tree.selector);
    }).runInContext(taskCtx);
}));

export const LoadParticlesVolume = StateAction.build({
    display: { name: 'Load Particles Volume', description: 'Load a volume and a particle list from URL or file and replicate the volume at each particle.' },
    from: PluginStateObject.Root,
    params(a, ctx: PluginContext) {
        const { options } = ctx.dataFormats;
        const volumeOptions = options.filter(o => o[2] === VolumeFormatCategory);
        const particlesOptions = options.filter(o => o[2] === ParticlesFormatCategory);

        const volumeExts: string[] = [];
        const particlesExts: string[] = [];
        for (const { provider } of ctx.dataFormats.list) {
            if (provider.category === VolumeFormatCategory) {
                if (provider.binaryExtensions) volumeExts.push(...provider.binaryExtensions);
                if (provider.stringExtensions) volumeExts.push(...provider.stringExtensions);
            } else if (provider.category === ParticlesFormatCategory) {
                if (provider.binaryExtensions) particlesExts.push(...provider.binaryExtensions);
                if (provider.stringExtensions) particlesExts.push(...provider.stringExtensions);
            }
        }

        return {
            source: PD.MappedStatic('file', {
                url: PD.Group({
                    volume: PD.Group({
                        url: PD.Url(''),
                        format: PD.Select(volumeOptions[0]?.[0] ?? '', volumeOptions),
                        isBinary: PD.Boolean(true),
                    }, { isExpanded: true }),
                    particles: PD.Group({
                        url: PD.Url(''),
                        format: PD.Select(particlesOptions[0]?.[0] ?? '', particlesOptions),
                        isBinary: PD.Boolean(false),
                    }, { isExpanded: true })
                }, { isFlat: true }),
                file: PD.Group({
                    volume: PD.File({ accept: volumeExts.map(e => `.${e}`).join(','), label: 'Volume' }),
                    particles: PD.File({ accept: particlesExts.map(e => `.${e}`).join(','), label: 'Particles' }),
                }, { isFlat: true }),
            }, { options: [['url', 'URL'], ['file', 'File']] })
        };
    }
})(({ params, state }, ctx: PluginContext) => Task.create('Load Particles Volume', taskCtx => {
    return state.transaction(async () => {
        const s = params.source;

        if (s.name === 'file' && (s.params.volume === null || s.params.particles === null)) {
            ctx.log.error('No file(s) selected');
            return;
        }

        if (s.name === 'url' && (!s.params.volume || !s.params.particles)) {
            ctx.log.error('No URL(s) given');
            return;
        }

        const processUrl = async (url: string | Asset.Url, format: string, isBinary: boolean) => {
            const data = await ctx.builders.data.download({ url, isBinary });
            const provider = ctx.dataFormats.get(format);

            if (!provider) {
                ctx.log.warn(`LoadParticlesVolume: could not find data provider for '${format}'`);
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
                ctx.log.warn(`LoadParticlesVolume: could not find data provider for '${info.ext}'`);
                await ctx.state.data.build().delete(data).commit();
                return;
            }

            return provider.parse(ctx, data);
        };

        try {
            const volumeParsed = s.name === 'url'
                ? await processUrl(s.params.volume.url, s.params.volume.format, s.params.volume.isBinary)
                : await processFile(s.params.volume);

            if (!volumeParsed || !('volume' in volumeParsed)) {
                ctx.log.error('Expected a volume format for the volume input');
                return;
            }

            //

            const particlesParsed = s.name === 'url'
                ? await processUrl(s.params.particles.url, s.params.particles.format, s.params.particles.isBinary)
                : await processFile(s.params.particles);

            if (!particlesParsed || !('list' in particlesParsed)) {
                ctx.log.error('Expected a particles format for the particles input');
                return;
            }

            //

            const dependsOn = [volumeParsed.volume.ref, particlesParsed.list.ref];
            const tree = state.build().toRoot()
                .apply(VolumeFromVolumeAndParticles, {
                    volumeRef: volumeParsed.volume.ref,
                    particlesRef: particlesParsed.list.ref,
                }, { dependsOn });

            await state.updateTree(tree).runInContext(taskCtx);
            await applyParticlesVolumeVisuals(ctx, tree.selector);
        } catch (e) {
            console.error(e);
            ctx.log.error(`Error loading particles volume`);
        }
    }).runInContext(taskCtx);
}));
