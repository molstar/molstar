/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginContext } from '../../mol-plugin/context';
import { StateAction, StateSelection, StateTransformer } from '../../mol-state';
import { Task } from '../../mol-task';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { PresetStructureRepresentations, StructureRepresentationPresetProvider } from '../builder/structure/representation-preset';
import { BuiltInTrajectoryFormat, BuiltInTrajectoryFormats, TrajectoryFormatCategory } from '../formats/trajectory';
import { RootStructureDefinition } from '../helpers/root-structure';
import { PluginStateObject } from '../objects';
import { StateTransforms } from '../transforms';
import { Download } from '../transforms/data';
import { CustomModelProperties, CustomStructureProperties, ModelFromTrajectory, TrajectoryFromModelAndCoordinates } from '../transforms/model';
import { Asset } from '../../mol-util/assets';
import { PluginConfig } from '../../mol-plugin/config';
import { getFileNameInfo } from '../../mol-util/file-info';
import { assertUnreachable } from '../../mol-util/type-helpers';
import { TopologyFormatCategory } from '../formats/topology';
import { CoordinatesFormatCategory } from '../formats/coordinates';

const DownloadModelRepresentationOptions = (plugin: PluginContext) => {
    const representationDefault = plugin.config.get(PluginConfig.Structure.DefaultRepresentationPreset) || PresetStructureRepresentations.auto.id;
    return PD.Group({
        type: RootStructureDefinition.getParams(void 0, 'auto').type,
        representation: PD.Select(representationDefault,
            plugin.builders.structure.representation.getPresets().map(p => [p.id, p.display.name, p.display.group] as any),
            { description: 'Which representation preset to use.' }),
        representationParams: PD.Group(StructureRepresentationPresetProvider.CommonParams, { isHidden: true }),
        asTrajectory: PD.Optional(PD.Boolean(false, { description: 'Load all entries into a single trajectory.' }))
    }, { isExpanded: false });
};

export const PdbDownloadProvider = {
    'rcsb': PD.Group({
        encoding: PD.Select('bcif', PD.arrayToOptions(['cif', 'bcif'] as const)),
    }, { label: 'RCSB PDB', isFlat: true }),
    'pdbe': PD.Group({
        variant: PD.Select('updated-bcif', [['updated-bcif', 'Updated (bcif)'], ['updated', 'Updated'], ['archival', 'Archival']] as ['updated' | 'updated-bcif' | 'archival', string][]),
    }, { label: 'PDBe', isFlat: true }),
    'pdbj': PD.EmptyGroup({ label: 'PDBj' }),
};
export type PdbDownloadProvider = keyof typeof PdbDownloadProvider;

export { DownloadStructure };
type DownloadStructure = typeof DownloadStructure
const DownloadStructure = StateAction.build({
    from: PluginStateObject.Root,
    display: { name: 'Download Structure', description: 'Load a structure from the provided source and create its representation.' },
    params: (_, plugin: PluginContext) => {
        const options = DownloadModelRepresentationOptions(plugin);
        const defaultPdbProvider = plugin.config.get(PluginConfig.Download.DefaultPdbProvider) || 'pdbe';
        return {
            source: PD.MappedStatic('pdb', {
                'pdb': PD.Group({
                    provider: PD.Group({
                        id: PD.Text('1tqn', { label: 'PDB Id(s)', description: 'One or more comma/space separated PDB ids.' }),
                        server: PD.MappedStatic(defaultPdbProvider, PdbDownloadProvider),
                    }, { pivot: 'id' }),
                    options
                }, { isFlat: true, label: 'PDB' }),
                'pdb-ihm': PD.Group({
                    provider: PD.Group({
                        id: PD.Text('8zzc', { label: 'PDB-IHM Id(s)', description: 'One or more comma/space separated ids.' }),
                        encoding: PD.Select('bcif', PD.arrayToOptions(['cif', 'bcif'] as const)),
                    }, { pivot: 'id' }),
                    options
                }, { isFlat: true, label: 'PDB-IHM' }),
                'swissmodel': PD.Group({
                    id: PD.Text('Q9Y2I8', { label: 'UniProtKB AC(s)', description: 'One or more comma/space separated ACs.' }),
                    options
                }, { isFlat: true, label: 'SWISS-MODEL', description: 'Loads the best homology model or experimental structure' }),
                'alphafolddb': PD.Group({
                    provider: PD.Group({
                        id: PD.Text('Q8W3K0', { label: 'UniProtKB AC(s)', description: 'One or more comma/space separated ACs.' }),
                        encoding: PD.Select('bcif', PD.arrayToOptions(['cif', 'bcif'] as const)),
                    }, { pivot: 'id' }),
                    options
                }, { isFlat: true, label: 'AlphaFold DB', description: 'Loads the predicted model if available' }),
                'modelarchive': PD.Group({
                    id: PD.Text('ma-bak-cepc-0003', { label: 'Accession Code(s)', description: 'One or more comma/space separated ACs.' }),
                    options
                }, { isFlat: true, label: 'Model Archive' }),
                'pubchem': PD.Group({
                    id: PD.Text('2244,2245', { label: 'PubChem ID', description: 'One or more comma/space separated IDs.' }),
                    options
                }, { isFlat: true, label: 'PubChem', description: 'Loads 3D conformer from PubChem.' }),
                'url': PD.Group({
                    url: PD.Url(''),
                    format: PD.Select<BuiltInTrajectoryFormat>('mmcif', PD.arrayToOptions(BuiltInTrajectoryFormats.map(f => f[0]), f => f)),
                    isBinary: PD.Boolean(false),
                    label: PD.Optional(PD.Text('')),
                    options
                }, { isFlat: true, label: 'URL' })
            })
        };
    }
})(({ params, state }, plugin: PluginContext) => Task.create('Download Structure', async ctx => {
    plugin.behaviors.layout.leftPanelTabName.next('data');

    const src = params.source;
    let downloadParams: StateTransformer.Params<Download>[];
    let asTrajectory = false;
    let format: BuiltInTrajectoryFormat = 'mmcif';

    switch (src.name) {
        case 'url':
            downloadParams = [{ url: src.params.url, isBinary: src.params.isBinary, label: src.params.label || undefined }];
            format = src.params.format;
            break;
        case 'pdb':
            downloadParams = await (
                src.params.provider.server.name === 'pdbe'
                    ? getPdbeDownloadParams(src)
                    : src.params.provider.server.name === 'pdbj'
                        ? getPdbjDownloadParams(src)
                        : src.params.provider.server.name === 'rcsb'
                            ? getRcsbDownloadParams(src)
                            : assertUnreachable(src as never)
            );
            asTrajectory = !!src.params.options.asTrajectory;
            break;
        case 'pdb-ihm':
            const map = (id: string) => id.startsWith('PDBDEV_') ? id : `PDBDEV_${id.padStart(8, '0')}`;
            downloadParams = await getDownloadParams(src.params.provider.id,
                id => {
                    // 4 character PDB id, TODO: support extended PDB ID
                    if (id.match(/^[1-9][A-Z0-9]{3}$/i) !== null) {
                        return src.params.provider.encoding === 'bcif'
                            ? `https://pdb-ihm.org/bcif/${id.toLowerCase()}.bcif`
                            : `https://pdb-ihm.org/cif/${id.toLowerCase()}.cif`;
                    }
                    const nId = map(id.toUpperCase());
                    return src.params.provider.encoding === 'bcif'
                        ? `https://pdb-ihm.org/bcif/${nId}.bcif`
                        : `https://pdb-ihm.org/cif/${nId}.cif`;
                },
                id => { const nId = id.toUpperCase(); return nId.match(/^[1-9][A-Z0-9]{3}$/) ? `PDB-IHM: ${nId}` : map(nId); },
                src.params.provider.encoding === 'bcif'
            );
            asTrajectory = !!src.params.options.asTrajectory;
            break;
        case 'swissmodel':
            downloadParams = await getDownloadParams(src.params.id, id => `https://swissmodel.expasy.org/repository/uniprot/${id.toUpperCase()}.pdb`, id => `SWISS-MODEL: ${id}`, false);
            asTrajectory = !!src.params.options.asTrajectory;
            format = 'pdb';
            break;
        case 'alphafolddb':
            downloadParams = await getDownloadParams(src.params.provider.id,
                async id => {
                    const url = `https://www.alphafold.ebi.ac.uk/api/prediction/${id.toUpperCase()}`;
                    const info = await plugin.runTask(plugin.fetch({ url, type: 'json' }));
                    if (Array.isArray(info) && info.length > 0) {
                        const prop = src.params.provider.encoding === 'bcif' ? 'bcifUrl' : 'cifUrl';
                        return info[0][prop];
                    }
                    throw new Error(`No AlphaFold DB entry for '${id}'`);
                },
                id => `AlphaFold DB: ${id}`,
                src.params.provider.encoding === 'bcif'
            );
            asTrajectory = !!src.params.options.asTrajectory;
            format = 'mmcif';
            break;
        case 'modelarchive':
            downloadParams = await getDownloadParams(src.params.id, id => `https://www.modelarchive.org/doi/10.5452/${id.toLowerCase()}.cif`, id => `Model Archive: ${id}`, false);
            asTrajectory = !!src.params.options.asTrajectory;
            format = 'mmcif';
            break;
        case 'pubchem':
            downloadParams = await getDownloadParams(src.params.id, id => `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/${id.trim()}/record/SDF/?record_type=3d`, id => `PubChem: ${id}`, false);
            asTrajectory = !!src.params.options.asTrajectory;
            format = 'mol';
            break;
        default: assertUnreachable(src);
    }

    const representationPreset: any = params.source.params.options.representation || plugin.config.get(PluginConfig.Structure.DefaultRepresentationPreset) || PresetStructureRepresentations.auto.id;
    const showUnitcell = representationPreset !== PresetStructureRepresentations.empty.id;

    const structure = src.params.options.type.name === 'auto' ? void 0 : src.params.options.type;
    await state.transaction(async () => {
        if (downloadParams.length > 0 && asTrajectory) {
            const blob = await plugin.builders.data.downloadBlob({
                sources: downloadParams.map((src, i) => ({ id: '' + i, url: src.url, isBinary: src.isBinary })),
                maxConcurrency: 6
            }, { state: { isGhost: true } });
            const trajectory = await plugin.builders.structure.parseTrajectory(blob, { formats: downloadParams.map((_, i) => ({ id: '' + i, format: 'cif' as const })) });

            await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default', {
                structure,
                showUnitcell,
                representationPreset,
                representationPresetParams: params.source.params.options.representationParams
            });
        } else {
            for (const download of downloadParams) {
                const data = await plugin.builders.data.download(download, { state: { isGhost: true } });
                const provider = plugin.dataFormats.get(format);
                if (!provider) throw new Error('unknown file format');
                const trajectory = await plugin.builders.structure.parseTrajectory(data, provider);

                await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default', {
                    structure,
                    showUnitcell,
                    representationPreset,
                    representationPresetParams: params.source.params.options.representationParams
                });
            }
        }
    }).runInContext(ctx);
}));

async function getDownloadParams(src: string, url: (id: string) => string | Promise<string>, label: (id: string) => string, isBinary: boolean): Promise<StateTransformer.Params<Download>[]> {
    const ids = src.split(/[,\s]/).map(id => id.trim()).filter(id => !!id && (id.length >= 4 || /^[1-9][0-9]*$/.test(id)));
    const ret: StateTransformer.Params<Download>[] = [];
    for (const id of ids) {
        ret.push({ url: Asset.Url(await url(id)), isBinary, label: label(id) });
    }
    return ret;
}

async function getPdbeDownloadParams(src: ReturnType<DownloadStructure['createDefaultParams']>['source']) {
    if (src.name !== 'pdb' || src.params.provider.server.name !== 'pdbe') throw new Error('expected pdbe');
    return src.params.provider.server.params.variant === 'updated'
        ? getDownloadParams(src.params.provider.id, id => `https://www.ebi.ac.uk/pdbe/static/entry/${id.toLowerCase()}_updated.cif`, id => `PDBe: ${id} (updated cif)`, false)
        : src.params.provider.server.params.variant === 'updated-bcif'
            ? getDownloadParams(src.params.provider.id, id => `https://www.ebi.ac.uk/pdbe/entry-files/download/${id.toLowerCase()}.bcif`, id => `PDBe: ${id} (updated cif)`, true)
            : getDownloadParams(src.params.provider.id, id => `https://www.ebi.ac.uk/pdbe/static/entry/${id.toLowerCase()}.cif`, id => `PDBe: ${id} (cif)`, false);
}

async function getPdbjDownloadParams(src: ReturnType<DownloadStructure['createDefaultParams']>['source']) {
    if (src.name !== 'pdb' || src.params.provider.server.name !== 'pdbj') throw new Error('expected pdbj');
    return getDownloadParams(src.params.provider.id, id => `https://data.pdbjlc1.pdbj.org/pub/pdb/data/structures/divided/mmCIF/${id.toLowerCase().substring(1, 3)}/${id.toLowerCase()}.cif`, id => `PDBj: ${id} (cif)`, false);
}

async function getRcsbDownloadParams(src: ReturnType<DownloadStructure['createDefaultParams']>['source']) {
    if (src.name !== 'pdb' || src.params.provider.server.name !== 'rcsb') throw new Error('expected rcsb');
    return src.params.provider.server.params.encoding === 'cif'
        ? getDownloadParams(src.params.provider.id, id => `https://files.rcsb.org/download/${id.toUpperCase()}.cif`, id => `RCSB PDB: ${id} (cif)`, false)
        : getDownloadParams(src.params.provider.id, id => `https://models.rcsb.org/${id.toUpperCase()}.bcif`, id => `RCSB PDB: ${id} (bcif)`, true);
}

export const UpdateTrajectory = StateAction.build({
    display: { name: 'Update Trajectory' },
    params: {
        action: PD.Select('advance', PD.arrayToOptions(['advance', 'reset'] as const)),
        by: PD.Optional(PD.Numeric(1, { min: -1, max: 1, step: 1 }))
    }
})(({ params, state }) => {
    const models = state.selectQ(q => q.ofTransformer(StateTransforms.Model.ModelFromTrajectory));

    const update = state.build();

    if (params.action === 'reset') {
        for (const m of models) {
            update.to(m).update({ modelIndex: 0 });
        }
    } else {
        for (const m of models) {
            const parent = StateSelection.findAncestorOfType(state.tree, state.cells, m.transform.ref, PluginStateObject.Molecule.Trajectory);
            if (!parent || !parent.obj) continue;
            const traj = parent.obj;
            update.to(m).update(old => {
                let modelIndex = (old.modelIndex + params.by!) % traj.data.frameCount;
                if (modelIndex < 0) modelIndex += traj.data.frameCount;
                return { modelIndex };
            });
        }
    }

    return state.updateTree(update);
});

export const EnableModelCustomProps = StateAction.build({
    display: { name: 'Custom Model Properties', description: 'Enable parameters for custom properties of the model.' },
    from: PluginStateObject.Molecule.Model,
    params(a, ctx: PluginContext) {
        return ctx.customModelProperties.getParams(a?.data);
    },
    isApplicable(a, t, ctx: PluginContext) {
        return t.transformer !== CustomModelProperties;
    }
})(({ ref, params }, ctx: PluginContext) => ctx.builders.structure.insertModelProperties(ref, params));

export const EnableStructureCustomProps = StateAction.build({
    display: { name: 'Custom Structure Properties', description: 'Enable parameters for custom properties of the structure.' },
    from: PluginStateObject.Molecule.Structure,
    params(a, ctx: PluginContext) {
        return ctx.customStructureProperties.getParams(a?.data);
    },
    isApplicable(a, t, ctx: PluginContext) {
        return t.transformer !== CustomStructureProperties;
    }
})(({ ref, params }, ctx: PluginContext) => ctx.builders.structure.insertStructureProperties(ref, params));

export const AddTrajectory = StateAction.build({
    display: { name: 'Add Trajectory', description: 'Add trajectory from existing model/topology and coordinates.' },
    from: PluginStateObject.Root,
    params(a, ctx: PluginContext) {
        const state = ctx.state.data;
        const models = [
            ...state.selectQ(q => q.rootsOfType(PluginStateObject.Molecule.Model)),
            ...state.selectQ(q => q.rootsOfType(PluginStateObject.Molecule.Topology)),
        ];
        const modelOptions = models.map(t => [t.transform.ref, t.obj!.label]) as [string, string][];
        const coords = state.selectQ(q => q.rootsOfType(PluginStateObject.Molecule.Coordinates));
        const coordOptions = coords.map(c => [c.transform.ref, c.obj!.label]) as [string, string][];
        return {
            model: PD.Select(modelOptions.length ? modelOptions[0][0] : '', modelOptions),
            coordinates: PD.Select(coordOptions.length ? coordOptions[0][0] : '', coordOptions)
        };
    }
})(({ params, state }, ctx: PluginContext) => Task.create('Add Trajectory', taskCtx => {
    return state.transaction(async () => {
        const dependsOn = [params.model, params.coordinates];
        const model = state.build().toRoot()
            .apply(TrajectoryFromModelAndCoordinates, {
                modelRef: params.model,
                coordinatesRef: params.coordinates
            }, { dependsOn })
            .apply(StateTransforms.Model.ModelFromTrajectory, { modelIndex: 0 });

        await state.updateTree(model).runInContext(taskCtx);
        const structure = await ctx.builders.structure.createStructure(model.selector);
        await ctx.builders.structure.representation.applyPreset(structure, 'auto');
    }).runInContext(taskCtx);
}));

export const LoadTrajectory = StateAction.build({
    display: { name: 'Load Trajectory', description: 'Load trajectory of model/topology and coordinates from URL or file.' },
    from: PluginStateObject.Root,
    params(a, ctx: PluginContext) {
        const { options } = ctx.dataFormats;
        const modelOptions = options.filter(o => o[2] === TrajectoryFormatCategory || o[2] === TopologyFormatCategory);
        const coordinatesOptions = options.filter(o => o[2] === CoordinatesFormatCategory);

        const modelExts: string[] = [];
        const coordinatesExts: string[] = [];
        for (const { provider } of ctx.dataFormats.list) {
            if (provider.category === TrajectoryFormatCategory || provider.category === TopologyFormatCategory) {
                if (provider.binaryExtensions) modelExts.push(...provider.binaryExtensions);
                if (provider.stringExtensions) modelExts.push(...provider.stringExtensions);
            } else if (provider.category === CoordinatesFormatCategory) {
                if (provider.binaryExtensions) coordinatesExts.push(...provider.binaryExtensions);
                if (provider.stringExtensions) coordinatesExts.push(...provider.stringExtensions);
            }
        }

        return {
            source: PD.MappedStatic('file', {
                url: PD.Group({
                    model: PD.Group({
                        url: PD.Url(''),
                        format: PD.Select(modelOptions[0][0], modelOptions),
                        isBinary: PD.Boolean(false),
                    }, { isExpanded: true }),
                    coordinates: PD.Group({
                        url: PD.Url(''),
                        format: PD.Select(coordinatesOptions[0][0], coordinatesOptions),
                    }, { isExpanded: true })
                }, { isFlat: true }),
                file: PD.Group({
                    model: PD.File({ accept: modelExts.map(e => `.${e}`).join(','), label: 'Model' }),
                    coordinates: PD.File({ accept: coordinatesExts.map(e => `.${e}`).join(','), label: 'Coordinates' }),
                }, { isFlat: true }),
            }, { options: [['url', 'URL'], ['file', 'File']] })
        };
    }
})(({ params, state }, ctx: PluginContext) => Task.create('Load Trajectory', taskCtx => {
    return state.transaction(async () => {
        const s = params.source;

        if (s.name === 'file' && (s.params.model === null || s.params.coordinates === null)) {
            ctx.log.error('No file(s) selected');
            return;
        }

        if (s.name === 'url' && (!s.params.model || !s.params.coordinates)) {
            ctx.log.error('No URL(s) given');
            return;
        }

        const processUrl = async (url: string | Asset.Url, format: string, isBinary: boolean) => {
            const data = await ctx.builders.data.download({ url, isBinary });
            const provider = ctx.dataFormats.get(format);

            if (!provider) {
                ctx.log.warn(`LoadTrajectory: could not find data provider for '${format}'`);
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
                ctx.log.warn(`LoadTrajectory: could not find data provider for '${info.ext}'`);
                await ctx.state.data.build().delete(data).commit();
                return;
            }

            return provider.parse(ctx, data);
        };

        try {
            const modelParsed = s.name === 'url'
                ? await processUrl(s.params.model.url, s.params.model.format, s.params.model.isBinary)
                : await processFile(s.params.model);

            let model;
            if ('trajectory' in modelParsed) {
                model = await state.build().to(modelParsed.trajectory)
                    .apply(ModelFromTrajectory, { modelIndex: 0 })
                    .commit();
            } else {
                model = modelParsed.topology;
            }

            //

            const coordinates = s.name === 'url'
                ? await processUrl(s.params.coordinates.url, s.params.coordinates.format, true)
                : await processFile(s.params.coordinates);

            //

            const dependsOn = [model.ref, coordinates.ref];
            const traj = state.build().toRoot()
                .apply(TrajectoryFromModelAndCoordinates, {
                    modelRef: model.ref,
                    coordinatesRef: coordinates.ref
                }, { dependsOn })
                .apply(StateTransforms.Model.ModelFromTrajectory, { modelIndex: 0 });

            await state.updateTree(traj).runInContext(taskCtx);
            const structure = await ctx.builders.structure.createStructure(traj.selector);
            await ctx.builders.structure.representation.applyPreset(structure, 'auto');
        } catch (e) {
            console.error(e);
            ctx.log.error(`Error loading trajectory`);
        }
    }).runInContext(taskCtx);
}));