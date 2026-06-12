/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StringLike } from '../../mol-io/common/string-like';
import { Volume } from '../../mol-model/volume';
import { OpenFiles } from '../../mol-plugin-state/actions/file';
import { DownloadStructure, PdbDownloadProvider } from '../../mol-plugin-state/actions/structure';
import { DownloadDensity } from '../../mol-plugin-state/actions/volume';
import { PresetTrajectoryHierarchy } from '../../mol-plugin-state/builder/structure/hierarchy-preset';
import { StructureRepresentationPresetProvider } from '../../mol-plugin-state/builder/structure/representation-preset';
import { BuiltInCoordinatesFormat } from '../../mol-plugin-state/formats/coordinates';
import { BuiltInTopologyFormat } from '../../mol-plugin-state/formats/topology';
import { BuiltInTrajectoryFormat } from '../../mol-plugin-state/formats/trajectory';
import { BuildInVolumeFormat } from '../../mol-plugin-state/formats/volume';
import { createVolumeRepresentationParams } from '../../mol-plugin-state/helpers/volume-representation-params';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { TrajectoryFromModelAndCoordinates } from '../../mol-plugin-state/transforms/model';
import { PluginCommands } from '../../mol-plugin/commands';
import { PluginConfig } from '../../mol-plugin/config';
import { PluginContext } from '../../mol-plugin/context';
import { PluginState } from '../../mol-plugin/state';
import { StateObjectSelector } from '../../mol-state';
import { Task } from '../../mol-task';
import { Asset } from '../../mol-util/assets';
import { Color } from '../../mol-util/color';
import { loadMVSData, loadMVSX } from '../mvs/components/formats';
import { loadMVS, MolstarLoadingExtension, MVSLoadOptions } from '../mvs/load';
import { MVSData } from '../mvs/mvs-data';

export async function setRemoteSnapshot(plugin: PluginContext, id: string) {
    await plugin.initialized;
    const url = `${plugin.config.get(PluginConfig.State.CurrentServer)}/get/${id}`;
    return PluginCommands.State.Snapshots.Fetch(plugin, { url });
}

export async function loadSnapshotFromUrl(plugin: PluginContext, url: string, type: PluginState.SnapshotType) {
    await plugin.initialized;
    return PluginCommands.State.Snapshots.OpenUrl(plugin, { url, type });
}

export async function loadStructureFromUrl(plugin: PluginContext, url: string, format: BuiltInTrajectoryFormat = 'mmcif', isBinary = false, options?: LoadStructureOptions & { label?: string }) {
    await plugin.initialized;
    const params = DownloadStructure.createDefaultParams(plugin.state.data.root.obj!, plugin);
    return plugin.runTask(plugin.state.data.applyAction(DownloadStructure, {
        source: {
            name: 'url',
            params: {
                url: Asset.Url(url),
                format: format as any,
                isBinary,
                label: options?.label,
                options: { ...params.source.params.options, representationParams: options?.representationParams as any },
            }
        }
    }));
}

export async function loadAllModelsOrAssemblyFromUrl(plugin: PluginContext, url: string, format: BuiltInTrajectoryFormat = 'mmcif', isBinary = false, options?: LoadStructureOptions) {
    await plugin.initialized;
    const data = await plugin.builders.data.download({ url, isBinary }, { state: { isGhost: true } });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, format);
    await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'all-models', { useDefaultIfSingleModel: true, representationPresetParams: options?.representationParams });
}

export async function loadStructureFromData(plugin: PluginContext, data: string | number[], format: BuiltInTrajectoryFormat, options?: { dataLabel?: string }) {
    await plugin.initialized;
    const _data = await plugin.builders.data.rawData({ data, label: options?.dataLabel });
    const trajectory = await plugin.builders.structure.parseTrajectory(_data, format);
    await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
}

export async function loadPdb(plugin: PluginContext, pdb: string, options?: LoadStructureOptions) {
    await plugin.initialized;
    const params = DownloadStructure.createDefaultParams(plugin.state.data.root.obj!, plugin);
    const provider = plugin.config.get(PluginConfig.Download.DefaultPdbProvider)!;
    return plugin.runTask(plugin.state.data.applyAction(DownloadStructure, {
        source: {
            name: 'pdb' as const,
            params: {
                provider: {
                    id: pdb,
                    server: {
                        name: provider,
                        params: PdbDownloadProvider[provider].defaultValue as any
                    }
                },
                options: { ...params.source.params.options, representationParams: options?.representationParams as any },
            }
        }
    }));
}

export async function loadPdbIhm(plugin: PluginContext, pdbIhm: string) {
    await plugin.initialized;
    const params = DownloadStructure.createDefaultParams(plugin.state.data.root.obj!, plugin);
    return plugin.runTask(plugin.state.data.applyAction(DownloadStructure, {
        source: {
            name: 'pdb-ihm' as const,
            params: {
                provider: {
                    id: pdbIhm,
                    encoding: 'bcif',
                },
                options: params.source.params.options,
            }
        }
    }));
}

export async function loadEmdb(plugin: PluginContext, emdb: string, options?: { detail?: number }) {
    await plugin.initialized;
    const provider = plugin.config.get(PluginConfig.Download.DefaultEmdbProvider)!;
    return plugin.runTask(plugin.state.data.applyAction(DownloadDensity, {
        source: {
            name: 'pdb-emd-ds' as const,
            params: {
                provider: {
                    id: emdb,
                    server: provider,
                },
                detail: options?.detail ?? 3,
            }
        }
    }));
}

export async function loadAlphaFoldDb(plugin: PluginContext, afdb: string) {
    await plugin.initialized;
    const params = DownloadStructure.createDefaultParams(plugin.state.data.root.obj!, plugin);
    return plugin.runTask(plugin.state.data.applyAction(DownloadStructure, {
        source: {
            name: 'alphafolddb' as const,
            params: {
                provider: {
                    id: afdb,
                    encoding: 'bcif'
                },
                options: {
                    ...params.source.params.options,
                    representation: 'preset-structure-representation-ma-quality-assessment-plddt'
                },
            }
        }
    }));
}

export async function loadModelArchive(plugin: PluginContext, id: string) {
    await plugin.initialized;
    const params = DownloadStructure.createDefaultParams(plugin.state.data.root.obj!, plugin);
    return plugin.runTask(plugin.state.data.applyAction(DownloadStructure, {
        source: {
            name: 'modelarchive' as const,
            params: {
                id,
                options: params.source.params.options,
            }
        }
    }));
}

/**
 * @example Load X-ray density from volume server
    viewer.loadVolumeFromUrl({
        url: 'https://www.ebi.ac.uk/pdbe/densities/x-ray/1tqn/cell?detail=3',
        format: 'dscif',
        isBinary: true
    }, [{
        type: 'relative',
        value: 1.5,
        color: 0x3362B2
    }, {
        type: 'relative',
        value: 3,
        color: 0x33BB33,
        volumeIndex: 1
    }, {
        type: 'relative',
        value: -3,
        color: 0xBB3333,
        volumeIndex: 1
    }], {
        entryId: ['2FO-FC', 'FO-FC'],
        isLazy: true
    });
 * *********************
 * @example Load EM density from volume server
    viewer.loadVolumeFromUrl({
        url: 'https://maps.rcsb.org/em/emd-30210/cell?detail=6',
        format: 'dscif',
        isBinary: true
    }, [{
        type: 'relative',
        value: 1,
        color: 0x3377aa
    }], {
        entryId: 'EMD-30210',
        isLazy: true
    });
 */
export async function loadVolumeFromUrl(plugin: PluginContext, { url, format, isBinary }: { url: string, format: BuildInVolumeFormat, isBinary: boolean }, isovalues: VolumeIsovalueInfo[], options?: { entryId?: string | string[], isLazy?: boolean }) {
    await plugin.initialized;

    if (!plugin.dataFormats.get(format)) {
        throw new Error(`Unknown density format: ${format}`);
    }

    if (options?.isLazy) {
        const update = plugin.build();
        update.toRoot().apply(StateTransforms.Data.LazyVolume, {
            url,
            format,
            entryId: options?.entryId,
            isBinary,
            isovalues: isovalues.map(v => ({ alpha: 1, volumeIndex: 0, ...v }))
        });
        return update.commit();
    }

    return plugin.dataTransaction(async () => {
        const data = await plugin.builders.data.download({ url, isBinary }, { state: { isGhost: true } });

        const parsed = await plugin.dataFormats.get(format)!.parse(plugin, data, { entryId: options?.entryId });
        const firstVolume = (parsed.volume || parsed.volumes[0]) as StateObjectSelector<PluginStateObject.Volume.Data>;
        if (!firstVolume?.isOk) throw new Error('Failed to parse any volume.');

        const repr = plugin.build();
        for (const iso of isovalues) {
            const volume: StateObjectSelector<PluginStateObject.Volume.Data> = parsed.volumes?.[iso.volumeIndex ?? 0] ?? parsed.volume;
            const volumeData = volume.cell!.obj!.data;
            repr
                .to(volume)
                .apply(StateTransforms.Representation.VolumeRepresentation3D, createVolumeRepresentationParams(plugin, firstVolume.data!, {
                    type: 'isosurface',
                    typeParams: { alpha: iso.alpha ?? 1, isoValue: Volume.adjustedIsoValue(volumeData, iso.value, iso.type) },
                    color: 'uniform',
                    colorParams: { value: iso.color }
                }));
        }

        await repr.commit();
    });
}

export async function loadFullResolutionEMDBMap(plugin: PluginContext, emdbId: string, options: { isoValue: Volume.IsoValue, color?: Color }) {
    await plugin.initialized;
    const numericId = parseInt(emdbId.toUpperCase().replace('EMD-', ''));
    const url = `https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-${numericId}/map/emd_${numericId}.map.gz`;

    return plugin.dataTransaction(async () => {
        const data = await plugin.build().toRoot()
            .apply(StateTransforms.Data.Download, { url, isBinary: true, label: emdbId }, { state: { isGhost: true } })
            .apply(StateTransforms.Data.DeflateData)
            .commit();

        const parsed = await plugin.dataFormats.get('ccp4')!.parse(plugin, data, { entryId: emdbId });
        const firstVolume = (parsed.volume || parsed.volumes[0]) as StateObjectSelector<PluginStateObject.Volume.Data>;
        if (!firstVolume?.isOk) throw new Error('Failed to parse any volume.');

        const volume: StateObjectSelector<PluginStateObject.Volume.Data> = parsed.volumes?.[0] ?? parsed.volume;
        await plugin.build()
            .to(volume)
            .apply(StateTransforms.Representation.VolumeRepresentation3D, createVolumeRepresentationParams(plugin, firstVolume.data!, {
                type: 'isosurface',
                typeParams: { alpha: 1, isoValue: options.isoValue },
                color: 'uniform',
                colorParams: { value: options.color ?? Color(0x33BB33) }
            }))
            .commit();
    });
}

/**
 * @example
 *  viewer.loadTrajectory({
 *      model: { kind: 'model-url', url: 'villin.gro', format: 'gro' },
 *      coordinates: { kind: 'coordinates-url', url: 'villin.xtc', format: 'xtc', isBinary: true },
 *      preset: 'all-models' // or 'default'
 *  });
 */
export async function loadTrajectory(plugin: PluginContext, params: LoadTrajectoryParams) {
    await plugin.initialized;

    let model: StateObjectSelector;

    if (params.model.kind === 'model-data' || params.model.kind === 'model-url') {
        const data = params.model.kind === 'model-data'
            ? await plugin.builders.data.rawData({ data: params.model.data, label: params.modelLabel })
            : await plugin.builders.data.download({ url: params.model.url, isBinary: params.model.isBinary, label: params.modelLabel });

        const trajectory = await plugin.builders.structure.parseTrajectory(data, params.model.format ?? 'mmcif');
        model = await plugin.builders.structure.createModel(trajectory);
    } else {
        const data = params.model.kind === 'topology-data'
            ? await plugin.builders.data.rawData({ data: params.model.data, label: params.modelLabel })
            : await plugin.builders.data.download({ url: params.model.url, isBinary: params.model.isBinary, label: params.modelLabel });

        const provider = plugin.dataFormats.get(params.model.format);
        const parsed = await provider!.parse(plugin, data);
        model = parsed.topology;
    }

    const data = params.coordinates.kind === 'coordinates-data'
        ? await plugin.builders.data.rawData({ data: params.coordinates.data, label: params.coordinatesLabel })
        : await plugin.builders.data.download({ url: params.coordinates.url, isBinary: params.coordinates.isBinary, label: params.coordinatesLabel });

    const provider = plugin.dataFormats.get(params.coordinates.format);
    const coords = await provider!.parse(plugin, data);

    const trajectory = await plugin.build().toRoot()
        .apply(TrajectoryFromModelAndCoordinates, {
            modelRef: model.ref,
            coordinatesRef: coords.ref
        }, { dependsOn: [model.ref, coords.ref] })
        .commit();

    const preset = await plugin.builders.structure.hierarchy.applyPreset(trajectory, params.preset ?? 'default');

    return { model, coords, preset };
}

export async function loadMVSFromUrl(plugin: PluginContext, url: string, format: 'mvsj' | 'mvsx', options?: { appendSnapshots?: boolean, keepCamera?: boolean, keepCameraOrientation?: boolean, extensions?: MolstarLoadingExtension<any>[] }) {
    await plugin.initialized;
    if (format === 'mvsj') {
        const data = await plugin.runTask(plugin.fetch({ url, type: 'string' }));
        const mvsData = MVSData.fromMVSJ(StringLike.toString(data));
        await loadMVS(plugin, mvsData, { sanityChecks: true, sourceUrl: url, ...options });
    } else if (format === 'mvsx') {
        const data = await plugin.runTask(plugin.fetch({ url, type: 'binary' }));
        await plugin.runTask(Task.create('Load MVSX file', async ctx => {
            const parsed = await loadMVSX(plugin, ctx, data, { doNotClearAssets: options?.appendSnapshots });
            await loadMVS(plugin, parsed.mvsData, { sanityChecks: true, sourceUrl: parsed.sourceUrl, ...options });
        }));
    } else {
        throw new Error(`Unknown MolViewSpec format: ${format}`);
    }
}

/** Load MolViewSpec from `data`.
 * If `format` is 'mvsj', `data` must be a string or a Uint8Array containing a UTF8-encoded string.
 * If `format` is 'mvsx', `data` must be a Uint8Array or a string containing base64-encoded binary data prefixed with 'base64,'. */
export async function loadMvsData(plugin: PluginContext, data: string | Uint8Array<ArrayBuffer>, format: 'mvsj' | 'mvsx', options?: { appendSnapshots?: boolean, keepCamera?: boolean, keepCameraOrientation?: boolean, extensions?: MolstarLoadingExtension<any>[] }) {
    await plugin.initialized;
    return loadMVSData(plugin, data, format, options);
}

export async function loadMvsState(plugin: PluginContext, mvsData: MVSData, options?: MVSLoadOptions) {
    await plugin.initialized;
    return loadMVS(plugin, mvsData, options);
}

export async function loadFiles(plugin: PluginContext, files: File[]) {
    await plugin.initialized;
    const sessions = files.filter(f => {
        const fn = f.name.toLowerCase();
        return fn.endsWith('.molx') || fn.endsWith('.molj');
    });

    if (sessions.length > 0) {
        return PluginCommands.State.Snapshots.OpenFile(plugin, { file: sessions[0] });
    } else {
        return plugin.runTask(plugin.state.data.applyAction(OpenFiles, {
            files: files.map(f => Asset.File(f)),
            format: { name: 'auto', params: {} },
            visuals: true
        }));
    }
}

export interface LoadStructureOptions {
    representationParams?: StructureRepresentationPresetProvider.CommonParams
}

export interface VolumeIsovalueInfo {
    type: 'absolute' | 'relative',
    value: number,
    color: Color,
    alpha?: number,
    volumeIndex?: number
}

export interface LoadTrajectoryParams {
    model: { kind: 'model-url', url: string, format?: BuiltInTrajectoryFormat /* mmcif */, isBinary?: boolean }
    | { kind: 'model-data', data: string | number[] | ArrayBuffer | Uint8Array<ArrayBuffer>, format?: BuiltInTrajectoryFormat /* mmcif */ }
    | { kind: 'topology-url', url: string, format: BuiltInTopologyFormat, isBinary?: boolean }
    | { kind: 'topology-data', data: string | number[] | ArrayBuffer | Uint8Array<ArrayBuffer>, format: BuiltInTopologyFormat },
    modelLabel?: string,
    coordinates: { kind: 'coordinates-url', url: string, format: BuiltInCoordinatesFormat, isBinary?: boolean }
    | { kind: 'coordinates-data', data: string | number[] | ArrayBuffer | Uint8Array<ArrayBuffer>, format: BuiltInCoordinatesFormat },
    coordinatesLabel?: string,
    preset?: keyof PresetTrajectoryHierarchy
}