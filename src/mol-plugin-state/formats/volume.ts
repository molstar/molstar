/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Aliaksei Chareshneu <chareshneu.tech@gmail.com>
 */

import { StateTransforms } from '../transforms';
import { DataFormatProvider, guessCifVariant } from './provider';
import { PluginContext } from '../../mol-plugin/context';
import { StateObjectSelector } from '../../mol-state';
import { PluginStateObject } from '../objects';
import { VolumeRepresentation3DHelpers } from '../transforms/representation';
import { ColorNames } from '../../mol-util/color/names';
import { Volume } from '../../mol-model/volume';
import { createVolumeRepresentationParams } from '../helpers/volume-representation-params';
import { objectForEach } from '../../mol-util/object';
import { RecommendedIsoValue } from '../../mol-model-formats/volume/property';
import { getContourLevelEmdb } from '../../mol-plugin/behavior/dynamic/volume-streaming/util';
import { Task } from '../../mol-task';

export const VolumeFormatCategory = 'Volume';
type Params = { entryId?: string };

async function tryObtainRecommendedIsoValue(plugin: PluginContext, volume?: Volume) {
    if (!volume) return;

    const { entryId } = volume;
    if (!entryId || !entryId.toLowerCase().startsWith('emd')) return;

    return plugin.runTask(Task.create('Try Set Recommended IsoValue', async ctx => {
        try {
            const absIsoLevel = await getContourLevelEmdb(plugin, ctx, entryId);
            RecommendedIsoValue.Provider.set(volume, Volume.IsoValue.absolute(absIsoLevel));
        } catch (e) {
            console.warn(e);
        }
    }));
}

function tryGetRecomendedIsoValue(volume: Volume) {
    const recommendedIsoValue = RecommendedIsoValue.Provider.get(volume);
    if (!recommendedIsoValue) return;
    if (recommendedIsoValue.kind === 'relative') return recommendedIsoValue;

    return Volume.adjustedIsoValue(volume, recommendedIsoValue.absoluteValue, 'absolute');
}

async function defaultVisuals(plugin: PluginContext, data: { volume: StateObjectSelector<PluginStateObject.Volume.Data> }) {

    const typeParams: { isoValue?: Volume.IsoValue } = {};
    const isoValue = data.volume.data && tryGetRecomendedIsoValue(data.volume.data);
    if (isoValue) typeParams.isoValue = isoValue;

    const visual = plugin.build().to(data.volume).apply(StateTransforms.Representation.VolumeRepresentation3D, createVolumeRepresentationParams(plugin, data.volume.data, {
        type: 'isosurface',
        typeParams,
    }));
    return [await visual.commit()];
}

export const Ccp4Provider = DataFormatProvider({
    label: 'CCP4/MRC/MAP',
    description: 'CCP4/MRC/MAP',
    category: VolumeFormatCategory,
    binaryExtensions: ['ccp4', 'mrc', 'map'],
    parse: async (plugin, data, params?: Params) => {
        const format = plugin.build()
            .to(data)
            .apply(StateTransforms.Data.ParseCcp4, {}, { state: { isGhost: true } });

        const volume = format.apply(StateTransforms.Volume.VolumeFromCcp4, { entryId: params?.entryId });

        await format.commit({ revertOnError: true });
        await tryObtainRecommendedIsoValue(plugin, volume.selector.data);

        return { format: format.selector, volume: volume.selector };
    },
    visuals: defaultVisuals
});

export const Dsn6Provider = DataFormatProvider({
    label: 'DSN6/BRIX',
    description: 'DSN6/BRIX',
    category: VolumeFormatCategory,
    binaryExtensions: ['dsn6', 'brix'],
    parse: async (plugin, data, params?: Params) => {
        const format = plugin.build()
            .to(data)
            .apply(StateTransforms.Data.ParseDsn6, {}, { state: { isGhost: true } });

        const volume = format.apply(StateTransforms.Volume.VolumeFromDsn6, { entryId: params?.entryId });

        await format.commit({ revertOnError: true });
        await tryObtainRecommendedIsoValue(plugin, volume.selector.data);

        return { format: format.selector, volume: volume.selector };
    },
    visuals: defaultVisuals
});

export const DxProvider = DataFormatProvider({
    label: 'DX',
    description: 'DX',
    category: VolumeFormatCategory,
    stringExtensions: ['dx'],
    binaryExtensions: ['dxbin'],
    parse: async (plugin, data, params?: Params) => {
        const format = plugin.build()
            .to(data)
            .apply(StateTransforms.Data.ParseDx, {}, { state: { isGhost: true } });

        const volume = format.apply(StateTransforms.Volume.VolumeFromDx, { entryId: params?.entryId });

        await volume.commit({ revertOnError: true });
        await tryObtainRecommendedIsoValue(plugin, volume.selector.data);

        return { volume: volume.selector };
    },
    visuals: defaultVisuals
});

export const CubeProvider = DataFormatProvider({
    label: 'Cube',
    description: 'Cube',
    category: VolumeFormatCategory,
    stringExtensions: ['cub', 'cube'],
    parse: async (plugin, data, params?: Params) => {
        const format = plugin.build()
            .to(data)
            .apply(StateTransforms.Data.ParseCube, {}, { state: { isGhost: true } });

        const volume = format.apply(StateTransforms.Volume.VolumeFromCube, { entryId: params?.entryId });
        const structure = format
            .apply(StateTransforms.Model.TrajectoryFromCube, void 0, { state: { isGhost: true } })
            .apply(StateTransforms.Model.ModelFromTrajectory)
            .apply(StateTransforms.Model.StructureFromModel);

        await format.commit({ revertOnError: true });
        await tryObtainRecommendedIsoValue(plugin, volume.selector.data);

        return { format: format.selector, volume: volume.selector, structure: structure.selector };
    },
    visuals: async (plugin: PluginContext, data: { volume: StateObjectSelector<PluginStateObject.Volume.Data>, structure: StateObjectSelector<PluginStateObject.Molecule.Structure> }) => {
        const surfaces = plugin.build();

        const volumeReprs: StateObjectSelector<PluginStateObject.Volume.Representation3D>[] = [];
        const volumeData = data.volume.cell?.obj?.data;
        if (volumeData && Volume.isOrbitals(volumeData)) {
            const volumePos = surfaces.to(data.volume).apply(StateTransforms.Representation.VolumeRepresentation3D, createVolumeRepresentationParams(plugin, volumeData, {
                type: 'isosurface',
                typeParams: { isoValue: Volume.IsoValue.relative(1), alpha: 0.4 },
                color: 'uniform',
                colorParams: { value: ColorNames.blue }
            }));
            const volumeNeg = surfaces.to(data.volume).apply(StateTransforms.Representation.VolumeRepresentation3D, createVolumeRepresentationParams(plugin, volumeData, {
                type: 'isosurface',
                typeParams: { isoValue: Volume.IsoValue.relative(-1), alpha: 0.4 },
                color: 'uniform',
                colorParams: { value: ColorNames.red }
            }));
            volumeReprs.push(volumePos.selector, volumeNeg.selector);
        } else {
            const volume = surfaces.to(data.volume).apply(StateTransforms.Representation.VolumeRepresentation3D, createVolumeRepresentationParams(plugin, volumeData, {
                type: 'isosurface',
                typeParams: { isoValue: Volume.IsoValue.relative(2), alpha: 0.4 },
                color: 'uniform',
                colorParams: { value: ColorNames.grey }
            }));
            volumeReprs.push(volume.selector);
        }

        const structure = await plugin.builders.structure.representation.applyPreset(data.structure, 'auto');
        await surfaces.commit();

        const structureReprs: StateObjectSelector<PluginStateObject.Molecule.Structure.Representation3D>[] = [];
        objectForEach(structure?.representations as any, (r: any) => {
            if (r) structureReprs.push(r);
        });

        return [...volumeReprs, ...structureReprs];
    }
});

type DsCifParams = { entryId?: string | string[] };

export const DscifProvider = DataFormatProvider({
    label: 'DensityServer CIF',
    description: 'DensityServer CIF',
    category: VolumeFormatCategory,
    stringExtensions: ['cif'],
    binaryExtensions: ['bcif'],
    isApplicable: (info, data) => {
        return guessCifVariant(info, data) === 'dscif' ? true : false;
    },
    parse: async (plugin, data, params?: DsCifParams) => {
        const cifCell = await plugin.build().to(data).apply(StateTransforms.Data.ParseCif).commit();
        const b = plugin.build().to(cifCell);
        const blocks = cifCell.obj!.data.blocks;

        if (blocks.length === 0) throw new Error('no data blocks');

        const volumes: StateObjectSelector<PluginStateObject.Volume.Data>[] = [];
        let i = 0;
        for (const block of blocks) {
            // Skip "server" data block.
            if (block.header.toUpperCase() === 'SERVER') continue;

            const entryId = Array.isArray(params?.entryId) ? params?.entryId[i] : params?.entryId;
            if (block.categories['volume_data_3d_info']?.rowCount > 0) {
                volumes.push(b.apply(StateTransforms.Volume.VolumeFromDensityServerCif, { blockHeader: block.header, entryId }).selector);
                i++;
            }
        }

        await b.commit();
        for (const v of volumes) await tryObtainRecommendedIsoValue(plugin, v.data);

        return { volumes };
    },
    visuals: async (plugin, data: { volumes: StateObjectSelector<PluginStateObject.Volume.Data>[] }) => {
        const { volumes } = data;
        const tree = plugin.build();
        const visuals: StateObjectSelector<PluginStateObject.Volume.Representation3D>[] = [];

        if (volumes.length > 0) {
            const isoValue = (volumes[0].data && tryGetRecomendedIsoValue(volumes[0].data)) || Volume.IsoValue.relative(1.5);

            visuals[0] = tree
                .to(volumes[0])
                .apply(StateTransforms.Representation.VolumeRepresentation3D, VolumeRepresentation3DHelpers.getDefaultParamsStatic(plugin, 'isosurface', { isoValue, alpha: 1 }, 'uniform', { value: ColorNames.teal }))
                .selector;
        }

        if (volumes.length > 1) {
            const posParams = VolumeRepresentation3DHelpers.getDefaultParamsStatic(plugin, 'isosurface', { isoValue: Volume.IsoValue.relative(3), alpha: 0.3 }, 'uniform', { value: ColorNames.green });
            const negParams = VolumeRepresentation3DHelpers.getDefaultParamsStatic(plugin, 'isosurface', { isoValue: Volume.IsoValue.relative(-3), alpha: 0.3 }, 'uniform', { value: ColorNames.red });
            visuals[visuals.length] = tree.to(volumes[1]).apply(StateTransforms.Representation.VolumeRepresentation3D, posParams).selector;
            visuals[visuals.length] = tree.to(volumes[1]).apply(StateTransforms.Representation.VolumeRepresentation3D, negParams).selector;
        }

        await tree.commit();

        return visuals;
    }
});

export const SegcifProvider = DataFormatProvider({
    label: 'Segmentation CIF',
    description: 'Segmentation CIF',
    category: VolumeFormatCategory,
    stringExtensions: ['cif'],
    binaryExtensions: ['bcif'],
    isApplicable: (info, data) => {
        return guessCifVariant(info, data) === 'segcif' ? true : false;
    },
    parse: async (plugin, data) => {
        const cifCell = await plugin.build().to(data).apply(StateTransforms.Data.ParseCif).commit();
        const b = plugin.build().to(cifCell);
        const blocks = cifCell.obj!.data.blocks;

        if (blocks.length === 0) throw new Error('no data blocks');

        const volumes: StateObjectSelector<PluginStateObject.Volume.Data>[] = [];
        for (const block of blocks) {
            // Skip "server" data block.
            if (block.header.toUpperCase() === 'SERVER') continue;

            if (block.categories['volume_data_3d_info']?.rowCount > 0) {
                volumes.push(b.apply(StateTransforms.Volume.VolumeFromSegmentationCif, { blockHeader: block.header }).selector);
            }
        }

        await b.commit();

        return { volumes };
    },
    visuals: async (plugin, data: { volumes: StateObjectSelector<PluginStateObject.Volume.Data>[] }) => {
        const { volumes } = data;
        const tree = plugin.build();
        const visuals: StateObjectSelector<PluginStateObject.Volume.Representation3D>[] = [];

        if (volumes.length > 0) {
            const segmentation = Volume.Segmentation.get(volumes[0].data!);
            if (segmentation) {
                visuals[visuals.length] = tree
                    .to(volumes[0])
                    .apply(StateTransforms.Representation.VolumeRepresentation3D, VolumeRepresentation3DHelpers.getDefaultParams(plugin, 'segment', volumes[0].data!, { alpha: 1, instanceGranularity: true }, 'volume-segment', { }))
                    .selector;
            }
        }

        await tree.commit();

        return visuals;
    }
});

export const BuiltInVolumeFormats = [
    ['ccp4', Ccp4Provider] as const,
    ['dsn6', Dsn6Provider] as const,
    ['cube', CubeProvider] as const,
    ['dx', DxProvider] as const,
    ['dscif', DscifProvider] as const,
    ['segcif', SegcifProvider] as const,
] as const;

export type BuildInVolumeFormat = (typeof BuiltInVolumeFormats)[number][0]