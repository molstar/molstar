/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */
import { CIF } from '../../mol-io/reader/cif';
import { volumeFromDensityServerData } from '../../mol-model-formats/volume/density-server';
import { volumeFromSegmentationData } from '../../mol-model-formats/volume/segmentation';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { PluginContext } from '../../mol-plugin/context';
import { StateTransformer } from '../../mol-state';
import { Task } from '../../mol-task';
import { ParamDefinition } from '../../mol-util/param-definition';
import { MeshData, VolsegMeshData, VolsegMeshDataParams, VolsegMeshSegmentation } from '../new-meshes/mesh-extension';

import { RawMeshSegmentData, VolsegEntry, VolsegEntryData, createVolsegEntryParams } from './entry-root';
import { VolsegState, VolsegStateParams, VOLSEG_STATE_FROM_ENTRY_TRANSFORMER_NAME } from './entry-state';
import { VolsegGlobalState, VolsegGlobalStateData, VolsegGlobalStateParams } from './global-state';
import { CreateTransformer } from './helpers';
import { VolsegGeometricSegmentation, VolsegShapePrimitivesData } from './shape_primitives';
import { GeometricSegmentationData, ShapePrimitiveData } from './volseg-api/data';

export const ProjectDataParams = {
    timeframeIndex: ParamDefinition.Numeric(0, { step: 1 }),
    channelId: ParamDefinition.Text('0'),
};

export const ProjectSegmentationDataParams = {
    timeframeIndex: ParamDefinition.Numeric(0, { step: 1 }),
    segmentationId: ParamDefinition.Text('0'),
    segmentLabels: ParamDefinition.ObjectList({ id: ParamDefinition.Numeric(-1), label: ParamDefinition.Text('') }, s => `${s.id} = ${s.label}`, { description: 'Mapping of segment IDs to segment labels' }),
    ownerId: ParamDefinition.Text('', { isHidden: true, description: 'Reference to the object which manages this volume' }),
};

export const ProjectMeshSegmentationDataParams = {
    timeframeIndex: ParamDefinition.Numeric(0, { step: 1 }),
    segmentationId: ParamDefinition.Text('0'),
    ...VolsegMeshDataParams
};

export type ProjectDataParamsValues = ParamDefinition.Values<typeof ProjectDataParams>;
export type ProjectSegmentationDataParamsValues = ParamDefinition.Values<typeof ProjectSegmentationDataParams>;
export type ProjectMeshSegmentationDataParamsValues = ParamDefinition.Values<typeof ProjectMeshSegmentationDataParams>;

export const ProjectVolumeData = CreateTransformer({
    name: 'project-volume-data',
    display: { name: 'Project Volume Data', description: 'Project Volume Data' },
    from: PluginStateObject.Root,
    to: PluginStateObject.Volume.Data,
    params: ProjectDataParams
})({
    apply({ a, params, spine }, plugin: PluginContext) {
        return Task.create('Project Volume Data', async ctx => {
            const { timeframeIndex, channelId } = params;
            const entry = spine.getAncestorOfType(VolsegEntry);
            const entryData = entry!.data;
            const rawData = await entryData.getData(timeframeIndex, channelId, 'volume') as Uint8Array | string;
            let label = entryData.metadata.getVolumeChannelLabel(channelId);
            if (!label) label = channelId.toString();

            const parsed = await CIF.parse(rawData).runInContext(ctx);
            if (parsed.isError) throw new Error(parsed.message);
            const cif = new PluginStateObject.Format.Cif(parsed.result);

            const header = cif.data.blocks[1].header; // zero block contains query meta-data
            const block = cif.data.blocks.find(b => b.header === header);
            if (!block) throw new Error(`Data block '${[header]}' not found.`);
            const volumeCif = CIF.schema.densityServer(block);
            const volume = await volumeFromDensityServerData(volumeCif).runInContext(ctx);
            const [x, y, z] = volume.grid.cells.space.dimensions;
            const props = { label: `Volume channel: ${label}`, description: `Volume ${x}\u00D7${y}\u00D7${z}` };
            return new PluginStateObject.Volume.Data(volume, props);
        });
    }
});

export const ProjectSegmentationData = CreateTransformer({
    name: 'project-segmentation-data',
    display: { name: 'Project Segmentation Data', description: 'Project Segmentation Data' },
    from: PluginStateObject.Root,
    to: PluginStateObject.Volume.Data,
    params: ProjectSegmentationDataParams
})({
    apply({ a, params, spine }, plugin: PluginContext) {
        return Task.create('Project Volume Data', async ctx => {
            const { timeframeIndex, segmentationId } = params;
            const entry = spine.getAncestorOfType(VolsegEntry);
            // const entry = a;
            const entryData = entry!.data;
            const rawData = await entryData.getData(timeframeIndex, segmentationId, 'segmentation') as Uint8Array | string;

            // NOTE: label - segmentationId;
            const label = segmentationId;

            const parsed = await CIF.parse(rawData).runInContext(ctx);
            // const parsed = await entryData.plugin.dataFormats.get('dscif')!.parse(entryData.plugin, data);
            if (parsed.isError) throw new Error(parsed.message);
            const cif = new PluginStateObject.Format.Cif(parsed.result);

            const header = cif.data.blocks[1].header; // zero block contains query meta-data
            const block = cif.data.blocks.find(b => b.header === header);
            if (!block) throw new Error(`Data block '${[header]}' not found.`);
            const segmentationCif = CIF.schema.segmentation(block);
            const segmentLabels: { [id: number]: string } = {};

            for (const segment of params.segmentLabels) segmentLabels[segment.id] = segment.label;
            // TODO: check here segment labels;
            // here params.segmentLabels has all 40 segment labels
            // should have only segment labels for that timeframe?
            const volume = await volumeFromSegmentationData(segmentationCif, { label: label, segmentLabels: segmentLabels, ownerId: params.ownerId }).runInContext(ctx);
            console.log(volume);
            const [x, y, z] = volume.grid.cells.space.dimensions;
            const props = { label: `ID: ${label}`, description: `Segmentation ID: ${label} ${x}\u00D7${y}\u00D7${z}` };
            return new PluginStateObject.Volume.Data(volume, props);
        });
    }
});

export const ProjectMeshData = CreateTransformer({
    name: 'project-mesh-data',
    display: { name: 'Project Mesh Data', description: 'Project Mesh Data' },
    from: VolsegEntry,
    to: VolsegMeshSegmentation,
    params: ProjectMeshSegmentationDataParams
})({
    apply({ a, params, spine }, plugin: PluginContext) {
        return Task.create('Project Mesh Data', async ctx => {
            const { timeframeIndex, segmentationId } = params;
            // TODO: alternatively to using a
            const entry = spine.getAncestorOfType(VolsegEntry);
            // const entry = a;
            const entryData = entry!.data;
            const segmentsToCreate = entryData.metadata.getMeshSegmentIdsForSegmentationIdAndTimeframe(segmentationId, timeframeIndex);

            const group = entryData.findNodesByTags('mesh-segmentation-group')[0]?.transform.ref;

            const totalVolume = entryData.metadata.gridTotalVolume;
            const meshData: MeshData[] = [];
            const segmentsParams = params.meshSegmentParams;
            const rawCifData = await entryData._loadRawMeshSegmentationData(timeframeIndex, segmentationId);
            const rawDataArray: RawMeshSegmentData[] = rawCifData.data as RawMeshSegmentData[];
            for (const segmentParam of segmentsParams) {
                const rawDataItem = rawDataArray.find(i => i.segmentId === segmentParam.id);
                if (!rawDataItem) throw new Error('no segment');
                const rawData = rawDataItem.data;
                const parsed = await CIF.parse(rawData).runInContext(ctx);
                if (parsed.isError) throw new Error(parsed.message);
                // const cif = new PluginStateObject.Format.Cif(parsed.result);
                const cif = parsed.result;
                const meshDataItem: MeshData = {
                    meshSegmentParams: segmentParam,
                    parsedCif: cif
                };
                meshData.push(meshDataItem);
            }
            return new VolsegMeshSegmentation(new VolsegMeshData(meshData));
        });
    }
});

// params should be similar to mesh segmentation data
export const ProjectGeometricSegmentationDataParams = {
    timeframeIndex: ParamDefinition.Numeric(0, { step: 1 }),
    segmentationId: ParamDefinition.Text('0'),
    // ...VolsegMeshDataParams
};

export type ProjectGeometricSegmentationDataParamsValues = ParamDefinition.Values<typeof ProjectGeometricSegmentationDataParams>;

// should have all geometric segmentation data
export const ProjectGeometricSegmentationData = CreateTransformer({
    name: 'project-geometric-segmentation',
    display: { name: 'Project Geometric Segmentation', description: 'Project Geometric Segmentation' },
    from: VolsegEntry,
    to: VolsegGeometricSegmentation,
    params: ProjectGeometricSegmentationDataParams
})({
    apply({ a, params, spine }, plugin: PluginContext) {
        return Task.create('Project Geometric Segmentation Data', async ctx => {
            const { timeframeIndex, segmentationId } = params;
            // TODO: alternatively to using a
            const entry = spine.getAncestorOfType(VolsegEntry);
            // const entry = a;
            const entryData = entry!.data;

            const shapePrimitiveData = await entryData._loadGeometricSegmentationData(timeframeIndex, segmentationId);
            return new VolsegGeometricSegmentation(new VolsegShapePrimitivesData(shapePrimitiveData));
        });
    }
});

export const VolsegEntryFromRoot = CreateTransformer({
    name: 'volseg-entry-from-root',
    display: { name: 'Vol & Seg Entry', description: 'Vol & Seg Entry' },
    from: PluginStateObject.Root,
    to: VolsegEntry,
    params: (a, plugin: PluginContext) => createVolsegEntryParams(plugin),
})({
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Load Vol & Seg Entry', async () => {
            const data = await VolsegEntryData.create(plugin, params);
            return new VolsegEntry(data, { label: data.entryId, description: 'Vol & Seg Entry' });
        });
    },
    update({ b, oldParams, newParams }) {
        Object.assign(newParams, oldParams);
        console.error('Changing params of existing VolsegEntry node is not allowed');
        return StateTransformer.UpdateResult.Unchanged;
    }
});


export const VolsegStateFromEntry = CreateTransformer({
    name: VOLSEG_STATE_FROM_ENTRY_TRANSFORMER_NAME,
    display: { name: 'Vol & Seg Entry State', description: 'Vol & Seg Entry State' },
    from: VolsegEntry,
    to: VolsegState,
    params: VolsegStateParams,
})({
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Create Vol & Seg Entry State', async () => {
            return new VolsegState(params, { label: 'State' });
        });
    }
});


export const VolsegGlobalStateFromRoot = CreateTransformer({
    name: 'volseg-global-state-from-root',
    display: { name: 'Vol & Seg Global State', description: 'Vol & Seg Global State' },
    from: PluginStateObject.Root,
    to: VolsegGlobalState,
    params: VolsegGlobalStateParams,
})({
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Create Vol & Seg Global State', async () => {
            const data = new VolsegGlobalStateData(plugin, params);
            return new VolsegGlobalState(data, { label: 'Global State', description: 'Vol & Seg Global State' });
        });
    },
    update({ b, oldParams, newParams }) {
        b.data.currentState.next(newParams);
        return StateTransformer.UpdateResult.Updated;
    }
});