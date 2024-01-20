/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { BehaviorSubject, distinctUntilChanged, Subject, throttleTime } from 'rxjs';
import { NewVolsegVolumeServerConfig } from '.';
import { Loci } from '../../mol-model/loci';

import { ShapeGroup } from '../../mol-model/shape';
import { Volume } from '../../mol-model/volume';
import { LociLabelProvider } from '../../mol-plugin-state/manager/loci-label';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { PluginBehavior } from '../../mol-plugin/behavior';
import { PluginCommands } from '../../mol-plugin/commands';
import { PluginContext } from '../../mol-plugin/context';
import { StateObjectCell, StateSelection, StateTransform } from '../../mol-state';
import { shallowEqualObjects } from '../../mol-util';
import { Choice } from '../../mol-util/param-choice';
import { ParamDefinition } from '../../mol-util/param-definition';
import { MeshlistData } from '../meshes/mesh-extension';

import { DEFAULT_VOLSEG_SERVER, VolumeApiV2 } from './volseg-api/api';
import { Segment, TimeInfo } from './volseg-api/data';
import { MetadataWrapper } from './volseg-api/utils';
import { DEFAULT_MESH_DETAIL, VolsegMeshSegmentationData } from './entry-meshes';
import { VolsegModelData } from './entry-models';
import { VolsegLatticeSegmentationData } from './entry-segmentation';
import { VolsegState, VolsegStateData, VolsegStateParams } from './entry-state';
import { VolsegVolumeData, SimpleVolumeParamValues, VOLUME_VISUAL_TAG } from './entry-volume';
import * as ExternalAPIs from './external-api';
import { VolsegGlobalStateData } from './global-state';
import { applyEllipsis, isDefined, lazyGetter, splitEntryId } from './helpers';
import { ProjectDataParamsValues, ProjectSegmentationDataParamsValues, type VolsegStateFromEntry } from './transformers';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { Asset } from '../../mol-util/assets';
import { PluginComponent } from '../../mol-plugin-state/component';
import { VolsegGeometricSegmentationData } from './entry-geometric-segmentation';


export const MAX_VOXELS = 10 ** 7;
// export const MAX_VOXELS = 10 ** 2; // DEBUG
export const BOX: [[number, number, number], [number, number, number]] | null = null;
// export const BOX: [[number, number, number], [number, number, number]] | null = [[-90, -90, -90], [90, 90, 90]]; // DEBUG

const MAX_ANNOTATIONS_IN_LABEL = 6;
export const VOLUME_NODE_TAG = 'volume-node-tag';
export const SEGMENTATION_NODE_TAG = 'segmenation-node-tag';

const SourceChoice = new Choice({ emdb: 'EMDB', empiar: 'EMPIAR', idr: 'IDR', pdbe: 'PDBe' }, 'emdb');
export type Source = Choice.Values<typeof SourceChoice>;


export function createLoadVolsegParams(plugin?: PluginContext, entrylists: { [source: string]: string[] } = {}) {
    const defaultVolumeServer = plugin?.config.get(NewVolsegVolumeServerConfig.DefaultServer) ?? DEFAULT_VOLSEG_SERVER;
    return {
        serverUrl: ParamDefinition.Text(defaultVolumeServer),
        source: ParamDefinition.Mapped(SourceChoice.values[0], SourceChoice.options, src => entryParam(entrylists[src])),
    };
}
function entryParam(entries: string[] = []) {
    const options: [string, string][] = entries.map(e => [e, e]);
    options.push(['__custom__', 'Custom']);
    return ParamDefinition.Group({
        entryId: ParamDefinition.Select(options[0][0], options, { description: 'Choose an entry from the list, or choose "Custom" and type any entry ID (useful when using other than default server).' }),
        customEntryId: ParamDefinition.Text('', { hideIf: p => p.entryId !== '__custom__', description: 'Entry identifier, including the source prefix, e.g. "emd-1832"' }),
    }, { isFlat: true });
}
type LoadVolsegParamValues = ParamDefinition.Values<ReturnType<typeof createLoadVolsegParams>>;
type RawDataKind = 'volume' | 'segmentation';

export function createVolsegEntryParams(plugin?: PluginContext) {
    const defaultVolumeServer = plugin?.config.get(NewVolsegVolumeServerConfig.DefaultServer) ?? DEFAULT_VOLSEG_SERVER;
    return {
        serverUrl: ParamDefinition.Text(defaultVolumeServer),
        source: SourceChoice.PDSelect(),
        entryId: ParamDefinition.Text('emd-1832', { description: 'Entry identifier, including the source prefix, e.g. "emd-1832"' }),
    };
}
type VolsegEntryParamValues = ParamDefinition.Values<ReturnType<typeof createVolsegEntryParams>>;

export namespace VolsegEntryParamValues {
    export function fromLoadVolsegParamValues(params: LoadVolsegParamValues): VolsegEntryParamValues {
        let entryId = (params.source.params as any).entryId;
        if (entryId === '__custom__') {
            entryId = (params.source.params as any).customEntryId;
        }
        return {
            serverUrl: params.serverUrl,
            source: params.source.name as Source,
            entryId: entryId
        };
    }
}


export class VolsegEntry extends PluginStateObject.CreateBehavior<VolsegEntryData>({ name: 'Vol & Seg Entry' }) { }

export type VolRepr3DT = typeof StateTransforms.Representation.VolumeRepresentation3D

export class RawMeshSegmentData {
    constructor(public segmentId: number, public data: Uint8Array | string) {
        // super();
    }
}

class RawChannelData extends PluginComponent {
    constructor(public timeframeIndex: number, public channelId: string, public data: Uint8Array | string | RawMeshSegmentData[]) {
        super();
    }
}

export interface StateHierarchyMirror {
    volumes: StateObjectCell<PluginStateObject.Volume.Data>[]
    segmentations: StateObjectCell<PluginStateObject.Volume.Data>[]
}

class RawTimeframesDataCache {
    private cache: Map<string, RawChannelData> = new Map<string, RawChannelData>();
    private maxEntries: number;
    private totalSizeLimitInBytes: number = 1_000_000_000;

    constructor(public entryData: VolsegEntryData, public kind: RawDataKind) {
    }

    private setCacheSizeInItems(data: RawChannelData) {
        if (!this.maxEntries) {
            const size = this._getEntrySize(data);
            this._setCacheSizeInItems(size);
        }
    }

    private _createKey(timeframeIndex: number, channelId: string) {
        return `${timeframeIndex.toString()}_${channelId}`;
    }

    private _getEntrySize(entry: RawChannelData | RawMeshSegmentData) {
        const data = entry.data;
        let bytes: number = 0;
        if (data instanceof Uint8Array) {
            bytes = data.length;
        } else if (data instanceof String) {
            // string
            bytes = new TextEncoder().encode(data as string).length;
        } else {
            // rawMeshSegmentData
            const arr: RawMeshSegmentData[] = data as RawMeshSegmentData[];
            for (const i of arr) {
                const b = this._getEntrySize(i);
                bytes = bytes + b;
            }
        }
        console.log('Entry size is: ', bytes, ' bytes');
        return bytes;
    }

    private _setCacheSizeInItems(entrySizeInBytes: number) {
        const limit = this.totalSizeLimitInBytes / entrySizeInBytes;
        this.maxEntries = Math.round(limit);
        console.log(`maxEntries of ${this.kind} cache set to: ${this.maxEntries}`);
    };

    add(data: RawChannelData) {
        if (this.cache.size >= this.maxEntries) {
            // least-recently used cache eviction strategy
            const keyToDelete = this.cache.keys().next().value;
            console.log('key to delete', keyToDelete);
            this.cache.delete(keyToDelete);
        }
        // check if exists
        const timeframeIndex = data.timeframeIndex;
        const channelId = data.channelId;
        const key = this._createKey(timeframeIndex, channelId);
        const hasKey = this.cache.has(key);
        if (hasKey) {
            return null;
        } else {
            // add
            this.cache.set(key, data);
            this.setCacheSizeInItems(data);
        }
    }

    async get(timeframeIndex: number, channelId: string) {
        const key = this._createKey(timeframeIndex, channelId);
        // check if exists
        // debugger;
        const hasKey = this.cache.has(key);
        if (hasKey) {
            // peek the entry, re-insert for LRU strategy
            console.log('found in cache');
            const entry = this.cache.get(key)!;
            this.cache.delete(key);
            this.cache.set(key, entry);
            // debugger;
            return entry;
        } else {
            if (this.cache.size >= this.maxEntries) {
                // least-recently used cache eviction strategy
                const keyToDelete = this.cache.keys().next().value;
                console.log('key to delete', keyToDelete);
                this.cache.delete(keyToDelete);
            }
            console.log('not in cache, loaded');
            const entry = await this.entryData._loadRawChannelData(timeframeIndex, channelId, this.kind);
            this.cache.delete(key);
            this.cache.set(key, entry);
            // debugger;
            this.setCacheSizeInItems(entry);
            return entry;
        }
    }
}

export class VolsegEntryData extends PluginBehavior.WithSubscribers<VolsegEntryParamValues> {
    plugin: PluginContext;
    ref: string = '';
    api: VolumeApiV2;
    source: Source;
    /** Number part of entry ID; e.g. '1832' */
    entryNumber: string;
    /** Full entry ID; e.g. 'emd-1832' */
    entryId: string;
    metadata: MetadataWrapper;
    pdbs: string[];

    public cachedVolumeTimeframesData = new RawTimeframesDataCache(this, 'volume');
    public cachedSegmentationTimeframesData = new RawTimeframesDataCache(this, 'segmentation');

    state = {
        hierarchy: new BehaviorSubject<StateHierarchyMirror | undefined>(undefined)
    };

    public readonly volumeData = new VolsegVolumeData(this);
    public readonly geometricSegmentationData = new VolsegGeometricSegmentationData(this);
    public readonly latticeSegmentationData = new VolsegLatticeSegmentationData(this);
    public readonly meshSegmentationData = new VolsegMeshSegmentationData(this);
    private readonly modelData = new VolsegModelData(this);
    private highlightRequest = new Subject<Segment | undefined>();

    private getStateNode = lazyGetter(() => this.plugin.state.data.selectQ(q => q.byRef(this.ref).subtree().ofType(VolsegState))[0] as StateObjectCell<VolsegState, StateTransform<typeof VolsegStateFromEntry>>, 'Missing VolsegState node. Must first create VolsegState for this VolsegEntry.');
    public currentState = new BehaviorSubject(ParamDefinition.getDefaultValues(VolsegStateParams));
    public currentVolume = new BehaviorSubject<StateTransform<VolRepr3DT>[]>([]);
    public currentTimeframe = new BehaviorSubject(0);

    private constructor(plugin: PluginContext, params: VolsegEntryParamValues) {
        super(plugin, params);
        this.plugin = plugin;
        this.api = new VolumeApiV2(params.serverUrl);
        this.source = params.source;
        this.entryId = params.entryId;
        this.entryNumber = splitEntryId(this.entryId).entryNumber;
    }

    private async initialize() {
        const metadata = await this.api.getMetadata(this.source, this.entryId);
        this.metadata = new MetadataWrapper(metadata);
        this.pdbs = await ExternalAPIs.getPdbIdsForEmdbEntry(this.metadata.raw.annotation?.entry_id.source_db_id ?? this.entryId);
        // TODO use Asset?
        await this.init();
    }

    static async create(plugin: PluginContext, params: VolsegEntryParamValues) {
        const result = new VolsegEntryData(plugin, params);
        await result.initialize();
        return result;
    }

    async getData(timeframeIndex: number, channelId: string, kind: RawDataKind) {
        // TODO: optimize if elif?
        if (kind === 'volume') {
            const channelData = await this.cachedVolumeTimeframesData.get(timeframeIndex, channelId);
            return channelData.data;
        } else {
            const channelData = await this.cachedSegmentationTimeframesData.get(timeframeIndex, channelId);
            return channelData.data;
        }
    }

    sync() {
        console.log('hierarchy synced');
        const volumes = this.findNodesByTags(VOLUME_NODE_TAG);
        const segmentations = this.findNodesByTags(SEGMENTATION_NODE_TAG);
        console.log('volumes, segmentations');
        console.log(volumes, segmentations);
        this.state.hierarchy.next({ volumes, segmentations });
    }

    private async init() {
        this.sync();
        this.subscribeObservable(this.plugin.state.data.events.changed, state => {
            console.log('data events changed emitted');
            this.sync();
        });
        const hasVolumes = this.metadata.raw.grid.volumes.volume_sampling_info.spatial_downsampling_levels.length > 0;
        if (hasVolumes) {
            await this.preloadVolumeTimeframesData();
        }
        // const hasLattices = this.metadata.raw.grid.segmentation_lattices;
        // if (hasLattices) {
        //     await this.preloadSegmentationTimeframesData();
        // }
    }


    async register(ref: string) {
        this.ref = ref;
        this.plugin.managers.lociLabels.addProvider(this.labelProvider);

        try {
            const params = this.getStateNode().obj?.data;
            if (params) {
                this.currentState.next(params);
            }
        } catch {
            // do nothing
        }

        const volumeVisual = this.findNodesByTags(VOLUME_VISUAL_TAG)[0];
        if (volumeVisual) this.addCurrentVolume(volumeVisual.transform);

        const volumeRefs = new Set<string>();
        this.subscribeObservable(this.plugin.state.data.events.cell.stateUpdated, e => {
            try { (this.getStateNode()); } catch { return; } // if state not does not exist yet
            if (e.cell.transform.ref === this.getStateNode().transform.ref) {
                const newState = this.getStateNode().obj?.data;
                if (newState && !shallowEqualObjects(newState, this.currentState.value)) { // avoid repeated update
                    this.currentState.next(newState);
                }
            } else if (e.cell.transform.tags?.includes(VOLUME_VISUAL_TAG)) {
                if (volumeRefs.has(e.ref)) {
                    this.addCurrentVolume(e.cell.transform);
                } else if (StateSelection.findAncestor(this.plugin.state.data.tree, this.plugin.state.data.cells, e.ref, a => a.transform.ref === ref)) {
                    volumeRefs.add(e.ref);
                    this.addCurrentVolume(e.cell.transform);
                }
            }
        });

        this.subscribeObservable(this.plugin.state.data.events.cell.removed, e => {
            if (volumeRefs.has(e.ref)) {
                volumeRefs.delete(e.ref);
                this.removeCurrentVolume(e.ref);
            }
        });

        this.subscribeObservable(this.plugin.behaviors.interaction.click, async e => {
            const loci = e.current.loci;
            const clickedSegment = this.getSegmentIdFromLoci(loci);
            if (clickedSegment === undefined) return;
            if (clickedSegment === this.currentState.value.selectedSegment) {
                this.actionSelectSegment(undefined);
            } else {
                this.actionSelectSegment(clickedSegment);
            }
        });

        this.subscribeObservable(
            this.highlightRequest.pipe(throttleTime(50, undefined, { leading: true, trailing: true })),
            async segment => await this.highlightSegment(segment)
        );

        this.subscribeObservable(
            this.currentState.pipe(distinctUntilChanged((a, b) => a.selectedSegment === b.selectedSegment)),
            async state => {
                if (VolsegGlobalStateData.getGlobalState(this.plugin)?.selectionMode) await this.selectSegment(state.selectedSegment);
            }
        );
    }

    async unregister() {
        this.plugin.managers.lociLabels.removeProvider(this.labelProvider);
    }

    async loadVolume() {
        // const result = await this.volumeData.loadVolume();
        // if (result) {
        //     const isovalue = result.isovalue.kind === 'relative' ? result.isovalue.relativeValue : result.isovalue.absoluteValue;
        //     await this.updateStateNode({ volumeIsovalueKind: result.isovalue.kind, volumeIsovalueValue: isovalue });
        // }
    }

    private async _resolveBinaryUrl(urlString: string) {
        const url = Asset.getUrlAsset(this.plugin.managers.asset, urlString);
        const asset = this.plugin.managers.asset.resolve(url, 'binary');
        const data = (await asset.run()).data;
        return data;
    }

    async _resolveStringUrl(urlString: string) {
        const url = Asset.getUrlAsset(this.plugin.managers.asset, urlString);
        const asset = this.plugin.managers.asset.resolve(url, 'string');
        const data = (await asset.run()).data;
        return data;
    }

    async _loadRawMeshChannelData(timeframe: number, channelId: string) {
        const segmentsData: RawMeshSegmentData[] = [];
        const segmentsToCreate = this.metadata.meshSegmentIds;
        for (const seg of segmentsToCreate) {
            const detail = this.metadata.getSufficientMeshDetail(seg, DEFAULT_MESH_DETAIL);
            const urlString = this.api.meshUrl_Bcif(this.source, this.entryId, seg, detail, timeframe, channelId);
            const data = await this._resolveBinaryUrl(urlString);
            segmentsData.push(
                new RawMeshSegmentData(
                    seg,
                    data
                )
            );
        }
        return new RawChannelData(
            timeframe,
            channelId,
            segmentsData
        );
    }

    async _loadRawChannelData(timeframe: number, channelId: string, kind: RawDataKind) {
        let urlString: string;
        if (kind === 'volume') {
            urlString = this.api.volumeUrl(this.source, this.entryId, timeframe, channelId, BOX, MAX_VOXELS);
        } else {
            urlString = this.api.latticeUrl(this.source, this.entryId, 0, timeframe, channelId, BOX, MAX_VOXELS);
        }
        // const url = Asset.getUrlAsset(this.plugin.managers.asset, urlString);
        // const asset = this.plugin.managers.asset.resolve(url, 'binary');
        const data = await this._resolveBinaryUrl(urlString);
        return new RawChannelData(
            timeframe,
            channelId,
            data
        );
    }

    private async loadRawChannelsData(timeInfo: TimeInfo, channelIds: number[], kind: RawDataKind) {
        // TODO: make it run in parallel without await
        // TODO: or to run it in ctx?

        const start = timeInfo.start;
        const end = timeInfo.end;
        for (let i = start; i <= end; i++) {
            const channelsData: RawChannelData[] = [];
            for (const channelId of channelIds) {
                const rawChannelData = await this._loadRawChannelData(i, channelId, kind);
                // channelsData.push(rawChannelData);
                if (kind === 'volume') {
                    this.cachedVolumeTimeframesData.add(
                        rawChannelData
                    );
                } else {
                    this.cachedSegmentationTimeframesData.add(
                        rawChannelData
                    );
                }
            }
        }
    }

    async preloadVolumeTimeframesData() {
        const timeInfo = this.metadata.raw.grid.volumes.time_info;
        const channelIds = this.metadata.raw.grid.volumes.channel_ids;
        this.loadRawChannelsData(timeInfo, channelIds, 'volume');
        console.log('cachedVolumeTimeframesData');
        console.log(this.cachedVolumeTimeframesData);
    }

    async preloadSegmentationTimeframesData() {
        // NOTE: for now single lattice
        const timeInfo = this.metadata.raw.grid.segmentation_lattices.time_info[0];
        const channelIds = this.metadata.raw.grid.segmentation_lattices.channel_ids[0];
        this.loadRawChannelsData(timeInfo, channelIds, 'segmentation');
        console.log('cachedSegmentationTimeframesData');
        console.log(this.cachedSegmentationTimeframesData);
    }

    async updateProjectData(timeframeIndex: number) {
        this.changeCurrentTimeframe(timeframeIndex);
        const volumes = this.state.hierarchy.value!.volumes;
        const segmenations = this.state.hierarchy.value!.segmentations;
        for (const v of volumes) {
            const projectDataTransform = v.transform.ref;
            const oldParams: ProjectDataParamsValues = v.transform.params;
            const params: ProjectDataParamsValues = {
                channelId: oldParams.channelId,
                timeframeIndex: timeframeIndex
            };
            await this.plugin.state.updateTransform(this.plugin.state.data, projectDataTransform, params, 'Project Data Transform');
        }

        for (const s of segmenations) {
            const projectDataTransform = s.transform.ref;
            const oldParams: ProjectSegmentationDataParamsValues = s.transform.params;
            const newParams: ProjectSegmentationDataParamsValues = {
                ...oldParams,
                channelId: oldParams.channelId,
                timeframeIndex: timeframeIndex
            };
            await this.plugin.state.updateTransform(this.plugin.state.data, projectDataTransform, newParams, 'Project Data Transform');
        }
    }

    async loadSegmentations() {
        await this.latticeSegmentationData.loadSegmentation();
        await this.meshSegmentationData.loadSegmentation();
        await this.actionShowSegments(this.metadata.allSegmentIds);
    }

    changeCurrentTimeframe(index: number) {
        this.currentTimeframe.next(index);
        console.log('Timeframe changed');
    }

    addCurrentVolume(t: StateTransform<VolRepr3DT>) {
        const current = this.currentVolume.value;
        const next: StateTransform<VolRepr3DT>[] = [];
        let added = false;
        for (const v of current) {
            if (v.ref === t.ref) {
                next.push(t);
                added = true;
            } else {
                next.push(v);
            }
        }
        if (!added) next.push(t);
        this.currentVolume.next(next);
        // console.log('volume added, currentVolume: ', this.currentVolume.value);
    }

    removeCurrentVolume(ref: string) {
        const current = this.currentVolume.value;
        const next: StateTransform<VolRepr3DT>[] = [];
        for (const v of current) {
            if (v.ref !== ref) {
                next.push(v);
            }
        }
        this.currentVolume.next(next);
        console.log('volume removed, currentVolume: ', this.currentVolume.value);
    }

    actionHighlightSegment(segment?: Segment) {
        this.highlightRequest.next(segment);
    }

    async actionToggleSegment(segment: number) {
        const current = this.currentState.value.visibleSegments.map(seg => seg.segmentId);
        if (current.includes(segment)) {
            await this.actionShowSegments(current.filter(s => s !== segment));
        } else {
            await this.actionShowSegments([...current, segment]);
        }
    }

    async actionToggleAllSegments() {
        const current = this.currentState.value.visibleSegments.map(seg => seg.segmentId);
        if (current.length !== this.metadata.allSegments.length) {
            await this.actionShowSegments(this.metadata.allSegmentIds);
        } else {
            await this.actionShowSegments([]);
        }
    }

    async actionSelectSegment(segment?: number) {
        if (segment !== undefined && this.currentState.value.visibleSegments.find(s => s.segmentId === segment) === undefined) {
            // first make the segment visible if it is not
            await this.actionToggleSegment(segment);
        }
        await this.updateStateNode({ selectedSegment: segment });
    }

    async actionSetOpacity(opacity: number) {
        if (opacity === this.getStateNode().obj?.data.segmentOpacity) return;
        this.latticeSegmentationData.updateOpacity(opacity);
        this.meshSegmentationData.updateOpacity(opacity);

        await this.updateStateNode({ segmentOpacity: opacity });
    }

    async actionShowFittedModel(pdbIds: string[]) {
        await this.modelData.showPdbs(pdbIds);
        await this.updateStateNode({ visibleModels: pdbIds.map(pdbId => ({ pdbId: pdbId })) });
    }

    async actionSetVolumeVisual(type: 'isosurface' | 'direct-volume' | 'off', channelId: string, transform: StateTransform) {
        await this.volumeData.setVolumeVisual(type, channelId, transform);
        const currentChannelsData = this.currentState.value.channelsData;
        // it needs to find object corresponding to that channel volume node and update only it!
        const channelToBeUpdated = currentChannelsData.filter(c => c.channelId === channelId)[0];
        channelToBeUpdated.volumeType = type;
        await this.updateStateNode({ channelsData: [...currentChannelsData] });
    }

    async actionUpdateVolumeVisual(params: SimpleVolumeParamValues, channelId: string, transform: StateTransform) {
        await this.volumeData.updateVolumeVisual(params, channelId, transform);
        const currentChannelsData = this.currentState.value.channelsData;
        // it needs to find object corresponding to that channel volume node and update only it!
        const channelToBeUpdated = currentChannelsData.filter(c => c.channelId === channelId)[0];
        channelToBeUpdated.volumeType = params.volumeType;
        channelToBeUpdated.volumeOpacity = params.opacity;
        await this.updateStateNode({ channelsData: [...currentChannelsData] });
    }



    async actionShowSegments(segments: number[]) {
        await this.latticeSegmentationData.showSegments(segments);
        await this.meshSegmentationData.showSegments(segments);
        await this.updateStateNode({ visibleSegments: segments.map(s => ({ segmentId: s })) });
    }

    private async highlightSegment(segment?: Segment) {
        await PluginCommands.Interactivity.ClearHighlights(this.plugin);
        if (segment) {
            await this.latticeSegmentationData.highlightSegment(segment);
            await this.meshSegmentationData.highlightSegment(segment);
        }
    }

    private async selectSegment(segment: number) {
        this.plugin.managers.interactivity.lociSelects.deselectAll();
        await this.latticeSegmentationData.selectSegment(segment);
        await this.meshSegmentationData.selectSegment(segment);
        await this.highlightSegment();
    }

    async updateStateNode(params: Partial<VolsegStateData>) {
        const oldParams = this.getStateNode().transform.params;
        const newParams = { ...oldParams, ...params };
        const state = this.plugin.state.data;
        const update = state.build().to(this.getStateNode().transform.ref).update(newParams);
        await PluginCommands.State.Update(this.plugin, { state, tree: update, options: { doNotUpdateCurrent: true } });
    }

    findNodesByRef(ref: string) {
        return this.plugin.state.data.selectQ(q => q.byRef(ref).subtree())[0];
    }

    /** Find the nodes under this entry root which have all of the given tags. */
    findNodesByTags(...tags: string[]) {
        return this.plugin.state.data.selectQ(q => {
            let builder = q.byRef(this.ref).subtree();
            for (const tag of tags) builder = builder.withTag(tag);
            return builder;
        });
    }

    newUpdate() {
        if (this.ref !== '') {
            return this.plugin.build().to(this.ref);
        } else {
            return this.plugin.build().toRoot();
        }
    }

    private readonly labelProvider: LociLabelProvider = {
        label: (loci: Loci): string | undefined => {
            const segmentId = this.getSegmentIdFromLoci(loci);
            if (segmentId === undefined) return;
            const segment = this.metadata.getSegment(segmentId);
            if (!segment) return;
            const annotLabels = segment.biological_annotation.external_references.map(annot => `${applyEllipsis(annot.label)} [${annot.resource}:${annot.accession}]`);
            if (annotLabels.length === 0) return;
            if (annotLabels.length > MAX_ANNOTATIONS_IN_LABEL + 1) {
                const nHidden = annotLabels.length - MAX_ANNOTATIONS_IN_LABEL;
                annotLabels.length = MAX_ANNOTATIONS_IN_LABEL;
                annotLabels.push(`(${nHidden} more annotations, click on the segment to see all)`);
            }
            return '<hr class="msp-highlight-info-hr"/>' + annotLabels.filter(isDefined).join('<br/>');
        }
    };

    private getSegmentIdFromLoci(loci: Loci): number | undefined {
        if (Volume.Segment.isLoci(loci) && loci.volume._propertyData.ownerId === this.ref) {
            if (loci.segments.length === 1) {
                return loci.segments[0];
            }
        }
        if (ShapeGroup.isLoci(loci)) {
            const meshData = (loci.shape.sourceData ?? {}) as MeshlistData;
            if (meshData.ownerId === this.ref && meshData.segmentId !== undefined) {
                return meshData.segmentId;
            }
        }
    }

    async setTryUseGpu(tryUseGpu: boolean) {
        await Promise.all([
            this.volumeData.setTryUseGpu(tryUseGpu),
            this.latticeSegmentationData.setTryUseGpu(tryUseGpu),
        ]);
    }
    async setSelectionMode(selectSegments: boolean) {
        if (selectSegments) {
            await this.selectSegment(this.currentState.value.selectedSegment);
        } else {
            this.plugin.managers.interactivity.lociSelects.deselectAll();
        }
    }

}



