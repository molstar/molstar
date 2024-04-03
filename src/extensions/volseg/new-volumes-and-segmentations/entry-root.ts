/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { BehaviorSubject, distinctUntilChanged, Subject, throttleTime } from 'rxjs';
import { NewVolsegVolumeServerConfig } from '.';
import { Loci } from '../../../mol-model/loci';

import { ShapeGroup } from '../../../mol-model/shape';
import { Volume } from '../../../mol-model/volume';
import { LociLabelProvider } from '../../../mol-plugin-state/manager/loci-label';
import { PluginStateObject } from '../../../mol-plugin-state/objects';
import { PluginBehavior } from '../../../mol-plugin/behavior';
import { PluginCommands } from '../../../mol-plugin/commands';
import { PluginContext } from '../../../mol-plugin/context';
import { StateObjectCell, StateSelection, StateTransform } from '../../../mol-state';
import { shallowEqualObjects } from '../../../mol-util';
import { Choice } from '../../../mol-util/param-choice';
import { ParamDefinition } from '../../../mol-util/param-definition';
import { isMeshlistData, MeshlistData, VolsegMeshSegmentation } from '../new-meshes/mesh-extension';

import { DEFAULT_VOLSEG_SERVER, VolumeApiV2 } from './volseg-api/api';
import { BoxPrimitive, Cylinder, DescriptionData, Ellipsoid, GridMetadata, Metadata, ParsedSegmentKey, PyramidPrimitive, SegmentAnnotationData, ShapePrimitiveData, Sphere, TimeInfo } from './volseg-api/data';
import { createSegmentKey, getCVSXGeometricSegmentationDataFromRaw, getCVSXMeshSegmentationDataFromRaw, getSegmentLabelsFromDescriptions, instanceOfShapePrimitiveData, MetadataWrapper, parseSegmentKey } from './volseg-api/utils';
import { DEFAULT_MESH_DETAIL, VolsegMeshSegmentationData } from './entry-meshes';
import { VolsegModelData } from './entry-models';
import { SEGMENT_VISUAL_TAG, VolsegLatticeSegmentationData } from './entry-segmentation';
import { VolsegState, VolsegStateData, VolsegStateParams } from './entry-state';
import { VolsegVolumeData, SimpleVolumeParamValues, VOLUME_VISUAL_TAG } from './entry-volume';
import * as ExternalAPIs from './external-api';
import { VolsegGlobalStateData } from './global-state';
import { applyEllipsis, isDefined, lazyGetter, splitEntryId } from './helpers';
import { ProjectDataParamsValues, ProjectGeometricSegmentationDataParamsValues, ProjectMeshSegmentationDataParamsValues, ProjectLatticeSegmentationDataParamsValues, type VolsegStateFromEntry } from './transformers';
import { StateTransforms } from '../../../mol-plugin-state/transforms';
import { Asset } from '../../../mol-util/assets';
import { PluginComponent } from '../../../mol-plugin-state/component';
import { VolsegGeometricSegmentationData } from './entry-geometric-segmentation';
import { createVolumeRepresentationParams } from '../../../mol-plugin-state/helpers/volume-representation-params';
import { CreateShapePrimitiveProviderParamsValues, isShapePrimitiveParamsValues, VolsegGeometricSegmentation } from './shape_primitives';
import { actionSelectSegment, parseCVSXJSON } from '../common';
import { RuntimeContext } from '../../../mol-task';
import { Unzip, unzip } from '../../../mol-util/zip/zip';
import { CVSXFilesData, CVSXLatticeSegmentationData, CVSXMeshSegmentationData, CVSXVolumeData, QueryArgs } from '../cvsx-extension/data';

export const GEOMETRIC_SEGMENTATION_NODE_TAG = 'geometric-segmentation-node';
export const MESH_SEGMENTATION_NODE_TAG = 'mesh-segmentation-node'

export const MAX_VOXELS = 10 ** 7;
// export const MAX_VOXELS = 10 ** 2; // DEBUG
export const BOX: [[number, number, number], [number, number, number]] | null = null;
// export const BOX: [[number, number, number], [number, number, number]] | null = [[-90, -90, -90], [90, 90, 90]]; // DEBUG

const MAX_ANNOTATIONS_IN_LABEL = 6;
export const VOLUME_NODE_TAG = 'volume-node-tag';
export const SEGMENTATION_NODE_TAG = 'segmenation-node-tag';

const SourceChoice = new Choice({ emdb: 'EMDB', empiar: 'EMPIAR', idr: 'IDR', pdbe: 'PDBe', custom: 'CUSTOM' }, 'emdb');
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
type RawDataKind = 'volume' | 'segmentation' | 'mesh' | 'primitive';

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

class RawSegmentationData extends PluginComponent {
    constructor(public timeframeIndex: number, public segmentationId: string, public data: Uint8Array | string | RawMeshSegmentData[] | ShapePrimitiveData) {
        super();
    }
}

export interface StateHierarchyMirror {
    volumes: StateObjectCell<PluginStateObject.Volume.Data>[]
    segmentations: StateObjectCell<PluginStateObject.Volume.Data>[]
    geometricSegmentations: StateObjectCell<VolsegGeometricSegmentation>[]
    meshSegmentations: StateObjectCell<VolsegMeshSegmentation>[]
}

class RawTimeframesDataCache {
    private cache: Map<string, RawChannelData | RawSegmentationData> = new Map<string, RawChannelData | RawSegmentationData>();
    private maxEntries: number;
    private totalSizeLimitInBytes: number = 1_000_000_000;

    constructor(public entryData: VolsegEntryData, public kind: RawDataKind) {
    }

    private setCacheSizeInItems(data: RawChannelData | RawSegmentationData) {
        if (!this.maxEntries) {
            const size = this._getEntrySize(data);
            this._setCacheSizeInItems(size);
        }
    }

    private _createKey(timeframeIndex: number, channelIdOrSegmentationId: string) {
        return `${timeframeIndex.toString()}_${channelIdOrSegmentationId}`;
    }

    private _getEntrySize(entry: RawChannelData | RawSegmentationData | RawMeshSegmentData) {
        const data = entry.data;
        let bytes: number = 0;
        if (data instanceof Uint8Array) {
            bytes = data.length;
        } else if (data instanceof String) {
            // string
            bytes = new TextEncoder().encode(data as string).length;
        } else if (instanceOfShapePrimitiveData(data)) {
            bytes = JSON.stringify(data).length;
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

    add(data: RawChannelData | RawSegmentationData) {
        if (this.cache.size >= this.maxEntries) {
            // least-recently used cache eviction strategy
            const keyToDelete = this.cache.keys().next().value;
            console.log('key to delete', keyToDelete);
            this.cache.delete(keyToDelete);
        }
        // check if exists
        const timeframeIndex = data.timeframeIndex;
        let key: string;
        if (data instanceof RawChannelData) {
            key = this._createKey(timeframeIndex, data.channelId);
        } else if (data instanceof RawSegmentationData) {
            key = this._createKey(timeframeIndex, data.segmentationId);
        } else {
            throw Error(`data type ${data}is not supported`);
        }
        const hasKey = this.cache.has(key);
        if (hasKey) {
            return null;
        } else {
            // add
            this.cache.set(key, data);
            this.setCacheSizeInItems(data);
        }
    }

    async get(timeframeIndex: number, channelIdOrSegmentationId: string, kind: RawDataKind) {
        const key = this._createKey(timeframeIndex, channelIdOrSegmentationId);
        // check if exists
        const hasKey = this.cache.has(key);
        if (hasKey) {
            // peek the entry, re-insert for LRU strategy
            console.log('found in cache');
            const entry = this.cache.get(key)!;
            this.cache.delete(key);
            this.cache.set(key, entry);
            return entry;
        } else {
            if (this.cache.size >= this.maxEntries) {
                // least-recently used cache eviction strategy
                const keyToDelete = this.cache.keys().next().value;
                console.log('key to delete', keyToDelete);
                this.cache.delete(keyToDelete);
            }
            console.log('not in cache, loaded');
            let entry;
            if (kind === 'volume') {
                entry = await this.entryData._loadRawChannelData(timeframeIndex, channelIdOrSegmentationId);
            } else if (kind === 'segmentation') {
                entry = await this.entryData._loadRawLatticeSegmentationData(timeframeIndex, channelIdOrSegmentationId);
            } else if (kind === 'mesh') {
                entry = await this.entryData._loadRawMeshSegmentationData(timeframeIndex, channelIdOrSegmentationId);
            } else if (kind === 'primitive') {
                entry = await this.entryData._loadGeometricSegmentationData(timeframeIndex, channelIdOrSegmentationId);
            } else {
                throw Error(`Data kind ${kind} is not supported`);
            }
            this.cache.delete(key);
            this.cache.set(key, entry);
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
    // metadata: MetadataWrapper;
    pdbs: string[];
    kind: 'api' | 'file' = 'api';
    filesData: CVSXFilesData | undefined = undefined;

    public metadata = new BehaviorSubject<MetadataWrapper | undefined>(undefined);

    public cachedVolumeTimeframesData = new RawTimeframesDataCache(this, 'volume');
    public cachedSegmentationTimeframesData = new RawTimeframesDataCache(this, 'segmentation');
    public cachedMeshesTimeframesData = new RawTimeframesDataCache(this, 'mesh');
    public cachedShapePrimitiveData = new RawTimeframesDataCache(this, 'primitive');
    state = {
        hierarchy: new BehaviorSubject<StateHierarchyMirror | undefined>(undefined)
    };

    public readonly volumeData = new VolsegVolumeData(this);
    public readonly geometricSegmentationData = new VolsegGeometricSegmentationData(this);
    public readonly latticeSegmentationData = new VolsegLatticeSegmentationData(this);
    public readonly meshSegmentationData = new VolsegMeshSegmentationData(this);
    private readonly modelData = new VolsegModelData(this);
    private highlightRequest = new Subject<string | undefined>();

    private getStateNode = lazyGetter(() => this.plugin.state.data.selectQ(q => q.byRef(this.ref).subtree().ofType(VolsegState))[0] as StateObjectCell<VolsegState, StateTransform<typeof VolsegStateFromEntry>>, 'Missing VolsegState node. Must first create VolsegState for this VolsegEntry.');
    public currentState = new BehaviorSubject(ParamDefinition.getDefaultValues(VolsegStateParams));
    public currentVolume = new BehaviorSubject<StateTransform<VolRepr3DT>[]>([]);
    public currentTimeframe = new BehaviorSubject(0);

    private constructor(plugin: PluginContext, params: VolsegEntryParamValues, filesData?: CVSXFilesData) {
        super(plugin, params);
        this.plugin = plugin;
        this.api = new VolumeApiV2(params.serverUrl);
        this.source = params.source;
        this.entryId = params.entryId;
        this.entryNumber = splitEntryId(this.entryId).entryNumber;
        this.filesData = filesData;
    }

    private async initialize() {
        const metadata = await this.api.getMetadata(this.source, this.entryId);
        this.metadata.next(new MetadataWrapper(metadata));
        this.pdbs = await ExternalAPIs.getPdbIdsForEmdbEntry(this.metadata.value!.raw.annotation?.entry_id.source_db_id ?? this.entryId);
        // TODO use Asset?
        await this.init();
    }

    private async initializeFromFile(metadata?: Metadata) {
        if (metadata) this.metadata.next(new MetadataWrapper(metadata));
        await this.initFromFile();
    }

    static async create(plugin: PluginContext, params: VolsegEntryParamValues) {
        const result = new VolsegEntryData(plugin, params);
        await result.initialize();
        return result;
    }

    static async createFromFile(plugin: PluginContext, file: Uint8Array, runtimeCtx: RuntimeContext) {
        const zippedFiles: { [path: string]: Uint8Array } = await unzip(runtimeCtx, file.buffer) as typeof zippedFiles;
        const zippedFilesEntries = Object.entries(zippedFiles);
        const annotationJSONEntry = zippedFilesEntries.find(z => z[0] === 'annotations.json');
        const metadataJSONEntry = zippedFilesEntries.find(z => z[0] === 'metadata.json');
        const rawQueryJSON = zippedFilesEntries.find(z => z[0] === 'query.json');
        const rawVolumes = zippedFilesEntries.filter(z => z[0].startsWith('volume'));
        const rawLatticeSegmentations = zippedFilesEntries.filter(z => z[0].startsWith('lattice'));
        const rawMeshSegmentations = zippedFilesEntries.filter(z => z[0].startsWith('mesh'));
        debugger;
        const geometricSegmentations = zippedFilesEntries.filter(z => z[0].startsWith('primitive'));
        const parsedQueryJSON: QueryArgs = await parseCVSXJSON(rawQueryJSON!, plugin);
        let params: VolsegEntryParamValues = {
            serverUrl: '',
            source: parsedQueryJSON.source_db,
            entryId: parsedQueryJSON.entry_id,
        };
        const parsedGridMetadata = await parseCVSXJSON(metadataJSONEntry!, plugin);
        const parsedAnnotationMetadata = await parseCVSXJSON(annotationJSONEntry!, plugin);
        const source: Source = parsedGridMetadata.entry_id.source_db_name;
        const entryId = parsedGridMetadata.entry_id.source_db_id;
        params = {
            serverUrl: '',
            source: source,
            entryId: entryId
        };
        const metadata: Metadata = {
            grid: parsedGridMetadata,
            annotation: parsedAnnotationMetadata
        };
        let gs = undefined;
        if (geometricSegmentations) {
            gs = await getCVSXGeometricSegmentationDataFromRaw(geometricSegmentations, parsedGridMetadata, parsedQueryJSON, plugin)
        }

        const meshSegmentationData = getCVSXMeshSegmentationDataFromRaw(rawMeshSegmentations, parsedGridMetadata, parsedQueryJSON)
        debugger;
        const filesData: CVSXFilesData = {
            // parsed everything
            volumes: rawVolumes ? rawVolumes.map(v => {
                const d: CVSXVolumeData = {
                    channelId: v[0].split('_')[1],
                    timeframeIndex: parseInt(v[0].split('_')[2]),
                    data: v[1]
                };
                return d;
            }) : undefined,
            latticeSegmentations: rawLatticeSegmentations ? rawLatticeSegmentations.map(v => {
                const d: CVSXLatticeSegmentationData = {
                    segmentationId: v[0].split('_')[1],
                    timeframeIndex: parseInt(v[0].split('_')[2]),
                    data: v[1]
                };
                return d;
            }) : undefined,
            geometricSegmentations: gs ? gs : undefined,
            meshSegmentations: meshSegmentationData ? meshSegmentationData : undefined,
            annotation: parsedAnnotationMetadata ? parsedAnnotationMetadata : undefined,
            metadata: parsedGridMetadata ? parsedGridMetadata : undefined,
            query: parsedQueryJSON
        };
        const result = new VolsegEntryData(plugin, params, filesData);
        await result.initializeFromFile(metadata);
        return result;
    }

    async getData(timeframeIndex: number, channelIdOrSegmentationId: string, kind: RawDataKind) {
        // TODO: optimize if elif?
        debugger;
        if (kind === 'volume') {
            const channelData = await this.cachedVolumeTimeframesData.get(timeframeIndex, channelIdOrSegmentationId, 'volume');
            return channelData.data;
        } else if (kind === 'segmentation') {
            const segmentationData = await this.cachedSegmentationTimeframesData.get(timeframeIndex, channelIdOrSegmentationId, 'segmentation');
            return segmentationData.data;
        } else if (kind === 'mesh') {
            const meshData = await this.cachedMeshesTimeframesData.get(timeframeIndex, channelIdOrSegmentationId, 'mesh');
            return meshData.data;
        } else if (kind === 'primitive') {
            const shapePrimitiveData = await this.cachedShapePrimitiveData.get(timeframeIndex, channelIdOrSegmentationId, 'primitive');
            return shapePrimitiveData.data;
        }
    }

    sync() {
        console.log('hierarchy synced');
        const volumes = this.findNodesByTags(VOLUME_NODE_TAG);
        const segmentations = this.findNodesByTags(SEGMENTATION_NODE_TAG);
        const geometricSegmentations = this.findNodesByTags(GEOMETRIC_SEGMENTATION_NODE_TAG);
        const meshSegmentations = this.findNodesByTags(MESH_SEGMENTATION_NODE_TAG);
        console.log('volumes, segmentations');
        console.log(volumes, segmentations);
        this.state.hierarchy.next({ volumes, segmentations, geometricSegmentations, meshSegmentations });
    }

    private async init() {
        this.sync();
        this.subscribeObservable(this.plugin.state.data.events.changed, state => {
            console.log('data events changed emitted');
            this.sync();
        });
        const hasVolumes = this.metadata.value!.raw.grid.volumes.volume_sampling_info.spatial_downsampling_levels.length > 0;
        if (hasVolumes) {
            await this.preloadVolumeTimeframesData();
        }
        const hasLattices = this.metadata.value!.raw.grid.segmentation_lattices;
        if (hasLattices) {
            await this.preloadSegmentationTimeframesData();
        }
        const hasMeshes = this.metadata.value!.raw.grid.segmentation_meshes;
        if (hasMeshes) {
            await this.preloadMeshesTimeframesData();
        }

        const hasShapePrimitives = this.metadata.value!.raw.grid.geometric_segmentation;
        if (hasShapePrimitives) {
            await this.preloadMeshesTimeframesData();
        }
    }

    private async initFromFile() {
        this.kind = 'file';
        this.sync();
        this.subscribeObservable(this.plugin.state.data.events.changed, state => {
            console.log('data events changed emitted');
            this.sync();
        });
        const hasVolumes = this.metadata.value!.raw.grid.volumes.volume_sampling_info.spatial_downsampling_levels.length > 0;
        if (hasVolumes) {
            this.preloadVolumeTimeframesDataFromFile();
        }
        const hasLattices = this.metadata.value!.raw.grid.segmentation_lattices;
        if (hasLattices) {
            this.preloadSegmentationTimeframesDataFromFile();
        }
        const hasGeometricSegmentation = this.metadata.value!.raw.grid.geometric_segmentation;
        if (hasGeometricSegmentation) {
            this.preloadShapePrimitivesTimeframesDataFromFile();
        }
        const hasMeshes = this.metadata.value!.raw.grid.segmentation_meshes;
        if (hasMeshes) {
            this.preloadMeshTimeframesDataFromFile();
        }
    }

    async updateMetadata() {
        const metadata = await this.api.getMetadata(this.source, this.entryId);
        this.metadata.next(new MetadataWrapper(metadata));
        // trigger update of state to re-render UI
        // const params = this.getStateNode().obj?.data;
        // if (params) {
        //     this.currentState.next(params);
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
            if (e.current.loci.kind === 'empty-loci') return;
            const loci = e.current.loci;
            const clickedSegmentId = this.getSegmentIdFromLoci(loci);
            const clickedSegmentSegmentationId = this.getSegmentationIdFromLoci(loci);
            const segmentationKind = this.getSegmentationKindFromLoci(loci);
            if (clickedSegmentSegmentationId === undefined) return;
            if (clickedSegmentId === undefined) return;
            if (segmentationKind === undefined) return;
            const clickedSegmentKey = createSegmentKey(clickedSegmentId, clickedSegmentSegmentationId, segmentationKind);
            if (clickedSegmentKey === this.currentState.value.selectedSegment) {
                actionSelectSegment(this, undefined);
            } else {
                actionSelectSegment(this, clickedSegmentKey);
            }
        });

        this.subscribeObservable(
            // TODO: get segment Id, segmentation id and segment kind from here
            this.highlightRequest.pipe(throttleTime(50, undefined, { leading: true, trailing: true })),
            async segmentKey => await this.highlightSegment(segmentKey)
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

    async removeSegmentAnnotation(segmentId: number, segmentationId: string, kind: 'lattice' | 'mesh' | 'primitive') {
        const targetAnnotation = this.metadata.value!.getSegmentAnnotation(segmentId, segmentationId, kind);
        this.api.removeSegmentAnnotationsUrl(this.source, this.entryId, [targetAnnotation!.id]);
        this.metadata.value!.removeSegmentAnnotation(targetAnnotation!.id);
    }

    async removeDescription(descriptionId: string) {
        this.api.removeDescriptionsUrl(this.source, this.entryId, [descriptionId]);
        this.metadata.value!.removeDescription(descriptionId);
    }

    async editDescriptions(descriptionData: DescriptionData[]) {
        await this.api.editDescriptionsUrl(this.source, this.entryId, descriptionData);
        // TODO: add to metadata
        // or fetch it again?
        // TODO: trigger re-rendering of UI
    }

    async editSegmentAnnotations(segmentAnnotationData: SegmentAnnotationData[]) {
        await this.api.editSegmentAnnotationsUrl(this.source, this.entryId, segmentAnnotationData);
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

    async _loadRawMeshSegmentationData(timeframe: number, segmentationId: string) {
        const segmentsData: RawMeshSegmentData[] = [];
        // need to have segmentsToCreate for given segmentationId and timeframe index
        const segmentsToCreate = this.metadata.value!.getMeshSegmentIdsForSegmentationIdAndTimeframe(segmentationId, timeframe);
        for (const seg of segmentsToCreate) {
            const detail = this.metadata.value!.getSufficientMeshDetail(segmentationId, timeframe, seg, DEFAULT_MESH_DETAIL);
            const urlString = this.api.meshUrl_Bcif(this.source, this.entryId, segmentationId, timeframe, seg, detail);
            const data = await this._resolveBinaryUrl(urlString);
            segmentsData.push(
                new RawMeshSegmentData(
                    seg,
                    data
                )
            );
        }
        return new RawSegmentationData(
            timeframe,
            segmentationId,
            segmentsData
        );
    }

    async _loadGeometricSegmentationData(timeframe: number, segmentationId: string) {
        // const primitivesData: ShapePrimitiveData
        const url = this.api.geometricSegmentationUrl(this.source, this.entryId, segmentationId, timeframe);
        const primitivesData = await this._resolveStringUrl(url);
        const parsedData: ShapePrimitiveData = JSON.parse(primitivesData);
        console.log('parsedData', parsedData);
        return new RawSegmentationData(
            timeframe,
            segmentationId,
            parsedData
        );
        // return parsedData;
    }

    async _loadRawChannelData(timeframe: number, channelId: string) {
        const urlString = this.api.volumeUrl(this.source, this.entryId, timeframe, channelId, BOX, MAX_VOXELS);
        // const url = Asset.getUrlAsset(this.plugin.managers.asset, urlString);
        // const asset = this.plugin.managers.asset.resolve(url, 'binary');
        const data = await this._resolveBinaryUrl(urlString);
        return new RawChannelData(
            timeframe,
            channelId,
            data
        );
    }

    async _loadRawLatticeSegmentationData(timeframe: number, segmentationId: string) {
        const urlString = this.api.latticeUrl(this.source, this.entryId, segmentationId, timeframe, BOX, MAX_VOXELS);
        const data = await this._resolveBinaryUrl(urlString);
        return new RawSegmentationData(
            timeframe,
            segmentationId,
            data
        );
    }

    private async loadRawChannelsData(timeInfo: TimeInfo, channelIds: string[]) {
        // TODO: make it run in parallel without await
        // TODO: or to run it in ctx?

        const start = timeInfo.start;
        const end = timeInfo.end;
        for (let i = start; i <= end; i++) {
            for (const channelId of channelIds) {
                const rawChannelData = await this._loadRawChannelData(i, channelId);
                this.cachedVolumeTimeframesData.add(
                    rawChannelData
                );
            }
        }
    }

    private async loadRawLatticeSegmentationData(timeInfoMapping: { [segmentation_id: string]: TimeInfo }, segmentationIds: string[]) {
        // TODO: make it run in parallel without await
        // TODO: or to run it in ctx?
        for (const segmentationId of segmentationIds) {
            const timeInfo = timeInfoMapping[segmentationId];
            const start = timeInfo.start;
            const end = timeInfo.end;
            for (let i = start; i <= end; i++) {
                const rawLatticeSegmentationData = await this._loadRawLatticeSegmentationData(i, segmentationId);
                // channelsData.push(rawChannelData);
                this.cachedSegmentationTimeframesData.add(
                    rawLatticeSegmentationData
                );
            }
        }
    }

    private async loadRawMeshSegmentationData(timeInfoMapping: { [segmentation_id: string]: TimeInfo }, segmentationIds: string[]) {
        for (const segmentationId of segmentationIds) {
            const timeInfo = timeInfoMapping[segmentationId];
            const start = timeInfo.start;
            const end = timeInfo.end;
            for (let i = start; i <= end; i++) {
                const rawMeshSegmentationData = await this._loadRawMeshSegmentationData(i, segmentationId);
                this.cachedMeshesTimeframesData.add(
                    rawMeshSegmentationData
                );
            }
        }
    }

    private async loadRawShapePrimitiveData(timeInfoMapping: { [segmentation_id: string]: TimeInfo }, segmentationIds: string[]) {
        for (const segmentationId of segmentationIds) {
            const timeInfo = timeInfoMapping[segmentationId];
            const start = timeInfo.start;
            const end = timeInfo.end;
            for (let i = start; i <= end; i++) {
                const shapePrimitiveData = await this._loadGeometricSegmentationData(i, segmentationId);
                this.cachedShapePrimitiveData.add(
                    shapePrimitiveData
                );
            }
        }
    }

    async preloadVolumeTimeframesData() {
        const timeInfo = this.metadata.value!.raw.grid.volumes.time_info;
        const channelIds = this.metadata.value!.raw.grid.volumes.channel_ids;
        this.loadRawChannelsData(timeInfo, channelIds);
        console.log('cachedVolumeTimeframesData');
        console.log(this.cachedVolumeTimeframesData);
    }

    preloadSegmentationTimeframesDataFromFile() {
        // const latticesMeta = this.metadata.value!.hasLatticeSegmentations();
        if (this.filesData!.latticeSegmentations) {
            for (const v of this.filesData!.latticeSegmentations) {
                const rawLatticeSegmentationData = new RawSegmentationData(
                    v.timeframeIndex, v.segmentationId, v.data
                );
                // channelsData.push(rawChannelData);
                this.cachedSegmentationTimeframesData.add(
                    rawLatticeSegmentationData
                );
                console.log('cachedSegmentationTimeframesData');
                console.log(this.cachedSegmentationTimeframesData);
            }
        } else {
            console.log('No segmentation data for this entry');
        }
    }

    preloadShapePrimitivesTimeframesDataFromFile() {
        // 
        if (this.filesData!.geometricSegmentations) {
            // if (this.metadata.value!.raw.grid.geometric_segmentation && this.metadata.value!.raw.grid.geometric_segmentation.segmentation_ids.length > 0) {
            for (const g of this.filesData!.geometricSegmentations) {
                const shapePrimitiveData = new RawSegmentationData(
                    g.timeframeIndex, g.segmentationId, g.data
                );
                debugger;
                // channelsData.push(rawChannelData);
                this.cachedShapePrimitiveData.add(
                    shapePrimitiveData
                );
            }
            console.log('cachedShapePrimitiveData');
            console.log(this.cachedShapePrimitiveData);
        } else {
            console.log('No shape primitive data for this entry');
        }
    }

    preloadMeshTimeframesDataFromFile() {
        debugger;
        if (this.filesData!.meshSegmentations) {
            for (const segmentationData of this.filesData!.meshSegmentations) {
                const segmentsData: RawMeshSegmentData[] = [];
                const rawData = segmentationData.data;
                for (const [filename, d] of rawData) {
                    const segmentId = parseInt((filename as string).split('_')[1]);
                    segmentsData.push(
                        new RawMeshSegmentData(
                            segmentId,
                            d
                        )
                    );
                }
                const data = new RawSegmentationData(
                    segmentationData.timeframeIndex,
                    segmentationData.segmentationId,
                    segmentsData
                );
                this.cachedMeshesTimeframesData.add(
                    data
                );
                console.log('cachedMeshesTimeframesData');
                console.log(this.cachedMeshesTimeframesData);
            }
        } else {
            console.log('No mesh segmentation data for this entry');
        }
    }

    preloadVolumeTimeframesDataFromFile() {
        // need to iterate over all volumes
        for (const v of this.filesData!.volumes!) {
            const rawChannelData = new RawChannelData(v.timeframeIndex, v.channelId, v.data);
            this.cachedVolumeTimeframesData.add(
                rawChannelData
            );
        }
        console.log('cachedVolumeTimeframesData');
        console.log(this.cachedVolumeTimeframesData);

    }

    async preloadSegmentationTimeframesData() {
        if (this.metadata.value!.raw.grid.segmentation_lattices) {
            const segmentationIds = this.metadata.value!.raw.grid.segmentation_lattices.segmentation_ids;
            const timeInfoMapping = this.metadata.value!.raw.grid.segmentation_lattices.time_info;
            this.loadRawLatticeSegmentationData(timeInfoMapping, segmentationIds);
            console.log('cachedSegmentationTimeframesData');
            console.log(this.cachedSegmentationTimeframesData);
        } else {
            console.log('No segmentation data for this entry');
        }
    }

    async preloadMeshesTimeframesData() {
        if (this.metadata.value!.raw.grid.segmentation_meshes) {
            const segmentationIds = this.metadata.value!.raw.grid.segmentation_meshes.segmentation_ids;
            const timeInfoMapping = this.metadata.value!.raw.grid.segmentation_meshes.time_info;
            this.loadRawMeshSegmentationData(timeInfoMapping, segmentationIds);
            console.log('cachedMeshesTimeframesData');
            console.log(this.cachedMeshesTimeframesData);
        } else {
            console.log('No mesh segmentation data for this entry');
        }
    }

    async preloadShapePrimitivesTimeframesData() {
        if (this.metadata.value!.raw.grid.geometric_segmentation) {
            const segmentationIds = this.metadata.value!.raw.grid.geometric_segmentation.segmentation_ids;
            const timeInfoMapping = this.metadata.value!.raw.grid.geometric_segmentation.time_info;
            this.loadRawShapePrimitiveData(timeInfoMapping, segmentationIds);
            console.log('cachedShapePrimitiveData');
            console.log(this.cachedShapePrimitiveData);
        } else {
            console.log('No shape primitive data for this entry');
        }
    }

    async updateProjectData(timeframeIndex: number) {
        // TODO: add meshes here as well and to state hierarchy mirror?
        this.changeCurrentTimeframe(timeframeIndex);
        const volumes = this.state.hierarchy.value!.volumes;
        const segmenations = this.state.hierarchy.value!.segmentations;
        const geometricSegmentations = this.state.hierarchy.value!.geometricSegmentations;
        const meshSegmentations = this.state.hierarchy.value!.meshSegmentations;
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
            const projectSegmentationDataTransform = s.transform.ref;
            const oldParams: ProjectLatticeSegmentationDataParamsValues = s.transform.params;
            // TODO: here get descriptions for segmentation and timeframe
            // and set segmentLabels
            const descriptionsForLattice = this.metadata.value!.getAllDescriptionsForSegmentationAndTimeframe(
                oldParams.segmentationId,
                'lattice',
                this.currentTimeframe.value
            );
            const segmentLabels = getSegmentLabelsFromDescriptions(descriptionsForLattice);
            const newParams: ProjectLatticeSegmentationDataParamsValues = {
                ...oldParams,
                segmentLabels: segmentLabels,
                timeframeIndex: timeframeIndex
            };
            await this.plugin.state.updateTransform(this.plugin.state.data, projectSegmentationDataTransform, newParams, 'Project Data Transform');
        }

        for (const s of geometricSegmentations) {
            const transform = s.transform.ref;
            const oldParams: ProjectGeometricSegmentationDataParamsValues = s.transform.params;
            // TODO: here get descriptions for segmentation and timeframe
            // and set segmentLabels
            // const descriptionsForLattice = this.metadata.value!.getAllDescriptionsForSegmentationAndTimeframe(
            //     oldParams.segmentationId,
            //     'lattice',
            //     this.currentTimeframe.value
            // );
            // const segmentLabels = getSegmentLabelsFromDescriptions(descriptionsForLattice);
            // ;
            const newParams: ProjectGeometricSegmentationDataParamsValues = {
                ...oldParams,
                // segmentLabels: segmentLabels,
                timeframeIndex: timeframeIndex
            };
            await this.plugin.state.updateTransform(this.plugin.state.data, transform, newParams, 'Project Data Transform');
        }

        for (const m of meshSegmentations) {
            const transform = m.transform.ref;
            const oldParams: ProjectMeshSegmentationDataParamsValues = m.transform.params;
            const newParams: ProjectMeshSegmentationDataParamsValues = {
                ...oldParams,
                // segmentLabels: segmentLabels,
                timeframeIndex: timeframeIndex
            };
            await this.plugin.state.updateTransform(this.plugin.state.data, transform, newParams, 'Project Data Transform');
        }
    }

    // async loadSegmentations() {
    //     await this.latticeSegmentationData.loadSegmentation();
    //     await this.meshSegmentationData.loadSegmentation();
    //     await this.actionShowSegments(this.metadata.value!.allSegmentIds);
    // }

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

    actionHighlightSegment(segmentKey?: string) {
        this.highlightRequest.next(segmentKey);
    }

    // async actionToggleSegment(segmentKey: string) {
    //     const current = this.currentState.value.visibleSegments.map(seg => seg.segmentKey);
    //     if (current.includes(segmentKey)) {
    //         await this.actionShowSegments(current.filter(s => s !== segmentKey));
    //     } else {
    //         await this.actionShowSegments([...current, segmentKey]);
    //     }
    // }

    // async actionToggleAllSegments() {
    //     const currentTimeframe = this.currentTimeframe.value;
    //     const current = this.currentState.value.visibleSegments.map(seg => seg.segmentKey);
    //     if (current.length !== this.metadata.value!.getAllAnnotationsForTimeframe(currentTimeframe).length) {
    //         const allSegmentKeys = this.metadata.value!.getAllAnnotationsForTimeframe(currentTimeframe).map(a =>
    //             createSegmentKey(a.segment_id, a.segmentation_id, a.segment_kind)
    //         );
    //         await this.actionShowSegments(allSegmentKeys);
    //     } else {
    //         await this.actionShowSegments([]);
    //     }
    // }

    // async actionSelectSegment(segment?: number) {
    // async actionSelectSegment(segmentKey?: string) {
    //     if (segmentKey !== undefined && this.currentState.value.visibleSegments.find(s => s.segmentKey === segmentKey) === undefined) {
    //         // first make the segment visible if it is not
    //         await this.actionToggleSegment(segmentKey);
    //     }
    //     await this.updateStateNode({ selectedSegment: segmentKey });
    // }

    async actionSetOpacity(opacity: number, segmentationId: string, kind: 'lattice' | 'mesh' | 'primitive') {
        debugger;
        // if (opacity === this.getStateNode().obj?.data.segmentOpacity) return;
        if (kind === 'lattice') this.latticeSegmentationData.updateOpacity(opacity, segmentationId);
        else if (kind === 'mesh') this.meshSegmentationData.updateOpacity(opacity, segmentationId);
        // else if TODO: geometric segmentation
        else if (kind === 'primitive') this.geometricSegmentationData.updateOpacity(opacity, segmentationId);
        // await this.updateStateNode({ segmentOpacity: opacity });
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

    // function to use in actionShowSegments for both lattice and mesh segments
    // async _actionShowSegments(parsedSegmentKeys: ParsedSegmentKey[], existingSegmentationIds: string[], kind: 'mesh' | 'lattice' | 'primitive') {
    //     const segmentKeys = parsedSegmentKeys.filter(k => k.kind === kind);
    //     const segmentIds = segmentKeys.map(k => k.segmentId);
    //     const SegmentationIdsToSegmentIds = new Map<string, number[]>;
    //     for (const key of segmentKeys) {
    //         if (!SegmentationIdsToSegmentIds.has(key.segmentationId)) {
    //             SegmentationIdsToSegmentIds.set(key.segmentationId, [key.segmentId]);
    //         } else {
    //             // have this key, add segmentId
    //             const currentSegmentationIds = SegmentationIdsToSegmentIds.get(key.segmentationId);
    //             SegmentationIdsToSegmentIds.set(key.segmentationId, [...currentSegmentationIds!, key.segmentId]);
    //         }
    //     }
    //     for (const id of existingSegmentationIds) {
    //         if (!SegmentationIdsToSegmentIds.has(id)) {
    //             SegmentationIdsToSegmentIds.set(id, []);
    //         }
    //     }
    //     console.log(SegmentationIdsToSegmentIds);
    //     const promises: Promise<void>[] = [];
    //     SegmentationIdsToSegmentIds.forEach((value, key) => {
    //         if (kind === 'lattice') promises.push(this.latticeSegmentationData.showSegments(value, key));
    //         else if (kind === 'mesh') promises.push(this.meshSegmentationData.showSegments(value, key));
    //         else if (kind === 'primitive') promises.push(this.geometricSegmentationData.showSegments(value, key));
    //     });
    //     await Promise.all(promises);
    // }

    // async actionShowSegments(segmentKeys: string[]) {
    //     const allExistingLatticeSegmentationIds = this.metadata.value!.raw.grid.segmentation_lattices!.segmentation_ids;
    //     const allExistingMeshSegmentationIds = this.metadata.value!.raw.grid.segmentation_meshes!.segmentation_ids;
    //     const allExistingGeometricSegmentationIds = this.metadata.value!.raw.grid.geometric_segmentation!.segmentation_ids;
    //     if (segmentKeys.length === 0) {
    //         for (const id of allExistingLatticeSegmentationIds) {
    //             await this.latticeSegmentationData.showSegments([], id);
    //         }
    //         for (const id of allExistingMeshSegmentationIds) {
    //             await this.meshSegmentationData.showSegments([], id);
    //         }
    //         for (const id of allExistingGeometricSegmentationIds) {
    //             await this.geometricSegmentationData.showSegments([], id);
    //         }
    //     }
    //     const parsedSegmentKeys = segmentKeys.map(
    //         k => parseSegmentKey(k)
    //     );
    //     // LATTICES PART
    //     this._actionShowSegments(parsedSegmentKeys, allExistingLatticeSegmentationIds, 'lattice');
    //     // MESHES PART
    //     this._actionShowSegments(parsedSegmentKeys, allExistingMeshSegmentationIds, 'mesh');
    //     // GEOMETRIC SEGMENTATION PAR
    //     this._actionShowSegments(parsedSegmentKeys, allExistingGeometricSegmentationIds, 'primitive');

    //     await this.updateStateNode({ visibleSegments: segmentKeys.map(s => ({ segmentKey: s })) });
    //     console.log('Current state');
    //     console.log(this.getStateNode());
    // }

    private async highlightSegment(segmentKey?: string) {
        await PluginCommands.Interactivity.ClearHighlights(this.plugin);
        if (segmentKey) {
            const parsedSegmentKey = parseSegmentKey(segmentKey);
            const { segmentId, segmentationId, kind } = parsedSegmentKey;
            if (kind === 'lattice') {
                await this.latticeSegmentationData.highlightSegment(segmentId, segmentationId);
            } else if (kind === 'mesh') {
                await this.meshSegmentationData.highlightSegment(segmentId, segmentationId);
            } else if (kind === 'primitive') {
                await this.geometricSegmentationData.highlightSegment(segmentId, segmentationId);
            }
            // TODO: support primitive
        }
    }

    private async selectSegment(segmentKey: string) {
        this.plugin.managers.interactivity.lociSelects.deselectAll();
        // TODO: parse segmentKey first
        const parsedSegmentKey = parseSegmentKey(segmentKey);
        if (parsedSegmentKey.kind === 'lattice') {
            await this.latticeSegmentationData.selectSegment(parsedSegmentKey.segmentId, parsedSegmentKey.segmentationId);
        } else if (parsedSegmentKey.kind === 'mesh') {
            await this.meshSegmentationData.selectSegment(parsedSegmentKey.segmentId, parsedSegmentKey.segmentationId);
        } else if (parsedSegmentKey.kind === 'primitive') {
            await this.geometricSegmentationData.selectSegment(parsedSegmentKey.segmentId, parsedSegmentKey.segmentationId);
        }

        // TODO: primitives
        // await this.latticeSegmentationData.selectSegment(segment);
        // await this.meshSegmentationData.selectSegment(segment);
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

    // Fix this
    private readonly labelProvider: LociLabelProvider = {
        label: (loci: Loci): string | undefined => {
            const segmentId = this.getSegmentIdFromLoci(loci);
            const segmentationId = this.getSegmentationIdFromLoci(loci);
            const segmentationKind = this.getSegmentationKindFromLoci(loci);
            if (segmentId === undefined || !segmentationId || !segmentationKind) return;
            // const segmentKey = createSegmentKey(segmentId, segmentationId, segmentationKind);
            const descriptions = this.metadata.value!.getSegmentDescription(segmentId, segmentationId, segmentationKind);
            if (!descriptions) return;

            // theoretically there can be multiple descriptions
            // for each segment
            // first try assuming there is a single one
            const annotLabels = descriptions[0].external_references?.map(e => `${applyEllipsis(e.label ? e.label : '')} [${e.resource}:${e.accession}]`);
            // TODO: try rendering multiple descriptions
            if (!annotLabels || annotLabels.length === 0 || descriptions[0].is_hidden === true) return;
            // if (annotLabels.length > MAX_ANNOTATIONS_IN_LABEL + 1) {
            //     const nHidden = annotLabels.length - MAX_ANNOTATIONS_IN_LABEL;
            //     annotLabels.length = MAX_ANNOTATIONS_IN_LABEL;
            //     annotLabels.push(`(${nHidden} more annotations, click on the segment to see all)`);
            // }
            return '<hr class="msp-highlight-info-hr"/>' + annotLabels.filter(isDefined).join('<br/>');
        }
    };

    private getSegmentationIdFromLoci(loci: Loci): string | undefined {
        if (Volume.Segment.isLoci(loci)) {
            return loci.volume.label;
        } else if (ShapeGroup.isLoci(loci)) {
            const sourceData = loci.shape.sourceData;
            debugger;
            // as any?
            if (isMeshlistData(sourceData as any)) {
                const meshData = (loci.shape.sourceData ?? {}) as MeshlistData;
                return meshData.segmentationId!;
            } else if (isShapePrimitiveParamsValues(sourceData as any)) {
                const shapePrimitiveParamsValues = (loci.shape.sourceData ?? {}) as CreateShapePrimitiveProviderParamsValues;
                return shapePrimitiveParamsValues.segmentationId;
            }
        }
    }

    private getSegmentationKindFromLoci(loci: Loci): 'lattice' | 'mesh' | 'primitive' | undefined {
        if (Volume.Segment.isLoci(loci)) {
            return 'lattice';
        } else if (ShapeGroup.isLoci(loci)) {
            // return 'mesh';
            const sourceData = loci.shape.sourceData;
            if (isMeshlistData(sourceData as any)) {
                const meshData = (loci.shape.sourceData ?? {}) as MeshlistData;
                return 'mesh';
            } else if (isShapePrimitiveParamsValues(sourceData as any)) {
                const shapePrimitiveParamsValues = (loci.shape.sourceData ?? {}) as CreateShapePrimitiveProviderParamsValues;
                return 'primitive';
            }
        } else {
            // TODO: fix in case of isosurface loci
            console.log(`Segmentation kind is not supported for ${loci}`);
        }
    }

    private getSegmentIdFromLoci(loci: Loci): number | undefined {
        if (Volume.Segment.isLoci(loci) && loci.volume._propertyData.ownerId === this.ref) {
            if (loci.segments.length === 1) {
                return loci.segments[0];
            }
        }
        if (ShapeGroup.isLoci(loci)) {
            const sourceData = loci.shape.sourceData;
            if (isMeshlistData(sourceData as any)) {
                const meshData = (loci.shape.sourceData ?? {}) as MeshlistData;
                if (meshData.ownerId === this.ref && meshData.segmentId !== undefined) {
                    return meshData.segmentId;
                }
                // TODO: check for ownerId? this would be entry root
            } else if (isShapePrimitiveParamsValues(sourceData as any)) {
                const shapePrimitiveParamsValues = (loci.shape.sourceData ?? {}) as CreateShapePrimitiveProviderParamsValues;
                return shapePrimitiveParamsValues.segmentId;
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



