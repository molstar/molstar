/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { BehaviorSubject, distinctUntilChanged, Subject, throttleTime } from 'rxjs';
import { VolsegVolumeServerConfig } from '.';
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
import { Segment } from './volseg-api/data';
import { MetadataWrapper } from './volseg-api/utils';
import { VolsegMeshSegmentationData } from './entry-meshes';
import { VolsegModelData } from './entry-models';
import { VolsegLatticeSegmentationData } from './entry-segmentation';
import { VolsegState, VolsegStateData, VolsegStateParams } from './entry-state';
import { VolsegVolumeData, SimpleVolumeParamValues, VOLUME_VISUAL_TAG } from './entry-volume';
import * as ExternalAPIs from './external-api';
import { VolsegGlobalStateData } from './global-state';
import { applyEllipsis, isDefined, lazyGetter, splitEntryId } from './helpers';
import { type VolsegStateFromEntry } from './transformers';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { OrderedSet } from '../../mol-data/int';


export const MAX_VOXELS = 10 ** 7;
// export const MAX_VOXELS = 10 ** 2; // DEBUG
export const BOX: [[number, number, number], [number, number, number]] | null = null;
// export const BOX: [[number, number, number], [number, number, number]] | null = [[-90, -90, -90], [90, 90, 90]]; // DEBUG

const MAX_ANNOTATIONS_IN_LABEL = 6;


const SourceChoice = new Choice({ emdb: 'EMDB', empiar: 'EMPIAR', idr: 'IDR' }, 'emdb');
export type Source = Choice.Values<typeof SourceChoice>;


export function createLoadVolsegParams(plugin?: PluginContext, entrylists: { [source: string]: string[] } = {}) {
    const defaultVolumeServer = plugin?.config.get(VolsegVolumeServerConfig.DefaultServer) ?? DEFAULT_VOLSEG_SERVER;
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

export function createVolsegEntryParams(plugin?: PluginContext) {
    const defaultVolumeServer = plugin?.config.get(VolsegVolumeServerConfig.DefaultServer) ?? DEFAULT_VOLSEG_SERVER;
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

type VolRepr3DT = typeof StateTransforms.Representation.VolumeRepresentation3D

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

    public readonly volumeData = new VolsegVolumeData(this);
    private readonly latticeSegmentationData = new VolsegLatticeSegmentationData(this);
    private readonly meshSegmentationData = new VolsegMeshSegmentationData(this);
    private readonly modelData = new VolsegModelData(this);
    private highlightRequest = new Subject<Segment | undefined>();

    private getStateNode = lazyGetter(() => this.plugin.state.data.selectQ(q => q.byRef(this.ref).subtree().ofType(VolsegState))[0] as StateObjectCell<VolsegState, StateTransform<typeof VolsegStateFromEntry>>, 'Missing VolsegState node. Must first create VolsegState for this VolsegEntry.');
    public currentState = new BehaviorSubject(ParamDefinition.getDefaultValues(VolsegStateParams));
    public currentVolume = new BehaviorSubject<StateTransform<VolRepr3DT> | undefined>(undefined);


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
        this.pdbs = await ExternalAPIs.getPdbIdsForEmdbEntry(this.metadata.raw.grid.general.source_db_id ?? this.entryId);
        // TODO use Asset?
    }

    static async create(plugin: PluginContext, params: VolsegEntryParamValues) {
        const result = new VolsegEntryData(plugin, params);
        await result.initialize();
        return result;
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
        if (volumeVisual) this.currentVolume.next(volumeVisual.transform);

        let volumeRef: string | undefined;
        this.subscribeObservable(this.plugin.state.data.events.cell.stateUpdated, e => {
            try { (this.getStateNode()); } catch { return; } // if state not does not exist yet
            if (e.cell.transform.ref === this.getStateNode().transform.ref) {
                const newState = this.getStateNode().obj?.data;
                if (newState && !shallowEqualObjects(newState, this.currentState.value)) { // avoid repeated update
                    this.currentState.next(newState);
                }
            } else if (e.cell.transform.tags?.includes(VOLUME_VISUAL_TAG)) {
                if (e.ref === volumeRef) {
                    this.currentVolume.next(e.cell.transform);
                } else if (StateSelection.findAncestor(this.plugin.state.data.tree, this.plugin.state.data.cells, e.ref, a => a.transform.ref === ref)) {
                    volumeRef = e.ref;
                    this.currentVolume.next(e.cell.transform);
                }
            }
        });

        this.subscribeObservable(this.plugin.state.data.events.cell.removed, e => {
            if (e.ref === volumeRef) {
                volumeRef = undefined;
                this.currentVolume.next(undefined);
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
        const result = await this.volumeData.loadVolume();
        if (result) {
            const isovalue = result.isovalue.kind === 'relative' ? result.isovalue.relativeValue : result.isovalue.absoluteValue;
            await this.updateStateNode({ volumeIsovalueKind: result.isovalue.kind, volumeIsovalueValue: isovalue });
        }
    }

    async loadSegmentations() {
        await this.latticeSegmentationData.loadSegmentation();
        await this.meshSegmentationData.loadSegmentation();
        await this.actionShowSegments(this.metadata.allSegmentIds);
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

    async actionSetVolumeVisual(type: 'isosurface' | 'direct-volume' | 'off') {
        await this.volumeData.setVolumeVisual(type);
        await this.updateStateNode({ volumeType: type });
    }

    async actionUpdateVolumeVisual(params: SimpleVolumeParamValues) {
        await this.volumeData.updateVolumeVisual(params);
        await this.updateStateNode({
            volumeType: params.volumeType,
            volumeOpacity: params.opacity,
        });
    }


    private async actionShowSegments(segments: number[]) {
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

    private async updateStateNode(params: Partial<VolsegStateData>) {
        const oldParams = this.getStateNode().transform.params;
        const newParams = { ...oldParams, ...params };
        const state = this.plugin.state.data;
        const update = state.build().to(this.getStateNode().transform.ref).update(newParams);
        await PluginCommands.State.Update(this.plugin, { state, tree: update, options: { doNotUpdateCurrent: true } });
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
            if (loci.elements.length === 1 && OrderedSet.size(loci.elements[0].segments) === 1) {
                return OrderedSet.start(loci.elements[0].segments);
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



