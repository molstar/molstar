import { BehaviorSubject, distinctUntilChanged, Subject, throttleTime } from 'rxjs';
import { CellstarVolumeServerConfig } from '.';
import { Loci } from '../../mol-model/loci';

import { ShapeGroup } from '../../mol-model/shape';
import { Volume } from '../../mol-model/volume';
import { LociLabelProvider } from '../../mol-plugin-state/manager/loci-label';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { PluginBehavior } from '../../mol-plugin/behavior';
import { PluginCommands } from '../../mol-plugin/commands';
import { PluginContext } from '../../mol-plugin/context';
import { StateObjectCell, StateTransform } from '../../mol-state';
import { shallowEqualObjects } from '../../mol-util';
import { ParamDefinition } from '../../mol-util/param-definition';
import { MeshlistData } from '../meshes/mesh-extension';

import { DEFAULT_VOLUME_SERVER_V2, VolumeApiV2 } from './cellstar-api/api';
import { Segment } from './cellstar-api/data';
import { MetadataWrapper } from './cellstar-api/utils';
import { CellstarMeshSegmentationData } from './entry-meshes';
import { CellstarModelData } from './entry-models';
import { CellstarLatticeSegmentationData } from './entry-segmentation';
import { CellstarState, CellstarStateData, CellstarStateParams } from './entry-state';
import { CellstarVolumeData } from './entry-volume';
import * as ExternalAPIs from './external-api';
import { applyEllipsis, Choice, isDefined, lazyGetter, splitEntryId } from './helpers';
import { type CellstarStateFromEntry } from './transformers';


export const MAX_VOXELS = 10 ** 7;
// export const MAX_VOXELS = 10 ** 2; // DEBUG
export const BOX: [[number, number, number], [number, number, number]] | null = null;
// export const BOX: [[number, number, number], [number, number, number]] | null = [[-90, -90, -90], [90, 90, 90]]; // DEBUG

const MAX_ANNOTATIONS_IN_LABEL = 6;


const SourceChoice = new Choice({ emdb: 'EMDB', empiar: 'EMPIAR', idr: 'IDR' }, 'emdb');
export type Source = Choice.Values<typeof SourceChoice>;


export function createLoadCellstarParams(plugin?: PluginContext, entrylists: { [source: string]: string[] } = {}) {
    console.log('createLoadCellstarParams', entrylists);
    const defaultVolumeServer = plugin?.config.get(CellstarVolumeServerConfig.DefaultServer) ?? DEFAULT_VOLUME_SERVER_V2;
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
type LoadCellstarParamValues = ParamDefinition.Values<ReturnType<typeof createLoadCellstarParams>>;

export function createCellstarEntryParams(plugin?: PluginContext) {
    const defaultVolumeServer = plugin?.config.get(CellstarVolumeServerConfig.DefaultServer) ?? DEFAULT_VOLUME_SERVER_V2;
    return {
        serverUrl: ParamDefinition.Text(defaultVolumeServer),
        source: SourceChoice.PDSelect(),
        entryId: ParamDefinition.Text('emd-1832', { description: 'Entry identifier, including the source prefix, e.g. "emd-1832"' }),
    };
}
type CellstarEntryParamValues = ParamDefinition.Values<ReturnType<typeof createCellstarEntryParams>>;

export namespace CellstarEntryParamValues {
    export function fromLoadCellstarParamValues(params: LoadCellstarParamValues): CellstarEntryParamValues {
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


export class CellstarEntry extends PluginStateObject.CreateBehavior<CellstarEntryData>({ name: 'Vol & Seg Entry' }) { }


export class CellstarEntryData extends PluginBehavior.WithSubscribers<CellstarEntryParamValues> {
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

    public readonly volumeData = new CellstarVolumeData(this);
    private readonly latticeSegmentationData = new CellstarLatticeSegmentationData(this);
    private readonly meshSegmentationData = new CellstarMeshSegmentationData(this);
    private readonly modelData = new CellstarModelData(this);
    private highlightRequest = new Subject<Segment | undefined>();

    private getStateNode = lazyGetter(() => this.plugin.state.data.selectQ(q => q.byRef(this.ref).subtree().ofType(CellstarState))[0] as StateObjectCell<CellstarState, StateTransform<typeof CellstarStateFromEntry>>, 'Missing CellstarState node. Must first create CellstarState for this CellstarEntry.');
    public currentState = new BehaviorSubject(ParamDefinition.getDefaultValues(CellstarStateParams));


    private constructor(plugin: PluginContext, params: CellstarEntryParamValues) {
        super(plugin, params);
        this.plugin = plugin;
        this.api = new VolumeApiV2(params.serverUrl);
        this.source = params.source;
        // this.entryNumber = params.entryNumber;
        // this.entryId = createEntryId(params.source, params.entryNumber);
        this.entryId = params.entryId;
        this.entryNumber = splitEntryId(this.entryId).entryNumber;
    }

    private async initialize() {
        const metadata = await this.api.getMetadata(this.source, this.entryId);
        this.metadata = new MetadataWrapper(metadata);
        this.pdbs = await ExternalAPIs.getPdbIdsForEmdbEntry(this.metadata.raw.grid.general.source_db_id ?? this.entryId);
        // TODO use Asset?
    }

    static async create(plugin: PluginContext, params: CellstarEntryParamValues) {
        const result = new CellstarEntryData(plugin, params);
        await result.initialize();
        return result;
    }

    async register(ref: string) {
        this.ref = ref;
        console.log('register', ref);
        this.plugin.managers.lociLabels.addProvider(this.labelProvider);

        try {
            const params = this.getStateNode().obj?.data;
            if (params) {
                this.currentState.next(params);
            }
        } catch {
            // do nothing
        }

        this.subscribeObservable(this.plugin.state.data.events.cell.stateUpdated, e => {
            try { (this.getStateNode()); } catch { return; } // if state not does not exist yet
            if (e.cell.transform.ref === this.getStateNode().transform.ref) {
                const newState = this.getStateNode().obj?.data;
                if (newState && !shallowEqualObjects(newState, this.currentState.value)) { // avoid repeated update
                    this.currentState.next(newState);
                }
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
            async state => await this.selectSegment(state.selectedSegment)
        );
    }

    async unregister() {
        console.log('unregister', this.ref);
        this.plugin.managers.lociLabels.removeProvider(this.labelProvider);
    }

    async loadVolume() {
        const result = await this.volumeData.loadVolume();
        if (result) {
            const isovalue = result.isovalue.kind === 'relative' ? result.isovalue.relativeValue : result.isovalue.absoluteValue;
            await this.updateStateNode({ volumeIsovalueKind: result.isovalue.kind, volumeIsovalueValue: isovalue })
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
        if (opacity === this.getStateNode().obj?.data.opacity) return;
        this.latticeSegmentationData.updateOpacity(opacity);
        this.meshSegmentationData.updateOpacity(opacity);

        await this.updateStateNode({ opacity: opacity });
    }

    async actionShowFittedModel(pdbIds: string[]) {
        await this.modelData.showPdbs(pdbIds);
        await this.updateStateNode({ visibleModels: pdbIds.map(pdbId => ({ pdbId: pdbId })) });
    }

    async actionSetVolumeVisual(type: 'isosurface' | 'direct-volume' | 'off') {
        await this.volumeData.setVolumeVisual(type);
        await this.updateStateNode({ volumeType: type });
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

    private async updateStateNode(params: Partial<CellstarStateData>) {
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
            const annotLabels = segment.biological_annotation.external_references.map(annot => `${applyEllipsis(annot.label)} [${annot.resource}]`);
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

}



