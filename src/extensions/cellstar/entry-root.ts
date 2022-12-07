import { BehaviorSubject, Subject, throttleTime } from 'rxjs';

import { ShapeGroup } from '../../mol-model/shape';
import { Volume } from '../../mol-model/volume';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { PluginBehavior } from '../../mol-plugin/behavior';
import { PluginCommands } from '../../mol-plugin/commands';
import { PluginContext } from '../../mol-plugin/context';
import { StateObjectCell, StateTransform } from '../../mol-state';
import { shallowEqualObjects } from '../../mol-util';
import { ParamDefinition } from '../../mol-util/param-definition';
import { sleep } from '../../mol-util/sleep';
import { MeshlistData } from '../meshes/mesh-extension';

import { DEFAULT_VOLUME_SERVER_V2, VolumeApiV2 } from './cellstar-api/api';
import { Metadata, Segment } from './cellstar-api/data';
import { CellStarMeshSegmentationData } from './entry-meshes';
import { CellStarModelData } from './entry-models';
import { CellStarLatticeSegmentationData } from './entry-segmentation';
import { CellStarState, CellStarStateData, CellStarStateParams } from './entry-state';
import { CellStarVolumeData } from './entry-volume';
import * as ExternalAPIs from './external-api';
import { Choice, createEntryId, lazyGetter, NodeManager } from './helpers';
import { type CellStarStateFromEntry } from './transformers';


export const MAX_VOXELS = 10 ** 7;
// export const MAX_VOXELS = 10 ** 2; // DEBUG
export const BOX: [[number, number, number], [number, number, number]] | null = null;
// export const BOX: [[number, number, number], [number, number, number]] | null = [[-90, -90, -90], [90, 90, 90]]; // DEBUG


const SourceChoice = new Choice({ emdb: 'EMDB', empiar: 'EMPIAR' }, 'emdb');
type Source = Choice.Values<typeof SourceChoice>


export const CellStarEntryParams = {
    serverUrl: ParamDefinition.Text(DEFAULT_VOLUME_SERVER_V2),
    source: SourceChoice.PDSelect(),
    entryNumber: ParamDefinition.Text('1832'),
};
type CellStarEntryParamValues = ParamDefinition.Values<typeof CellStarEntryParams>;


export class CellStarEntry extends PluginStateObject.CreateBehavior<CellStarEntryData>({ name: 'CellStar Entry' }) { }


export class CellStarEntryData extends PluginBehavior.WithSubscribers<CellStarEntryParamValues> {
    plugin: PluginContext;
    ref: string = '';
    api: VolumeApiV2;
    source: Source;
    /** Number part of entry ID; e.g. '1832' */
    entryNumber: string;
    /** Full entry ID; e.g. 'emd-1832' */
    entryId: string;
    metadata: Metadata;
    pdbs: string[];

    public readonly groupNodeMgr = new NodeManager();
    public readonly volumeData = new CellStarVolumeData(this);
    public readonly latticeSegmentationData = new CellStarLatticeSegmentationData(this);
    public readonly meshSegmentationData = new CellStarMeshSegmentationData(this);
    public readonly modelData = new CellStarModelData(this);
    currentSegment = new BehaviorSubject<Segment | undefined>(undefined);
    visibleSegments = new BehaviorSubject<Segment[]>([]);
    private highlightRequest = new Subject<Segment | undefined>();

    private getStateNode = lazyGetter(() => this.plugin.state.data.selectQ(q => q.byRef(this.ref).subtree().ofType(CellStarState))[0] as StateObjectCell<CellStarState, StateTransform<typeof CellStarStateFromEntry>>, 'Missing CellStarState node. Must first create CellStarState for this CellStarEntry.');
    public currentState = new BehaviorSubject(ParamDefinition.getDefaultValues(CellStarStateParams));


    private constructor(plugin: PluginContext, params: CellStarEntryParamValues) {
        super(plugin, params);

        this.plugin = plugin;
        this.api = new VolumeApiV2(params.serverUrl);
        this.source = params.source;
        this.entryNumber = params.entryNumber;
        this.entryId = createEntryId(params.source, params.entryNumber);
    }


    async register(ref: string) {
        this.ref = ref;
        console.log('register', ref);
        try {
            const params = this.getStateNode().obj?.data;
            if (params) {
                this.currentState.next(params);
            }
        } catch {
            // do nothing
        }

        this.subscribeObservable(this.plugin.state.data.events.cell.stateUpdated, e => {
            try { (this.getStateNode()) } catch { return } // if state not does not exist yet
            if (e.cell.transform.ref === this.getStateNode().transform.ref) {
                const newState = this.getStateNode().obj?.data;
                if (newState && !shallowEqualObjects(newState, this.currentState.value)) { // avoid repeated update
                    console.log('stateUpdated', this.getStateNode().obj?.id, e.cell.obj?.id, this.getStateNode().obj?.data.opacity, e.cell.obj?.data.opacity);
                    this.currentState.next(newState);
                }
            }
        })

        this.subscribeObservable(this.plugin.behaviors.interaction.click, e => {
            const loci = e.current.loci;
            console.log('click', this.ref, e.current.loci, e.current.repr, e.current.repr?.state); // DEBUG

            if (Volume.Segment.isLoci(loci)) console.log('ownerId', loci.volume._propertyData.ownerId); // DEBUG
            if (Volume.Segment.isLoci(loci) && loci.volume._propertyData.ownerId === this.ref) {
                const clickedSegmentId = loci.segments.length === 1 ? loci.segments[0] : undefined;
                const clickedSegment = this.metadata.annotation?.segment_list.find(seg => seg.id === clickedSegmentId);
                this.currentSegment.next(clickedSegment);
            }
            if (ShapeGroup.isLoci(loci)) {
                const meshData = (loci.shape.sourceData ?? {}) as MeshlistData;
                if (meshData.ownerId === this.ref && meshData.segmentId !== undefined) {
                    this.currentSegment.next(this.metadata.annotation?.segment_list.find(segment => segment.id === meshData.segmentId));
                }
            }
        });
        this.subscribeObservable(this.highlightRequest.pipe(throttleTime(50, undefined, { leading: true, trailing: true })),
            async segment => {
                await PluginCommands.Interactivity.ClearHighlights(this.plugin);
                if (segment) {
                    await this.latticeSegmentationData.highlightSegment(segment);
                    await this.meshSegmentationData.highlightSegment(segment);
                }
            }
        );
    }

    async unregister() {
        console.log('unregister', this.ref);
    }

    private async initialize() {
        this.metadata = await this.api.getMetadata(this.source, this.entryId);
        this.pdbs = await ExternalAPIs.getPdbIdsForEmdbEntry(this.metadata.grid.general.source_db_id ?? this.entryId);
        // TODO use Asset?
    }

    static async create(plugin: PluginContext, params: CellStarEntryParamValues) {
        const result = new CellStarEntryData(plugin, params);
        await result.initialize();
        return result;
    }

    async updateOpacity(opacity: number) {
        const stateNode = this.getStateNode().obj?.data;
        console.log('updateOpacity', stateNode?.opacity, '->', opacity);
        if (opacity === stateNode?.opacity) return;

        await sleep(1000); // TODO remove
        this.latticeSegmentationData.updateOpacity(opacity);
        this.meshSegmentationData.updateOpacity(opacity);

        this.updateStateNode({ opacity: opacity });
    }

    private async updateStateNode(params: Partial<CellStarStateData>) {
        const state = this.plugin.state.data;
        const update = state.build().to(this.getStateNode().transform.ref).update(params);
        await PluginCommands.State.Update(this.plugin, { state, tree: update, options: { doNotUpdateCurrent: true } });
    }

    public newUpdate() {
        if (this.ref !== '') {
            return this.plugin.build().to(this.ref);
        } else {
            return this.plugin.build().toRoot();
        }
    }

    public async showSegmentations() {
        await this.latticeSegmentationData.showSegmentation();
        await this.meshSegmentationData.showSegmentation();
        this.visibleSegments.next(this.metadata.annotation?.segment_list ?? []);
    }

    highlightSegment(segment?: Segment) {
        this.highlightRequest.next(segment);
    }

    async toggleSegment(segment: Segment) {
        const current = this.visibleSegments.value;
        if (current.includes(segment)) {
            this.showSegments(current.filter(s => s !== segment));
        } else {
            this.showSegments([...current, segment]);
        }
    }

    toggleAllSegments() {
        const current = this.visibleSegments.value;
        if (current.length !== (this.metadata.annotation?.segment_list.length ?? 0)) {
            this.showSegments(this.metadata.annotation?.segment_list ?? []);
        } else {
            this.showSegments([]);
        }
    }

    showAnnotation(segment: Segment) {
        this.currentSegment.next(segment);
    }
    public async showSegments(segments: Segment[]) {
        console.log('showSegments', segments.map(seg => seg.id), this.ref);
        await this.latticeSegmentationData.showSegments(segments, { opacity: this.currentState.value.opacity });
        await this.meshSegmentationData.showSegments(segments);
        this.visibleSegments.next(segments);
    }
}



