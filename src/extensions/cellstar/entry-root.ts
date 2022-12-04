import { BehaviorSubject, Subject, throttleTime } from 'rxjs';

import { ShapeGroup } from '../../mol-model/shape';
import { Volume } from '../../mol-model/volume';
import { PluginComponent } from '../../mol-plugin-state/component';
import { PluginStateObject as SO, PluginStateTransform } from '../../mol-plugin-state/objects';
import { PluginCommands } from '../../mol-plugin/commands';
import { PluginContext } from '../../mol-plugin/context';
import { StateObjectSelector } from '../../mol-state';
import { Task } from '../../mol-task';
import { ParamDefinition } from '../../mol-util/param-definition';
import { MeshlistData } from '../meshes/mesh-extension';

import { DEFAULT_VOLUME_SERVER_V2, VolumeApiV2 } from './cellstar-api/api';
import { Metadata, Segment } from './cellstar-api/data';
import { CellStarMeshSegmentationData } from './entry-meshes';
import { CellStarModelData } from './entry-models';
import { CellStarLatticeSegmentationData } from './entry-segmentation';
import { CellStarLatticeSegmentationData2 } from './entry-segmentation-2';
import { CellStarVolumeData } from './entry-volume';
import * as ExternalAPIs from './external-api';
import { Choice, createEntryId, NodeManager } from './helpers';


export const MAX_VOXELS = 10 ** 7;
// export const MAX_VOXELS = 10**2; // DEBUG


const SourceChoice = new Choice({ emdb: 'EMDB', empiar: 'EMPIAR' }, 'emdb');
type Source = Choice.Values<typeof SourceChoice>


export const CellStarEntryParams = {
    serverUrl: ParamDefinition.Text(DEFAULT_VOLUME_SERVER_V2),
    source: SourceChoice.PDSelect(),
    entryNumber: ParamDefinition.Text('1832'),
};

export class CellStarEntryData extends PluginComponent {
    plugin: PluginContext;
    entryRoot?: StateObjectSelector;
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
    public readonly latticeSegmentationData2 = new CellStarLatticeSegmentationData2(this);
    public readonly meshSegmentationData = new CellStarMeshSegmentationData(this);
    public readonly modelData = new CellStarModelData(this);
    currentSegment = new BehaviorSubject<Segment | undefined>(undefined);
    visibleSegments = new BehaviorSubject<Segment[]>([]);
    opacity = new BehaviorSubject(1);
    private highlightRequest = new Subject<Segment | undefined>();


    private constructor(plugin: PluginContext, serverUrl: string, source: Source, entryNumber: string) {
        super();

        this.plugin = plugin;
        this.api = new VolumeApiV2(serverUrl);
        this.source = source;
        this.entryNumber = entryNumber;
        this.entryId = createEntryId(source, entryNumber);

        this.subscribe(plugin.behaviors.interaction.click, e => {
            const loci = e.current.loci;
            if (Volume.isLoci(loci) && loci.volume._propertyData.ownerId === this.entryRoot?.ref && loci.volume._propertyData.segment) {
                this.currentSegment.next(loci.volume._propertyData.segment); // TODO remove with old latticeSegmentationData
            }
            if (Volume.Segment.isLoci(loci) && loci.volume._propertyData.ownerId === this.entryRoot?.ref) {
                const clickedSegmentId = loci.segments.length === 1 ? loci.segments[0] : undefined;
                const clickedSegment = this.metadata.annotation?.segment_list.find(seg => seg.id === clickedSegmentId);
                this.currentSegment.next(clickedSegment);
            }
            if (ShapeGroup.isLoci(loci)) {
                const meshData = (loci.shape.sourceData ?? {}) as MeshlistData;
                if (meshData.ownerId === this.entryRoot?.ref && meshData.segmentId !== undefined) {
                    this.currentSegment.next(this.metadata.annotation?.segment_list.find(segment => segment.id === meshData.segmentId));
                }
            }
        });
        this.subscribe(this.opacity, opacity => {
            this.latticeSegmentationData.updateOpacity(opacity);
            this.latticeSegmentationData2.updateOpacity(opacity);
            this.meshSegmentationData.updateOpacity(opacity);
        });
        this.subscribe(this.highlightRequest.pipe(throttleTime(50, undefined, { leading: true, trailing: true })),
            async segment => {
                await PluginCommands.Interactivity.ClearHighlights(this.plugin);
                if (segment) {
                    await this.latticeSegmentationData.highlightSegment(segment);
                    await this.latticeSegmentationData2.highlightSegment(segment);
                    await this.meshSegmentationData.highlightSegment(segment);
                }
            }
        );
    }

    private async initialize() {
        this.metadata = await this.api.getMetadata(this.source, this.entryId);
        this.pdbs = await ExternalAPIs.getPdbIdsForEmdbEntry(this.metadata.grid.general.source_db_id ?? this.entryId);
    }

    static async create(plugin: PluginContext, serverUrl: string, source: Source, entryNumber: string) {
        const result = new CellStarEntryData(plugin, serverUrl, source, entryNumber);
        await result.initialize();
        return result;
    }

    public newUpdate() {
        if (this.entryRoot) {
            return this.plugin.build().to(this.entryRoot);
        } else {
            return this.plugin.build().toRoot();
        }
    }

    public async showSegmentations() {
        await this.latticeSegmentationData.showSegmentation();
        await this.latticeSegmentationData2.showSegmentation();
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
        await this.latticeSegmentationData.showSegments(segments, { opacity: this.opacity.value });
        await this.latticeSegmentationData2.showSegments(segments, { opacity: this.opacity.value });
        await this.meshSegmentationData.showSegments(segments);
        this.visibleSegments.next(segments);
    }
}


export class CellStarEntry extends SO.Create<CellStarEntryData>({ name: 'CellStar Entry', typeClass: 'Object' }) { }


export const CellStarEntryFromRoot = PluginStateTransform.BuiltIn({
    name: 'cellstar-entry-from-root',
    display: { name: 'CellStar Entry', description: 'Load CellStar Entry' },
    from: SO.Root,
    to: CellStarEntry,
    params: CellStarEntryParams,
})({
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Load CellStar Entry', async () => {
            const data = await CellStarEntryData.create(plugin, params.serverUrl, params.source, params.entryNumber);
            return new CellStarEntry(data, { label: data.entryId, description: 'CellStar Entry' });
        });
    }
});