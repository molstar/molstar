import { PluginStateObject as SO, PluginStateTransform } from '../../mol-plugin-state/objects';
import { PluginContext } from '../../mol-plugin/context';
import { StateObjectSelector } from '../../mol-state';
import { Task } from '../../mol-task';
import { ParamDefinition } from '../../mol-util/param-definition';

import { DEFAULT_VOLUME_SERVER_V2, VolumeApiV2 } from './cellstar-api/api';
import { Metadata, Segment } from './cellstar-api/data';
import { Choice, createEntryId, NodeManager } from './helpers';
import { CellStarVolumeData } from './entry-volume';
import { CellStarLatticeSegmentationData } from './entry-segmentation';
import { CellStarModelData } from './entry-models';
import * as ExternalAPIs from './external-api';
import { CellStarMeshSegmentationData } from './entry-meshes';


export const MAX_VOXELS = 10**7;


const SourceChoice = new Choice({ emdb: 'EMDB', empiar: 'EMPIAR' }, 'emdb');
type Source = Choice.Values<typeof SourceChoice>


export const CellStarEntryParams = {
    serverUrl: ParamDefinition.Text(DEFAULT_VOLUME_SERVER_V2),
    source: SourceChoice.PDSelect(),
    entryNumber: ParamDefinition.Text('1832'),
}

export class CellStarEntryData {
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
    public readonly meshSegmentationData = new CellStarMeshSegmentationData(this);
    public readonly modelData = new CellStarModelData(this);


    private constructor(plugin: PluginContext, serverUrl: string, source: Source, entryNumber: string) {
        this.plugin = plugin
        this.api = new VolumeApiV2(serverUrl);
        this.source = source;
        this.entryNumber = entryNumber;
        this.entryId = createEntryId(source, entryNumber);
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
        await this.meshSegmentationData.showSegmentation();
    }
    public async showSegments(segments: Segment[]) {
        await this.latticeSegmentationData.showSegments(segments);
        await this.meshSegmentationData.showSegments(segments);
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