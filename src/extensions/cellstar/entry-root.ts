import { PluginStateObject as SO, PluginStateTransform } from '../../mol-plugin-state/objects';
import { PluginContext } from '../../mol-plugin/context';
import { Task } from '../../mol-task';
import { ParamDefinition } from '../../mol-util/param-definition';

import { DEFAULT_VOLUME_SERVER_V2, VolumeApiV2 } from './cellstar-api/api';
import { Metadata } from './cellstar-api/data';
import { Choice, createEntryId } from './helpers';


const SourceChoice = new Choice({emdb: 'EMDB', empiar: 'EMPIAR'}, 'emdb');
type Source = Choice.Values<typeof SourceChoice>

export const CellStarEntryParams = {
    serverUrl: ParamDefinition.Text(DEFAULT_VOLUME_SERVER_V2),
    source: SourceChoice.PDSelect(),
    entryNumber: ParamDefinition.Text('1832'),
}

interface CellStarEntryData {
    api: VolumeApiV2,
    source: Source,
    /** Number part of entry ID, e.g. '1832' */
    entryNumber: string,
    /** Full entry ID, e.g. 'emd-1832' */
    entryId: string,
    metadata: Metadata,
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
            const entryId = createEntryId(params.source, params.entryNumber);
            const api = new VolumeApiV2(params.serverUrl);
            const metadata = await api.getMetadata(params.source, entryId);
            const data: CellStarEntryData = {api: api, source: params.source, entryNumber: params.entryNumber, entryId: entryId, metadata: metadata};
            return new CellStarEntry(data, {label: entryId, description: 'CellStar Entry'});
        });
    }
});