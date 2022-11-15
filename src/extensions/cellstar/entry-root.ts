import { PluginStateObject as SO, PluginStateTransform } from '../../mol-plugin-state/objects';
import { PluginContext } from '../../mol-plugin/context';
import { Task } from '../../mol-task';
import { ParamDefinition } from '../../mol-util/param-definition';


export const CellStarEntryParams = {
    serverUrl: ParamDefinition.Url('http://localhost:9000/v2'),
    entryId: ParamDefinition.Text('emd-1832') 
}

interface CellStarEntryData {
    // TODO Cell* entry metadata
}


export class CellStarEntry extends SO.Create<CellStarEntryData>({ name: 'CellStar Entry', typeClass: 'Object' }) { }


export const CellStarEntryFromRoot = PluginStateTransform.BuiltIn({
    name: 'cellstar-entry-from-root',
    display: { name: 'CellStar Entry', description: 'Load CellStar Entry' },
    from: SO.Root,
    to: CellStarEntry,
    params: CellStarEntryParams
})({
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Load CellStar Entry', async () => {
            const data = {};
            return new CellStarEntry(data, {label: params.entryId, description: 'CellStar Entry'});
        });
    }
});