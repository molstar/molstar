import { PluginStateObject, PluginStateTransform } from '../../mol-plugin-state/objects';
import { PluginContext } from '../../mol-plugin/context';
import { StateTransformer } from '../../mol-state';
import { Task } from '../../mol-task';

import { CellStarEntry, CellStarEntryData, createCellStarEntryParams } from './entry-root';
import { CellStarState, CellStarStateParams, CELLSTAR_STATE_FROM_ENTRY_TRANSFORMER_NAME } from './entry-state';


export const CellStarEntryFromRoot = PluginStateTransform.BuiltIn({
    name: 'cellstar-entry-from-root',
    display: { name: 'CellStar Entry', description: 'CellStar Entry' },
    from: PluginStateObject.Root,
    to: CellStarEntry,
    params: (a, plugin: PluginContext) => createCellStarEntryParams(plugin),
})({
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Load CellStar Entry', async () => {
            const data = await CellStarEntryData.create(plugin, params);
            return new CellStarEntry(data, { label: data.entryId, description: 'CellStar Entry' });
        });
    },
    update({ b, oldParams, newParams }) {
        Object.assign(newParams, oldParams);
        console.error('Changing params of existing CellStarEntry node is not allowed');
        return StateTransformer.UpdateResult.Unchanged;
    }
});


export const CellStarStateFromEntry = PluginStateTransform.BuiltIn({
    name: CELLSTAR_STATE_FROM_ENTRY_TRANSFORMER_NAME,
    display: { name: 'CellStar Entry State', description: 'CellStar Entry State' },
    from: CellStarEntry,
    to: CellStarState,
    params: CellStarStateParams,
})({
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Create CellStar Entry State', async () => {
            return new CellStarState(params, { label: 'State' });
        });
    }
});