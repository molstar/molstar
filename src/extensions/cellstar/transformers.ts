/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { PluginStateObject, PluginStateTransform } from '../../mol-plugin-state/objects';
import { PluginContext } from '../../mol-plugin/context';
import { StateTransformer } from '../../mol-state';
import { Task } from '../../mol-task';

import { CellstarEntry, CellstarEntryData, createCellstarEntryParams } from './entry-root';
import { CellstarState, CellstarStateParams, CELLSTAR_STATE_FROM_ENTRY_TRANSFORMER_NAME } from './entry-state';
import { CellstarGlobalState, CellstarGlobalStateData, CellstarGlobalStateParams } from './global-state';


export const CellstarEntryFromRoot = PluginStateTransform.BuiltIn({
    name: 'cellstar-entry-from-root',
    display: { name: 'Vol & Seg Entry', description: 'Vol & Seg Entry' },
    from: PluginStateObject.Root,
    to: CellstarEntry,
    params: (a, plugin: PluginContext) => createCellstarEntryParams(plugin),
})({
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Load Vol & Seg Entry', async () => {
            const data = await CellstarEntryData.create(plugin, params);
            return new CellstarEntry(data, { label: data.entryId, description: 'Vol & Seg Entry' });
        });
    },
    update({ b, oldParams, newParams }) {
        Object.assign(newParams, oldParams);
        console.error('Changing params of existing CellstarEntry node is not allowed');
        return StateTransformer.UpdateResult.Unchanged;
    }
});


export const CellstarStateFromEntry = PluginStateTransform.BuiltIn({
    name: CELLSTAR_STATE_FROM_ENTRY_TRANSFORMER_NAME,
    display: { name: 'Vol & Seg Entry State', description: 'Vol & Seg Entry State' },
    from: CellstarEntry,
    to: CellstarState,
    params: CellstarStateParams,
})({
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Create Vol & Seg Entry State', async () => {
            return new CellstarState(params, { label: 'State' });
        });
    }
});


export const CellstarGlobalStateFromRoot = PluginStateTransform.BuiltIn({
    name: 'cellstar-global-state-from-root',
    display: { name: 'Vol & Seg Global State', description: 'Vol & Seg Global State' },
    from: PluginStateObject.Root,
    to: CellstarGlobalState,
    params: CellstarGlobalStateParams,
})({
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Create Vol & Seg Global State', async () => {
            const data = new CellstarGlobalStateData(plugin, params);
            return new CellstarGlobalState(data, { label: 'Global State', description: 'Vol & Seg Global State' });
        });
    },
    update({ b, oldParams, newParams }) {
        b.data.currentState.next(newParams);
        return StateTransformer.UpdateResult.Updated;
    }
});