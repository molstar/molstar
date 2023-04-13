/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { PluginStateObject as SO } from '../../mol-plugin-state/objects';
import { PluginBehavior } from '../../mol-plugin/behavior';
import { PluginConfigItem } from '../../mol-plugin/config';
import { PluginContext } from '../../mol-plugin/context';
import { StateAction } from '../../mol-state';
import { Task } from '../../mol-task';
import { DEFAULT_VOLSEG_SERVER, VolumeApiV2 } from './volseg-api/api';

import { VolsegEntryData, VolsegEntryParamValues, createLoadVolsegParams } from './entry-root';
import { VolsegGlobalState } from './global-state';
import { createEntryId } from './helpers';
import { VolsegEntryFromRoot, VolsegGlobalStateFromRoot, VolsegStateFromEntry } from './transformers';
import { VolsegUI } from './ui';


const DEBUGGING = typeof window !== 'undefined' ? window?.location?.hostname === 'localhost' : false;

export const VolsegVolumeServerConfig = {
    // DefaultServer: new PluginConfigItem('volseg-volume-server', DEFAULT_VOLUME_SERVER_V2),
    DefaultServer: new PluginConfigItem('volseg-volume-server', DEBUGGING ? 'http://localhost:9000/v2' : DEFAULT_VOLSEG_SERVER),
};


export const Volseg = PluginBehavior.create<{ autoAttach: boolean, showTooltip: boolean }>({
    name: 'volseg',
    category: 'misc',
    display: {
        name: 'Volseg',
        description: 'Volseg'
    },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean, showTooltip: boolean }> {
        register() {
            this.ctx.state.data.actions.add(LoadVolseg);
            this.ctx.customStructureControls.set('volseg', VolsegUI as any);
            this.initializeEntryLists(); // do not await

            const entries = new Map<string, VolsegEntryData>();
            this.subscribeObservable(this.ctx.state.data.events.cell.created, o => {
                if (o.cell.obj instanceof VolsegEntryData) entries.set(o.ref, o.cell.obj);
            });

            this.subscribeObservable(this.ctx.state.data.events.cell.removed, o => {
                if (entries.has(o.ref)) {
                    entries.get(o.ref)!.dispose();
                    entries.delete(o.ref);
                }
            });
        }
        unregister() {
            this.ctx.state.data.actions.remove(LoadVolseg);
            this.ctx.customStructureControls.delete('volseg');
        }
        private async initializeEntryLists() {
            const apiUrl = this.ctx.config.get(VolsegVolumeServerConfig.DefaultServer) ?? DEFAULT_VOLSEG_SERVER;
            const api = new VolumeApiV2(apiUrl);
            const entryLists = await api.getEntryList(10 ** 6);
            Object.values(entryLists).forEach(l => l.sort());
            (this.ctx.customState as any).volsegAvailableEntries = entryLists;
        }
    }
});


export const LoadVolseg = StateAction.build({
    display: { name: 'Load Volume & Segmentation' },
    from: SO.Root,
    params: (a, plugin: PluginContext) => {
        const res = createLoadVolsegParams(plugin, (plugin.customState as any).volsegAvailableEntries);
        return res;
    },
})(({ params, state }, ctx: PluginContext) => Task.create('Loading Volume & Segmentation', taskCtx => {
    return state.transaction(async () => {
        const entryParams = VolsegEntryParamValues.fromLoadVolsegParamValues(params);
        if (entryParams.entryId.trim().length === 0) {
            alert('Must specify Entry Id!');
            throw new Error('Specify Entry Id');
        }
        if (!entryParams.entryId.includes('-')) {
            // add source prefix if the user omitted it (e.g. 1832 -> emd-1832)
            entryParams.entryId = createEntryId(entryParams.source, entryParams.entryId);
        }
        ctx.behaviors.layout.leftPanelTabName.next('data');

        const globalStateNode = ctx.state.data.selectQ(q => q.ofType(VolsegGlobalState))[0];
        if (!globalStateNode) {
            await state.build().toRoot().apply(VolsegGlobalStateFromRoot, {}, { state: { isGhost: !DEBUGGING } }).commit();
        }

        const entryNode = await state.build().toRoot().apply(VolsegEntryFromRoot, entryParams).commit();
        await state.build().to(entryNode).apply(VolsegStateFromEntry, {}, { state: { isGhost: !DEBUGGING } }).commit();
        if (entryNode.data) {
            await entryNode.data.loadVolume();
            await entryNode.data.loadSegmentations();
        }
    }).runInContext(taskCtx);
}));