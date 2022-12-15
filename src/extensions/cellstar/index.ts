import { PluginStateObject as SO } from '../../mol-plugin-state/objects';
import { PluginBehavior } from '../../mol-plugin/behavior';
import { PluginConfigItem } from '../../mol-plugin/config';
import { PluginContext } from '../../mol-plugin/context';
import { StateAction } from '../../mol-state';
import { Task } from '../../mol-task';
import { DEFAULT_VOLUME_SERVER_V2, VolumeApiV2 } from './cellstar-api/api';

import { CellstarEntryData, CellstarEntryParamValues, createLoadCellstarParams } from './entry-root';
import { createEntryId } from './helpers';
import { CellstarEntryFromRoot, CellstarStateFromEntry } from './transformers';
import { CellstarUI } from './ui';


export const CellstarVolumeServerConfig = {
    // DefaultServer: new PluginConfigItem('cellstar-volume-server', DEFAULT_VOLUME_SERVER_V2),
    DefaultServer: new PluginConfigItem('cellstar-volume-server', window.location.hostname !== 'localhost' ? DEFAULT_VOLUME_SERVER_V2 : 'http://localhost:9000/v2'), // DEBUG, TODO remove
};


export const Cellstar = PluginBehavior.create<{ autoAttach: boolean, showTooltip: boolean }>({
    name: 'cellstar',
    category: 'misc',
    display: {
        name: 'Cellstar',
        description: 'Cellstar'
    },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean, showTooltip: boolean }> {
        register() {
            this.ctx.state.data.actions.add(LoadCellstar);
            this.ctx.customStructureControls.set('cellstar', CellstarUI as any);
            this.initializeEntryLists(); // do not await

            const entries = new Map<string, CellstarEntryData>();
            this.subscribeObservable(this.ctx.state.data.events.cell.created, o => {
                if (o.cell.obj instanceof CellstarEntryData) entries.set(o.ref, o.cell.obj);
            });

            this.subscribeObservable(this.ctx.state.data.events.cell.removed, o => {
                if (entries.has(o.ref)) {
                    entries.get(o.ref)!.dispose();
                    entries.delete(o.ref);
                }
            });
        }
        unregister() {
            this.ctx.state.data.actions.remove(LoadCellstar);
            this.ctx.customStructureControls.delete('cellstar');
        }
        private async initializeEntryLists() {
            const apiUrl = this.ctx.config.get(CellstarVolumeServerConfig.DefaultServer) ?? DEFAULT_VOLUME_SERVER_V2;
            const api = new VolumeApiV2(apiUrl);
            const entryLists = await api.getEntryList(10 ** 6);
            Object.values(entryLists).forEach(l => l.sort());
            (this.ctx.customState as any).cellstarAvailableEntries = entryLists;
        }
    }
});


export const LoadCellstar = StateAction.build({
    display: { name: 'Load Volume & Segmentation' },
    from: SO.Root,
    params: (a, plugin: PluginContext) => {
        const res = createLoadCellstarParams(plugin, (plugin.customState as any).cellstarAvailableEntries);
        console.log('default entry id for EMDB:', res.source.map('emdb').defaultValue.entryId);
        return res;
    },
})(({ params, state }, ctx: PluginContext) => Task.create('Loading Volume & Segmentation', taskCtx => {
    return state.transaction(async () => {
        const entryParams = CellstarEntryParamValues.fromLoadCellstarParamValues(params);
        if (entryParams.entryId.trim().length === 0) {
            alert('Must specify Entry Id!');
            throw new Error('Specify Entry Id');
        }
        if (!entryParams.entryId.includes('-')) {
            // add source prefix if the user omitted it (e.g. 1832 -> emd-1832)
            entryParams.entryId = createEntryId(entryParams.source, entryParams.entryId);
        }
        ctx.behaviors.layout.leftPanelTabName.next('data');

        const entryNode = await state.build().toRoot().apply(CellstarEntryFromRoot, entryParams).commit();
        await state.build().to(entryNode).apply(CellstarStateFromEntry, {}).commit(); // TODO isGhost
        if (entryNode.data) {
            await entryNode.data.loadVolume();
            await entryNode.data.loadSegmentations();
        }
    }).runInContext(taskCtx);
}));