import { PluginStateObject as SO } from '../../mol-plugin-state/objects';
import { PluginBehavior } from '../../mol-plugin/behavior';
import { PluginConfigItem } from '../../mol-plugin/config';
import { PluginContext } from '../../mol-plugin/context';
import { StateAction } from '../../mol-state';
import { Task } from '../../mol-task';
import { DEFAULT_VOLUME_SERVER_V2 } from './cellstar-api/api';

import { CellStarEntryData, createCellStarEntryParams } from './entry-root';
import { CellStarEntryFromRoot, CellStarStateFromEntry } from './transformers';
import { CellStarUI } from './ui';


export const CellStarVolumeServerConfig = {
    // DefaultServer: new PluginConfigItem('cellstar-volume-server', DEFAULT_VOLUME_SERVER_V2),
    DefaultServer: new PluginConfigItem('cellstar-volume-server', window.location.hostname !== 'localhost' ? DEFAULT_VOLUME_SERVER_V2 : 'http://localhost:9000/v2'), // DEBUG, TODO remove
};


export const CellStar = PluginBehavior.create<{ autoAttach: boolean, showTooltip: boolean }>({
    name: 'cellstar',
    category: 'misc',
    display: {
        name: 'CellStar',
        description: 'CellStar'
    },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean, showTooltip: boolean }> {
        register() {
            this.ctx.state.data.actions.add(LoadCellStar);
            this.ctx.customStructureControls.set('cellstar', CellStarUI as any);

            const entries = new Map<string, CellStarEntryData>();
            this.subscribeObservable(this.ctx.state.data.events.cell.created, o => {
                if (o.cell.obj instanceof CellStarEntryData) entries.set(o.ref, o.cell.obj);
            });

            this.subscribeObservable(this.ctx.state.data.events.cell.removed, o => {
                if (entries.has(o.ref)) {
                    entries.get(o.ref)!.dispose();
                    entries.delete(o.ref);
                }
            });
        }
        unregister() {
            this.ctx.state.data.actions.remove(LoadCellStar);
            this.ctx.customStructureControls.delete('cellstar');
        }
    }
});


export const LoadCellStar = StateAction.build({
    display: { name: 'Load Volume & Segmentation' },
    from: SO.Root,
    params: (a, plugin: PluginContext) => createCellStarEntryParams(plugin),
})(({ params, state }, ctx: PluginContext) => Task.create('Loading Volume & Segmentation', taskCtx => {
    return state.transaction(async () => {
        if (params.entryNumber.trim().length === 0) {
            throw new Error('Specify Entry ID');
        }
        ctx.behaviors.layout.leftPanelTabName.next('data');

        const entryNode = await state.build().toRoot().apply(CellStarEntryFromRoot, params).commit();
        await state.build().to(entryNode).apply(CellStarStateFromEntry, {}).commit(); // TODO isGhost
        if (entryNode.data) {
            await entryNode.data.volumeData.loadVolume();
            await entryNode.data.loadSegmentations();
        }
    }).runInContext(taskCtx);
}));