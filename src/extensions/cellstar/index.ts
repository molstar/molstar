import { PluginBehavior } from '../../mol-plugin/behavior';
import { PluginContext } from '../../mol-plugin/context';
import { StateAction } from '../../mol-state';
import { PluginStateObject as SO } from '../../mol-plugin-state/objects';
import { Task } from '../../mol-task';

import { CellStarEntryData, CellStarEntryFromRoot, CellStarEntryParams } from './entry-root';
import { createEntryId } from './helpers';
import { CellStarUI } from './ui';


export const CellStar = PluginBehavior.create<{ autoAttach: boolean, showTooltip: boolean }>({
    name: 'cellstar',
    category: 'misc',
    display: {
        name: 'CellStar',
        description: 'CellStar'
    },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean, showTooltip: boolean }> {
        register() {
            console.log('Registering CellStar extension behavior');
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
            console.log('Unregistering CellStar extension behavior');
            this.ctx.state.data.actions.remove(LoadCellStar);
            this.ctx.customStructureControls.delete('cellstar');
        }
    }
});


export const LoadCellStar = StateAction.build({
    display: { name: 'Load Cell*', description: 'Load entry from Cell* volume.' },
    from: SO.Root,
    params: CellStarEntryParams,
})(({ params, state }, ctx: PluginContext) => Task.create('CellStar', taskCtx => {
    return state.transaction(async () => {
        if (params.entryNumber.trim().length === 0) {
            throw new Error('Specify Entry ID');
        }
        console.log('Cell* loading', createEntryId(params.source, params.entryNumber));

        ctx.behaviors.layout.leftPanelTabName.next('data');

        const entryNode = await state.build().toRoot().apply(CellStarEntryFromRoot, params).commit();
        console.log('entryNode', entryNode.ref);
        if (entryNode.data) {
            await entryNode.data.volumeData.showVolume();
            await entryNode.data.showSegmentations();
        }
    }).runInContext(taskCtx);
}));