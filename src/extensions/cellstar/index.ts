import { PluginBehavior } from '../../mol-plugin/behavior';
import { PluginContext } from '../../mol-plugin/context';
import { StateAction, StateObjectRef } from '../../mol-state';
import { PluginStateObject as SO, PluginStateTransform } from '../../mol-plugin-state/objects';
import { ParamDefinition } from '../../mol-util/param-definition';
import { Task } from '../../mol-task';

import { CellStarEntry, CellStarEntryFromRoot, CellStarEntryParams } from './entry-root';


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
        }
        unregister() {
            console.log('Unregistering CellStar extension behavior');
            this.ctx.state.data.actions.remove(LoadCellStar);
        }
    }
});


export const LoadCellStar = StateAction.build({
    display: { name: 'Load Cell*', description: 'Load entry from Cell* volume.' },
    from: SO.Root,
    params: CellStarEntryParams,
})(({ params, state }, ctx: PluginContext) => Task.create('CellStar', taskCtx => {
    return state.transaction(async () => {
        if (params.entryId.trim().length === 0) {
            throw new Error('Specify Entry ID');
        }
        console.log('Cell* loading', params.entryId);

        ctx.behaviors.layout.leftPanelTabName.next('data');

        // const trajectory = await state.build().toRoot()
        //     .apply(G3DHeaderFromUrl, { url: params.url })
        //     .apply(G3DTrajectory)
        //     .commit();

        // await defaultStructure(ctx, { trajectory });
        state.build().toRoot().apply(CellStarEntryFromRoot, params).commit();
    }).runInContext(taskCtx);
}));