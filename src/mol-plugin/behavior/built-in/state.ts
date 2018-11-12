/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginCommands } from '../../command';
import { PluginContext } from '../../context';

export function registerAll(ctx: PluginContext) {
    SetCurrentObject(ctx);
    Update(ctx);
    ApplyAction(ctx);
    RemoveObject(ctx);
    ToggleExpanded(ctx);
}

export function SetCurrentObject(ctx: PluginContext) {
    PluginCommands.State.SetCurrentObject.subscribe(ctx, ({ state, ref }) => state.setCurrent(ref));
}

export function Update(ctx: PluginContext) {
    PluginCommands.State.Update.subscribe(ctx, ({ state, tree }) => ctx.runTask(state.update(tree)));
}

export function ApplyAction(ctx: PluginContext) {
    PluginCommands.State.ApplyAction.subscribe(ctx, ({ state, action, ref }) => ctx.runTask(state.apply(action.action, action.params, ref)));
}

export function RemoveObject(ctx: PluginContext) {
    PluginCommands.State.RemoveObject.subscribe(ctx, ({ state, ref }) => {
        const tree = state.tree.build().delete(ref).getTree();
        return ctx.runTask(state.update(tree));
    });
}

export function ToggleExpanded(ctx: PluginContext) {
    PluginCommands.State.ToggleExpanded.subscribe(ctx, ({ state, ref }) => state.updateCellState(ref, ({ isCollapsed }) => ({ isCollapsed: !isCollapsed })));
}

// export const SetCurrentObject = PluginBehavior.create({
//     name: 'set-current-data-object-behavior',
//     ctor: PluginBehavior.simpleCommandHandler(PluginCommands.State.SetCurrentObject, ({ state, ref }, ctx) => state.setCurrent(ref)),
//     display: { name: 'Set Current Handler', group: 'Data' }
// });

// export const Update = PluginBehavior.create({
//     name: 'update-data-behavior',
//     ctor: PluginBehavior.simpleCommandHandler(PluginCommands.State.Update, ({ state, tree }, ctx) => ctx.runTask(state.update(tree))),
//     display: { name: 'Update Data Handler', group: 'Data' }
// });

// export const ApplyAction = PluginBehavior.create({
//     name: 'update-data-behavior',
//     ctor: PluginBehavior.simpleCommandHandler(PluginCommands.State.Update, ({ state, tree }, ctx) => ctx.runTask(state.update(tree))),
//     display: { name: 'Update Data Handler', group: 'Data' }
// });

// export const RemoveObject = PluginBehavior.create({
//     name: 'remove-object-data-behavior',
//     ctor: PluginBehavior.simpleCommandHandler(PluginCommands.State.RemoveObject, ({ state, ref }, ctx) => {
//         const tree = state.tree.build().delete(ref).getTree();
//         return ctx.runTask(state.update(tree));
//     }),
//     display: { name: 'Remove Object Handler', group: 'Data' }
// });