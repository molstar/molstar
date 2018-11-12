/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginCommands } from '../../command';
import { PluginContext } from '../../context';
import { StateTree, Transform, State } from 'mol-state';

export function registerDefault(ctx: PluginContext) {
    SetCurrentObject(ctx);
    Update(ctx);
    ApplyAction(ctx);
    RemoveObject(ctx);
    ToggleExpanded(ctx);
    ToggleVisibility(ctx);
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

export function ToggleVisibility(ctx: PluginContext) {
    PluginCommands.State.ToggleVisibility.subscribe(ctx, ({ state, ref }) => setVisibility(state, ref, !state.tree.cellStates.get(ref).isHidden));
}

function setVisibility(state: State, root: Transform.Ref, value: boolean) {
    StateTree.doPreOrder(state.tree, state.tree.transforms.get(root), { state, value }, setVisibilityVisitor);
}

function setVisibilityVisitor(t: Transform, tree: StateTree, ctx: { state: State, value: boolean }) {
    ctx.state.updateCellState(t.ref, { isHidden: ctx.value });
}