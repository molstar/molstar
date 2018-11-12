/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginCommands } from '../../command';
import { PluginContext } from '../../context';
import { StateTree, Transform, State } from 'mol-state';
import { PluginStateSnapshotManager } from 'mol-plugin/state/snapshots';

export function registerDefault(ctx: PluginContext) {
    SetCurrentObject(ctx);
    Update(ctx);
    ApplyAction(ctx);
    RemoveObject(ctx);
    ToggleExpanded(ctx);
    ToggleVisibility(ctx);
    Snapshots(ctx);
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

export function Snapshots(ctx: PluginContext) {
    PluginCommands.State.Snapshots.Clear.subscribe(ctx, () => {
        ctx.state.snapshots.clear();
    });

    PluginCommands.State.Snapshots.Remove.subscribe(ctx, ({ id }) => {
        ctx.state.snapshots.remove(id);
    });

    PluginCommands.State.Snapshots.Add.subscribe(ctx, ({ name, description }) => {
        const entry = PluginStateSnapshotManager.Entry(name || new Date().toLocaleTimeString(), ctx.state.getSnapshot(), description);
        ctx.state.snapshots.add(entry);
    });

    PluginCommands.State.Snapshots.Apply.subscribe(ctx, ({ id }) => {
        const e = ctx.state.snapshots.getEntry(id);
        return ctx.state.setSnapshot(e.snapshot);
    });
}