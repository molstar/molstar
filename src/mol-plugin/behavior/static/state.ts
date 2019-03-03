/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginCommands } from '../../command';
import { PluginContext } from '../../context';
import { StateTree, StateTransform, State } from 'mol-state';
import { PluginStateSnapshotManager } from 'mol-plugin/state/snapshots';
import { PluginStateObject as SO } from '../../state/objects';
import { getFormattedTime } from 'mol-util/date';
import { readFromFile } from 'mol-util/data-source';

export function registerDefault(ctx: PluginContext) {
    SyncBehaviors(ctx);
    SetCurrentObject(ctx);
    Update(ctx);
    ApplyAction(ctx);
    RemoveObject(ctx);
    ToggleExpanded(ctx);
    ToggleVisibility(ctx);
    Highlight(ctx);
    ClearHighlight(ctx);
    Snapshots(ctx);
}

export function SyncBehaviors(ctx: PluginContext) {
    ctx.events.state.object.created.subscribe(o => {
        if (!SO.isBehavior(o.obj)) return;
        o.obj.data.register(o.ref);
    });

    ctx.events.state.object.removed.subscribe(o => {
        if (!SO.isBehavior(o.obj)) return;
        o.obj.data.unregister();
    });

    ctx.events.state.object.updated.subscribe(o => {
        if (o.action === 'recreate') {
            if (o.oldObj && SO.isBehavior(o.oldObj)) o.oldObj.data.unregister();
            if (o.obj && SO.isBehavior(o.obj)) o.obj.data.register(o.ref);
        }
    });
}

export function SetCurrentObject(ctx: PluginContext) {
    PluginCommands.State.SetCurrentObject.subscribe(ctx, ({ state, ref }) => state.setCurrent(ref));
}

export function Update(ctx: PluginContext) {
    PluginCommands.State.Update.subscribe(ctx, ({ state, tree, options }) => ctx.runTask(state.updateTree(tree, options)));
}

export function ApplyAction(ctx: PluginContext) {
    PluginCommands.State.ApplyAction.subscribe(ctx, ({ state, action, ref }) => ctx.runTask(state.applyAction(action.action, action.params, ref)));
}

export function RemoveObject(ctx: PluginContext) {
    function remove(state: State, ref: string) {
        const tree = state.build().delete(ref).getTree();
        return ctx.runTask(state.updateTree(tree));
    }

    PluginCommands.State.RemoveObject.subscribe(ctx, ({ state, ref, removeParentGhosts }) => {
        if (removeParentGhosts) {
            const tree = state.tree;
            let curr = tree.transforms.get(ref);
            if (curr.parent === ref) return remove(state, ref);

            while (true) {
                const children = tree.children.get(curr.parent);
                if (curr.parent === curr.ref || children.size > 1) return remove(state, curr.ref);
                const parent = tree.transforms.get(curr.parent);
                if (!parent.props || !parent.props.isGhost) return remove(state, curr.ref);
                curr = parent;
            }
        } else {
            remove(state, ref);
        }
    });
}

export function ToggleExpanded(ctx: PluginContext) {
    PluginCommands.State.ToggleExpanded.subscribe(ctx, ({ state, ref }) => state.updateCellState(ref, ({ isCollapsed }) => ({ isCollapsed: !isCollapsed })));
}

export function ToggleVisibility(ctx: PluginContext) {
    PluginCommands.State.ToggleVisibility.subscribe(ctx, ({ state, ref }) => setVisibility(state, ref, !state.cellStates.get(ref).isHidden));
}

function setVisibility(state: State, root: StateTransform.Ref, value: boolean) {
    StateTree.doPreOrder(state.tree, state.transforms.get(root), { state, value }, setVisibilityVisitor);
}

function setVisibilityVisitor(t: StateTransform, tree: StateTree, ctx: { state: State, value: boolean }) {
    ctx.state.updateCellState(t.ref, { isHidden: ctx.value });
}

// TODO make isHighlighted and isSelect part of StateObjectCell.State and subscribe from there???
// TODO select structures of subtree
// TODO should also work for volumes and shapes
export function Highlight(ctx: PluginContext) {
    PluginCommands.State.Highlight.subscribe(ctx, ({ state, ref }) => {
        // const cell = state.select(ref)[0]
        // const repr = cell && SO.isRepresentation3D(cell.obj) ? cell.obj.data : undefined
        // if (cell && cell.obj && cell.obj.type === PluginStateObject.Molecule.Structure.type) {
        //     ctx.events.canvas3d.highlight.next({ current: { loci: Structure.Loci(cell.obj.data) } });
        // } else if (repr) {
        //     ctx.events.canvas3d.highlight.next({ current: { loci: EveryLoci, repr } });
        // }
    });
}

export function ClearHighlight(ctx: PluginContext) {
    PluginCommands.State.ClearHighlight.subscribe(ctx, ({ state, ref }) => {
        // ctx.events.canvas3d.highlight.next({ current: { loci: EmptyLoci } });
    });
}

export function Snapshots(ctx: PluginContext) {
    PluginCommands.State.Snapshots.Clear.subscribe(ctx, () => {
        ctx.state.snapshots.clear();
    });

    PluginCommands.State.Snapshots.Remove.subscribe(ctx, ({ id }) => {
        ctx.state.snapshots.remove(id);
    });

    PluginCommands.State.Snapshots.Add.subscribe(ctx, ({ name, description }) => {
        const entry = PluginStateSnapshotManager.Entry(ctx.state.getSnapshot(), name, description);
        ctx.state.snapshots.add(entry);
    });

    PluginCommands.State.Snapshots.Apply.subscribe(ctx, ({ id }) => {
        const e = ctx.state.snapshots.getEntry(id);
        return ctx.state.setSnapshot(e.snapshot);
    });

    PluginCommands.State.Snapshots.Upload.subscribe(ctx, ({ name, description, serverUrl }) => {
        return fetch(`${serverUrl}/set?name=${encodeURIComponent(name || '')}&description=${encodeURIComponent(description || '')}`, {
            method: 'POST',
            mode: 'cors',
            referrer: 'no-referrer',
            headers: { 'Content-Type': 'application/json; charset=utf-8' },
            body: JSON.stringify(ctx.state.getSnapshot())
        }) as any as Promise<void>;
    });

    PluginCommands.State.Snapshots.Fetch.subscribe(ctx, async ({ url }) => {
        const req = await fetch(url, { referrer: 'no-referrer' });
        const json = await req.json();
        return ctx.state.setSnapshot(json.data);
    });

    PluginCommands.State.Snapshots.DownloadToFile.subscribe(ctx, ({ name }) => {
        const element = document.createElement('a');
        const json = JSON.stringify(ctx.state.getSnapshot(), null, 2);
        element.setAttribute('href', 'data:application/json;charset=utf-8,' + encodeURIComponent(json));
        element.setAttribute('download', `mol-star_state_${(name || getFormattedTime())}.json`);
        element.style.display = 'none';
        document.body.appendChild(element);
        element.click();
        document.body.removeChild(element);
    });

    PluginCommands.State.Snapshots.OpenFile.subscribe(ctx, async ({ file }) => {
        try {
            const data = await readFromFile(file, 'string').run();
            const json = JSON.parse(data as string);
            await ctx.state.setSnapshot(json);
        } catch (e) {
            ctx.log.error(`Reading JSON state: ${e}`);
        }
    });
}