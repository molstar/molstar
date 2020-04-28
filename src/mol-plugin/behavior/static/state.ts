/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure } from '../../../mol-model/structure';
import { PluginStateSnapshotManager } from '../../../mol-plugin-state/manager/snapshots';
import { PluginStateObject as SO } from '../../../mol-plugin-state/objects';
import { State, StateTransform, StateTree } from '../../../mol-state';
import { getFormattedTime } from '../../../mol-util/date';
import { download } from '../../../mol-util/download';
import { urlCombine } from '../../../mol-util/url';
import { PluginCommands } from '../../commands';
import { PluginConfig } from '../../config';
import { PluginContext } from '../../context';

export function registerDefault(ctx: PluginContext) {
    SyncBehaviors(ctx);
    SetCurrentObject(ctx);
    Update(ctx);
    ApplyAction(ctx);
    RemoveObject(ctx);
    ToggleExpanded(ctx);
    ToggleVisibility(ctx);
    Highlight(ctx);
    ClearHighlights(ctx);
    Snapshots(ctx);
}

export function SyncBehaviors(ctx: PluginContext) {
    ctx.state.events.object.created.subscribe(o => {
        if (!SO.isBehavior(o.obj)) return;
        o.obj.data.register(o.ref);
    });

    ctx.state.events.object.removed.subscribe(o => {
        if (!SO.isBehavior(o.obj)) return;
        o.obj.data.unregister();
    });

    ctx.state.events.object.updated.subscribe(o => {
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
        const tree = state.build().delete(ref);
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
                // TODO: should this use "cell state" instead?
                if (!parent.state.isGhost) return remove(state, curr.ref);
                curr = parent;
            }
        } else {
            return remove(state, ref);
        }
    });
}

export function ToggleExpanded(ctx: PluginContext) {
    PluginCommands.State.ToggleExpanded.subscribe(ctx, ({ state, ref }) => state.updateCellState(ref, ({ isCollapsed }) => ({ isCollapsed: !isCollapsed })));
}

export function ToggleVisibility(ctx: PluginContext) {
    PluginCommands.State.ToggleVisibility.subscribe(ctx, ({ state, ref }) => setSubtreeVisibility(state, ref, !state.cells.get(ref)!.state.isHidden));
}

export function setSubtreeVisibility(state: State, root: StateTransform.Ref, value: boolean) {
    StateTree.doPreOrder(state.tree, state.transforms.get(root), { state, value }, setVisibilityVisitor);
}

function setVisibilityVisitor(t: StateTransform, tree: StateTree, ctx: { state: State, value: boolean }) {
    ctx.state.updateCellState(t.ref, { isHidden: ctx.value });
}

export function Highlight(ctx: PluginContext) {
    PluginCommands.Interactivity.Object.Highlight.subscribe(ctx, ({ state, ref }) => {
        ctx.managers.interactivity.lociHighlights.clearHighlights();

        const refs = typeof ref === 'string' ? [ref] : ref;
        for (const r of refs) {
            const cell = state.cells.get(r);
            if (!cell) continue;
            if (SO.Molecule.Structure.is(cell.obj)) {
                ctx.managers.interactivity.lociHighlights.highlight({ loci: Structure.Loci(cell.obj.data) }, false);
            } else if (cell && SO.isRepresentation3D(cell.obj)) {
                const { repr } = cell.obj.data;
                ctx.managers.interactivity.lociHighlights.highlight({ loci: repr.getLoci(), repr }, false);
            } else if (SO.Molecule.Structure.Selections.is(cell.obj)) {
                for (const entry of cell.obj.data) {
                    ctx.managers.interactivity.lociHighlights.highlight({ loci: entry.loci }, false);
                }
            }
        }

        // TODO: highlight volumes?
        // TODO: select structures of subtree?
    });
}

export function ClearHighlights(ctx: PluginContext) {
    PluginCommands.Interactivity.ClearHighlights.subscribe(ctx, () => {
        ctx.managers.interactivity.lociHighlights.clearHighlights();
    });
}

export function Snapshots(ctx: PluginContext) {
    ctx.config.set(PluginConfig.State.CurrentServer, ctx.config.get(PluginConfig.State.DefaultServer));

    PluginCommands.State.Snapshots.Clear.subscribe(ctx, () => {
        ctx.managers.snapshot.clear();
    });

    PluginCommands.State.Snapshots.Remove.subscribe(ctx, ({ id }) => {
        ctx.managers.snapshot.remove(id);
    });

    PluginCommands.State.Snapshots.Add.subscribe(ctx, ({ name, description, params }) => {
        const entry = PluginStateSnapshotManager.Entry(ctx.state.getSnapshot(params), { name, description });
        ctx.managers.snapshot.add(entry);
    });

    PluginCommands.State.Snapshots.Replace.subscribe(ctx, ({ id, params }) => {
        ctx.managers.snapshot.replace(id, ctx.state.getSnapshot(params));
    });

    PluginCommands.State.Snapshots.Move.subscribe(ctx, ({ id, dir }) => {
        ctx.managers.snapshot.move(id, dir);
    });

    PluginCommands.State.Snapshots.Apply.subscribe(ctx, ({ id }) => {
        const snapshot = ctx.managers.snapshot.setCurrent(id);
        if (!snapshot) return;
        return ctx.state.setSnapshot(snapshot);
    });

    PluginCommands.State.Snapshots.Upload.subscribe(ctx, ({ name, description, playOnLoad, serverUrl, params }) => {
        return fetch(urlCombine(serverUrl, `set?name=${encodeURIComponent(name || '')}&description=${encodeURIComponent(description || '')}`), {
            method: 'POST',
            mode: 'cors',
            referrer: 'no-referrer',
            headers: { 'Content-Type': 'application/json; charset=utf-8' },
            body: JSON.stringify(ctx.managers.snapshot.getStateSnapshot({ name, description, playOnLoad }))
        }) as any as Promise<void>;
    });

    PluginCommands.State.Snapshots.Fetch.subscribe(ctx, async ({ url }) => {
        const json = await ctx.runTask(ctx.fetch({ url, type: 'json' })); //  fetch(url, { referrer: 'no-referrer' });
        await ctx.managers.snapshot.setStateSnapshot(json.data);
    });

    PluginCommands.State.Snapshots.DownloadToFile.subscribe(ctx, async ({ name, type, params }) => {
        const filename = `mol-star_state_${(name || getFormattedTime())}.${type === 'json' ? 'molj' : 'molx'}`;
        const data = await ctx.managers.snapshot.serialize({ type, params });
        download(data, `${filename}`);
    });

    PluginCommands.State.Snapshots.OpenFile.subscribe(ctx, ({ file }) => {
        return ctx.managers.snapshot.open(file);
    });

    PluginCommands.State.Snapshots.OpenUrl.subscribe(ctx, async ({ url, type }) => {
        const data = await ctx.runTask(ctx.fetch({ url, type: 'binary' }));
        return ctx.managers.snapshot.open(new File([data], `state.${type}`));
    });
}