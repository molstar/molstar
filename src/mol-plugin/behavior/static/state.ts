/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginCommands } from '../../command';
import { PluginContext } from '../../context';
import { StateTree, StateTransform, State } from '../../../mol-state';
import { PluginStateSnapshotManager } from '../../../mol-plugin/state/snapshots';
import { PluginStateObject as SO } from '../../state/objects';
import { getFormattedTime } from '../../../mol-util/date';
import { readFromFile } from '../../../mol-util/data-source';
import { download } from '../../../mol-util/download';
import { Structure } from '../../../mol-model/structure';
import { urlCombine } from '../../../mol-util/url';
import { PluginConfig } from '../../config';

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
    PluginCommands.State.ToggleVisibility.subscribe(ctx, ({ state, ref }) => setVisibility(state, ref, !state.cells.get(ref)!.state.isHidden));
}

function setVisibility(state: State, root: StateTransform.Ref, value: boolean) {
    StateTree.doPreOrder(state.tree, state.transforms.get(root), { state, value }, setVisibilityVisitor);
}

function setVisibilityVisitor(t: StateTransform, tree: StateTree, ctx: { state: State, value: boolean }) {
    ctx.state.updateCellState(t.ref, { isHidden: ctx.value });
}

export function Highlight(ctx: PluginContext) {
    PluginCommands.State.Highlight.subscribe(ctx, ({ state, ref }) => {
        const cell = state.select(ref)[0];
        if (!cell) return;
        if (SO.Molecule.Structure.is(cell.obj)) {
            ctx.interactivity.lociHighlights.highlightOnly({ loci: Structure.Loci(cell.obj.data) }, false);
        } else if (cell && SO.isRepresentation3D(cell.obj)) {
            const { repr } = cell.obj.data
            ctx.interactivity.lociHighlights.highlightOnly({ loci: repr.getLoci(), repr }, false);
        } else if (SO.Molecule.Structure.Selections.is(cell.obj)) {
            ctx.interactivity.lociHighlights.clearHighlights();
            for (const entry of cell.obj.data) {
                ctx.interactivity.lociHighlights.highlight({ loci: entry.loci }, false);
            }
        }

        // TODO: highlight volumes?
        // TODO: select structures of subtree?
    });
}

export function ClearHighlight(ctx: PluginContext) {
    PluginCommands.State.ClearHighlight.subscribe(ctx, ({ state, ref }) => {
        ctx.interactivity.lociHighlights.clearHighlights();
    });
}

export function Snapshots(ctx: PluginContext) {
    ctx.config.set(PluginConfig.State.CurrentServer, ctx.config.get(PluginConfig.State.DefaultServer));

    PluginCommands.State.Snapshots.Clear.subscribe(ctx, () => {
        ctx.state.snapshots.clear();
    });

    PluginCommands.State.Snapshots.Remove.subscribe(ctx, ({ id }) => {
        ctx.state.snapshots.remove(id);
    });

    PluginCommands.State.Snapshots.Add.subscribe(ctx, ({ name, description, params }) => {
        const entry = PluginStateSnapshotManager.Entry(ctx.state.getSnapshot(params), { name, description });
        ctx.state.snapshots.add(entry);
    });

    PluginCommands.State.Snapshots.Replace.subscribe(ctx, ({ id, params }) => {
        ctx.state.snapshots.replace(id, ctx.state.getSnapshot(params));
    });

    PluginCommands.State.Snapshots.Move.subscribe(ctx, ({ id, dir }) => {
        ctx.state.snapshots.move(id, dir);
    });

    PluginCommands.State.Snapshots.Apply.subscribe(ctx, ({ id }) => {
        const snapshot = ctx.state.snapshots.setCurrent(id);
        if (!snapshot) return;
        return ctx.state.setSnapshot(snapshot);
    });

    PluginCommands.State.Snapshots.Upload.subscribe(ctx, ({ name, description, playOnLoad, serverUrl }) => {
        return fetch(urlCombine(serverUrl, `set?name=${encodeURIComponent(name || '')}&description=${encodeURIComponent(description || '')}`), {
            method: 'POST',
            mode: 'cors',
            referrer: 'no-referrer',
            headers: { 'Content-Type': 'application/json; charset=utf-8' },
            body: JSON.stringify(ctx.state.snapshots.getRemoteSnapshot({ name, description, playOnLoad }))
        }) as any as Promise<void>;
    });

    PluginCommands.State.Snapshots.Fetch.subscribe(ctx, async ({ url }) => {
        const json = await ctx.runTask(ctx.fetch({ url, type: 'json' })); //  fetch(url, { referrer: 'no-referrer' });
        await ctx.state.snapshots.setRemoteSnapshot(json.data);
    });

    PluginCommands.State.Snapshots.DownloadToFile.subscribe(ctx, ({ name }) => {
        const json = JSON.stringify(ctx.state.getSnapshot(), null, 2);
        const blob = new Blob([json], {type : 'application/json;charset=utf-8'});
        download(blob, `mol-star_state_${(name || getFormattedTime())}.json`)
    });

    PluginCommands.State.Snapshots.OpenFile.subscribe(ctx, async ({ file }) => {
        try {
            const data = await readFromFile(file, 'string').run();
            const snapshot = JSON.parse(data as string);
            return ctx.state.setSnapshot(snapshot);
        } catch (e) {
            ctx.log.error(`Reading JSON state: ${e}`);
        }
    });
}