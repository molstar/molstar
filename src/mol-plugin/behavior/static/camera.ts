/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from 'mol-plugin/context';
import { PluginCommands } from 'mol-plugin/command';
import { PluginStateObject as SO } from '../../state/objects';
import { CameraSnapshotManager } from 'mol-plugin/state/camera';

export function registerDefault(ctx: PluginContext) {
    Reset(ctx);
    SetSnapshot(ctx);
    Snapshots(ctx);
}

export function Reset(ctx: PluginContext) {
    PluginCommands.Camera.Reset.subscribe(ctx, () => {
        const sel = ctx.state.dataState.select(q => q.root.subtree().ofType(SO.Molecule.Structure));
        if (!sel.length) return;

        const center = (sel[0].obj! as SO.Molecule.Structure).data.boundary.sphere.center;
        ctx.canvas3d.camera.setState({ target: center });
        ctx.canvas3d.requestDraw(true);

        // TODO
        // ctx.canvas3d.resetCamera();
    })
}

export function SetSnapshot(ctx: PluginContext) {
    PluginCommands.Camera.SetSnapshot.subscribe(ctx, ({ snapshot }) => {
        ctx.canvas3d.camera.setState(snapshot);
        ctx.canvas3d.requestDraw();
    })
}

export function Snapshots(ctx: PluginContext) {
    PluginCommands.Camera.Snapshots.Clear.subscribe(ctx, () => {
        ctx.state.cameraSnapshots.clear();
    });

    PluginCommands.Camera.Snapshots.Remove.subscribe(ctx, ({ id }) => {
        ctx.state.cameraSnapshots.remove(id);
    });

    PluginCommands.Camera.Snapshots.Add.subscribe(ctx, ({ name, description }) => {
        const entry = CameraSnapshotManager.Entry(name || new Date().toLocaleTimeString(), ctx.canvas3d.camera.getSnapshot(), description);
        ctx.state.cameraSnapshots.add(entry);
    });

    PluginCommands.Camera.Snapshots.Apply.subscribe(ctx, ({ id }) => {
        const e = ctx.state.cameraSnapshots.getEntry(id);
        return PluginCommands.Camera.SetSnapshot.dispatch(ctx, { snapshot: e.snapshot });
    });
}