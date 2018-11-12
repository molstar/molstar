/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from 'mol-plugin/context';
import { PluginCommands } from 'mol-plugin/command';
import { PluginStateObject as SO } from '../../state/objects';

export function registerDefault(ctx: PluginContext) {
    Reset(ctx);
}

export function Reset(ctx: PluginContext) {
    PluginCommands.Camera.Reset.subscribe(ctx, () => {
        const sel = ctx.state.data.select(q => q.root.subtree().ofType(SO.Molecule.Structure));
        if (!sel.length) return;

        const center = (sel[0].obj! as SO.Molecule.Structure).data.boundary.sphere.center;
        ctx.canvas3d.camera.setState({ target: center });
        ctx.canvas3d.requestDraw(true);

        // TODO
        // ctx.canvas3d.resetCamera();
    })
}