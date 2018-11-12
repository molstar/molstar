/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateObject as SO } from '../../state/objects';
import { PluginContext } from 'mol-plugin/context';

export function registerDefault(ctx: PluginContext) {
    SyncRepresentationToCanvas(ctx);
}

export function SyncRepresentationToCanvas(ctx: PluginContext) {
    ctx.events.state.data.object.created.subscribe(e => {
        if (!SO.isRepresentation3D(e.obj)) return;
        ctx.canvas3d.add(e.obj.data);
        ctx.canvas3d.requestDraw(true);

        // TODO: update visiblity
    });
    ctx.events.state.data.object.updated.subscribe(e => {
        if (e.oldObj && SO.isRepresentation3D(e.oldObj)) {
            ctx.canvas3d.remove(e.oldObj.data);
            ctx.canvas3d.requestDraw(true);
            e.oldObj.data.destroy();
        }

        if (!SO.isRepresentation3D(e.obj)) return;

        // TODO: update visiblity
        ctx.canvas3d.add(e.obj.data);
        ctx.canvas3d.requestDraw(true);
    });
    ctx.events.state.data.object.removed.subscribe(e => {
        const oo = e.obj;
        if (!SO.isRepresentation3D(oo)) return;
        ctx.canvas3d.remove(oo.data);
        ctx.canvas3d.requestDraw(true);
        oo.data.destroy();
    });
}

export function UpdateRepresentationVisibility(ctx: PluginContext) {
    ctx.events.state.data.object.cellState.subscribe(e => {
        const cell = e.state.cells.get(e.ref)!;
        if (!SO.isRepresentation3D(cell.obj)) return;

        // TODO: update visiblity
    })
}