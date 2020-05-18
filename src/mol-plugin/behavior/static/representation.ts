/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateObject as SO } from '../../../mol-plugin-state/objects';
import { PluginContext } from '../../../mol-plugin/context';
import { Representation } from '../../../mol-repr/representation';
import { StateObjectCell } from '../../../mol-state';

export function registerDefault(ctx: PluginContext) {
    SyncRepresentationToCanvas(ctx);
    SyncStructureRepresentation3DState(ctx); // should be AFTER SyncRepresentationToCanvas
    UpdateRepresentationVisibility(ctx);
}

export function SyncRepresentationToCanvas(ctx: PluginContext) {
    const events = ctx.state.data.events;
    events.object.created.subscribe(e => {
        if (!SO.isRepresentation3D(e.obj)) return;
        updateVisibility(e.state.cells.get(e.ref)!, e.obj.data.repr);
        e.obj.data.repr.setState({ syncManually: true });
        ctx.canvas3d?.add(e.obj.data.repr);
    });
    events.object.updated.subscribe(e => {
        if (e.oldObj && SO.isRepresentation3D(e.oldObj)) {
            ctx.canvas3d?.remove(e.oldObj.data.repr);
            e.oldObj.data.repr.destroy();
        }

        if (!SO.isRepresentation3D(e.obj)) {
            return;
        }

        updateVisibility(e.state.cells.get(e.ref)!, e.obj.data.repr);
        if (e.action === 'recreate') {
            e.obj.data.repr.setState({ syncManually: true });
        }
        ctx.canvas3d?.add(e.obj.data.repr);
    });
    events.object.removed.subscribe(e => {
        if (!SO.isRepresentation3D(e.obj)) return;
        ctx.canvas3d?.remove(e.obj.data.repr);

        e.obj.data.repr.destroy();
    });
}


export function SyncStructureRepresentation3DState(ctx: PluginContext) {
    // TODO: figure out how to do transform composition here?
    const events = ctx.state.data.events;
    events.object.created.subscribe(e => {
        if (!SO.Molecule.Structure.Representation3DState.is(e.obj)) return;
        const data = e.obj.data as SO.Molecule.Structure.Representation3DStateData;
        data.source.data.repr.setState(data.state);
        ctx.canvas3d?.update(data.source.data.repr);
    });
    events.object.updated.subscribe(e => {
        if (!SO.Molecule.Structure.Representation3DState.is(e.obj)) return;
        const data = e.obj.data as SO.Molecule.Structure.Representation3DStateData;
        data.source.data.repr.setState(data.state);
        ctx.canvas3d?.update(data.source.data.repr);
    });
    events.object.removed.subscribe(e => {
        if (!SO.Molecule.Structure.Representation3DState.is(e.obj)) return;
        const data = e.obj.data as SO.Molecule.Structure.Representation3DStateData;
        data.source.data.repr.setState(data.initialState);
        ctx.canvas3d?.update(data.source.data.repr);
    });
}


export function UpdateRepresentationVisibility(ctx: PluginContext) {
    ctx.state.data.events.cell.stateUpdated.subscribe(e => {
        const cell = e.state.cells.get(e.ref)!;
        if (!SO.isRepresentation3D(cell.obj)) return;
        updateVisibility(cell, cell.obj.data.repr);
        ctx.canvas3d?.syncVisibility();
    });
}

function updateVisibility(cell: StateObjectCell, r: Representation<any>) {
    r.setState({ visible: !cell.state.isHidden });
}