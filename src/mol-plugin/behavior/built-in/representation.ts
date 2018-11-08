/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginBehavior } from '../behavior';
import { PluginStateObject as SO } from '../../state/base';

class _AddRepresentationToCanvas extends PluginBehavior.Handler {
    register(): void {
        this.subscribeObservable(this.ctx.events.state.data.object.created, o => {
            if (!SO.isRepresentation3D(o.obj)) return;
            this.ctx.canvas3d.add(o.obj.data);
            this.ctx.canvas3d.requestDraw(true);
        });
        this.subscribeObservable(this.ctx.events.state.data.object.updated, o => {
            const oo = o.obj;
            if (!SO.isRepresentation3D(oo)) return;
            this.ctx.canvas3d.add(oo.data);
            this.ctx.canvas3d.requestDraw(true);
        });
        this.subscribeObservable(this.ctx.events.state.data.object.removed, o => {
            const oo = o.obj;
            if (!SO.isRepresentation3D(oo)) return;
            this.ctx.canvas3d.remove(oo.data);
            this.ctx.canvas3d.requestDraw(true);
            oo.data.destroy();
        });
        this.subscribeObservable(this.ctx.events.state.data.object.replaced, o => {
            if (o.oldObj && SO.isRepresentation3D(o.oldObj)) {
                this.ctx.canvas3d.remove(o.oldObj.data);
                this.ctx.canvas3d.requestDraw(true);
                o.oldObj.data.destroy();
            }
            if (o.newObj && SO.isRepresentation3D(o.newObj)) {
                this.ctx.canvas3d.add(o.newObj.data);
                this.ctx.canvas3d.requestDraw(true);
            }
        });
    }
}

export const AddRepresentationToCanvas = PluginBehavior.create({
    name: 'add-representation-to-canvas',
    ctor: _AddRepresentationToCanvas,
    display: { name: 'Add Representation To Canvas' }
});