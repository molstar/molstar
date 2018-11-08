/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginBehavior } from '../behavior';
import { PluginStateObject as SO } from '../../state/base';
import { EmptyLoci, Loci, areLociEqual } from 'mol-model/loci';
import { MarkerAction } from 'mol-geo/geometry/marker-data';

export const AddRepresentationToCanvas = PluginBehavior.create({
    name: 'add-representation-to-canvas',
    ctor: class extends PluginBehavior.Handler {
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
    },
    display: { name: 'Add Representation To Canvas' }
});

export const HighlightLoci = PluginBehavior.create({
    name: 'representation-highlight-loci',
    ctor: class extends PluginBehavior.Handler {
        register(): void {
            let prevLoci: Loci = EmptyLoci, prevRepr: any = void 0;
            this.subscribeObservable(this.ctx.behaviors.canvas.highlightLoci, current => {
                if (!this.ctx.canvas3d) return;

                if (current.repr !== prevRepr || !areLociEqual(current.loci, prevLoci)) {
                    this.ctx.canvas3d.mark(prevLoci, MarkerAction.RemoveHighlight);
                    this.ctx.canvas3d.mark(current.loci, MarkerAction.Highlight);
                    prevLoci = current.loci;
                    prevRepr = current.repr;
                }
            });
        }
    },
    display: { name: 'Highlight Loci on Canvas' }
});

export const SelectLoci = PluginBehavior.create({
    name: 'representation-select-loci',
    ctor: class extends PluginBehavior.Handler {
        register(): void {
            this.subscribeObservable(this.ctx.behaviors.canvas.selectLoci, ({ loci }) => {
                if (!this.ctx.canvas3d) return;
                this.ctx.canvas3d.mark(loci, MarkerAction.Toggle);
            });
        }
    },
    display: { name: 'Select Loci on Canvas' }
});