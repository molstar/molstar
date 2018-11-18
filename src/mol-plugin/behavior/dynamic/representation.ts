/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginBehavior } from '../behavior';
import { EmptyLoci, Loci, areLociEqual } from 'mol-model/loci';
import { MarkerAction } from 'mol-geo/geometry/marker-data';
import { labelFirst } from 'mol-theme/label';
import { PluginContext } from 'mol-plugin/context';

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
    display: { name: 'Highlight Loci on Canvas', group: 'Representation' }
});

export const SelectLoci = PluginBehavior.create({
    name: 'representation-select-loci',
    ctor: class extends PluginBehavior.Handler {
        register(): void {
            let prevLoci: Loci = EmptyLoci, prevRepr: any = void 0;
            this.subscribeObservable(this.ctx.behaviors.canvas.selectLoci, current => {
                if (!this.ctx.canvas3d) return;
                if (current.repr !== prevRepr || !areLociEqual(current.loci, prevLoci)) {
                    this.ctx.canvas3d.mark(prevLoci, MarkerAction.Deselect);
                    this.ctx.canvas3d.mark(current.loci, MarkerAction.Select);
                    prevLoci = current.loci;
                    prevRepr = current.repr;
                } else {
                    this.ctx.canvas3d.mark(current.loci, MarkerAction.Toggle);
                }
                // this.ctx.canvas3d.mark(loci, MarkerAction.Toggle);
            });
        }
    },
    display: { name: 'Select Loci on Canvas', group: 'Representation' }
});

export const DefaultLociLabelProvider = PluginBehavior.create({
    name: 'default-loci-label-provider',
    ctor: class implements PluginBehavior<undefined> {
        private f = labelFirst;
        register(): void { this.ctx.lociLabels.addProvider(this.f); }
        unregister() { this.ctx.lociLabels.removeProvider(this.f); }
        constructor(protected ctx: PluginContext) { }
    },
    display: { name: 'Provide Default Loci Label', group: 'Representation' }
});