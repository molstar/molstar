/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginBehavior } from '../behavior';
import { EmptyLoci, Loci, areLociEqual } from 'mol-model/loci';
import { MarkerAction } from 'mol-geo/geometry/marker-data';

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
    display: { name: 'Highlight Loci on Canvas', group: 'Data' },
    params: () => ({})
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
    display: { name: 'Select Loci on Canvas', group: 'Data' },
    params: () => ({})
});