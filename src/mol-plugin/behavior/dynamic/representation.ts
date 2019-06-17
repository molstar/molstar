/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { MarkerAction } from '../../../mol-geo/geometry/marker-data';
import { PluginContext } from '../../../mol-plugin/context';
import { labelFirst } from '../../../mol-theme/label';
import { PluginBehavior } from '../behavior';
import { Interaction } from '../../util/interaction';

export const HighlightLoci = PluginBehavior.create({
    name: 'representation-highlight-loci',
    category: 'interaction',
    ctor: class extends PluginBehavior.Handler {
        private lociMarkProvider = (loci: Interaction.Loci, action: MarkerAction) => {
            if (!this.ctx.canvas3d) return;
            this.ctx.canvas3d.mark(loci, action)
        }
        register() {
            this.ctx.lociHighlights.addProvider(this.lociMarkProvider)
        }
        unregister() {
            this.ctx.lociHighlights.removeProvider(this.lociMarkProvider)
        }
    },
    display: { name: 'Highlight Loci on Canvas' }
});

export const SelectLoci = PluginBehavior.create({
    name: 'representation-select-loci',
    category: 'interaction',
    ctor: class extends PluginBehavior.Handler {
        private lociMarkProvider = (loci: Interaction.Loci, action: MarkerAction) => {
            if (!this.ctx.canvas3d) return;
            this.ctx.canvas3d.mark(loci, action)
        }
        register() {
            this.ctx.lociSelections.addProvider(this.lociMarkProvider)
        }
        unregister() {
            this.ctx.lociSelections.removeProvider(this.lociMarkProvider)
        }
    },
    display: { name: 'Select Loci on Canvas' }
});

export const DefaultLociLabelProvider = PluginBehavior.create({
    name: 'default-loci-label-provider',
    category: 'interaction',
    ctor: class implements PluginBehavior<undefined> {
        private f = labelFirst;
        register() { this.ctx.lociLabels.addProvider(this.f); }
        unregister() { this.ctx.lociLabels.removeProvider(this.f); }
        constructor(protected ctx: PluginContext) { }
    },
    display: { name: 'Provide Default Loci Label' }
});
