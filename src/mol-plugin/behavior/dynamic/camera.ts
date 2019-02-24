/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Loci } from 'mol-model/loci';
import { ParamDefinition } from 'mol-util/param-definition';
import { PluginBehavior } from '../behavior';

export const FocusLociOnSelect = PluginBehavior.create<{ minRadius: number, extraRadius: number }>({
    name: 'focus-loci-on-select',
    category: 'interaction',
    ctor: class extends PluginBehavior.Handler<{ minRadius: number, extraRadius: number }> {
        register(): void {
            this.subscribeObservable(this.ctx.behaviors.canvas.selectLoci, current => {
                if (!this.ctx.canvas3d) return;
                const sphere = Loci.getBoundingSphere(current.loci);
                if (!sphere) return;
                this.ctx.canvas3d.camera.focus(sphere.center, Math.max(sphere.radius + this.params.extraRadius, this.params.minRadius));
            });
        }
    },
    params: () => ({
        minRadius: ParamDefinition.Numeric(10, { min: 1, max: 50, step: 1 }),
        extraRadius: ParamDefinition.Numeric(4, { min: 1, max: 50, step: 1 }, { description: 'Value added to the boundning sphere radius of the Loci.' })
    }),
    display: { name: 'Focus Loci on Select', group: 'Camera' }
});