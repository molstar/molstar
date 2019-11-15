/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Loci } from '../../../mol-model/loci';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { PluginBehavior } from '../behavior';
import { ButtonsType, ModifiersKeys } from '../../../mol-util/input/input-observer';
import { Binding } from '../../../mol-util/binding';

const B = ButtonsType
const M = ModifiersKeys
const Trigger = Binding.Trigger

const DefaultFocusLociBindings = {
    clickCenterFocus: Binding([Trigger(B.Flag.Auxilary, M.create())], 'Center and focus the clicked element using ${triggers}.'),
}
const FocusLociParams = {
    minRadius: PD.Numeric(8, { min: 1, max: 50, step: 1 }),
    extraRadius: PD.Numeric(4, { min: 1, max: 50, step: 1 }, { description: 'Value added to the bounding-sphere radius of the Loci.' }),
    durationMs: PD.Numeric(250, { min: 0, max: 1000, step: 1 }, { description: 'Camera transition duration.' }),

    bindings: PD.Value(DefaultFocusLociBindings, { isHidden: true }),
}
type FocusLociProps = PD.Values<typeof FocusLociParams>

export const FocusLoci = PluginBehavior.create<FocusLociProps>({
    name: 'camera-focus-loci',
    category: 'interaction',
    ctor: class extends PluginBehavior.Handler<FocusLociProps> {
        register(): void {
            this.subscribeObservable(this.ctx.behaviors.interaction.click, ({ current, buttons, modifiers }) => {
                if (!this.ctx.canvas3d) return;

                const p = this.params;
                if (Binding.match(this.params.bindings.clickCenterFocus, buttons, modifiers)) {
                    const sphere = Loci.getBoundingSphere(current.loci);
                    if (sphere) {
                        const radius = Math.max(sphere.radius + p.extraRadius, p.minRadius);
                        this.ctx.canvas3d.camera.focus(sphere.center, radius, p.durationMs);
                    }
                }
            });
        }
    },
    params: () => FocusLociParams,
    display: { name: 'Focus Loci on Canvas' }
});