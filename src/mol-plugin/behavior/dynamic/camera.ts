/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Loci } from '../../../mol-model/loci';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { PluginBehavior } from '../behavior';
import { ButtonsType, ModifiersKeys } from '../../../mol-util/input/input-observer';
import { Binding } from '../../../mol-util/binding';
import { PluginCommands } from '../../commands';

const B = ButtonsType;
const M = ModifiersKeys;
const Trigger = Binding.Trigger;

const DefaultFocusLociBindings = {
    clickCenterFocus: Binding([
        Trigger(B.Flag.Primary, M.create()),
        Trigger(B.Flag.Secondary, M.create()),
        Trigger(B.Flag.Primary, M.create({ control: true }))
    ], 'Camera center and focus', 'Click element using ${triggers}'),
    clickCenterFocusSelectMode: Binding([
        Trigger(B.Flag.Secondary, M.create()),
        Trigger(B.Flag.Primary, M.create({ control: true }))
    ], 'Camera center and focus', 'Click element using ${triggers}'),
};
const FocusLociParams = {
    minRadius: PD.Numeric(8, { min: 1, max: 50, step: 1 }),
    extraRadius: PD.Numeric(4, { min: 1, max: 50, step: 1 }, { description: 'Value added to the bounding-sphere radius of the Loci' }),
    durationMs: PD.Numeric(250, { min: 0, max: 1000, step: 1 }, { description: 'Camera transition duration' }),

    bindings: PD.Value(DefaultFocusLociBindings, { isHidden: true }),
};
type FocusLociProps = PD.Values<typeof FocusLociParams>

export const FocusLoci = PluginBehavior.create<FocusLociProps>({
    name: 'camera-focus-loci',
    category: 'interaction',
    ctor: class extends PluginBehavior.Handler<FocusLociProps> {
        register(): void {
            this.subscribeObservable(this.ctx.behaviors.interaction.click, ({ current, button, modifiers }) => {
                if (!this.ctx.canvas3d) return;

                const binding = this.ctx.selectionMode
                    ? this.params.bindings.clickCenterFocusSelectMode
                    : this.params.bindings.clickCenterFocus;

                if (Binding.match(binding, button, modifiers)) {
                    if (Loci.isEmpty(current.loci)) {
                        PluginCommands.Camera.Reset(this.ctx, { });
                        return;
                    }

                    const loci = Loci.normalize(current.loci, this.ctx.managers.interactivity.props.granularity);
                    this.ctx.managers.camera.focusLoci(loci, this.params);
                }
            });
        }
    },
    params: () => FocusLociParams,
    display: { name: 'Camera Focus Loci on Canvas' }
});