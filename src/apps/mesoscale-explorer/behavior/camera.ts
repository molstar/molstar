/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Loci } from '../../../mol-model/loci';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { PluginBehavior } from '../../../mol-plugin/behavior';
import { ButtonsType, ModifiersKeys } from '../../../mol-util/input/input-observer';
import { Binding } from '../../../mol-util/binding';
import { PluginCommands } from '../../../mol-plugin/commands';
import { Sphere3D } from '../../../mol-math/geometry';

const B = ButtonsType;
const M = ModifiersKeys;
const Trigger = Binding.Trigger;

const DefaultMesoFocusLociBindings = {
    clickCenter: Binding([
        Trigger(B.Flag.Primary, M.create()),
    ], 'Camera center', 'Click element using ${triggers}'),
    clickCenterFocus: Binding([
        Trigger(B.Flag.Secondary, M.create()),
    ], 'Camera center and focus', 'Click element using ${triggers}'),
};
const MesoFocusLociParams = {
    minRadius: PD.Numeric(8, { min: 1, max: 50, step: 1 }),
    extraRadius: PD.Numeric(4, { min: 1, max: 50, step: 1 }, { description: 'Value added to the bounding-sphere radius of the Loci' }),
    durationMs: PD.Numeric(250, { min: 0, max: 1000, step: 1 }, { description: 'Camera transition duration' }),
    centerOnly: PD.Boolean(true, { description: 'Keep current camera distance' }),

    bindings: PD.Value(DefaultMesoFocusLociBindings, { isHidden: true }),
};
type MesoFocusLociProps = PD.Values<typeof MesoFocusLociParams>

export const MesoFocusLoci = PluginBehavior.create<MesoFocusLociProps>({
    name: 'camera-meso-focus-loci',
    category: 'interaction',
    ctor: class extends PluginBehavior.Handler<MesoFocusLociProps> {
        register(): void {
            this.subscribeObservable(this.ctx.behaviors.interaction.click, ({ current, button, modifiers }) => {
                const { canvas3d } = this.ctx;
                if (!canvas3d) return;

                const loci = Loci.normalize(current.loci, this.ctx.managers.interactivity.props.granularity);
                const sphere = Loci.getBoundingSphere(loci) || Sphere3D();

                const { clickCenter, clickCenterFocus } = this.params.bindings;
                const { durationMs, extraRadius, minRadius } = this.params;
                const radius = Math.max(sphere.radius + extraRadius, minRadius);

                if (Binding.match(clickCenter, button, modifiers)) {
                    if (Loci.isEmpty(current.loci)) {
                        PluginCommands.Camera.Reset(this.ctx, { });
                        return;
                    }

                    const snapshot = canvas3d.camera.getCenter(sphere.center);
                    canvas3d.requestCameraReset({ durationMs, snapshot });
                } else if (Binding.match(clickCenterFocus, button, modifiers)) {
                    if (Loci.isEmpty(current.loci)) {
                        PluginCommands.Camera.Reset(this.ctx, { });
                        return;
                    }

                    const snapshot = canvas3d.camera.getCenter(sphere.center, radius);
                    canvas3d.requestCameraReset({ durationMs, snapshot });
                }
            });
        }
    },
    params: () => MesoFocusLociParams,
    display: { name: 'Camera Meso Focus Loci on Canvas' }
});
