/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
import { CameraHelperAxis, isCameraAxesLoci } from '../../../mol-canvas3d/helper/camera-helper';
import { Vec3 } from '../../../mol-math/linear-algebra';

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

export const CameraAxisHelper = PluginBehavior.create<{}>({
    name: 'camera-axis-helper',
    category: 'interaction',
    ctor: class extends PluginBehavior.Handler<{}> {
        register(): void {

            let lastPlane = CameraHelperAxis.None;
            let state = 0;

            this.subscribeObservable(this.ctx.behaviors.interaction.click, ({ current }) => {
                if (!this.ctx.canvas3d || !isCameraAxesLoci(current.loci)) return;

                const axis = current.loci.elements[0].groupId;
                if (axis === CameraHelperAxis.None) {
                    lastPlane = CameraHelperAxis.None;
                    state = 0;
                    return;
                }

                const { camera } = this.ctx.canvas3d;
                let dir: Vec3, up: Vec3;

                if (axis >= CameraHelperAxis.X && axis <= CameraHelperAxis.Z) {
                    lastPlane = CameraHelperAxis.None;
                    state = 0;

                    const d = Vec3.sub(Vec3(), camera.target, camera.position);
                    const c = Vec3.cross(Vec3(), d, camera.up);

                    up = Vec3();
                    up[axis - 1] = 1;
                    dir = Vec3.cross(Vec3(), up, c);
                    if (Vec3.magnitude(dir) === 0) dir = d;
                } else {
                    if (lastPlane === axis) {
                        state = (state + 1) % 2;
                    } else {
                        lastPlane = axis;
                        state = 0;
                    }

                    if (axis === CameraHelperAxis.XY) {
                        up = state ? Vec3.unitX : Vec3.unitY;
                        dir = Vec3.negUnitZ;
                    } else if (axis === CameraHelperAxis.XZ) {
                        up = state ? Vec3.unitX : Vec3.unitZ;
                        dir = Vec3.negUnitY;
                    } else {
                        up = state ? Vec3.unitY : Vec3.unitZ;
                        dir = Vec3.negUnitX;
                    }
                }

                this.ctx.canvas3d.requestCameraReset({
                    snapshot: (scene, camera) => camera.getInvariantFocus(scene.boundingSphereVisible.center, scene.boundingSphereVisible.radius, up, dir)
                });
            });
        }
    },
    params: () => ({}),
    display: { name: 'Camera Axis Helper' }
});