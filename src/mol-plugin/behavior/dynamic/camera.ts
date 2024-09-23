/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Jason Pattle <jpattle.exscientia.co.uk>
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
const Key = Binding.TriggerKey;

export const DefaultClickResetCameraOnEmpty = Binding([
    Trigger(B.Flag.Primary, M.create()),
    Trigger(B.Flag.Secondary, M.create()),
    Trigger(B.Flag.Primary, M.create({ control: true }))
], 'Reset camera focus', 'Click on nothing using ${triggers}');
export const DefaultClickResetCameraOnEmptySelectMode = Binding([
    Trigger(B.Flag.Secondary, M.create()),
    Trigger(B.Flag.Primary, M.create({ control: true }))
], 'Reset camera focus', 'Click on nothing using ${triggers}');

type FocusLociBindings = {
    clickCenterFocus: Binding
    clickCenterFocusSelectMode: Binding
    clickResetCameraOnEmpty?: Binding
    clickResetCameraOnEmptySelectMode?: Binding
}
export const DefaultFocusLociBindings: FocusLociBindings = {
    clickCenterFocus: Binding([
        Trigger(B.Flag.Primary, M.create()),
        Trigger(B.Flag.Secondary, M.create()),
        Trigger(B.Flag.Primary, M.create({ control: true }))
    ], 'Camera center and focus', 'Click element using ${triggers}'),
    clickCenterFocusSelectMode: Binding([
        Trigger(B.Flag.Secondary, M.create()),
        Trigger(B.Flag.Primary, M.create({ control: true }))
    ], 'Camera center and focus', 'Click element using ${triggers}'),
    clickResetCameraOnEmpty: DefaultClickResetCameraOnEmpty,
    clickResetCameraOnEmptySelectMode: DefaultClickResetCameraOnEmptySelectMode,
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

                const resetBinding = this.ctx.selectionMode
                    ? (this.params.bindings.clickResetCameraOnEmptySelectMode ?? DefaultClickResetCameraOnEmptySelectMode)
                    : (this.params.bindings.clickResetCameraOnEmpty ?? DefaultClickResetCameraOnEmpty);

                if (Loci.isEmpty(current.loci) && Binding.match(resetBinding, button, modifiers)) {
                    PluginCommands.Camera.Reset(this.ctx, { });
                    return;
                }

                if (Binding.match(binding, button, modifiers)) {
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

const DefaultCameraControlsBindings = {
    keySpinAnimation: Binding([Key('I')], 'Spin Animation', 'Press ${triggers}'),
    keyRockAnimation: Binding([Key('O')], 'Rock Animation', 'Press ${triggers}'),
    keyToggleFlyMode: Binding([Key('Space', M.create({ shift: true }))], 'Toggle Fly Mode', 'Press ${triggers}'),
    keyResetView: Binding([Key('T')], 'Reset View', 'Press ${triggers}'),
    keyGlobalIllumination: Binding([Key('G')], 'Gobal Illumination', 'Press ${triggers}'),
};
const CameraControlsParams = {
    bindings: PD.Value(DefaultCameraControlsBindings, { isHidden: true }),
};
type CameraControlsProps = PD.Values<typeof CameraControlsParams>

export const CameraControls = PluginBehavior.create<CameraControlsProps>({
    name: 'camera-controls',
    category: 'interaction',
    ctor: class extends PluginBehavior.Handler<CameraControlsProps> {
        register(): void {
            this.subscribeObservable(this.ctx.behaviors.interaction.key, ({ code, key, modifiers }) => {
                if (!this.ctx.canvas3d) return;

                // include defaults for backwards state compatibility
                const b = { ...DefaultCameraControlsBindings, ...this.params.bindings };
                const tp = this.ctx.canvas3d.props.trackball;
                const ip = this.ctx.canvas3d.props.illumination;

                if (Binding.matchKey(b.keySpinAnimation, code, modifiers, key)) {
                    const name = tp.animate.name !== 'spin' ? 'spin' : 'off';
                    if (name === 'off') {
                        this.ctx.canvas3d.setProps({
                            trackball: { animate: { name, params: {} } }
                        });
                    } else {
                        this.ctx.canvas3d.setProps({
                            trackball: { animate: {
                                name, params: { speed: 1 } }
                            }
                        });
                    }
                }

                if (Binding.matchKey(b.keyRockAnimation, code, modifiers, key)) {
                    const name = tp.animate.name !== 'rock' ? 'rock' : 'off';
                    if (name === 'off') {
                        this.ctx.canvas3d.setProps({
                            trackball: { animate: { name, params: {} } }
                        });
                    } else {
                        this.ctx.canvas3d.setProps({
                            trackball: { animate: {
                                name, params: { speed: 0.3, angle: 10 } }
                            }
                        });
                    }
                }

                if (Binding.matchKey(b.keyToggleFlyMode, code, modifiers, key)) {
                    const flyMode = !tp.flyMode;

                    this.ctx.canvas3d.setProps({
                        trackball: { flyMode }
                    });

                    if (this.ctx.canvas3dContext?.canvas) {
                        this.ctx.canvas3dContext.canvas.style.cursor = flyMode ? 'crosshair' : 'unset';
                    }
                }

                if (Binding.matchKey(b.keyResetView, code, modifiers, key)) {
                    PluginCommands.Camera.Reset(this.ctx, {});
                }

                if (Binding.matchKey(b.keyGlobalIllumination, code, modifiers, key)) {
                    this.ctx.canvas3d.setProps({
                        illumination: {
                            ...ip,
                            enabled: !ip.enabled,
                        }
                    });
                }
            });
        }
    },
    params: () => CameraControlsParams,
    display: { name: 'Camera Controls on Canvas' }
});