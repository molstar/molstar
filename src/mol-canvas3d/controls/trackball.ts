/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 *
 * This code has been modified from https://github.com/mrdoob/three.js/,
 * copyright (c) 2010-2018 three.js authors. MIT License
 */

import { Quat, Vec2, Vec3, EPSILON } from '../../mol-math/linear-algebra';
import { Viewport } from '../camera/util';
import InputObserver, { DragInput, WheelInput, PinchInput, ButtonsType, ModifiersKeys } from '../../mol-util/input/input-observer';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Camera } from '../camera';
import { absMax } from '../../mol-math/misc';
import { Binding } from '../../mol-util/binding';

const B = ButtonsType;
const M = ModifiersKeys;
const Trigger = Binding.Trigger;

export const DefaultTrackballBindings = {
    dragRotate: Binding([Trigger(B.Flag.Primary, M.create())], 'Rotate', 'Drag using ${triggers}'),
    dragRotateZ: Binding([Trigger(B.Flag.Primary, M.create({ shift: true }))], 'Rotate around z-axis', 'Drag using ${triggers}'),
    dragPan: Binding([Trigger(B.Flag.Secondary, M.create()), Trigger(B.Flag.Primary, M.create({ control: true }))], 'Pan', 'Drag using ${triggers}'),
    dragZoom: Binding.Empty,
    dragFocus: Binding([Trigger(B.Flag.Forth, M.create())], 'Focus', 'Drag using ${triggers}'),
    dragFocusZoom: Binding([Trigger(B.Flag.Auxilary, M.create())], 'Focus and zoom', 'Drag using ${triggers}'),

    scrollZoom: Binding([Trigger(B.Flag.Auxilary, M.create())], 'Zoom', 'Scroll using ${triggers}'),
    scrollFocus: Binding([Trigger(B.Flag.Auxilary, M.create({ shift: true }))], 'Clip', 'Scroll using ${triggers}'),
    scrollFocusZoom: Binding.Empty,
};

export const TrackballControlsParams = {
    noScroll: PD.Boolean(true, { isHidden: true }),

    rotateSpeed: PD.Numeric(3.0, { min: 0.1, max: 10, step: 0.1 }),
    zoomSpeed: PD.Numeric(6.0, { min: 0.1, max: 10, step: 0.1 }),
    panSpeed: PD.Numeric(0.8, { min: 0.1, max: 5, step: 0.1 }),

    spin: PD.Boolean(false, { description: 'Spin the 3D scene around the x-axis in view space' }),
    spinSpeed: PD.Numeric(1, { min: -20, max: 20, step: 1 }),

    staticMoving: PD.Boolean(true, { isHidden: true }),
    dynamicDampingFactor: PD.Numeric(0.2, {}, { isHidden: true }),

    minDistance: PD.Numeric(0.01, {}, { isHidden: true }),
    maxDistance: PD.Numeric(1e150, {}, { isHidden: true }),

    bindings: PD.Value(DefaultTrackballBindings, { isHidden: true })
};
export type TrackballControlsProps = PD.Values<typeof TrackballControlsParams>

export { TrackballControls };
interface TrackballControls {
    viewport: Viewport

    readonly props: Readonly<TrackballControlsProps>
    setProps: (props: Partial<TrackballControlsProps>) => void

    update: (t: number) => void
    reset: () => void
    dispose: () => void
}
namespace TrackballControls {
    export function create(input: InputObserver, camera: Camera, props: Partial<TrackballControlsProps> = {}): TrackballControls {
        const p = { ...PD.getDefaultValues(TrackballControlsParams), ...props };

        const viewport = Viewport();

        let disposed = false;

        const dragSub = input.drag.subscribe(onDrag);
        const interactionEndSub = input.interactionEnd.subscribe(onInteractionEnd);
        const wheelSub = input.wheel.subscribe(onWheel);
        const pinchSub = input.pinch.subscribe(onPinch);

        let _isInteracting = false;

        // For internal use
        const lastPosition = Vec3();

        const _eye = Vec3();

        const _rotPrev = Vec2();
        const _rotCurr = Vec2();
        const _rotLastAxis = Vec3();
        let _rotLastAngle = 0;

        const _zRotPrev = Vec2();
        const _zRotCurr = Vec2();
        let _zRotLastAngle = 0;

        const _zoomStart = Vec2();
        const _zoomEnd = Vec2();

        const _focusStart = Vec2();
        const _focusEnd = Vec2();

        const _panStart = Vec2();
        const _panEnd = Vec2();

        // Initial values for reseting
        const target0 = Vec3.clone(camera.target);
        const position0 = Vec3.clone(camera.position);
        const up0 = Vec3.clone(camera.up);

        const mouseOnScreenVec2 = Vec2();
        function getMouseOnScreen(pageX: number, pageY: number) {
            return Vec2.set(
                mouseOnScreenVec2,
                (pageX - viewport.x) / viewport.width,
                (pageY - viewport.y) / viewport.height
            );
        }

        const mouseOnCircleVec2 = Vec2();
        function getMouseOnCircle(pageX: number, pageY: number) {
            return Vec2.set(
                mouseOnCircleVec2,
                (pageX - viewport.width * 0.5 - viewport.x) / (viewport.width * 0.5),
                (viewport.height + 2 * (viewport.y - pageY)) / viewport.width // screen.width intentional
            );
        }

        const rotAxis = Vec3();
        const rotQuat = Quat();
        const rotEyeDir = Vec3();
        const rotObjUpDir = Vec3();
        const rotObjSideDir = Vec3();
        const rotMoveDir = Vec3();

        function rotateCamera() {
            const dx = _rotCurr[0] - _rotPrev[0];
            const dy = _rotCurr[1] - _rotPrev[1];
            Vec3.set(rotMoveDir, dx, dy, 0);

            const angle = Vec3.magnitude(rotMoveDir) * p.rotateSpeed;

            if (angle) {
                Vec3.sub(_eye, camera.position, camera.target);

                Vec3.normalize(rotEyeDir, _eye);
                Vec3.normalize(rotObjUpDir, camera.up);
                Vec3.normalize(rotObjSideDir, Vec3.cross(rotObjSideDir, rotObjUpDir, rotEyeDir));

                Vec3.setMagnitude(rotObjUpDir, rotObjUpDir, dy);
                Vec3.setMagnitude(rotObjSideDir, rotObjSideDir, dx);

                Vec3.add(rotMoveDir, rotObjUpDir, rotObjSideDir);
                Vec3.normalize(rotAxis, Vec3.cross(rotAxis, rotMoveDir, _eye));
                Quat.setAxisAngle(rotQuat, rotAxis, angle);

                Vec3.transformQuat(_eye, _eye, rotQuat);
                Vec3.transformQuat(camera.up, camera.up, rotQuat);

                Vec3.copy(_rotLastAxis, rotAxis);
                _rotLastAngle = angle;
            } else if (!p.staticMoving && _rotLastAngle) {
                _rotLastAngle *= Math.sqrt(1.0 - p.dynamicDampingFactor);
                Vec3.sub(_eye, camera.position, camera.target);
                Quat.setAxisAngle(rotQuat, _rotLastAxis, _rotLastAngle);

                Vec3.transformQuat(_eye, _eye, rotQuat);
                Vec3.transformQuat(camera.up, camera.up, rotQuat);
            }

            Vec2.copy(_rotPrev, _rotCurr);
        }

        const zRotQuat = Quat();

        function zRotateCamera() {
            const dx = _zRotCurr[0] - _zRotPrev[0];
            const dy = _zRotCurr[1] - _zRotPrev[1];
            const angle = p.rotateSpeed * (-dx + dy) * -0.05;

            if (angle) {
                Vec3.sub(_eye, camera.position, camera.target);
                Quat.setAxisAngle(zRotQuat, _eye, angle);
                Vec3.transformQuat(camera.up, camera.up, zRotQuat);
                _zRotLastAngle = angle;
            } else if (!p.staticMoving && _zRotLastAngle) {
                _zRotLastAngle *= Math.sqrt(1.0 - p.dynamicDampingFactor);
                Vec3.sub(_eye, camera.position, camera.target);
                Quat.setAxisAngle(zRotQuat, _eye, _zRotLastAngle);
                Vec3.transformQuat(camera.up, camera.up, zRotQuat);
            }

            Vec2.copy(_zRotPrev, _zRotCurr);
        }

        function zoomCamera() {
            const factor = 1.0 + (_zoomEnd[1] - _zoomStart[1]) * p.zoomSpeed;
            if (factor !== 1.0 && factor > 0.0) {
                Vec3.scale(_eye, _eye, factor);
            }

            if (p.staticMoving) {
                Vec2.copy(_zoomStart, _zoomEnd);
            } else {
                _zoomStart[1] += (_zoomEnd[1] - _zoomStart[1]) * p.dynamicDampingFactor;
            }
        }

        function focusCamera() {
            const factor = (_focusEnd[1] - _focusStart[1]) * p.zoomSpeed;
            if (factor !== 0.0) {
                const radius = Math.max(1, camera.state.radius + camera.state.radius * factor);
                camera.setState({ radius });
            }

            if (p.staticMoving) {
                Vec2.copy(_focusStart, _focusEnd);
            } else {
                _focusStart[1] += (_focusEnd[1] - _focusStart[1]) * p.dynamicDampingFactor;
            }
        }

        const panMouseChange = Vec2();
        const panObjUp = Vec3();
        const panOffset = Vec3();

        function panCamera() {
            Vec2.sub(panMouseChange, Vec2.copy(panMouseChange, _panEnd), _panStart);

            if (Vec2.squaredMagnitude(panMouseChange)) {
                Vec2.scale(panMouseChange, panMouseChange, Vec3.magnitude(_eye) * p.panSpeed);

                Vec3.cross(panOffset, Vec3.copy(panOffset, _eye), camera.up);
                Vec3.setMagnitude(panOffset, panOffset, panMouseChange[0]);

                Vec3.setMagnitude(panObjUp, camera.up, panMouseChange[1]);
                Vec3.add(panOffset, panOffset, panObjUp);

                Vec3.add(camera.position, camera.position, panOffset);
                Vec3.add(camera.target, camera.target, panOffset);

                if (p.staticMoving) {
                    Vec2.copy(_panStart, _panEnd);
                } else {
                    Vec2.sub(panMouseChange, _panEnd, _panStart);
                    Vec2.scale(panMouseChange, panMouseChange, p.dynamicDampingFactor);
                    Vec2.add(_panStart, _panStart, panMouseChange);
                }
            }
        }

        /**
         * Ensure the distance between object and target is within the min/max distance
         * and not too large compared to `camera.state.radiusMax`
         */
        function checkDistances() {
            const maxDistance = Math.min(Math.max(camera.state.radiusMax * 1000, 0.01), p.maxDistance);
            if (Vec3.squaredMagnitude(_eye) > maxDistance * maxDistance) {
                Vec3.setMagnitude(_eye, _eye, maxDistance);
                Vec3.add(camera.position, camera.target, _eye);
                Vec2.copy(_zoomStart, _zoomEnd);
                Vec2.copy(_focusStart, _focusEnd);
            }

            if (Vec3.squaredMagnitude(_eye) < p.minDistance * p.minDistance) {
                Vec3.setMagnitude(_eye, _eye, p.minDistance);
                Vec3.add(camera.position, camera.target, _eye);
                Vec2.copy(_zoomStart, _zoomEnd);
                Vec2.copy(_focusStart, _focusEnd);
            }
        }

        let lastUpdated = -1;
        /** Update the object's position, direction and up vectors */
        function update(t: number) {
            if (lastUpdated === t) return;
            if (p.spin) spin(t - lastUpdated);

            Vec3.sub(_eye, camera.position, camera.target);

            rotateCamera();
            zRotateCamera();
            zoomCamera();
            focusCamera();
            panCamera();

            Vec3.add(camera.position, camera.target, _eye);
            checkDistances();

            if (Vec3.squaredDistance(lastPosition, camera.position) > EPSILON) {
                Vec3.copy(lastPosition, camera.position);
            }

            lastUpdated = t;
        }

        /** Reset object's vectors and the target vector to their initial values */
        function reset() {
            Vec3.copy(camera.target, target0);
            Vec3.copy(camera.position, position0);
            Vec3.copy(camera.up, up0);

            Vec3.sub(_eye, camera.position, camera.target);
            Vec3.copy(lastPosition, camera.position);
        }

        // listeners

        function onDrag({ pageX, pageY, buttons, modifiers, isStart }: DragInput) {
            _isInteracting = true;

            const dragRotate = Binding.match(p.bindings.dragRotate, buttons, modifiers);
            const dragRotateZ = Binding.match(p.bindings.dragRotateZ, buttons, modifiers);
            const dragPan = Binding.match(p.bindings.dragPan, buttons, modifiers);
            const dragZoom = Binding.match(p.bindings.dragZoom, buttons, modifiers);
            const dragFocus = Binding.match(p.bindings.dragFocus, buttons, modifiers);
            const dragFocusZoom = Binding.match(p.bindings.dragFocusZoom, buttons, modifiers);

            getMouseOnCircle(pageX, pageY);
            getMouseOnScreen(pageX, pageY);

            if (isStart) {
                if (dragRotate) {
                    Vec2.copy(_rotCurr, mouseOnCircleVec2);
                    Vec2.copy(_rotPrev, _rotCurr);
                }
                if (dragRotateZ) {
                    Vec2.copy(_zRotCurr, mouseOnCircleVec2);
                    Vec2.copy(_zRotPrev, _zRotCurr);
                }
                if (dragZoom || dragFocusZoom) {
                    Vec2.copy(_zoomStart, mouseOnScreenVec2);
                    Vec2.copy(_zoomEnd, _zoomStart);
                }
                if (dragFocus) {
                    Vec2.copy(_focusStart, mouseOnScreenVec2);
                    Vec2.copy(_focusEnd, _focusStart);
                }
                if (dragPan) {
                    Vec2.copy(_panStart, mouseOnScreenVec2);
                    Vec2.copy(_panEnd, _panStart);
                }
            }

            if (dragRotate) Vec2.copy(_rotCurr, mouseOnCircleVec2);
            if (dragRotateZ) Vec2.copy(_zRotCurr, mouseOnCircleVec2);
            if (dragZoom || dragFocusZoom) Vec2.copy(_zoomEnd, mouseOnScreenVec2);
            if (dragFocus) Vec2.copy(_focusEnd, mouseOnScreenVec2);
            if (dragFocusZoom) {
                const dist = Vec3.distance(camera.state.position, camera.state.target);
                camera.setState({ radius: dist / 5 });
            }
            if (dragPan) Vec2.copy(_panEnd, mouseOnScreenVec2);
        }

        function onInteractionEnd() {
            _isInteracting = false;
        }

        function onWheel({ dx, dy, dz, buttons, modifiers }: WheelInput) {
            const delta = absMax(dx, dy, dz);
            if (Binding.match(p.bindings.scrollZoom, buttons, modifiers)) {
                _zoomEnd[1] += delta * 0.0001;
            }
            if (Binding.match(p.bindings.scrollFocus, buttons, modifiers)) {
                _focusEnd[1] += delta * 0.0001;
            }
        }

        function onPinch({ fraction, buttons, modifiers }: PinchInput) {
            if (Binding.match(p.bindings.scrollZoom, buttons, modifiers)) {
                _isInteracting = true;
                _zoomEnd[1] += (fraction - 1) * 0.1;
            }
        }

        function dispose() {
            if (disposed) return;
            disposed = true;

            dragSub.unsubscribe();
            wheelSub.unsubscribe();
            pinchSub.unsubscribe();
            interactionEndSub.unsubscribe();
        }

        const _spinSpeed = Vec2.create(0.005, 0);
        function spin(deltaT: number) {
            const frameSpeed = (p.spinSpeed || 0) / 1000;
            _spinSpeed[0] = 60 * Math.min(Math.abs(deltaT), 1000 / 8) / 1000 * frameSpeed;
            if (!_isInteracting) Vec2.add(_rotCurr, _rotPrev, _spinSpeed);
        }

        // force an update at start
        update(0);

        return {
            viewport,

            get props() { return p as Readonly<TrackballControlsProps>; },
            setProps: (props: Partial<TrackballControlsProps>) => {
                Object.assign(p, props);
            },

            update,
            reset,
            dispose
        };
    }
}