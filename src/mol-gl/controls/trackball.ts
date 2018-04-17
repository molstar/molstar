/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/*
 * This code has been modified from https://github.com/mrdoob/three.js/,
 * copyright (c) 2010-2018 three.js authors. MIT License
 */

import { Quat, Vec2, Vec3, EPSILON } from 'mol-math/linear-algebra';
import { cameraLookAt, Viewport } from '../camera/util';
import InputObserver, { DragInput, WheelInput, MouseButtonsFlag, PinchInput } from 'mol-util/input/input-observer';

export const DefaultTrackballControlsProps = {
    noScroll: true,

    rotateSpeed: 3.0,
    zoomSpeed: 2.0,
    panSpeed: 0.1,

    staticMoving: false,
    dynamicDampingFactor: 0.2,

    minDistance: 0,
    maxDistance: Infinity
}
export type TrackballControlsProps = Partial<typeof DefaultTrackballControlsProps>

interface Object {
    position: Vec3,
    direction: Vec3,
    up: Vec3,
}

interface TrackballControls {
    viewport: Viewport

    dynamicDampingFactor: number
    rotateSpeed: number
    zoomSpeed: number
    panSpeed: number

    update: () => void
    reset: () => void
    dispose: () => void
}

namespace TrackballControls {
    export function create (element: Element, object: Object, props: TrackballControlsProps = {}): TrackballControls {
        const p = { ...DefaultTrackballControlsProps, ...props }

        const viewport: Viewport = { x: 0, y: 0, width: 0, height: 0 }

        let { rotateSpeed, zoomSpeed, panSpeed } = p
        let { staticMoving, dynamicDampingFactor } = p
        let { minDistance, maxDistance } = p

        let disposed = false

        const input = InputObserver.create(element)
        input.drag.subscribe(onDrag)
        input.wheel.subscribe(onWheel)
        input.pinch.subscribe(onPinch)

        // internals
        const target = Vec3.zero()
        const lastPosition = Vec3.zero()

        const _eye = Vec3.zero()

        const _movePrev = Vec2.zero()
        const _moveCurr = Vec2.zero()

        const _lastAxis = Vec3.zero()
        let _lastAngle = 0

        const _zoomStart = Vec2.zero()
        const _zoomEnd = Vec2.zero()

        // let _touchZoomDistanceStart = 0
        // let _touchZoomDistanceEnd = 0

        const _panStart = Vec2.zero()
        const _panEnd = Vec2.zero()

        // for reset
        const target0 = Vec3.clone(target)
        const position0 = Vec3.clone(object.position)
        const up0 = Vec3.clone(object.up)

        const mouseOnScreenVec2 = Vec2.zero()
        function getMouseOnScreen(pageX: number, pageY: number) {
            Vec2.set(
                mouseOnScreenVec2,
                (pageX - viewport.x) / viewport.width,
                (pageY - viewport.y) / viewport.height
            );
            return mouseOnScreenVec2;
        }

        const mouseOnCircleVec2 = Vec2.zero()
        function getMouseOnCircle(pageX: number, pageY: number) {
            Vec2.set(
                mouseOnCircleVec2,
                ((pageX - viewport.width * 0.5 - viewport.x) / (viewport.width * 0.5)),
                ((viewport.height + 2 * (viewport.y - pageY)) / viewport.width) // screen.width intentional
            );
            return mouseOnCircleVec2;
        }

        const rotAxis = Vec3.zero()
        const rotQuat = Quat.zero()
        const rotEyeDir = Vec3.zero()
        const rotObjUpDir = Vec3.zero()
        const rotObjSideDir = Vec3.zero()
        const rotMoveDir = Vec3.zero()

        function rotateCamera() {
            Vec3.set(rotMoveDir, _moveCurr[0] - _movePrev[0], _moveCurr[1] - _movePrev[1], 0);
            let angle = Vec3.magnitude(rotMoveDir);

            if (angle) {
                Vec3.copy(_eye, object.position)
                Vec3.sub(_eye, _eye, target)

                Vec3.normalize(rotEyeDir, Vec3.copy(rotEyeDir, _eye))
                Vec3.normalize(rotObjUpDir, Vec3.copy(rotObjUpDir, object.up))
                Vec3.normalize(rotObjSideDir, Vec3.cross(rotObjSideDir, rotObjUpDir, rotEyeDir))

                Vec3.setMagnitude(rotObjUpDir, rotObjUpDir, _moveCurr[1] - _movePrev[1])
                Vec3.setMagnitude(rotObjSideDir, rotObjSideDir, _moveCurr[0] - _movePrev[0])

                Vec3.add(rotMoveDir, Vec3.copy(rotMoveDir, rotObjUpDir), rotObjSideDir)

                Vec3.normalize(rotAxis, Vec3.cross(rotAxis, rotMoveDir, _eye))

                angle *= rotateSpeed;
                Quat.setAxisAngle(rotQuat, rotAxis, angle )

                Vec3.transformQuat(_eye, _eye, rotQuat)
                Vec3.transformQuat(object.up, object.up, rotQuat)

                Vec3.copy(_lastAxis, rotAxis)
                _lastAngle = angle;
            } else if (!staticMoving && _lastAngle) {
                _lastAngle *= Math.sqrt(1.0 - dynamicDampingFactor);
                Vec3.sub(_eye, Vec3.copy(_eye, object.position), target)
                Quat.setAxisAngle(rotQuat, _lastAxis, _lastAngle)

                Vec3.transformQuat(_eye, _eye, rotQuat)
                Vec3.transformQuat(object.up, object.up, rotQuat)
            }

            Vec2.copy(_movePrev, _moveCurr)
        }


        function zoomCamera () {
            // if (_state === STATE.TOUCH_ZOOM_PAN) {
            //     const factor = _touchZoomDistanceStart / _touchZoomDistanceEnd
            //     _touchZoomDistanceStart = _touchZoomDistanceEnd;
            //     Vec3.scale(_eye, _eye, factor)
            // } else {
            //     const factor = 1.0 + ( _zoomEnd[1] - _zoomStart[1] ) * zoomSpeed
            //     if (factor !== 1.0 && factor > 0.0) {
            //         Vec3.scale(_eye, _eye, factor)
            //     }

            //     if (staticMoving) {
            //         Vec2.copy(_zoomStart, _zoomEnd)
            //     } else {
            //         _zoomStart[1] += ( _zoomEnd[1] - _zoomStart[1] ) * dynamicDampingFactor
            //     }
            // }

            const factor = 1.0 + ( _zoomEnd[1] - _zoomStart[1] ) * zoomSpeed
            if (factor !== 1.0 && factor > 0.0) {
                Vec3.scale(_eye, _eye, factor)
            }

            if (staticMoving) {
                Vec2.copy(_zoomStart, _zoomEnd)
            } else {
                _zoomStart[1] += ( _zoomEnd[1] - _zoomStart[1] ) * dynamicDampingFactor
            }
        }

        const panMouseChange = Vec2.zero()
        const panObjUp = Vec3.zero()
        const panOffset = Vec3.zero()

        function panCamera() {
            Vec2.sub(panMouseChange, Vec2.copy(panMouseChange, _panEnd), _panStart)

            if (Vec2.squaredMagnitude(panMouseChange)) {
                Vec2.scale(panMouseChange, panMouseChange, Vec3.magnitude(_eye) * panSpeed)

                Vec3.cross(panOffset, Vec3.copy(panOffset, _eye), object.up)
                Vec3.setMagnitude(panOffset, panOffset, panMouseChange[0])

                Vec3.setMagnitude(panObjUp, object.up, panMouseChange[1])
                Vec3.add(panOffset, panOffset, panObjUp)

                Vec3.add(object.position, object.position, panOffset)
                Vec3.add(target, target, panOffset)

                if (staticMoving) {
                    Vec2.copy(_panStart, _panEnd)
                } else {
                    Vec2.sub(panMouseChange, _panEnd, _panStart)
                    Vec2.scale(panMouseChange, panMouseChange, dynamicDampingFactor)
                    Vec2.add(_panStart, _panStart, panMouseChange)
                }
            }
        }

        function checkDistances() {
            if (Vec3.squaredMagnitude(_eye) > maxDistance * maxDistance) {
                Vec3.setMagnitude(_eye, _eye, maxDistance)
                Vec3.add(object.position, target, _eye)
                Vec2.copy(_zoomStart, _zoomEnd)
            }

            if (Vec3.squaredMagnitude(_eye) < minDistance * minDistance) {
                Vec3.setMagnitude(_eye, _eye, minDistance)
                Vec3.add(object.position, target, _eye)
                Vec2.copy(_zoomStart, _zoomEnd)
            }
        }

        function update() {
            Vec3.sub( _eye, object.position, target)

            rotateCamera()
            zoomCamera()
            panCamera()

            Vec3.add(object.position, target, _eye)

            checkDistances()

            cameraLookAt(object.position, object.up, object.direction, target)

            if (Vec3.squaredDistance(lastPosition, object.position) > EPSILON.Value) {
                Vec3.copy(lastPosition, object.position)
            }
        }

        function reset() {
            Vec3.copy(target, target0)
            Vec3.copy(object.position, position0)
            Vec3.copy(object.up, up0)

            Vec3.sub(_eye, object.position, target)
            cameraLookAt(object.position, object.up, object.direction, target)
            Vec3.copy(lastPosition, object.position)
        }

        // listeners

        function onDrag({ pageX, pageY, buttons, modifiers, started }: DragInput) {
            if (started) {
                if (buttons === MouseButtonsFlag.Primary) {
                    Vec2.copy(_moveCurr, getMouseOnCircle(pageX, pageY))
                    Vec2.copy(_movePrev, _moveCurr)
                } else if (buttons === MouseButtonsFlag.Auxilary) {
                    Vec2.copy(_zoomStart, getMouseOnScreen(pageX, pageY))
                    Vec2.copy(_zoomEnd, _zoomStart)
                } else if (buttons === MouseButtonsFlag.Secondary) {
                    Vec2.copy(_panStart, getMouseOnScreen(pageX, pageY))
                    Vec2.copy(_panEnd, _panStart)
                }
            }

            if (buttons === MouseButtonsFlag.Primary) {
                Vec2.copy(_movePrev, _moveCurr)
                Vec2.copy(_moveCurr, getMouseOnCircle(pageX, pageY))
            } else if (buttons === MouseButtonsFlag.Auxilary) {
                Vec2.copy(_zoomEnd, getMouseOnScreen(pageX, pageY))
            } else if (buttons === MouseButtonsFlag.Secondary) {
                Vec2.copy(_panEnd, getMouseOnScreen(pageX, pageY))
            }
        }

        function onWheel({ dy }: WheelInput) {
            _zoomStart[1] -= dy
        }

        function onPinch({ delta }: PinchInput) {
            console.log(delta)
            _zoomStart[1] -= delta
        }

        // function touchstart(event: TouchEvent ) {
        //     switch ( event.touches.length ) {
        //         case 1:
        //             // _state = STATE.TOUCH_ROTATE;
        //             Vec2.copy(_moveCurr, getMouseOnCircle(event.touches[0].pageX, event.touches[0].pageY))
        //             Vec2.copy(_movePrev, _moveCurr)
        //             break;
        //         default: // 2 or more
        //             // _state = STATE.TOUCH_ZOOM_PAN;
        //             const dx = event.touches[0].pageX - event.touches[1].pageX;
        //             const dy = event.touches[0].pageY - event.touches[1].pageY;
        //             _touchZoomDistanceEnd = _touchZoomDistanceStart = Math.sqrt(dx * dx + dy * dy);

        //             const x = ( event.touches[0].pageX + event.touches[1].pageX) / 2;
        //             const y = ( event.touches[0].pageY + event.touches[1].pageY) / 2;
        //             Vec2.copy(_panStart, getMouseOnScreen(x, y))
        //             Vec2.copy(_panEnd, _panStart)
        //             break;
        //     }
        // }

        // function touchmove(event: TouchEvent) {
        //     event.preventDefault();
        //     event.stopPropagation();

        //     switch ( event.touches.length ) {
        //         case 1:
        //             Vec2.copy(_movePrev, _moveCurr)
        //             Vec2.copy(_moveCurr, getMouseOnCircle(event.touches[0].pageX, event.touches[0].pageY))
        //             break;
        //         default: // 2 or more
        //             const dx = event.touches[0].pageX - event.touches[1].pageX;
        //             const dy = event.touches[0].pageY - event.touches[1].pageY;
        //             _touchZoomDistanceEnd = Math.sqrt(dx * dx + dy * dy);

        //             const x = (event.touches[0].pageX + event.touches[1].pageX) / 2;
        //             const y = (event.touches[0].pageY + event.touches[1].pageY) / 2;
        //             Vec2.copy(_panEnd, getMouseOnScreen(x, y))
        //             break;
        //     }
        // }

        // function touchend(event: TouchEvent) {
        //     switch ( event.touches.length ) {
        //         case 0:
        //             // _state = STATE.NONE;
        //             break;
        //         case 1:
        //             // _state = STATE.TOUCH_ROTATE;
        //             Vec2.copy(_moveCurr, getMouseOnCircle(event.touches[0].pageX, event.touches[0].pageY))
        //             Vec2.copy(_movePrev, _moveCurr)
        //             break;

        //     }
        // }

        function dispose() {
            if (disposed) return
            disposed = true
            input.dispose()

            // element.removeEventListener( 'touchstart', touchstart as any, false );
            // element.removeEventListener( 'touchend', touchend as any, false );
            // element.removeEventListener( 'touchmove', touchmove as any, false );
        }

        // element.addEventListener( 'touchstart', touchstart as any, false );
        // element.addEventListener( 'touchend', touchend as any, false );
        // element.addEventListener( 'touchmove', touchmove as any, false );

        // force an update at start
        update();

        return {
            viewport,

            get dynamicDampingFactor() { return dynamicDampingFactor },
            set dynamicDampingFactor(value: number ) { dynamicDampingFactor = value },
            get rotateSpeed() { return rotateSpeed },
            set rotateSpeed(value: number ) { rotateSpeed = value },
            get zoomSpeed() { return zoomSpeed },
            set zoomSpeed(value: number ) { zoomSpeed = value },
            get panSpeed() { return panSpeed },
            set panSpeed(value: number ) { panSpeed = value },

            update,
            reset,
            dispose
        }
    }
}

export default TrackballControls