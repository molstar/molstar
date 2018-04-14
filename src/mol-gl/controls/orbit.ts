/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { map, filter, scan } from 'rxjs/operators';

import { Quat, Vec2, Vec3, EPSILON } from 'mol-math/linear-algebra';
import { clamp } from 'mol-math/interpolate';
import InputObserver from 'mol-util/input/input-observer';

const Y_UP = Vec3.create(0, 1, 0)
const tmpVec3 = Vec3.zero()

function cameraLookAt (direction: Vec3, up: Vec3, position: Vec3, target: Vec3) {
    Vec3.copy(direction, target)
    Vec3.sub(direction, direction, position)
    Vec3.normalize(direction, direction)
}

export const DefaultOrbitControlsProps = {
    parent: window as Window | Element,
    noScroll: true,

    phi: Math.PI / 2,
    theta: 0,

    position: Vec3.zero(),
    up: Vec3.create(0, 1, 0),
    target: Vec3.zero(),

    distance: undefined as (number|undefined),
    damping: 0.25,
    rotateSpeed: 0.28,
    zoomSpeed: 0.0075,
    pinchSpeed: 0.0075,
    translateSpeed: 1.0,
}
export type OrbitControlsProps = Partial<typeof DefaultOrbitControlsProps>

interface OrbitControls {
    update: () => void
    copyInto: (positionOut: Vec3, directionOut: Vec3, upOut: Vec3) => void

    position: Vec3
    direction: Vec3
    up: Vec3
    target: Vec3

    distance: number
    damping: number
    rotateSpeed: number
    zoomSpeed: number
    pinchSpeed: number
    translateSpeed: number

    phi: number
    theta: number
}

namespace OrbitControls {
    export function create (element: Element, props: OrbitControlsProps = {}): OrbitControls {
        const p = { ...DefaultOrbitControlsProps, ...props }

        const inputDelta = Vec3.zero() // x, y, zoom
        const offset = Vec3.zero()

        const upQuat = Quat.identity()
        const upQuatInverse = Quat.identity()
        const translateVec3 = Vec3.zero()

        const position = Vec3.clone(p.position)
        const direction = Vec3.zero()
        const up = Vec3.clone(p.up)
        const target = Vec3.clone(p.target)

        // const phiBounds = Vec2.create(0, Math.PI)
        const phiBounds = Vec2.create(-Infinity, Infinity)
        const thetaBounds = Vec2.create(-Infinity, Infinity)
        const distanceBounds = Vec2.create(0, Infinity)

        let { damping, rotateSpeed, zoomSpeed, pinchSpeed, translateSpeed, phi, theta } = p
        let distance = 0

        // Compute distance if not defined in user options
        if (p.distance === undefined) {
            Vec3.sub(tmpVec3, position, target)
            distance = Vec3.magnitude(tmpVec3)
        }

        const input = InputObserver.create(element, {
            parent: p.parent,
            noScroll: p.noScroll
        })
        input.drag.pipe(filter(v => v.buttons === 1)).subscribe(inputRotate)
        input.drag.pipe(filter(v => v.buttons === 4)).subscribe(inputTranslate)
        input.wheel.subscribe(inputZoom)
        input.pinch.subscribe(inputPinch)

        // Apply an initial phi and theta
        applyPhiTheta()

        return {
            update,
            copyInto,

            position,
            direction,
            up,
            target,

            get distance() { return distance },
            set distance(value: number ) { distance = value },
            get damping() { return damping },
            set damping(value: number ) { damping = value },
            get rotateSpeed() { return rotateSpeed },
            set rotateSpeed(value: number ) { rotateSpeed = value },
            get zoomSpeed() { return zoomSpeed },
            set zoomSpeed(value: number ) { zoomSpeed = value },
            get pinchSpeed() { return pinchSpeed },
            set pinchSpeed(value: number ) { pinchSpeed = value },
            get translateSpeed() { return translateSpeed },
            set translateSpeed(value: number ) { translateSpeed = value },

            get phi() { return phi },
            set phi(value: number ) { phi = value; applyPhiTheta() },
            get theta() { return theta },
            set theta(value: number ) { theta = value; applyPhiTheta() },
        }

        function copyInto(positionOut: Vec3, directionOut: Vec3, upOut: Vec3) {
            Vec3.copy(positionOut, position)
            Vec3.copy(directionOut, direction)
            Vec3.copy(upOut, up)
        }

        function inputRotate ({ dx, dy }: { dx: number, dy: number }) {
            const PI2 = Math.PI * 2
            inputDelta[0] -= PI2 * dx * rotateSpeed
            inputDelta[1] -= PI2 * dy * rotateSpeed
        }

        function inputZoom ({ dy }: { dy: number }) {
            inputDelta[2] += dy * zoomSpeed
        }

        function inputPinch (delta: number) {
            inputDelta[2] -= delta * pinchSpeed
        }

        function inputTranslate ({ dx, dy }: { dx: number, dy: number }) {
            // TODO
            console.log('translate', { dx, dy })
            const x = dx * translateSpeed * distance
            const y = dy * translateSpeed * distance
            // Vec3.set(translateVec3, x, y, 0)
            // Vec3.transformQuat(translateVec3, translateVec3, upQuat)

            // pan.copy( _eye ).cross( _this.object.up ).setLength( mouseChange.x );
            // pan.add( objectUp.copy( _this.object.up ).setLength( mouseChange.y ) );

            Vec3.copy(translateVec3, position)
            Vec3.cross(translateVec3, translateVec3, up)
            Vec3.normalize(translateVec3, translateVec3)
            Vec3.scale(translateVec3, translateVec3, x )

            const up2 = Vec3.clone(up)
            Vec3.normalize(up2, up2)
            Vec3.scale(up2, up2, y )
            Vec3.add(translateVec3, translateVec3, up2)

            Vec3.add(target, target, translateVec3)
            Vec3.add(position, position, translateVec3)
        }

        function updateDirection () {
            Quat.fromUnitVec3(upQuat, up, Y_UP)
            Quat.invert(upQuatInverse, upQuat)

            Vec3.sub(offset, position, target)
            Vec3.transformQuat(offset, offset, upQuat)

            let _distance = distance
            let _theta = Math.atan2(offset[0], offset[2])
            let _phi = Math.atan2(Math.sqrt(offset[0] * offset[0] + offset[2] * offset[2]), offset[1])

            _theta += inputDelta[0]
            _phi += inputDelta[1]

            _theta = clamp(_theta, thetaBounds[0], thetaBounds[1])
            _phi = clamp(_phi, phiBounds[0], phiBounds[1])
            _phi = clamp(_phi, EPSILON.Value, Math.PI - EPSILON.Value)

            _distance += inputDelta[2]
            _distance = clamp(_distance, distanceBounds[0], distanceBounds[1])

            const radius = Math.abs(_distance) <= EPSILON.Value ? EPSILON.Value : _distance
            offset[0] = radius * Math.sin(_phi) * Math.sin(_theta)
            offset[1] = radius * Math.cos(_phi)
            offset[2] = radius * Math.sin(_phi) * Math.cos(_theta)

            phi = _phi
            theta = _theta
            distance = _distance

            Vec3.transformQuat(offset, offset, upQuatInverse)
            Vec3.add(position, target, offset)
            cameraLookAt(direction, up, position, target)
        }

        function update () {
            updateDirection()
            for (let i = 0; i < inputDelta.length; i++) {
                inputDelta[i] *= 1 - damping
            }
        }

        function applyPhiTheta () {
            let _phi = phi
            let _theta = theta
            _theta = clamp(_theta, thetaBounds[0], thetaBounds[1])
            _phi = clamp(_phi, phiBounds[0], phiBounds[1])
            _phi = clamp(_phi, EPSILON.Value, Math.PI - EPSILON.Value)

            const dist = Math.max(EPSILON.Value, distance)
            position[0] = dist * Math.sin(_phi) * Math.sin(_theta)
            position[1] = dist * Math.cos(_phi)
            position[2] = dist * Math.sin(_phi) * Math.cos(_theta)
            Vec3.add(position, position, target)

            updateDirection()
        }
    }
}

export default OrbitControls