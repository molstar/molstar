/**
 * Copyright (c) 2020-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec2, Vec3 } from '../../mol-math/linear-algebra';
import { InputObserver } from '../../mol-util/input/input-observer';
import { Camera } from '../camera';
import { Viewport } from '../camera/util';

export namespace ObjectControls {
    function mouseOnScreen(out: Vec2, page: Vec2, viewport: Viewport) {
        return Vec2.set(
            out,
            (page[0] - viewport.x) / viewport.width,
            (page[1] - viewport.y) / viewport.height
        );
    }

    const panMouseChange = Vec2();
    const panObjUp = Vec3();
    const panOffset = Vec3();
    const eye = Vec3();
    const panStart = Vec2();
    const panEnd = Vec2();

    const target = Vec3();

    /**
     * Get vector for movement in camera projection plane:
     * `pageStart` and `pageEnd` are 2d window coordinates;
     * `ref` defines the plane depth, if not given `camera.target` is used
     */
    export function panDirection(out: Vec3, pageStart: Vec2, pageEnd: Vec2, input: InputObserver, camera: Camera, ref?: Vec3) {
        mouseOnScreen(panStart, pageStart, camera.viewport);
        mouseOnScreen(panEnd, pageEnd, camera.viewport);
        Vec2.sub(panMouseChange, Vec2.copy(panMouseChange, panEnd), panStart);
        Vec3.sub(eye, camera.position, camera.target);

        if (!ref || camera.state.mode === 'orthographic') Vec3.copy(target, camera.target);
        else Vec3.projectPointOnVector(target, ref, eye, camera.position);

        const dist = Vec3.distance(camera.position, target);
        const height = 2 * Math.tan(camera.state.fov / 2) * dist;
        const zoom = camera.viewport.height / height;

        panMouseChange[0] *= (1 / zoom) * camera.viewport.width * input.pixelRatio;
        panMouseChange[1] *= (1 / zoom) * camera.viewport.height * input.pixelRatio;

        Vec3.cross(panOffset, Vec3.copy(panOffset, eye), camera.up);
        Vec3.setMagnitude(panOffset, panOffset, panMouseChange[0]);

        Vec3.setMagnitude(panObjUp, camera.up, panMouseChange[1]);
        Vec3.add(panOffset, panOffset, panObjUp);

        return Vec3.negate(out, panOffset);
    }
}