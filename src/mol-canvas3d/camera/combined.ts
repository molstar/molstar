/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PerspectiveCamera,  } from './perspective';
import { OrthographicCamera } from './orthographic';
import { Camera, DefaultCameraProps } from './base';
import { Vec3 } from 'mol-math/linear-algebra';

export type CombinedCameraMode = 'perspective' | 'orthographic'

export interface CombinedCamera extends Camera {
    target: Vec3
    fov: number
    zoom: number
    mode: CombinedCameraMode
}

export const DefaultCombinedCameraProps = {
    ...DefaultCameraProps,
    target: Vec3.zero(),
    fov: Math.PI / 4,
    zoom: 1,
    mode: 'perspective' as CombinedCameraMode
}
export type CombinedCameraProps = Partial<typeof DefaultCombinedCameraProps>

export namespace CombinedCamera {
    export function create(props: CombinedCameraProps = {}): CombinedCamera {
        const { zoom, fov, mode, target: t } = { ...DefaultCombinedCameraProps, ...props };
        const target = Vec3.create(t[0], t[1], t[2])
        const camera = { ...Camera.create(props), zoom, fov, mode, target }
        update(camera)

        return camera
    }

    export function update(camera: CombinedCamera) {
        const height = 2 * Math.tan(camera.fov / 2) * Vec3.distance(camera.position, camera.target)
        camera.zoom = camera.viewport.height / height

        switch (camera.mode) {
            case 'orthographic': OrthographicCamera.update(camera); break
            case 'perspective': PerspectiveCamera.update(camera); break
        }
    }
}