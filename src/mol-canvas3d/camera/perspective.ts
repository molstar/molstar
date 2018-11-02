/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Vec3 } from 'mol-math/linear-algebra'
import { DefaultCameraProps, Camera } from './base'

export interface PerspectiveCamera extends Camera {
    fov: number
}

export const DefaultPerspectiveCameraProps = {
    ...DefaultCameraProps,
    fov: Math.PI / 4,
}
export type PerspectiveCameraProps = Partial<typeof DefaultPerspectiveCameraProps>

export namespace PerspectiveCamera {
    export function create(props: PerspectiveCameraProps = {}): PerspectiveCamera {
        const { fov } = { ...DefaultPerspectiveCameraProps, ...props }
        const camera = { ...Camera.create(props), fov }
        update(camera)

        return camera
    }

    const center = Vec3.zero()
    export function update(camera: PerspectiveCamera) {
        const aspect = camera.viewport.width / camera.viewport.height

        // build projection matrix
        Mat4.perspective(camera.projection, camera.fov, aspect, Math.abs(camera.near), Math.abs(camera.far))

        // build view matrix
        Vec3.add(center, camera.position, camera.direction)
        Mat4.lookAt(camera.view, camera.position, center, camera.up)

        // update projection * view and invert
        Camera.update(camera)
    }
}