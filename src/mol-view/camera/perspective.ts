/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Vec3 } from 'mol-math/linear-algebra'
import { DefaultCameraProps, Camera } from './base'

export interface PerspectiveCamera extends Camera {
    fov: number,
    near: number,
    far: number
}

export const DefaultPerspectiveCameraProps = {
    fov: Math.PI / 4,
    near: 0.1,
    far: 10000,
    ...DefaultCameraProps
}
export type PerspectiveCameraProps = Partial<typeof DefaultPerspectiveCameraProps>

export namespace PerspectiveCamera {
    export function create(props: PerspectiveCameraProps = {}): PerspectiveCamera {
        let { fov, near, far } = { ...DefaultPerspectiveCameraProps, ...props };

        const camera = Camera.create(props)
        const center = Vec3.zero()

        function update () {
            const aspect = camera.viewport.width / camera.viewport.height

            // build projection matrix
            Mat4.perspective(camera.projection, fov, aspect, Math.abs(near), Math.abs(far))

            // build view matrix
            Vec3.add(center, camera.position, camera.direction)
            Mat4.lookAt(camera.view, camera.position, center, camera.up)

            // update projection * view and invert
            camera.update()
        }

        update()

        return {
            ...camera,
            update,

            get far() { return far },
            set far(value: number) { far = value },

            get near() { return near },
            set near(value: number) { near = value },

            get fov() { return fov },
            set fov(value: number) { fov = value },
        }
    }
}