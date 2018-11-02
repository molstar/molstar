/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Vec3 } from 'mol-math/linear-algebra'
import { DefaultCameraProps, Camera } from './base'

export interface OrthographicCamera extends Camera {
    zoom: number
}

export const DefaultOrthographicCameraProps = {
    ...DefaultCameraProps,
    zoom: 1
}
export type OrthographicCameraProps = Partial<typeof DefaultOrthographicCameraProps>

export namespace OrthographicCamera {
    export function create(props: OrthographicCameraProps = {}): OrthographicCamera {
        const { zoom } = { ...DefaultOrthographicCameraProps, ...props };
        const camera = { ...Camera.create(props), zoom }
        update(camera)

        return camera
    }

    const center = Vec3.zero()
    export function update(camera: OrthographicCamera) {
        const { viewport, zoom } = camera

        const fullLeft = (viewport.width - viewport.x) / -2
        const fullRight = (viewport.width - viewport.x) / 2
        const fullTop = (viewport.height - viewport.y) / 2
        const fullBottom = (viewport.height - viewport.y) / -2

        const dx = (fullRight - fullLeft) / (2 * zoom)
        const dy = (fullTop - fullBottom) / (2 * zoom)
        const cx = (fullRight + fullLeft) / 2
        const cy = (fullTop + fullBottom) / 2

        const left = cx - dx
        const right = cx + dx
        const top = cy + dy
        const bottom = cy - dy

        // build projection matrix
        Mat4.ortho(camera.projection, left, right, bottom, top, Math.abs(camera.near), Math.abs(camera.far))

        // build view matrix
        Vec3.add(center, camera.position, camera.direction)
        Mat4.lookAt(camera.view, camera.position, center, camera.up)

        // update projection * view and invert
        Camera.update(camera)
    }
}