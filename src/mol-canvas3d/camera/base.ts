/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Vec3, Vec4 } from 'mol-math/linear-algebra'
import { cameraProject, cameraUnproject, cameraLookAt, Viewport } from './util';
import { Object3D } from 'mol-gl/object3d';

export interface Camera extends Object3D {
    readonly projection: Mat4,
    readonly projectionView: Mat4,
    readonly inverseProjectionView: Mat4,
    readonly viewport: Viewport,

    near: number,
    far: number,
    fogNear: number,
    fogFar: number,
}

export const DefaultCameraProps = {
    position: Vec3.zero(),
    direction: Vec3.create(0, 0, -1),
    up: Vec3.create(0, 1, 0),
    viewport: Viewport.create(-1, -1, 1, 1),
    target: Vec3.create(0, 0, 0),

    near: 0.1,
    far: 10000,
    fogNear: 0.1,
    fogFar: 10000,
}
export type CameraProps = typeof DefaultCameraProps

export namespace Camera {
    export function create(props?: Partial<CameraProps>): Camera {
        const p = { ...DefaultCameraProps, ...props };

        const { view, position, direction, up } = Object3D.create()
        Vec3.copy(position, p.position)
        Vec3.copy(direction, p.direction)
        Vec3.copy(up, p.up)

        const projection = Mat4.identity()
        const viewport = Viewport.clone(p.viewport)
        const projectionView = Mat4.identity()
        const inverseProjectionView = Mat4.identity()

        return {
            projection,
            projectionView,
            inverseProjectionView,
            viewport,

            view,
            position,
            direction,
            up,

            near: p.near,
            far: p.far,
            fogNear: p.fogNear,
            fogFar: p.fogFar,
        }
    }

    export function update (camera: Camera) {
        Mat4.mul(camera.projectionView, camera.projection, camera.view)
        Mat4.invert(camera.inverseProjectionView, camera.projectionView)
        return camera
    }

    export function lookAt (camera: Camera, target: Vec3) {
        cameraLookAt(camera.direction, camera.up, camera.position, target)
    }

    export function reset (camera: Camera, props: CameraProps) {
        Vec3.copy(camera.position, props.position)
        Vec3.copy(camera.direction, props.direction)
        Vec3.copy(camera.up, props.up)
        Mat4.setIdentity(camera.view)
        Mat4.setIdentity(camera.projection)
        Mat4.setIdentity(camera.projectionView)
        Mat4.setIdentity(camera.inverseProjectionView)
    }

    export function translate (camera: Camera, v: Vec3) {
        Vec3.add(camera.position, camera.position, v)
    }

    export function project (camera: Camera, out: Vec4, point: Vec3) {
        return cameraProject(out, point, camera.viewport, camera.projectionView)
    }

    export function unproject (camera: Camera, out: Vec3, point: Vec3) {
        return cameraUnproject(out, point, camera.viewport, camera.inverseProjectionView)
    }
}