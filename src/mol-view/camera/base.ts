/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Vec3, Vec4 } from 'mol-math/linear-algebra'
import { cameraProject, cameraUnproject, cameraLookAt, Viewport } from './util';
import { Object3D, createObject3D } from 'mol-gl/object3d';

export interface Camera extends Object3D {
    readonly projection: Mat4,
    readonly projectionView: Mat4,
    readonly inverseProjectionView: Mat4,
    readonly viewport: Viewport,

    near: number,
    far: number,
    fogNear: number,
    fogFar: number,

    translate: (v: Vec3) => void,
    reset: () => void,
    lookAt: (target: Vec3) => void,
    project: (out: Vec4, point: Vec3) => Vec4,
    unproject: (out: Vec3, point: Vec3) => Vec3
}

export const DefaultCameraProps = {
    position: Vec3.zero(),
    direction: Vec3.create(0, 0, -1),
    up: Vec3.create(0, 1, 0),
    viewport: Viewport.create(-1, -1, 1, 1),
    near: 0.1,
    far: 10000,
    fogNear: 0.1,
    fogFar: 10000,
}
export type CameraProps = Partial<typeof DefaultCameraProps>

export namespace Camera {
    export function create(props?: CameraProps): Camera {
        const p = { ...DefaultCameraProps, ...props };

        const { view, position, direction, up } = createObject3D()
        Vec3.copy(position, p.position)
        Vec3.copy(direction, p.direction)
        Vec3.copy(up, p.up)

        const projection = Mat4.identity()
        const viewport = Viewport.clone(p.viewport)
        const projectionView = Mat4.identity()
        const inverseProjectionView = Mat4.identity()

        function update () {
            Mat4.mul(projectionView, projection, view)
            Mat4.invert(inverseProjectionView, projectionView)
        }

        function lookAt (target: Vec3) {
            cameraLookAt(direction, up, position, target)
        }

        function reset () {
            Vec3.copy(position, p.position)
            Vec3.copy(direction, p.direction)
            Vec3.copy(up, p.up)
            Mat4.setIdentity(view)
            Mat4.setIdentity(projection)
            Mat4.setIdentity(projectionView)
            Mat4.setIdentity(inverseProjectionView)
        }

        function translate (v: Vec3) {
            Vec3.add(position, position, v)
        }

        function project (out: Vec4, point: Vec3) {
            return cameraProject(out, point, viewport, projectionView)
        }

        function unproject (out: Vec3, point: Vec3) {
            return cameraUnproject(out, point, viewport, inverseProjectionView)
        }

        return {
            view,
            projection,
            projectionView,
            inverseProjectionView,

            viewport,
            position,
            direction,
            up,

            get near() { return p.near },
            set near(value: number) { p.near = value },
            get far() { return p.far },
            set far(value: number) { p.far = value },
            get fogNear() { return p.fogNear },
            set fogNear(value: number) { p.fogNear = value },
            get fogFar() { return p.fogFar },
            set fogFar(value: number) { p.fogFar = value },

            translate,
            reset,
            lookAt,
            update,
            project,
            unproject
        }
    }
}