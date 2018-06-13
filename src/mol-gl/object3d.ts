/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3, Mat4 } from 'mol-math/linear-algebra';

export interface Object3D {
    readonly view: Mat4,
    readonly position: Vec3,
    readonly direction: Vec3,
    readonly up: Vec3,

    update: () => void
}

export function createObject3D(): Object3D {
    const view = Mat4.identity()
    const position = Vec3.create(0, 0, 0)
    const direction = Vec3.create(0, 0, -1)
    const up = Vec3.create(0, 1, 0)

    const center = Vec3.zero()

    return {
        view,
        position,
        direction,
        up,

        update() {
            // console.log('position', position)
            // console.log('direction', direction)
            // console.log('up', up)
            Vec3.add(center, position, direction)
            Mat4.lookAt(view, position, center, up)
            // Mat4.lookAt(view, center, position, up)
            // console.log('view', view)
        }
    }
}