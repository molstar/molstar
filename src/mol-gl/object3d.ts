/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3, Mat4 } from '../mol-math/linear-algebra';

export interface Object3D {
    readonly view: Mat4
    readonly position: Vec3
    readonly direction: Vec3
    readonly up: Vec3
}

export namespace Object3D {
    export function create(): Object3D {
        return {
            view: Mat4.identity(),
            position: Vec3.create(0, 0, 0),
            direction: Vec3.create(0, 0, -1),
            up: Vec3.create(0, 1, 0),
        };
    }
}