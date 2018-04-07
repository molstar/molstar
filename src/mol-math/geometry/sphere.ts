/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '../linear-algebra'

export interface Sphere {
    center: Vec3
    radius: number
}

export namespace Sphere {
    export function create(center: Vec3, radius: number): Sphere {
        return { center, radius }
    }
}

export default Sphere