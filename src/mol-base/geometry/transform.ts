/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Mat4 } from '../math/linear-algebra-3d'

interface Transform extends Readonly<{
    transform: Mat4,
    // cache the inverse of the transform
    inverse: Mat4,
    // optimize the identity case
    isIdentity: boolean
}> { }

namespace Transform {
    export function isIdentity(m: Mat4) {
        throw 'nyi'
    }

    export function isRotationAndTranslation(m: Mat4) {
        throw 'nyi'
    }
}

export default Transform